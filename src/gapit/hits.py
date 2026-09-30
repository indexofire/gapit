"""Hit processing: the SPEC.md §4 algorithm, in order, nothing more.

No interval merging of any kind — upstream reports overlapping genes at
different query spans, and so do we. Subject ids decode through
:mod:`gapit.dbcodec`: legacy ``~~~`` ids parse by the frozen db.py rules,
``gapit|``-tagged ids by the strict native codec — a malformed native header
raises DatabaseError ``HEADER_MALFORMED`` (exit 4), never a silent fallback.
"""

import re
from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING, Literal

from pydantic import BaseModel

from gapit.db import IDSEP
from gapit.dbcodec import decode_seqid, is_gapit_header
from gapit.minimap import minimap

# Perl: $product =~ s/^\S+\s+// if $product =~ m/~~~/  — strips the leading
# makeblastdb id token only when followed by whitespace; a bare ~~~id stays.
# Native gapit| ids get the same substitution: makeblastdb prefixes stitle
# with the id, and there the prefix is the whitespace-free tagged seqid.
_LEADING_TOKEN_RE = re.compile(r"^\S+\s+")

if TYPE_CHECKING:
    # TYPE_CHECKING-only: blast -> report -> hits would otherwise be a cycle.
    from gapit.blast import BlastRow


class Fragment(BaseModel, frozen=True):
    """One contributing fragment of a merged hit (``--merge-fragments``)."""

    contig: str
    start: int
    end: int
    strand: Literal["+", "-"]
    identity_pct: float
    coverage_pct: float


class Hit(BaseModel, frozen=True):
    """One surviving hit — raw semantic values only (display formatting is Phase 3)."""

    sequence: str
    start: int
    end: int
    strand: Literal["+", "-"]
    gene: str
    database: str
    accession: str
    product: str
    function: str
    s_start: int
    s_end: int
    s_len: int
    coverage_map: str
    gap_openings: int
    gaps: int
    identity_pct: float
    coverage_pct: float
    # Ungapped aligned columns (BLAST length - gaps); weight for fragment
    # identity means and anchor selection (internal, never rendered).
    aligned_len: int
    # --merge-fragments extension: merged hits carry their contributing
    # fragments; default-path hits never do (parity output is unchanged).
    merged: bool = False
    fragments: tuple[Fragment, ...] = ()


def _swap_dedup(rows: Iterable["BlastRow"]) -> Iterator[tuple["BlastRow", int, int]]:
    """SPEC.md §4 steps 1-2: minus-strand swap (subject coords only) and dedup
    on ``(qseqid, qstart, qend)`` — first row wins, the key ignores strand,
    and the key is claimed even if the row is later coverage-filtered."""
    seen: set[tuple[str, int, int]] = set()
    for row in rows:
        if row.sstrand == "minus":
            s_start, s_end = row.send, row.sstart
        else:
            s_start, s_end = row.sstart, row.send
        key = (row.qseqid, row.qstart, row.qend)
        if key in seen:
            continue
        seen.add(key)
        yield row, s_start, s_end


def _hit(row: "BlastRow", s_start: int, s_end: int, coverage_pct: float, default_db: str) -> Hit:
    """SPEC.md §4 steps 4-5: subject-id decode (native ``gapit|`` codec or
    legacy ``~~~`` rules) and product cleanup, assembled into a Hit."""
    header = decode_seqid(row.sseqid, default_db)
    product = row.stitle or "n/a"
    product = product.replace(",", "").replace("\t", "")
    if IDSEP in product or is_gapit_header(row.sseqid):
        product = _LEADING_TOKEN_RE.sub("", product, count=1)
    return Hit(
        sequence=row.qseqid,
        start=row.qstart,
        end=row.qend,
        strand="-" if row.sstrand == "minus" else "+",
        gene=header.gene,
        database=header.database,
        accession=header.accession,
        function=header.function,
        product=product,
        s_start=s_start,
        s_end=s_end,
        s_len=row.slen,
        coverage_map=minimap(s_start, s_end, row.slen, row.gapopen),
        gap_openings=row.gapopen,
        gaps=row.gaps,
        identity_pct=row.pident,
        coverage_pct=coverage_pct,
        aligned_len=row.length - row.gaps,
    )


def process_rows(rows: Iterable["BlastRow"], *, mincov: float, default_db: str) -> list[Hit]:
    """Turn BLAST rows into hits, strictly in SPEC.md §4 order:

    1. minus-strand swap (subject coords only), 2. dedup on
    (qseqid, qstart, qend) — first row wins, key ignores strand, and the key
    is claimed even if the row is later coverage-filtered, 3. coverage filter
    on the unrounded float, 4. header decode (native ``gapit|`` codec or
    legacy ``~~~`` rules), 5. product cleanup.
    """
    hits: list[Hit] = []
    for row, s_start, s_end in _swap_dedup(rows):
        # 3. coverage filter on the unrounded float (decode stays after the
        # filter: a malformed native header on a filtered row never raises)
        coverage_pct = 100.0 * (row.length - row.gaps) / row.slen
        if coverage_pct < mincov:
            continue
        hits.append(_hit(row, s_start, s_end, coverage_pct, default_db))
    return hits


def process_rows_with_partials(
    rows: Iterable["BlastRow"], *, mincov: float, default_db: str
) -> tuple[list[Hit], list[Hit]]:
    """Merge-mode variant of :func:`process_rows`: returns ``(hits, partials)``
    where ``partials`` are the sub-mincov dedup survivors, decoded so the
    fragment-merge pass (:mod:`gapit.fragments`) can reason about them.

    Opt-in ``--merge-fragments`` path only — every dedup survivor is decoded
    (step 4 runs before the split), unlike the parity path. ``process_rows``
    stays the default and byte-identical to abricate.
    """
    hits: list[Hit] = []
    partials: list[Hit] = []
    for row, s_start, s_end in _swap_dedup(rows):
        coverage_pct = 100.0 * (row.length - row.gaps) / row.slen
        hit = _hit(row, s_start, s_end, coverage_pct, default_db)
        (hits if coverage_pct >= mincov else partials).append(hit)
    return hits, partials
