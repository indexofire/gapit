"""Hit processing: the SPEC.md §4 algorithm, in order, nothing more.

No interval merging of any kind — upstream reports overlapping genes at
different query spans, and so do we.
"""

import re
from collections.abc import Iterable
from typing import TYPE_CHECKING, Literal

from pydantic import BaseModel

from gapit.db import IDSEP, parse_db_header
from gapit.minimap import minimap

# Perl: $product =~ s/^\S+\s+// if $product =~ m/~~~/  — strips the leading
# makeblastdb id token only when followed by whitespace; a bare ~~~id stays.
_LEADING_TOKEN_RE = re.compile(r"^\S+\s+")

if TYPE_CHECKING:
    # TYPE_CHECKING-only: blast -> report -> hits would otherwise be a cycle.
    from gapit.blast import BlastRow


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
    resistance: str
    s_start: int
    s_end: int
    s_len: int
    coverage_map: str
    gap_openings: int
    gaps: int
    identity_pct: float
    coverage_pct: float


def process_rows(rows: Iterable["BlastRow"], *, mincov: float, default_db: str) -> list[Hit]:
    """Turn BLAST rows into hits, strictly in SPEC.md §4 order:

    1. minus-strand swap (subject coords only), 2. dedup on
    (qseqid, qstart, qend) — first row wins, key ignores strand, and the key
    is claimed even if the row is later coverage-filtered, 3. coverage filter
    on the unrounded float, 4. ~~~ header parse, 5. product cleanup.
    """
    hits: list[Hit] = []
    seen: set[tuple[str, int, int]] = set()
    for row in rows:
        # 1. minus-strand normalize: swap sstart/send (query coords untouched)
        if row.sstrand == "minus":
            s_start, s_end = row.send, row.sstart
        else:
            s_start, s_end = row.sstart, row.send
        # 2. dedup on the query span
        key = (row.qseqid, row.qstart, row.qend)
        if key in seen:
            continue
        seen.add(key)
        # 3. coverage filter on the unrounded float
        coverage_pct = 100.0 * (row.length - row.gaps) / row.slen
        if coverage_pct < mincov:
            continue
        # 4. subject id parse with fallbacks
        header = parse_db_header(row.sseqid, default_db)
        # 5. product cleanup: n/a fallback, strip ',' and tab, drop leading id token
        product = row.stitle or "n/a"
        product = product.replace(",", "").replace("\t", "")
        if IDSEP in product:
            product = _LEADING_TOKEN_RE.sub("", product, count=1)
        hits.append(
            Hit(
                sequence=row.qseqid,
                start=row.qstart,
                end=row.qend,
                strand="-" if row.sstrand == "minus" else "+",
                gene=header.gene,
                database=header.database,
                accession=header.accession,
                resistance=header.resistance,
                product=product,
                s_start=s_start,
                s_end=s_end,
                s_len=row.slen,
                coverage_map=minimap(s_start, s_end, row.slen, row.gapopen),
                gap_openings=row.gapopen,
                gaps=row.gaps,
                identity_pct=row.pident,
                coverage_pct=coverage_pct,
            )
        )
    return hits
