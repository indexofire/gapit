"""Cross-contig gene-fragment merging — the opt-in ``--merge-fragments`` path.

Draft assemblies split genes across contig boundaries; each partial hit fails
mincov individually, so the gene goes unreported. This module merges the
sub-mincov dedup survivors of one gene whose UNION subject coverage reaches
mincov into a single reported hit. Pure functions over :class:`gapit.hits.Hit`
lists; the default (flag-off) path never calls in here — abricate parity is
untouched (SPEC.md §4 explicitly forbids merging on the default path).
"""

from collections.abc import Iterable

from gapit.hits import Fragment, Hit

Span = tuple[int, int]


def union_spans(spans: Iterable[Span]) -> list[Span]:
    """Merge overlapping 1-based-inclusive intervals; abutting pairs stay
    separate (union length is identical either way). Output is start-sorted."""
    merged: list[Span] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            previous = merged[-1]
            merged[-1] = (previous[0], max(previous[1], end))
        else:
            merged.append((start, end))
    return merged


def union_coverage_pct(spans: Iterable[Span], slen: int) -> float:
    """``100 * union_length / slen`` on the unrounded float — the same shape
    as the SPEC.md §4 coverage formula, over the union of subject intervals."""
    return 100.0 * sum(end - start + 1 for start, end in union_spans(spans)) / slen


def union_coverage_map(spans: Iterable[Span], slen: int, *, broken: bool) -> str:
    """15-char coverage bar over the interval union, mirroring minimap.py's
    binning arithmetic exactly (``int(x/scale)`` truncation, 1-based quirks,
    and the 14-box + '/' broken form when any fragment opened a gap)."""
    width = 15 - (1 if broken else 0)
    scale = slen / width
    covered: set[int] = set()
    for start, end in union_spans(spans):
        covered.update(range(int(start / scale), int(end / scale) + 1))
    chars: list[str] = []
    for i in range(width):
        chars.append("=" if i in covered else ".")
        if broken and i == width // 2:
            chars.append("/")
    return "".join(chars)


def _merged_hit(fragments: list[Hit], *, mincov: float) -> Hit | None:
    """Assemble one merged hit from >= 2 same-gene partials, or None.

    Field rules (the gapit extension contract):
    - coverage_pct: union subject coverage, unrounded; kept iff >= mincov.
    - identity_pct: aligned-length-weighted mean of fragment identities.
    - coverage_map: minimap binning over the union of subject intervals.
    - gaps/gap_openings: summed totals.
    - sequence: comma-joined contig names, sorted (tab-safe TSV cell).
    - start/end/strand and the header fields: taken from the ANCHOR — the
      fragment with the largest aligned length (ties: first in stable
      (contig, start) order).
    - COVERAGE column span: the union's bounding [min(s_start), max(s_end)].
    """
    if len(fragments) < 2:
        return None
    ordered = sorted(fragments, key=lambda hit: (hit.sequence, hit.start))
    spans = [(hit.s_start, hit.s_end) for hit in ordered]
    coverage_pct = union_coverage_pct(spans, ordered[0].s_len)
    if coverage_pct < mincov:
        return None
    anchor = ordered[0]
    for hit in ordered[1:]:
        if hit.aligned_len > anchor.aligned_len:
            anchor = hit
    gap_openings = sum(hit.gap_openings for hit in ordered)
    return Hit(
        sequence=",".join(sorted({hit.sequence for hit in ordered})),
        start=anchor.start,
        end=anchor.end,
        strand=anchor.strand,
        gene=anchor.gene,
        database=anchor.database,
        accession=anchor.accession,
        product=anchor.product,
        function=anchor.function,
        s_start=min(hit.s_start for hit in ordered),
        s_end=max(hit.s_end for hit in ordered),
        s_len=anchor.s_len,
        coverage_map=union_coverage_map(spans, anchor.s_len, broken=gap_openings > 0),
        gap_openings=gap_openings,
        gaps=sum(hit.gaps for hit in ordered),
        identity_pct=sum(hit.identity_pct * hit.aligned_len for hit in ordered)
        / sum(hit.aligned_len for hit in ordered),
        coverage_pct=coverage_pct,
        aligned_len=sum(hit.aligned_len for hit in ordered),
        merged=True,
        fragments=tuple(
            Fragment(
                contig=hit.sequence,
                start=hit.start,
                end=hit.end,
                strand=hit.strand,
                identity_pct=hit.identity_pct,
                coverage_pct=hit.coverage_pct,
            )
            for hit in ordered
        ),
    )


def merge_fragment_hits(hits: list[Hit], partials: list[Hit], *, mincov: float) -> list[Hit]:
    """Fold mergeable partial genes into ``hits`` (opt-in extension).

    ``hits`` are the >= mincov survivors (returned untouched, in order);
    ``partials`` are the sub-mincov dedup survivors. A gene that already has
    a hit is never merged (its stray partials are ignored) — merging only
    rescues genes that would otherwise go unreported. Merged hits are
    appended after ``hits``; the caller's (sequence, start) report sort
    produces the final order.
    """
    reported_genes = {hit.gene for hit in hits}
    groups: dict[str, list[Hit]] = {}
    for partial in partials:
        if partial.gene not in reported_genes:
            groups.setdefault(partial.gene, []).append(partial)
    merged = (_merged_hit(group, mincov=mincov) for group in groups.values())
    return hits + [hit for hit in merged if hit is not None]
