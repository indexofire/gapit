"""gapit.cluster/1 Markdown renderer.

Split from formats/cluster.py at the 250-LOC ceiling (stage 3): the per-file
summary row (with a ``Phenotype`` column on typed databases), the best
locus's gene table, and YAML frontmatter mirroring the engine params. The
frontmatter is STATIC metadata only and each file's section carries its own
summary row + gene table, so the document streams one chunk per file exactly
like the TSV — no all-files summary table, no run totals.
"""

from collections.abc import Iterable
from datetime import UTC, datetime

from gapit import __version__
from gapit.cluster import ClusterParams, ClusterReport, LocusCall
from gapit.formats.md import md_cell

_SUMMARY_COLUMNS = (
    "File",
    "Best locus",
    "Type",
    "Coverage%",
    "Identity%",
    "Present",
    "Partial",
    "Missing",
)
_SUMMARY_COLUMNS_TYPED = (
    "File",
    "Best locus",
    "Type",
    "Phenotype",
    "Coverage%",
    "Identity%",
    "Present",
    "Partial",
    "Missing",
)
_GENE_COLUMNS = ("Gene", "Start", "End", "Strand", "Coverage%", "Identity%", "Verdict")


def _md_table(columns: tuple[str, ...], rows: Iterable[tuple[str, ...]]) -> list[str]:
    """Header + separator + one escaped row per tuple (md.py table shape)."""
    table = ["| " + " | ".join(columns) + " |", "|" + "---|" * len(columns)]
    table += ["| " + " | ".join(md_cell(cell) for cell in row) + " |" for row in rows]
    return table


def _phenotype_cell(report: ClusterReport) -> str:
    """The called phenotype, ``-`` when untyped/uncalled/ambiguous."""
    if report.best is None or report.best.phenotype is None:
        return "-"
    return report.best.phenotype


def _summary_row(report: ClusterReport, typed: bool) -> tuple[str, ...]:
    """One file's summary row; a null best renders the untyped dashes."""
    if report.best is None:
        return (report.file, *("-",) * (8 if typed else 7))
    best = report.best
    phenotype = (_phenotype_cell(report),) if typed else ()
    return (
        report.file,
        best.locus,
        best.type,
        *phenotype,
        f"{best.coverage_pct:.2f}",
        f"{best.identity_pct:.2f}",
        str(best.genes_present),
        str(best.genes_partial),
        ";".join(
            gene_id
            for locus in report.loci
            if locus.locus == best.locus
            for gene_id in locus.missing
        )
        or "-",
    )


def _gene_rows(top: LocusCall) -> list[tuple[str, ...]]:
    return [
        (
            gene.gene_id,
            str(gene.start),
            str(gene.end),
            gene.strand,
            f"{gene.coverage_pct:.2f}",
            f"{gene.identity_pct:.2f}",
            gene.verdict,
        )
        for gene in top.genes
    ]


def cluster_md_head(params: ClusterParams, *, now: datetime) -> str:
    """The STATIC frontmatter + title chunk (schema, tool, timestamp,
    thresholds — everything knowable before file 1), so a streaming caller
    emits it first; the per-file sections follow one chunk per file."""
    lines: list[str] = [
        "---",
        "schema: gapit.cluster/1",
        f"tool: gapit {__version__}",
        f"created_at: {now.astimezone(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"db: {params.db}",
        f"preset: {params.preset}",
        f"min_gene_cov: {params.min_gene_cov}",
        f"min_gene_id: {params.min_gene_id}",
        f"min_cluster_cov: {params.min_cluster_cov}",
        f"threads: {params.threads}",
        "---",
        "",
        "# gapit cluster screening report",
        "",
    ]
    return "\n".join(lines) + "\n"


def cluster_md_file_chunk(report: ClusterReport, *, typed: bool) -> str:
    """One file's section, self-contained so it streams the moment the file
    completes: heading, the file's summary row (one-row table with the same
    columns the former all-files summary carried; a ``Phenotype`` column on
    typed databases), phenotype detail line on typed databases, and the best
    locus's gene table."""
    lines: list[str] = [f"## `{report.file}`", ""]
    lines.extend(
        _md_table(
            _SUMMARY_COLUMNS_TYPED if typed else _SUMMARY_COLUMNS, [_summary_row(report, typed)]
        )
    )
    lines.append("")
    top = next(
        (locus for locus in report.loci if report.best and locus.locus == report.best.locus),
        None,
    )
    if top is None:
        lines.append("_No locus detected._")
        lines.append("")
        return "\n".join(lines) + "\n"
    if typed and report.best is not None and report.best.phenotype_detail is not None:
        detail = report.best.phenotype_detail
        lines.append(
            f"Phenotype `{_phenotype_cell(report)}` (score {detail.score:.4f},"
            f" {detail.confidence} confidence)"
        )
        lines.append("")
    lines.append(f"Best locus `{top.locus}` ({top.label}, type {top.type}):")
    lines.append("")
    lines.extend(_md_table(_GENE_COLUMNS, _gene_rows(top)))
    lines.append("")
    return "\n".join(lines) + "\n"


def render_cluster_md(
    reports: Iterable[ClusterReport], params: ClusterParams, *, now: datetime, typed: bool = False
) -> str:
    """Render cluster screening as deterministic Markdown: static YAML
    frontmatter, then one section per file (its summary row with the
    phenotype column on typed databases, then the best locus's gene table).
    The head + file chunks concatenate to exactly this document."""
    files = list(reports)
    return cluster_md_head(params, now=now) + "".join(
        cluster_md_file_chunk(report, typed=typed) for report in files
    )
