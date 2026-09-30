"""gapit.cluster/1 Markdown renderer.

Split from formats/cluster.py at the 250-LOC ceiling (stage 3): the summary
table (one row per file, a ``Phenotype`` column on typed databases), the
best locus's gene table, and YAML frontmatter mirroring the engine params.
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


def _summary_rows(files: list[ClusterReport], typed: bool) -> list[tuple[str, ...]]:
    """One summary row per file; a null best renders the untyped dashes."""
    rows: list[tuple[str, ...]] = []
    for report in files:
        if report.best is None:
            rows.append((report.file, *("-",) * (8 if typed else 7)))
            continue
        best = report.best
        phenotype = (_phenotype_cell(report),) if typed else ()
        rows.append(
            (
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
        )
    return rows


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


def render_cluster_md(
    reports: Iterable[ClusterReport], params: ClusterParams, *, now: datetime, typed: bool = False
) -> str:
    """Render cluster screening as deterministic Markdown: YAML frontmatter,
    a per-file summary table (with the phenotype column on typed databases),
    then the best locus's gene table per file."""
    files = list(reports)
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
        f"files: {len(files)}",
        "---",
        "",
        "# gapit cluster screening report",
        "",
        *_md_table(
            _SUMMARY_COLUMNS_TYPED if typed else _SUMMARY_COLUMNS, _summary_rows(files, typed)
        ),
        "",
    ]
    for report in files:
        lines.append(f"## `{report.file}`")
        lines.append("")
        top = next(
            (locus for locus in report.loci if report.best and locus.locus == report.best.locus),
            None,
        )
        if top is None:
            lines.append("_No locus detected._")
            lines.append("")
            continue
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
