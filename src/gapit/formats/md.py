"""Human- and agent-readable Markdown reports (gapit.report/1, gapit.reads/1)."""

from collections.abc import Iterable
from datetime import UTC, datetime

from gapit import __version__
from gapit.reads import ReadsParams, ReadsReport
from gapit.report import Report, ScreeningParams

_COLUMNS = (
    "Sequence",
    "Start",
    "End",
    "Strand",
    "Gene",
    "Coverage",
    "Map",
    "Gaps",
    "%Coverage",
    "%Identity",
    "Database",
    "Accession",
    "Product",
    "Resistance",
)


def render_markdown(reports: Iterable[Report], params: ScreeningParams, *, now: datetime) -> str:
    """Render screening results as deterministic Markdown: YAML frontmatter,
    one section per file, stable table columns."""
    files = list(reports)
    total_hits = sum(len(report.hits) for report in files)
    lines: list[str] = [
        "---",
        "schema: gapit.report/1",
        f"tool: gapit {__version__}",
        f"created_at: {now.astimezone(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"db: {params.db}",
        f"minid: {params.minid}",
        f"mincov: {params.mincov}",
        f"threads: {params.threads}",
        f"files: {len(files)}",
        f"hits: {total_hits}",
        "---",
        "",
        "# gapit screening report",
        "",
    ]
    for report in files:
        lines.append(f"## `{report.file}`")
        lines.append("")
        if not report.hits:
            lines.append("_No hits._")
            lines.append("")
            continue
        lines.append("| " + " | ".join(_COLUMNS) + " |")
        lines.append("|" + "---|" * len(_COLUMNS))
        for hit in report.hits:
            cells = (
                hit.sequence,
                str(hit.start),
                str(hit.end),
                hit.strand,
                hit.gene,
                f"{hit.s_start}-{hit.s_end}/{hit.s_len}",
                hit.coverage_map,
                f"{hit.gap_openings}/{hit.gaps}",
                f"{hit.coverage_pct:.2f}",
                f"{hit.identity_pct:.2f}",
                hit.database,
                hit.accession,
                hit.product,
                hit.function,
            )
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"


_READS_COLUMNS = (
    "Gene",
    "Breadth%",
    "Depth",
    "Reads",
    "Present",
    "Database",
    "Accession",
    "Product",
    "Resistance",
)


def render_reads_markdown(
    reports: Iterable[ReadsReport], params: ReadsParams, *, now: datetime
) -> str:
    """Render read-screening results as deterministic Markdown."""
    files = list(reports)
    genes_found = sum(1 for report in files for gene in report.genes if gene.present)
    lines: list[str] = [
        "---",
        "schema: gapit.reads/1",
        f"tool: gapit {__version__}",
        f"created_at: {now.astimezone(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"db: {params.db}",
        f"read_type: {params.read_type}",
        f"min_breadth: {params.min_breadth}",
        f"threads: {params.threads}",
        f"files: {len(files)}",
        f"genes_found: {genes_found}",
        "---",
        "",
        "# gapit read screening report",
        "",
    ]
    for report in files:
        lines.append(f"## `{', '.join(report.reads)}`")
        lines.append("")
        if not report.genes:
            lines.append("_No genes detected._")
            lines.append("")
            continue
        lines.append("| " + " | ".join(_READS_COLUMNS) + " |")
        lines.append("|" + "---|" * len(_READS_COLUMNS))
        for gene in report.genes:
            cells = (
                gene.gene,
                f"{gene.breadth_pct:.2f}",
                f"{gene.mean_depth:.2f}",
                str(gene.reads_mapped),
                "yes" if gene.present else "no",
                gene.database,
                gene.accession,
                gene.product,
                gene.function,
            )
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"
