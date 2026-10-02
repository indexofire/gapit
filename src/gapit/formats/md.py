"""Human- and agent-readable Markdown reports (gapit.report/1, gapit.reads/1).

The report is pure gene detection (typed and untyped gene databases render
byte-identically; designation renders from formats/typing_result.py).
"""

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


def md_cell(value: str) -> str:
    """Escape Markdown table metacharacters so cells cannot break the table
    (victors gene ids like ``gi|115534241:2616-3152`` contain pipes)."""
    return value.replace("\\", "\\\\").replace("|", "\\|")


def md_report_preamble(params: ScreeningParams, *, now: datetime) -> str:
    """The frontmatter + title chunk: STATIC metadata only (schema, tool,
    timestamp, thresholds — everything knowable before file 1), so a
    streaming caller emits it first and then one section per file as it
    completes. Run totals (files/hits) live in the JSON document, not here."""
    lines: list[str] = [
        "---",
        "schema: gapit.report/1",
        f"tool: gapit {__version__}",
        f"created_at: {now.astimezone(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"db: {params.db}",
        f"minid: {params.minid}",
        f"mincov: {params.mincov}",
        f"threads: {params.threads}",
        "---",
        "",
        "# gapit screening report",
        "",
    ]
    return "\n".join(lines) + "\n"


def md_report_chunk(report: Report) -> str:
    """One file's section: hits table and merged-fragment lines — the unit a
    streaming caller computes as the file completes."""
    lines: list[str] = [f"## `{report.file}`", ""]
    if not report.hits:
        lines.append("_No hits._")
        lines.append("")
        return "\n".join(lines) + "\n"
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
        lines.append("| " + " | ".join(md_cell(cell) for cell in cells) + " |")
    for hit in report.hits:
        if not hit.merged:
            continue
        spans = ", ".join(
            f"`{fragment.contig}:{fragment.start}-{fragment.end}({fragment.strand})`"
            for fragment in hit.fragments
        )
        lines.append(f"- `{hit.gene}` merged from {len(hit.fragments)} fragments: {spans}")
    lines.append("")
    return "\n".join(lines) + "\n"


def render_markdown(reports: Iterable[Report], params: ScreeningParams, *, now: datetime) -> str:
    """Render screening results as deterministic Markdown: static YAML
    frontmatter, then one section per file, stable table columns. The
    preamble + section chunks concatenate to exactly this document."""
    reports = list(reports)
    return md_report_preamble(params, now=now) + "".join(
        md_report_chunk(report) for report in reports
    )


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

_READS2_COLUMNS = (*_READS_COLUMNS, "Identity%")


def _render_reads_markdown(
    reports: Iterable[ReadsReport],
    params: ReadsParams,
    *,
    now: datetime,
    schema: str,
    columns: tuple[str, ...],
    with_identity: bool,
) -> str:
    """Shared reads-report Markdown core; gapit.reads/2 adds the two filter
    lines to the frontmatter and an Identity% cell per gene."""
    files = list(reports)
    genes_found = sum(1 for report in files for gene in report.genes if gene.present)
    lines: list[str] = [
        "---",
        f"schema: {schema}",
        f"tool: gapit {__version__}",
        f"created_at: {now.astimezone(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"db: {params.db}",
        f"read_type: {params.read_type}",
        f"min_breadth: {params.min_breadth}",
        f"threads: {params.threads}",
    ]
    if with_identity:
        lines.append(f"min_identity: {params.min_identity}")
        lines.append(f"min_mapq: {params.min_mapq}")
    lines += [
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
        lines.append("| " + " | ".join(columns) + " |")
        lines.append("|" + "---|" * len(columns))
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
            if with_identity:
                cells += (f"{gene.mean_identity_pct:.2f}",)
            lines.append("| " + " | ".join(md_cell(cell) for cell in cells) + " |")
        lines.append("")
    return "\n".join(lines) + "\n"


def render_reads_markdown(
    reports: Iterable[ReadsReport], params: ReadsParams, *, now: datetime
) -> str:
    """Render read-screening results as deterministic Markdown."""
    return _render_reads_markdown(
        reports,
        params,
        now=now,
        schema="gapit.reads/1",
        columns=_READS_COLUMNS,
        with_identity=False,
    )


def render_reads2_markdown(
    reports: Iterable[ReadsReport], params: ReadsParams, *, now: datetime
) -> str:
    """Render filtered read-screening results as deterministic Markdown
    (gapit.reads/2): the reads/1 shape plus the filter thresholds in
    frontmatter and an Identity% column per gene."""
    return _render_reads_markdown(
        reports,
        params,
        now=now,
        schema="gapit.reads/2",
        columns=_READS2_COLUMNS,
        with_identity=True,
    )
