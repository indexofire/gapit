"""abricate-byte-compatible TSV/CSV rendering of Reports (SPEC.md §5)."""

from collections.abc import Iterable
from pathlib import PurePath

from gapit.report import Report

HEADER = (
    "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t"
    "%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE"
)


def format_tsv(reports: Iterable[Report], *, csv: bool, noheader: bool, nopath: bool) -> str:
    """Render reports as abricate-format lines: one header (unless noheader),
    then one line per hit per report in report order. Line = sep.join(fields)
    + "\\n"; fields never contain the separator (product cleanup strips commas).
    """
    sep = "," if csv else "\t"
    lines: list[str] = [] if noheader else [HEADER.replace("\t", sep)]
    for report in reports:
        file_column = PurePath(report.file).name if nopath else report.file
        for hit in report.hits:
            lines.append(
                sep.join(
                    (
                        file_column,
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
                        hit.resistance,
                    )
                )
            )
    return "".join(line + "\n" for line in lines)
