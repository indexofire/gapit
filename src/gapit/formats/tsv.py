"""abricate-byte-compatible TSV/CSV rendering of Reports (SPEC.md §5).

The 15-column abricate table is the frozen parity surface and the only
report table: typed gene databases screen byte-identically to untyped
ones (designation is the ``gapit typing`` command's job). The cluster
TSV's typed-only PHENOTYPE column (formats.cluster) is the cluster-side
exception — locus databases have no re-readable screen table.
"""

from collections.abc import Iterable
from pathlib import PurePath

from gapit.formats.csvutil import delimited_line
from gapit.report import Report

HEADER = (
    "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t"
    "%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE"
)
_HEADER_CELLS = tuple(HEADER.split("\t"))


def tsv_preamble(*, csv: bool, noheader: bool) -> str:
    """The header chunk: the ``#FILE`` line, or "" under noheader (sinks
    skip empty chunks, so a noheader stream starts straight at file rows)."""
    if noheader:
        return ""
    return delimited_line(_HEADER_CELLS, as_csv=csv)


def tsv_file_chunk(report: Report, *, csv: bool, nopath: bool) -> str:
    """One file's chunk: one line per hit in report order, the streaming
    unit the use-case emits as soon as the file's report exists."""
    file_column = PurePath(report.file).name if nopath else report.file
    return "".join(
        delimited_line(
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
                hit.function,
            ),
            as_csv=csv,
        )
        for hit in report.hits
    )


def format_tsv(reports: Iterable[Report], *, csv: bool, noheader: bool, nopath: bool) -> str:
    """Render reports as abricate-format lines: one header (unless noheader),
    then one line per hit per report in report order.

    CSV is RFC 4180-quoted by :mod:`gapit.formats.csvutil`; TSV keeps the
    frozen raw-tab shape, so control characters in DB fields are rejected at
    the database-build boundary (BUILD_INVALID) instead of being escaped.
    """
    return tsv_preamble(csv=csv, noheader=noheader) + "".join(
        tsv_file_chunk(report, csv=csv, nopath=nopath) for report in reports
    )
