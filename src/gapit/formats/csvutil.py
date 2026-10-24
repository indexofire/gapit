"""One CSV/TSV line at a time, with RFC 4180 quoting for CSV.

The table renderers stream row-by-row, so they cannot hand a whole document
to ``csv.writer``. This leaf helper writes one row to a small buffer with the
standard library writer (``QUOTE_MINIMAL``, LF terminator) and returns the
exact chunk. TSV keeps the abricate-compatible raw tab join: TSV has no
quoting convention, so TSV field validation happens at the database-build
boundary instead.
"""

import csv
import io
from collections.abc import Sequence


def csv_line(cells: Sequence[str]) -> str:
    """One RFC 4180 row with LF endings (no CRLF)."""
    buffer = io.StringIO()
    csv.writer(buffer, lineterminator="\n").writerow(cells)
    return buffer.getvalue()


def delimited_line(cells: Sequence[str], *, as_csv: bool) -> str:
    """One row in the selected table dialect: quoted CSV or raw TSV."""
    if as_csv:
        return csv_line(cells)
    return "\t".join(cells) + "\n"
