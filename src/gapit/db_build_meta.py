"""The `gapit db build --tsv` metadata sidecar: parse + merge.

``read_metadata`` parses the optional TSV (header row with a mandatory
``gene`` column; ``accession``/``function`` optional, extra columns
ignored) and ``merge_metadata`` overlays the carried columns onto the
FASTA-parsed records. :mod:`gapit.db_build_ops` orchestrates; this module
owns the sidecar's rules only.
"""

import csv
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from gapit.errors import InputError
from gapit.records import Record


@dataclass(frozen=True, slots=True)
class Metadata:
    """Parsed --tsv: which merge columns the header carries, and the row
    values per gene (accession, ';'-joined function; '' for absent columns)."""

    columns: frozenset[str]
    rows: dict[str, tuple[str, str]]


def split_function(joined: str) -> tuple[str, ...]:
    """';'-joined classes -> tuple, empty pieces dropped."""
    return tuple(part for part in joined.split(";") if part)


def read_metadata(path: Path, warn: Callable[[str], None]) -> Metadata:
    """Parse the metadata TSV: header row mandatory and must carry ``gene``
    (InputError METADATA_MALFORMED otherwise); ``accession``/``function``
    are optional per-file, extra columns are ignored. Duplicate genes keep
    the FIRST row (warn); short rows read their missing cells as empty.
    """
    if not path.is_file():
        raise InputError(
            f"metadata file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next((row for row in rows if row), None)
        if header is None or "gene" not in header:
            raise InputError(
                f"metadata TSV needs a header row with a 'gene' column: {path}",
                code="METADATA_MALFORMED",
                context={"file": str(path), "needed": "gene"},
            )
        columns = {column: index for index, column in enumerate(header)}

        def cell(padded: list[str], column: str) -> str:
            return padded[columns[column]]

        metadata: dict[str, tuple[str, str]] = {}
        for row in rows:
            if not any(value.strip() for value in row):
                continue
            padded = row + [""] * (len(header) - len(row))
            gene = cell(padded, "gene")
            if not gene:
                continue
            if gene in metadata:
                warn(f"duplicate gene {gene!r} in metadata TSV: keeping the first row")
                continue
            metadata[gene] = (
                cell(padded, "accession") if "accession" in columns else "",
                cell(padded, "function") if "function" in columns else "",
            )
    return Metadata(columns=frozenset(columns), rows=metadata)


def merge_metadata(
    records: list[Record], metadata: Metadata, warn: Callable[[str], None]
) -> list[Record]:
    """Overwrite accession/function on records whose gene has a metadata row
    (only the columns the TSV actually carries); genes only present in the
    TSV are warned about and skipped."""
    merged: list[Record] = []
    for record in records:
        row = metadata.rows.get(record.gene)
        if row is None:
            merged.append(record)
            continue
        update: dict[str, object] = {}
        if "accession" in metadata.columns:
            update["accession"] = row[0]
        if "function" in metadata.columns:
            update["function"] = split_function(row[1])
        merged.append(record.model_copy(update=update))
    fasta_genes = {record.gene for record in records}
    for gene in sorted(set(metadata.rows) - fasta_genes):
        warn(f"gene {gene!r} in metadata TSV not found in FASTA: skipped")
    return merged
