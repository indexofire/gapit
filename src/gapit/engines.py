"""Engine and format enums shared by the screening use-cases.

A leaf module: the gene path (screening.py) and the reads path
(screening_reads.py) both branch on these; they previously lived in
screening.py.
"""

import enum


class OutputFormat(enum.Enum):
    """Screen output formats."""

    tsv = "tsv"
    csv = "csv"
    json = "json"
    md = "md"


class AlignerEnum(enum.Enum):
    """Alignment engines for screen (SPEC.md §1/§10)."""

    blastn = "blastn"
    minimap2 = "minimap2"
