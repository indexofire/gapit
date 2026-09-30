"""Engine and format enums shared by the screening use-cases.

A leaf module (stage 3): the gene path (screening.py), the reads path
(screening_reads.py), and the cluster path (screening_cluster.py) all branch
on these; they previously lived in screening.py, which the cluster split
would have made a cycle.
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
