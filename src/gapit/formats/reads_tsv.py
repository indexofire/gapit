"""Streaming reads TSV/CSV (the positional wildcard table format).

The positional all-FASTQ wildcard defaults to tsv so a terminal sees each
sample's rows the moment its screening completes; json remains the opt-in
single document. Chunk concatenation is byte-identical to the buffered
``format_reads_tsv`` render (the reads_md_preamble/chunk precedent).
"""

from collections.abc import Iterable

from gapit.reads import GeneCoverage, ReadsReport

_READS_HEADER = "#SAMPLE\tGENE\tBREADTH%\tDEPTH\tREADS\tPRESENT\tDATABASE\tACCESSION\tPRODUCT"
_READS2_HEADER = f"{_READS_HEADER}\tIDENTITY%"


def _gene_row(gene: GeneCoverage, sample: str, *, reads2: bool) -> list[str]:
    """One table row's cells: the sample key leading, then the gene call
    (reads/2 appends the Identity% cell) — the md `_reads_md_section`
    columns in abricate-table shape; PRODUCT closes the row (the RESISTANCE
    cell rides the contig table only)."""
    cells = [
        sample,
        gene.gene,
        f"{gene.breadth_pct:.2f}",
        f"{gene.mean_depth:.2f}",
        str(gene.reads_mapped),
        "yes" if gene.present else "no",
        gene.database,
        gene.accession,
        gene.product,
    ]
    if reads2:
        cells.append(f"{gene.mean_identity_pct:.2f}")
    return cells


def reads_tsv_preamble(*, reads2: bool = False, csv: bool = False) -> str:
    """The header chunk (abricate's ``#``-prefixed first line): static,
    knowable before sample 1, so a streaming caller emits it first."""
    header = _READS2_HEADER if reads2 else _READS_HEADER
    return header.replace("\t", "," if csv else "\t") + "\n"


def reads_tsv_chunk(
    report: ReadsReport, *, reads2: bool = False, csv: bool = False, all_genes: bool = False
) -> str:
    """One sample's chunk: one line per PRESENT gene by default (absent
    calls join with ``all_genes``), the streaming unit emitted as that
    sample's screening completes (a gene-less sample contributes no
    rows, like the contig table's hit-less report)."""
    sep = "," if csv else "\t"
    sample = ", ".join(report.reads)
    genes = report.genes if all_genes else [g for g in report.genes if g.present]
    return "".join(sep.join(_gene_row(gene, sample, reads2=reads2)) + "\n" for gene in genes)


def format_reads_tsv(
    reports: Iterable[ReadsReport],
    *,
    reads2: bool = False,
    csv: bool = False,
    all_genes: bool = False,
) -> str:
    """Buffered render: the header, then one line per (present) gene per
    report in report order; preamble + chunks concatenate to exactly this."""
    return reads_tsv_preamble(reads2=reads2, csv=csv) + "".join(
        reads_tsv_chunk(report, reads2=reads2, csv=csv, all_genes=all_genes) for report in reports
    )
