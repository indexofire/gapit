"""ncbi database provider — NCBI AMRFinderPlus curated AMR (Wave B1).

Transform-only module mirroring abricate-get_db ``get_ncbi``: pair
``AMR_CDS.fa`` records with ``ReferenceGeneCatalog.txt`` rows keyed by
``refseq_nucleotide_accession`` (column 10, 0-based), keeping only plain
(non-fusion) genes whose catalog row is scope core / type AMR / subtype AMR.
``PROVIDER`` wires the pinned metadata into the frozen B0
:class:`gapit.providers.common.Provider` contract.
"""

import re
from collections.abc import Iterator
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "ncbi"
AMR_CDS_FILE = "AMR_CDS.fa"
CATALOG_FILE = "ReferenceGeneCatalog.txt"
_ACCESSION_COLUMN = 10
_VERSIONED_ACCESSION = re.compile(r"\.\d+$")
_LATEST = (
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
)


def _load_catalog(path: Path) -> dict[str, dict[str, str]]:
    """ReferenceGeneCatalog rows as ``{accession: {header name: value}}``.

    The first line is the header; every later row is keyed by column 10
    (``refseq_nucleotide_accession``). Duplicate accessions keep the FIRST
    row (upstream's ``||=``), and rows are padded to the header width so
    lookups of unused trailing columns cannot fail on short rows.
    """
    rows: dict[str, dict[str, str]] = {}
    header: list[str] | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            columns = line.rstrip("\r\n").split("\t")
            if header is None:
                header = columns
                continue
            width = max(len(header), _ACCESSION_COLUMN + 1)
            padded = columns + [""] * (width - len(columns))
            if (acc := padded[_ACCESSION_COLUMN]) and acc not in rows:
                rows[acc] = dict(zip(header, padded, strict=False))  # perl zip: shorter wins
    return rows


def transform(workdir: Path) -> Iterator[Record]:
    """Parse ``workdir/AMR_CDS.fa`` + ``ReferenceGeneCatalog.txt`` (db = ncbi).

    AMRFinderPlus fasta ids are ``pi|acc|fp|fn|gene|fam|prod`` — perl's
    7-variable ``split`` assignment imposes LIMIT 7, so everything past the
    sixth pipe belongs to ``prod``. Records are skipped — passively, like
    upstream — unless the id yields 7 non-empty fields, fp/fn are both "1"
    (fusion filter), and the accession (``.1`` appended unless already
    versioned, e.g. ``NG_050200`` -> ``NG_050200.1``) matches a catalog row
    with scope ``core``, type ``AMR``, and subtype ``AMR``. Sequences and
    function categories stay raw: fetch_provider owns normalization.
    """
    catalog = _load_catalog(workdir / CATALOG_FILE)
    for fasta in iter_fasta(workdir / AMR_CDS_FILE):
        fields = fasta.id.split("|")
        if len(fields) < 7:
            continue
        pi, acc, fp, fn, gene, fam = fields[:6]
        prod = "|".join(fields[6:])
        if not all((pi, acc, fp, fn, gene, fam, prod)):
            continue
        if fp != "1" or fn != "1":
            continue
        if not _VERSIONED_ACCESSION.search(acc):
            acc += ".1"
        row = catalog.get(acc)
        if row is None:
            continue
        if row["scope"] != "core" or row["type"] != "AMR" or row["subtype"] != "AMR":
            continue
        yield Record(
            db=NAME,
            gene=gene,
            accession=row["refseq_nucleotide_accession"],
            function=tuple(part for part in row["subclass"].split("/") if part),
            product=prod.replace("_", " "),
            sequence=fasta.sequence,
            source_id=pi,
        )


PROVIDER = Provider(
    name=NAME,
    description="NCBI AMRFinderPlus (reference finder) curated AMR",
    source_urls=(
        f"{_LATEST}/AMR_CDS.fa",
        f"{_LATEST}/ReferenceGeneCatalog.txt",
    ),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
