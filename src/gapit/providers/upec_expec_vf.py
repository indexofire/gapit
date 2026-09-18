"""upec_expec_vf provider — UPEC/ExPEC virulence genes (FordeGenomics).

Transform-only module (Wave B12): turns the downloaded
``UPEC_ExPEC_VF.tsv`` into typed ``Record``s exactly like upstream
``abricate-get_db get_upec_expec_vf`` — ``load_tabular`` keys every row by
its ``Gene name`` column with first occurrence winning (perl ``||=``), then
ID=Gene name, ACC=``Accession/Source:Begin-End``, DESC=Description,
SEQ=Sequence. Function is the locked ``virulence`` constant: the TSV Class
column is source metadata, NOT the function slot (Wave F2c). Sequences
stay RAW here: normalization, sequence dedupe, and sorting belong to the
generic :func:`gapit.providers.common.fetch_provider` pipeline.
"""

import csv
from collections.abc import Iterator
from pathlib import Path

from gapit.errors import DatabaseError
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "upec_expec_vf"

_FUNCTION = ("virulence",)  # locked gapit/v1 func vocabulary (Wave F2c)

_FILENAME = "UPEC_ExPEC_VF.tsv"
_REQUIRED_COLUMNS = (
    "Gene name",
    "Accession/Source",
    "Begin",
    "End",
    "Description",
    "Class",
    "Sequence",
)


def transform(workdir: Path) -> Iterator[Record]:
    """Yield one Record per well-formed row of ``workdir/UPEC_ExPEC_VF.tsv``.

    First row is the header; every column is looked up by NAME. Rows whose
    Gene name or Sequence is empty are skipped (upstream would splice undef
    pieces into the record — deviation noted in the Phase 7 notepad). A
    header missing required columns raises DatabaseError PROVIDER_MALFORMED.
    """
    path = workdir / _FILENAME
    with path.open(encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t", restval="")
        missing = [
            column for column in _REQUIRED_COLUMNS if column not in set(rows.fieldnames or ())
        ]
        if missing:
            raise DatabaseError(
                f"upec_expec_vf source is missing column(s): {', '.join(missing)}",
                code="PROVIDER_MALFORMED",
                context={"file": str(path), "missing_columns": ",".join(missing)},
            )
        seen: set[str] = set()
        for row in rows:
            gene = row["Gene name"]
            sequence = row["Sequence"]
            if not gene or not sequence or gene in seen:
                continue
            seen.add(gene)
            yield Record(
                db=NAME,
                gene=gene,
                accession=f"{row['Accession/Source']}:{row['Begin']}-{row['End']}",
                function=_FUNCTION,
                product=row["Description"],
                sequence=sequence,
                source_id=gene,
            )


PROVIDER = Provider(
    name=NAME,
    description="UPEC/ExPEC virulence genes (FordeGenomics)",
    source_urls=(
        "https://raw.githubusercontent.com/FordeGenomics/ST167_Code/refs/heads/main/UPEC-ExPEC_VF/UPEC_ExPEC_VF.tsv",
    ),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
