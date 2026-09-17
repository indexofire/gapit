"""Tests for the upec_expec_vf provider transform (Wave B12).

Offline and transform-only: the committed synthetic TSV fixture exercises
upstream ``abricate-get_db get_upec_expec_vf`` semantics — rows keyed by
'Gene name' with first occurrence winning (perl ``||=``), accession built as
``Accession/Source:Begin-End``, and function the locked ``virulence``
constant (the TSV Class column is source metadata, NOT the function slot —
Wave F2c). Sequences stay RAW here: normalization, sequence dedupe, and
sorting belong to the generic :func:`gapit.providers.common.fetch_provider`
pipeline.
"""

import shutil
from pathlib import Path

import pytest

from gapit.errors import DatabaseError
from gapit.records import Record

FIXTURE = Path(__file__).parent / "data" / "providers" / "upec_expec_vf"

# Sequences as they appear in the fixture (the transform must NOT normalize).
CHUA_SEQUENCE = "ACGT" * 15
VAT_SEQUENCE = "AGCT" * 15


def workdir_with_fixture(tmp_path: Path) -> Path:
    """A download workdir holding a copy of the synthetic source TSV."""
    workdir = tmp_path / "work"
    workdir.mkdir()
    shutil.copy2(FIXTURE / "UPEC_ExPEC_VF.tsv", workdir / "UPEC_ExPEC_VF.tsv")
    return workdir


def test_transform_yields_exact_records_when_fixture_is_well_formed(
    tmp_path: Path,
) -> None:
    """Given the synthetic fixture, When transformed, Then exactly the
    well-formed rows become Records in file order — db/gene/source_id set,
    accession joined as Accession/Source:Begin-End, product from
    Description, function the locked virulence constant."""
    from gapit.providers.upec_expec_vf import transform

    records = tuple(transform(workdir_with_fixture(tmp_path)))

    assert records == (
        Record(
            db="upec_expec_vf",
            gene="chuA",
            sequence=CHUA_SEQUENCE,
            accession="NC_000913:4342631-4343534",
            function=("virulence",),
            product="outer membrane heme receptor",
            source_id="chuA",
        ),
        Record(
            db="upec_expec_vf",
            gene="vat",
            sequence=VAT_SEQUENCE,
            accession="NC_000913:1165553-1166518",
            function=("virulence",),
            product="vacuolating autotransporter toxin",
            source_id="vat",
        ),
    )


def test_transform_keeps_first_occurrence_when_gene_name_duplicates(
    tmp_path: Path,
) -> None:
    """Given two rows sharing Gene name 'chuA' with entirely different data,
    When transformed, Then only the FIRST row survives (upstream ||=) — the
    later duplicate's accession, sequence, and product are nowhere."""
    from gapit.providers.upec_expec_vf import transform

    records = tuple(transform(workdir_with_fixture(tmp_path)))

    assert [record.gene for record in records] == ["chuA", "vat"]
    chu_a = records[0]
    assert chu_a.accession == "NC_000913:4342631-4343534"
    assert chu_a.sequence == CHUA_SEQUENCE
    assert chu_a.product == "outer membrane heme receptor"
    assert all(record.sequence != "T" * 60 for record in records)


def test_transform_skips_row_when_sequence_is_empty(tmp_path: Path) -> None:
    """Given a row whose Sequence column is empty, When transformed, Then no
    Record is emitted for that gene (upstream would splice undef pieces; the
    deviation is noted in the Phase 7 notepad)."""
    from gapit.providers.upec_expec_vf import transform

    records = tuple(transform(workdir_with_fixture(tmp_path)))

    assert "iroN" not in [record.gene for record in records]


def test_transform_function_is_locked_constant_when_class_is_multi_word(
    tmp_path: Path,
) -> None:
    """Given a row whose Class is multi-word ('Virulence assoc'), When
    transformed, Then function is exactly the locked ('virulence',)
    constant — the TSV Class column is source metadata and must NOT leak
    into the function slot (Wave F2c)."""
    from gapit.providers.upec_expec_vf import transform

    records = tuple(transform(workdir_with_fixture(tmp_path)))

    assert records[1].function == ("virulence",)
    assert all("Virulence assoc" not in value for record in records for value in record.function)


def test_transform_raises_when_required_column_missing_from_header(
    tmp_path: Path,
) -> None:
    """Given a source TSV whose header lacks required columns, When
    transformed, Then DatabaseError PROVIDER_MALFORMED (exit 4) names the
    file and the missing columns in required-order (parse, don't validate:
    upstream would silently build all-undef records instead)."""
    from gapit.providers.upec_expec_vf import transform

    workdir = tmp_path / "work"
    workdir.mkdir()
    (workdir / "UPEC_ExPEC_VF.tsv").write_text(
        "Gene name\tDescription\tSequence\nchuA\tsome product\tACGT\n",
        encoding="utf-8",
    )

    with pytest.raises(DatabaseError) as excinfo:
        tuple(transform(workdir))
    error = excinfo.value
    assert error.code == "PROVIDER_MALFORMED"
    assert error.exit_code == 4
    assert error.context["missing_columns"] == "Accession/Source,Begin,End,Class"
    assert error.context["file"].endswith("UPEC_ExPEC_VF.tsv")


def test_provider_metadata_when_inspected() -> None:
    """Given the module's PROVIDER value, When inspected, Then name,
    description, source URLs, dbtype, and the transform wiring match the
    frozen Wave B contract exactly."""
    from gapit.providers.upec_expec_vf import PROVIDER, transform

    assert PROVIDER.name == "upec_expec_vf"
    assert PROVIDER.description == "UPEC/ExPEC virulence genes (FordeGenomics)"
    assert PROVIDER.source_urls == (
        "https://raw.githubusercontent.com/FordeGenomics/ST167_Code/refs/heads/main/UPEC-ExPEC_VF/UPEC_ExPEC_VF.tsv",
    )
    assert PROVIDER.dbtype == "nucl"
    assert PROVIDER.transform is transform
