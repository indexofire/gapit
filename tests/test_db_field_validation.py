"""Database-build boundary: output fields must not break the table dialect."""

from pathlib import Path

import pytest

from gapit.db_build_ops import perform_build
from gapit.dbbuild import validate_records
from gapit.errors import DatabaseError
from gapit.records import Record

SEQ = "ACGT" * 20


def _record(**overrides: object) -> Record:
    values: dict[str, object] = {
        "db": "db",
        "gene": "g1",
        "sequence": SEQ,
        "accession": "A1",
        "function": ("class",),
        "product": "product",
    }
    values.update(overrides)
    return Record.model_validate(values)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("gene", "g\t1"),
        ("gene", "g\n1"),
        ("accession", "A\t1"),
        ("db", "d\tb"),
        ("function", ("class\tone",)),
    ],
)
def test_validate_records_rejects_control_characters(field: str, value: object) -> None:
    """Given a gene/accession/db/function value with a table delimiter, When
    validated, Then BUILD_INVALID names the gene and field before any write."""
    with pytest.raises(DatabaseError) as excinfo:
        validate_records([_record(**{field: value})])

    assert excinfo.value.code == "BUILD_INVALID"
    assert excinfo.value.exit_code == 4
    assert excinfo.value.context["gene"] is not None
    assert excinfo.value.context["field"].startswith(field)


def test_build_rejects_control_char_before_creating_the_target(tmp_path: Path) -> None:
    """Given a custom FASTA whose gapit| gene decodes to a tab, When built,
    Then the typed error fires before the target directory/manifest exists."""
    fasta = tmp_path / "db.fa"
    fasta.write_text(
        f">gapit|db=x|gene=bad%09gene|acc=A|func=f\n{SEQ}\n",
        encoding="utf-8",
    )
    datadir = tmp_path / "dd"

    with pytest.raises(DatabaseError) as excinfo:
        perform_build(
            "bad",
            fasta,
            None,
            None,
            "",
            datadir,
            False,
            warn=lambda message: None,
            quiet=True,
        )

    assert excinfo.value.code == "BUILD_INVALID"
    assert not (datadir / "bad" / "gapit-manifest.json").exists()
