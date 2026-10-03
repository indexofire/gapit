"""CLI integration tests: `gapit db fetch NAME` (provider path) + `db list`.

Offline by construction: ``gapit.db_ops.REGISTRY`` is monkeypatched down to
two synthetic Providers whose source_urls point at a ``file://`` fasta fixture
in tmp_path, so the full download -> transform -> build pipeline runs against
the pixi env's real makeblastdb/minimap2 without touching the network (the
test_providers_common.py pattern, driven through the CLI surface instead).
"""

import io
import json
import sys
from collections.abc import Iterable
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from rich.console import Console
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.cmd_db import db_list_command
from gapit.db_ops import DbListEntry, db_list_table
from gapit.errors import ErrorEnvelope
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

SYN = "synamr"
OTHER = "synavail"
SYN_DESCRIPTION = "synthetic provider for the db fetch CLI tests"

# 60 bp sequences (Wave A3 pinned minimap2 indexing 60-70 bp synthetic seqs).
SEQ_A = "ACGTAG" * 10
SEQ_B = "ACGGCTAGAT" * 6
UPSTREAM_FASTA = f">syn_a first synthetic gene\n{SEQ_A}\n>syn_b second synthetic gene\n{SEQ_B}\n"

runner = CliRunner()
envelope_adapter = TypeAdapter(ErrorEnvelope)


def last_envelope(stderr: str) -> ErrorEnvelope:
    """Parse the last non-empty stderr line as the error envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return envelope_adapter.validate_json(lines[-1])


def syn_transform(workdir: Path) -> Iterable[Record]:
    """Load the downloaded upstream fasta and yield Records (db = SYN)."""
    for fasta in iter_fasta(workdir / "upstream.fa"):
        yield Record(db=SYN, gene=fasta.id, sequence=fasta.sequence)


def patch_registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point ``gapit.db_ops.REGISTRY`` at two synthetic file:// providers
    (SYN fetchable, OTHER never fetched) and create an empty datadir under
    tmp_path; returns that datadir."""
    source = tmp_path / "upstream.fa"
    source.write_text(UPSTREAM_FASTA, encoding="utf-8")
    url = source.as_uri()
    monkeypatch.setattr(
        "gapit.db_ops.REGISTRY",
        {
            SYN: Provider(
                name=SYN,
                description=SYN_DESCRIPTION,
                vendor="Synthetica",
                source_urls=(url,),
                dbtype="nucl",
                transform=syn_transform,
            ),
            OTHER: Provider(
                name=OTHER,
                description="available but never fetched in these tests",
                vendor="Example Org",
                source_urls=(url,),
                dbtype="nucl",
                transform=syn_transform,
            ),
        },
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    return datadir


def fetch(*extra: str) -> Result:
    return runner.invoke(app, ["db", "fetch", *extra])


def test_db_fetch_builds_provider_database(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a registry provider backed by a file:// fasta, When
    `db fetch NAME --datadir D`, Then exit 0, stdout is a one-line JSON
    receipt (db/records/dbtype/destination) naming D, and the database
    artifacts exist under D (the --datadir flag is respected)."""
    datadir = patch_registry(tmp_path, monkeypatch)
    result = fetch(SYN, "--datadir", str(datadir))
    assert result.exit_code == 0
    db_dir = datadir / SYN
    assert json.loads(result.stdout) == {
        "db": SYN,
        "records": 2,
        "dbtype": "nucl",
        "destination": str(db_dir),
    }
    for name in ("sequences", "sequences.nin", "gapit-manifest.json"):
        assert (db_dir / name).is_file(), name


def test_db_fetch_unknown_provider_is_usage_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a patched registry, When `db fetch` names a provider that is not
    in it, Then exit 2 with a USAGE_ERROR envelope whose message lists the
    available provider names."""
    patch_registry(tmp_path, monkeypatch)
    result = fetch("ncbi")
    assert result.exit_code == 2
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "USAGE_ERROR"
    assert "ncbi" in envelope.message
    assert SYN in envelope.message
    assert OTHER in envelope.message


def test_db_fetch_refuses_overwrite_without_force(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an already-fetched database, When `db fetch NAME` runs again
    without --force, Then exit 4 with a DB_ALREADY_EXISTS envelope and nothing
    on stdout."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = fetch(SYN, "--datadir", str(datadir))
    assert result.exit_code == 4
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope.code == "DB_ALREADY_EXISTS"
    assert envelope.context["db"] == SYN


def test_db_fetch_force_rebuilds_existing_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an already-fetched database, When `db fetch NAME --force` runs,
    Then exit 0 with a fresh receipt for the same record count."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = fetch(SYN, "--datadir", str(datadir), "--force")
    assert result.exit_code == 0
    assert json.loads(result.stdout)["records"] == 2


def test_db_fetch_creates_missing_nested_datadir(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a --datadir whose nested path does not exist yet (fresh
    machine), When `db fetch NAME`, Then exit 0 with the receipt naming the
    created datadir and the database artifacts land under it — not a
    DATADIR_NOT_FOUND refusal (Wave E bootstrap gap)."""
    datadir = patch_registry(tmp_path, monkeypatch)
    nested = datadir / "fresh" / "machine" / "db"
    result = fetch(SYN, "--datadir", str(nested))
    assert result.exit_code == 0
    assert json.loads(result.stdout)["destination"] == str(nested / SYN)
    assert (nested / SYN / "gapit-manifest.json").is_file()


def test_db_fetch_debug_echoes_index_argv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a registry provider backed by a file:// fasta, When
    `db fetch NAME --debug`, Then exit 0 with the SAME one-line stdout
    receipt as a plain fetch, the makeblastdb argv line is echoed to stderr
    (`gapit: run:`), and no `minimap2 -d` line appears (no .mmi is built)."""
    datadir = patch_registry(tmp_path, monkeypatch)
    plain = fetch(SYN, "--datadir", str(datadir), "--force")
    assert plain.exit_code == 0
    assert "gapit: run:" not in plain.stderr
    result = fetch(SYN, "--datadir", str(datadir), "--force", "--debug")
    assert result.exit_code == 0
    assert result.stdout == plain.stdout
    run_lines = [line for line in result.stderr.splitlines() if line.startswith("gapit: run:")]
    assert any(line.startswith("gapit: run: makeblastdb -in ") for line in run_lines)
    assert not any("minimap2 -d " in line for line in run_lines)


def test_db_list_text_shows_installed_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir holding one fetched provider, When `db list`, Then the
    text table shows that provider as installed with its record count and the
    unfetched one as available."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir)])
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0] == "NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION"
    assert f"{SYN}\tSynthetica\tinstalled (2)\tnucl\t{SYN_DESCRIPTION}" in lines
    other_row = f"{OTHER}\tExample Org\tavailable\tnucl\tavailable but never fetched"
    assert f"{other_row} in these tests" in lines


def test_db_list_json_emits_gapit_dblist_document(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a datadir holding one fetched provider, When `db list --json`,
    Then stdout parses as a gapit.dblist/1 document whose installed entry
    carries records and whose available entry omits the field."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])
    assert result.exit_code == 0
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.dblist/1"
    by_name = {entry["name"]: entry for entry in document["providers"]}
    assert by_name[SYN]["installed"] is True
    assert by_name[SYN]["records"] == 2
    assert by_name[SYN]["dbtype"] == "nucl"
    assert by_name[SYN]["vendor"] == "Synthetica"
    assert by_name[OTHER]["installed"] is False
    assert "records" not in by_name[OTHER]


class FakeStdout(io.StringIO):
    """A stdout stand-in whose TTY-ness is pinned — the isatty seam that
    decides rich table vs TSV in `db list`."""

    def __init__(self, tty: bool) -> None:
        super().__init__()
        self._tty = tty

    def isatty(self) -> bool:
        return self._tty


def test_db_list_table_builds_rich_table() -> None:
    """Given a mixed installed/available entry list, When db_list_table
    renders at width 100, Then the title, all five column headers, both
    STATUS texts, and the vendor strings appear — and the captured text
    carries no ANSI escapes."""
    console = Console(record=True, width=100)
    console.print(
        db_list_table(
            (
                DbListEntry(
                    name="card",
                    vendor="McMaster University",
                    description="CARD protein homolog resistance models",
                    dbtype="nucl",
                    installed=True,
                    records=203,
                ),
                DbListEntry(
                    name="vfdb",
                    vendor="USTC (VFDB)",
                    description="VFDB virulence factors (set A, nucleotide)",
                    dbtype="nucl",
                    installed=False,
                ),
            )
        )
    )
    text = console.export_text()
    assert "Databases" in text
    for header in ("Name", "Provider", "Status", "DBTYPE", "Description"):
        assert header in text
    assert "card" in text
    assert "McMaster University" in text
    assert "USTC (VFDB)" in text
    assert "installed (203)" in text
    assert "available" in text
    assert "\x1b" not in text


def test_db_list_non_tty_emits_exact_tsv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a datadir with one fetched provider and a NON-TTY stdout
    (isatty False), When `db list` runs, Then stdout is byte-identical TSV —
    header plus one row per provider — never a rich table."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    fake = FakeStdout(tty=False)
    monkeypatch.setattr(sys, "stdout", fake)
    db_list_command(datadir=datadir)
    assert fake.getvalue() == (
        "NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION\n"
        f"{SYN}\tSynthetica\tinstalled (2)\tnucl\t{SYN_DESCRIPTION}\n"
        f"{OTHER}\tExample Org\tavailable\tnucl\tavailable but never fetched in these tests\n"
        "ecoh\tHolt lab (srst2)\tbundled\tnucl\tE. coli O and H antigens (srst2 EcOH)\n"
        "ecoli_dec\tgapit-curated (public-domain sources)\tbundled\tnucl\t"
        "Diarrheagenic E. coli marker panel (GB 4789.6 + risk-monitoring designation)\n"
        "lm_doumith\tgapit-curated (public-domain INSDC sources)\tbundled\tnucl\t"
        "Listeria monocytogenes serogrouping (Doumith 2004)\n"
        "ncbi\tNCBI\tbundled\tnucl\tNCBI AMRFinderPlus (reference finder) curated AMR\n"
        "resfinder\tDTU CGE\tbundled\tnucl\tCGE ResFinder acquired resistance genes\n"
        "upec_expec_vf\tFordeGenomics\tbundled\tnucl\tUPEC/ExPEC virulence genes (FordeGenomics)\n"
    )


def test_db_list_tty_renders_rich_table_not_tsv(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the same datadir but a TTY stdout (isatty True), When `db list`
    runs, Then stdout carries the rich table (title + status texts) and
    never the TSV header line."""
    datadir = patch_registry(tmp_path, monkeypatch)
    assert fetch(SYN, "--datadir", str(datadir)).exit_code == 0
    fake = FakeStdout(tty=True)
    monkeypatch.setattr(sys, "stdout", fake)
    db_list_command(datadir=datadir)
    out = fake.getvalue()
    assert "Databases" in out
    assert "installed (2)" in out
    assert "available" in out
    assert "NAME\tPROVIDER\tSTATUS\tDBTYPE\tDESCRIPTION" not in out


def test_db_fetch_all_installs_default_dbs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a registry whose card+vfdb providers point at a file:// fasta,
    When `db fetch all`, Then both defaults install in DEFAULT_DBS order
    over the download path with one JSON receipt line per db on stdout."""
    source = tmp_path / "upstream.fa"
    source.write_text(UPSTREAM_FASTA, encoding="utf-8")
    url = source.as_uri()
    monkeypatch.setattr(
        "gapit.db_ops.REGISTRY",
        {
            name: Provider(
                name=name,
                description=f"synthetic {name} provider for the default-set test",
                vendor="Synthetica",
                source_urls=(url,),
                dbtype="nucl",
                transform=syn_transform,
            )
            for name in ("card", "vfdb")
        },
    )
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = fetch("all", "--datadir", str(datadir))

    assert result.exit_code == 0
    receipts = [json.loads(line) for line in result.stdout.splitlines()]
    assert [receipt["db"] for receipt in receipts] == ["card", "vfdb"]
    assert [receipt["records"] for receipt in receipts] == [2, 2]
    for name in ("card", "vfdb"):
        assert (datadir / name / "gapit-manifest.json").is_file()


def test_db_fetch_from_source_flag_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given any registry, When `db fetch NAME --from-source`, Then the
    command is a usage error (exit 2, nothing on stdout) — the flag was
    removed with the bundled snapshots it used to bypass."""
    datadir = patch_registry(tmp_path, monkeypatch)

    result = fetch(SYN, "--datadir", str(datadir), "--from-source")

    assert result.exit_code == 2
    assert result.stdout == ""


def test_db_list_json_surfaces_provider_licenses(tmp_path: Path) -> None:
    """Given the REAL registry and an empty datadir, When `db list --json`,
    Then the gapit.dblist/1 entries carry a license field exactly for the
    providers that pin one — card (McMaster terms), vfdb + ecoli_vf
    (CC BY-NC), the seven kaptive cluster providers (GPL) — and omit it for
    the rest; every entry also names its upstream maintainer in vendor."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = runner.invoke(app, ["db", "list", "--datadir", str(datadir), "--json"])

    assert result.exit_code == 0
    by_name = {entry["name"]: entry for entry in json.loads(result.stdout)["providers"]}
    assert by_name["card"]["license"].startswith("Custom (McMaster University)")
    assert by_name["vfdb"]["license"] == "CC BY-NC 4.0 (non-commercial)"
    assert by_name["ecoli_vf"]["license"] == "CC BY-NC 4.0 (VFDB-derived content)"
    assert by_name["kpsc_k"]["license"] == "GPL-3.0 (database content)"
    assert "license" not in by_name["ncbi"]
    assert "license" not in by_name["argannot"]
    assert by_name["kpsc_k"]["vendor"] == "Kaptive (klebgenomics)"
    assert by_name["ab_k"]["vendor"] == "Kaptive (Kenyon lab)"
    assert by_name["ecoli_kps"]["vendor"] == "Kaptive (Gladstone lab)"
    assert by_name["ncbi"]["vendor"] == "NCBI"
    assert all(entry["vendor"] for entry in by_name.values())
