"""CLI integration tests: `gapit db outdated` and `gapit db search`.

Offline by construction: fixtures are hand-written manifest + records.jsonl
databases (no binaries run — neither command builds anything). The
snapshot-update tests read the REAL bundled card snapshot in-memory, so the
packaged tar doubles as fixture: its manifest fetched_at (2026-09-17) is
frozen, which makes "installed older/newer than the bundle" deterministic
for any test run date.
"""

import json
from datetime import UTC, datetime, timedelta
from pathlib import Path

from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.records import Manifest, Record, write_manifest, write_records

# The real bundled card snapshot's manifest fetched_at (frozen in the archive).
CARD_SNAPSHOT_FETCHED_AT = "2026-09-17T23:14:25Z"
SEQ = "ACGTAG" * 10  # 60 bp, matches the synthetic-fixture convention

runner = CliRunner()


def make_db(
    datadir: Path,
    name: str,
    fetched_at: str,
    records: tuple[Record, ...] | None = None,
    *,
    upstream_version: str = "",
) -> None:
    """Hand-write an installed database: manifest always, records.jsonl only
    when ``records`` is given (None leaves the dir incomplete on purpose)."""
    db_dir = datadir / name
    db_dir.mkdir()
    if records is not None:
        write_records(records, db_dir / "records.jsonl")
    write_manifest(
        Manifest(
            name=name,
            source_urls=(),
            fetched_at=fetched_at,
            sha256="0" * 64,
            n_records=len(records) if records is not None else 0,
            dbtype="nucl",
            upstream_version=upstream_version,
        ),
        db_dir / "gapit-manifest.json",
    )


def rec(
    db: str,
    gene: str,
    *,
    accession: str = "",
    function: tuple[str, ...] = (),
    product: str = "n/a",
) -> Record:
    """One synthetic record (60 bp sequence)."""
    return Record(
        db=db, gene=gene, sequence=SEQ, accession=accession, function=function, product=product
    )


def rows(stdout: str) -> dict[str, list[str]]:
    """Data rows of a TSV table keyed by the first column (header skipped)."""
    return {line.split("\t")[0]: line.split("\t") for line in stdout.splitlines()[1:]}


def last_envelope(stderr: str) -> dict[str, object]:
    """The last non-empty stderr line parsed as the gapit.error/1 envelope."""
    lines = [line for line in stderr.splitlines() if line.strip()]
    assert lines, "expected an error envelope on stderr"
    return json.loads(lines[-1])


def outdated(*extra: str) -> Result:
    return runner.invoke(app, ["db", "outdated", *extra])


def search(term: str, *extra: str) -> Result:
    return runner.invoke(app, ["db", "search", term, *extra])


# ---------------------------------------------------------------- outdated --


def test_db_outdated_table_marks_ok_and_stale(tmp_path: Path) -> None:
    """Given one fresh db and old dbs (one provider, one not), When
    `db outdated`, Then exit 0 with a NAME/FETCHED_AT/AGE_DAYS/STATUS table:
    fresh -> ok, old -> stale (ncbi proves a provider without a bundled
    snapshot never reports snapshot-update), and AGE_DAYS agrees with UTC
    now-minus-fetched_at math."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fresh_at = (datetime.now(UTC) - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    make_db(datadir, "tinyamr", fresh_at, ())  # custom db, no provider
    make_db(datadir, "zold", "2020-01-01T00:00:00Z", ())  # not a provider name
    make_db(datadir, "ncbi", "2020-01-01T00:00:00Z", ())  # provider, no snapshot

    result = outdated("--datadir", str(datadir))

    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0] == "NAME\tFETCHED_AT\tAGE_DAYS\tSTATUS"
    table = rows(result.stdout)
    assert table["tinyamr"][1] == fresh_at
    assert table["tinyamr"][3] == "ok"
    assert 0.99 <= float(table["tinyamr"][2]) <= 2.0
    assert table["zold"][1] == "2020-01-01T00:00:00Z"
    assert table["zold"][3] == "stale"
    assert table["ncbi"][3] == "stale"
    expected = (
        datetime.now(UTC) - datetime.fromisoformat("2020-01-01T00:00:00Z")
    ).total_seconds() / 86400
    assert abs(float(table["zold"][2]) - round(expected, 2)) < 0.1


def test_db_outdated_json_document(tmp_path: Path) -> None:
    """Given two installed dbs, When `db outdated --json`, Then stdout parses
    as ONE gapit.dboutdated/1 document whose per-db entries carry
    db/fetched_at/age_days/status/upstream_version."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fresh_at = (datetime.now(UTC) - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    make_db(datadir, "tinyamr", fresh_at, (), upstream_version="v1")
    make_db(datadir, "zold", "2020-01-01T00:00:00Z", ())

    result = outdated("--datadir", str(datadir), "--json")

    assert result.exit_code == 0
    document = json.loads(result.stdout)
    assert document["schema"] == "gapit.dboutdated/1"
    by_db = {entry["db"]: entry for entry in document["databases"]}
    assert by_db["tinyamr"]["status"] == "ok"
    assert by_db["tinyamr"]["upstream_version"] == "v1"
    assert by_db["tinyamr"]["fetched_at"] == fresh_at
    assert abs(by_db["tinyamr"]["age_days"] - 1.0) < 0.1
    assert by_db["zold"]["status"] == "stale"
    assert by_db["zold"]["upstream_version"] == ""


def test_db_outdated_days_zero_marks_everything_stale(tmp_path: Path) -> None:
    """Given a one-day-old database, When `db outdated --days 0`, Then its
    status is stale — any positive age exceeds a zero threshold."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fresh_at = (datetime.now(UTC) - timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%SZ")
    make_db(datadir, "tinyamr", fresh_at, ())

    result = outdated("--datadir", str(datadir), "--days", "0")

    assert result.exit_code == 0
    assert rows(result.stdout)["tinyamr"][3] == "stale"


def test_db_outdated_days_threshold_boundary(tmp_path: Path) -> None:
    """Given a 10-day-old database, When `db outdated --days 11/--days 9`,
    Then the statuses are ok / stale — stale means age strictly past the
    threshold, so clearly separated ages avoid timing skew."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    fetched_at = (datetime.now(UTC) - timedelta(days=10)).strftime("%Y-%m-%dT%H:%M:%SZ")
    make_db(datadir, "tinyamr", fetched_at, ())

    ok_result = outdated("--datadir", str(datadir), "--days", "11")
    stale_result = outdated("--datadir", str(datadir), "--days", "9")

    assert rows(ok_result.stdout)["tinyamr"][3] == "ok"
    assert rows(stale_result.stdout)["tinyamr"][3] == "stale"


def test_db_outdated_snapshot_update_when_bundled_is_newer(tmp_path: Path) -> None:
    """Given an installed card older than the bundled snapshot (real REGISTRY,
    real packaged tar), When `db outdated`, Then card reports
    stale+snapshot-update while a custom db of the same age reports plain
    stale (no provider -> no bundle to compare)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    make_db(datadir, "card", "2020-01-01T00:00:00Z", ())
    make_db(datadir, "tinyamr", "2020-01-01T00:00:00Z", ())

    result = outdated("--datadir", str(datadir))

    assert result.exit_code == 0
    table = rows(result.stdout)
    assert table["card"][3] == "stale+snapshot-update"
    assert table["tinyamr"][3] == "stale"


def test_db_outdated_no_snapshot_update_when_installed_is_newer(tmp_path: Path) -> None:
    """Given an installed card fetched AFTER the bundled snapshot date, When
    `db outdated`, Then the status is plain ok — the bundle is not newer, so
    no snapshot-update flag even for a snapshot-backed provider."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    newer_than_bundle = (
        datetime.now(UTC) - timedelta(days=1)
    ).strftime(  # any run date > 2026-09-18 keeps this newer than the bundle
        "%Y-%m-%dT%H:%M:%SZ"
    )
    assert newer_than_bundle > CARD_SNAPSHOT_FETCHED_AT  # fixture premise
    make_db(datadir, "card", newer_than_bundle, ())

    result = outdated("--datadir", str(datadir))

    assert result.exit_code == 0
    assert rows(result.stdout)["card"][3] == "ok"


def test_db_outdated_snapshot_update_only_when_stale_threshold_not_hit(
    tmp_path: Path,
) -> None:
    """Given an installed card older than the bundle but within --days, When
    `db outdated --days 36500` (a century), Then the status is
    snapshot-update alone — the two flags combine but stay independent."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    # 2026-09-17T00:00:00Z: older than the bundle (23:14 same day), ~3 days
    # old at test time — inside the 36500-day threshold.
    make_db(datadir, "card", "2026-09-17T00:00:00Z", ())

    result = outdated("--datadir", str(datadir), "--days", "36500")

    assert result.exit_code == 0
    assert rows(result.stdout)["card"][3] == "snapshot-update"


def test_db_outdated_empty_datadir_is_database_error(tmp_path: Path) -> None:
    """Given a datadir with no manifest-carrying subdirectory, When
    `db outdated`, Then exit 4 with a DATADIR_EMPTY envelope and no stdout."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()

    result = outdated("--datadir", str(datadir))

    assert result.exit_code == 4
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope["code"] == "DATADIR_EMPTY"


def test_db_outdated_missing_datadir_is_database_error(tmp_path: Path) -> None:
    """Given a --datadir that does not exist, When `db outdated`, Then exit 4
    with DATADIR_NOT_FOUND (the shared resolve_datadir semantics)."""
    result = outdated("--datadir", str(tmp_path / "nope"))

    assert result.exit_code == 4
    assert last_envelope(result.stderr)["code"] == "DATADIR_NOT_FOUND"


def test_db_outdated_malformed_fetched_at_is_manifest_malformed(
    tmp_path: Path,
) -> None:
    """Given a manifest whose fetched_at is not an ISO-8601 timestamp, When
    `db outdated`, Then exit 5 with MANIFEST_MALFORMED (timestamp parsing is
    part of reading the manifest)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    make_db(datadir, "tinyamr", "not-a-timestamp", ())

    result = outdated("--datadir", str(datadir))

    assert result.exit_code == 5
    envelope = last_envelope(result.stderr)
    assert envelope["code"] == "MANIFEST_MALFORMED"


# ------------------------------------------------------------------ search --


def search_datadir(tmp_path: Path) -> Path:
    """Two complete databases (alpha, beta) + one incomplete (gamma:
    manifest but no records.jsonl) + a legacy abricate-style dir (delta:
    sequences only, invisible to search)."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    make_db(
        datadir,
        "alpha",
        "2020-01-01T00:00:00Z",
        (
            rec("alpha", "sul1", accession="U12338.4:1-940", product="sulfonamide resistance"),
            rec(
                "alpha",
                "ermB",
                function=("macrolide", "lincosamide"),
                product="ERYthromycin resistance",
            ),
            rec("alpha", "blaCTX-M", product="class A beta-lactamase"),
        ),
    )
    make_db(
        datadir,
        "beta",
        "2020-01-01T00:00:00Z",
        (rec("beta", "SUL2", accession="X99045.1:1-816"),),
    )
    make_db(datadir, "gamma", "2020-01-01T00:00:00Z", None)
    legacy = datadir / "delta"
    legacy.mkdir()
    (legacy / "sequences").write_text(">old~~~sul3\nACGT\n", encoding="utf-8")
    return datadir


def test_db_search_substring_matches_across_dbs(tmp_path: Path) -> None:
    """Given two complete databases, When `db search sul`, Then case- and
    db-insensitive substring rows stream on stdout in db-then-file order with
    DB/GENE/ACCESSION/FUNCTION/PRODUCT/LENGTH columns."""
    datadir = search_datadir(tmp_path)

    result = search("sul", "--datadir", str(datadir))

    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        f"alpha\tsul1\tU12338.4:1-940\t\tsulfonamide resistance\t{len(SEQ)}",
        f"beta\tSUL2\tX99045.1:1-816\t\tn/a\t{len(SEQ)}",
    ]


def test_db_search_exact_full_field_equality(tmp_path: Path) -> None:
    """Given a gene named sul1, When `db search SUL1 --exact` / `db search sul
    --exact`, Then only the full-field (case-insensitive) equality hits —
    substring near-misses print nothing and still exit 0."""
    datadir = search_datadir(tmp_path)

    hit = search("SUL1", "--datadir", str(datadir), "--exact")
    miss = search("sul", "--datadir", str(datadir), "--exact")

    assert hit.exit_code == 0
    assert hit.stdout.splitlines() == [
        f"alpha\tsul1\tU12338.4:1-940\t\tsulfonamide resistance\t{len(SEQ)}"
    ]
    assert miss.exit_code == 0
    assert miss.stdout == ""


def test_db_search_field_gene_vs_product(tmp_path: Path) -> None:
    """Given a term present only in a PRODUCT, When `db search` with
    --field gene vs --field product, Then only the product search hits."""
    datadir = search_datadir(tmp_path)

    gene = search("lactamase", "--datadir", str(datadir), "--field", "gene")
    product = search("lactamase", "--datadir", str(datadir), "--field", "product")

    assert gene.exit_code == 0
    assert gene.stdout == ""
    assert product.exit_code == 0
    assert product.stdout.startswith("alpha\tblaCTX-M\t")


def test_db_search_field_function_items(tmp_path: Path) -> None:
    """Given a record with function (macrolide, lincosamide), When
    `db search Macrolide --field function --exact`, Then it hits — matching
    any ONE function item with full-item equality, case-insensitively."""
    datadir = search_datadir(tmp_path)

    result = search("Macrolide", "--datadir", str(datadir), "--field", "function", "--exact")

    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        f"alpha\termB\t\tmacrolide;lincosamide\tERYthromycin resistance\t{len(SEQ)}"
    ]


def test_db_search_field_accession(tmp_path: Path) -> None:
    """Given an accession-bearing record, When `db search u12338 --field
    accession`, Then the substring hit's row carries the accession."""
    datadir = search_datadir(tmp_path)

    result = search("u12338", "--datadir", str(datadir), "--field", "accession")

    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        f"alpha\tsul1\tU12338.4:1-940\t\tsulfonamide resistance\t{len(SEQ)}"
    ]


def test_db_search_field_any_includes_all_fields(tmp_path: Path) -> None:
    """Given terms living in different fields, When `db search` with default
    --field any, Then gene, accession, function-item, and product terms all
    hit."""
    datadir = search_datadir(tmp_path)

    for term in ("ctx", "x99045", "lincosamide", "ERYTHRO"):
        result = search(term, "--datadir", str(datadir))
        assert result.exit_code == 0, term
        assert result.stdout != "", term


def test_db_search_db_subset_and_unknown_db(tmp_path: Path) -> None:
    """Given two complete databases, When `db search --db beta`, Then only
    beta streams; When `--db nosuch`, Then exit 2 USAGE_ERROR whose message
    lists the installed database names."""
    datadir = search_datadir(tmp_path)

    subset = search("sul", "--datadir", str(datadir), "--db", "beta")
    unknown = search("sul", "--datadir", str(datadir), "--db", "nosuch")

    assert subset.exit_code == 0
    assert subset.stdout.startswith("beta\tSUL2\t")
    assert unknown.exit_code == 2
    assert unknown.stdout == ""
    envelope = last_envelope(unknown.stderr)
    assert envelope["code"] == "USAGE_ERROR"
    message = str(envelope["message"])
    assert "alpha" in message and "beta" in message and "gamma" in message


def test_db_search_limit_and_truncation_note(tmp_path: Path) -> None:
    """Given three matching records in one complete database, When
    `db search --limit 2`, Then two rows stream and stderr notes the
    truncation; When `--limit 3` (exactly the hit count) and `--limit 0`
    (unlimited), Then all rows stream and no note appears."""
    datadir = tmp_path / "datadir"
    datadir.mkdir()
    make_db(
        datadir,
        "alpha",
        "2020-01-01T00:00:00Z",
        (rec("alpha", "sul1"), rec("alpha", "sul2"), rec("alpha", "sul3")),
    )

    truncated = search("sul", "--datadir", str(datadir), "--limit", "2")
    exact = search("sul", "--datadir", str(datadir), "--limit", "3")
    unlimited = search("sul", "--datadir", str(datadir), "--limit", "0")

    assert truncated.exit_code == 0
    assert len(truncated.stdout.splitlines()) == 2
    note_lines = [line for line in truncated.stderr.splitlines() if line.startswith("gapit:")]
    assert any("truncated" in line for line in note_lines)
    assert exact.exit_code == 0
    assert len(exact.stdout.splitlines()) == 3
    assert not any("truncated" in line for line in exact.stderr.splitlines())
    assert unlimited.exit_code == 0
    assert len(unlimited.stdout.splitlines()) == 3
    assert not any("truncated" in line for line in unlimited.stderr.splitlines())


def test_db_search_zero_hits_empty_stdout(tmp_path: Path) -> None:
    """Given a term matching nothing, When `db search`, Then exit 0 with a
    completely empty stdout."""
    datadir = search_datadir(tmp_path)

    result = search("zzznope", "--datadir", str(datadir))

    assert result.exit_code == 0
    assert result.stdout == ""


def test_db_search_skips_incomplete_db_with_warning(tmp_path: Path) -> None:
    """Given an installed db without records.jsonl (gamma) beside a complete
    one, When `db search` scans all, Then gamma is skipped with a stderr
    warning naming it and alpha's hits still stream; When `--quiet`, Then the
    warning is silenced."""
    datadir = search_datadir(tmp_path)

    result = search("sul", "--datadir", str(datadir))
    quiet = search("sul", "--datadir", str(datadir), "--quiet")

    assert result.exit_code == 0
    assert result.stdout.startswith("alpha\tsul1\t")
    warnings = [line for line in result.stderr.splitlines() if line.startswith("gapit:")]
    assert any("gamma" in line for line in warnings)
    assert quiet.exit_code == 0
    assert "gapit:" not in quiet.stderr


def test_db_search_explicit_incomplete_db_errors(tmp_path: Path) -> None:
    """Given gamma (manifest but no records.jsonl), When `db search --db
    gamma`, Then exit 4 with DB_INCOMPLETE and context {db: gamma}."""
    datadir = search_datadir(tmp_path)

    result = search("sul", "--datadir", str(datadir), "--db", "gamma")

    assert result.exit_code == 4
    assert result.stdout == ""
    envelope = last_envelope(result.stderr)
    assert envelope["code"] == "DB_INCOMPLETE"
    assert envelope["context"] == {"db": "gamma"}


def test_db_search_json_emits_jsonl_rows(tmp_path: Path) -> None:
    """Given matching records, When `db search --json`, Then stdout is one
    JSON object per line with snake_case db/gene/accession/function (array)/
    product/length fields."""
    datadir = search_datadir(tmp_path)

    result = search("ermB", "--datadir", str(datadir), "--json")

    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert len(lines) == 1
    entry = json.loads(lines[0])
    assert entry == {
        "db": "alpha",
        "gene": "ermB",
        "accession": "",
        "function": ["macrolide", "lincosamide"],
        "product": "ERYthromycin resistance",
        "length": len(SEQ),
    }
