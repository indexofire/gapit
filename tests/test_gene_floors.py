"""Unit tests for gapit.floors/1 per-gene identity floors (reads-mode gate).

The gate contract (gapit.gene_floors): alignments to floored genes drop
sub-floor rows BEFORE aggregation, so breadth is recomputed from the
survivors; unfloored genes are untouched; the gapit.reads/2
``--min-identity`` filter and floors are conjunctive (a row must clear
both — max-of-both).
"""

import json
from pathlib import Path

import pytest

from gapit.db import Database
from gapit.errors import DatabaseError, InputError
from gapit.gene_floors import (
    GeneFloors,
    apply_gene_floors,
    gene_floors,
    read_floors,
    validate_floors,
)
from gapit.paf import PafRecord, filter_alignments
from gapit.reads import aggregate_coverage

GENE_A = "db~~~geneA~~~ACC~~~RES"
GENE_B = "db~~~geneB~~~ACC~~~RES"

ECOLI_DEC_FLOORS = '{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}'


def row(
    tname: str = GENE_A,
    *,
    identity_pct: float = 100.0,
    tstart: int = 0,
    tend: int = 50,
    tlen: int = 100,
    qname: str = "r1",
    nm: int | None = None,
) -> PafRecord:
    """A primary alignment on `tname` spanning [tstart, tend); identity is
    set through the NM tag (100*(alen-nm)/alen), None = no NM tag."""
    alen = tend - tstart
    mismatches = nm if nm is not None else round(alen * (1.0 - identity_pct / 100.0))
    return PafRecord(
        qname=qname,
        qlen=alen,
        qstart=0,
        qend=alen,
        strand="+",
        tname=tname,
        tlen=tlen,
        tstart=tstart,
        tend=tend,
        nmatch=alen - mismatches,
        alen=alen,
        mapq=60,
        nm=mismatches,
    )


def db_at(tmp_path: Path, floors_json: str | None) -> Database:
    """A gene Database whose directory optionally carries floors.json."""
    (tmp_path / "sequences").write_text(">x\nACGT\n", encoding="utf-8")
    if floors_json is not None:
        (tmp_path / "floors.json").write_text(floors_json, encoding="utf-8")
    return Database(name="db", path=tmp_path, sequences_path=tmp_path / "sequences")


# ------------------------------------------------------------ model + load --


def test_model_roundtrips_the_ecoli_dec_document() -> None:
    """Given the bundled ecoli_dec floors bytes, When parsed, Then the
    schema alias, null default, and the pic floor decode (and re-serialize
    by alias)."""
    floors = GeneFloors.model_validate_json(ECOLI_DEC_FLOORS)
    assert floors.schema_name == "gapit.floors/1"
    assert floors.default is None
    assert floors.genes == {"pic": 90.0}
    assert json.loads(floors.model_dump_json(by_alias=True)) == json.loads(ECOLI_DEC_FLOORS)


def test_read_floors_missing_file_is_input_not_found(tmp_path: Path) -> None:
    """Given a floors path that does not exist, When read, Then the typed
    INPUT_NOT_FOUND input error names the file (the --typing precedent)."""
    with pytest.raises(InputError) as exc_info:
        read_floors(tmp_path / "nope.json")
    assert exc_info.value.code == "INPUT_NOT_FOUND"
    assert exc_info.value.context["file"] == str(tmp_path / "nope.json")


@pytest.mark.parametrize(
    "text",
    [
        "{nope",  # broken JSON
        '{"schema": "gapit.floors/2", "genes": {}}',  # wrong tag
        '{"schema": "gapit.floors/1", "genes": {"pic": 100.5}}',  # above range
        '{"schema": "gapit.floors/1", "genes": {"pic": -1}}',  # below range
        '{"schema": "gapit.floors/1", "default": 101.0}',  # default above range
        '{"schema": "gapit.floors/1", "genes": {"pic": "high"}}',  # wrong type
    ],
)
def test_read_floors_bad_documents_are_floors_malformed(tmp_path: Path, text: str) -> None:
    """Given structurally or range-invalid floors content, When read, Then
    the FLOORS_MALFORMED DatabaseError carries the file in context."""
    path = tmp_path / "floors.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(DatabaseError) as exc_info:
        read_floors(path)
    assert exc_info.value.code == "FLOORS_MALFORMED"
    assert exc_info.value.context["file"] == str(path)


def test_gene_floors_loader_is_none_without_sidecar(tmp_path: Path) -> None:
    """Given a database directory without floors.json, When loaded, Then
    the loader yields None — the floorless majority keeps today's path."""
    assert gene_floors(db_at(tmp_path, floors_json=None)) is None


def test_gene_floors_loader_parses_the_sidecar(tmp_path: Path) -> None:
    """Given a database directory carrying floors.json, When loaded, Then
    the parsed document comes back (also for an explicit no-op document)."""
    floors = gene_floors(db_at(tmp_path, floors_json=ECOLI_DEC_FLOORS))
    assert floors is not None
    assert floors.genes == {"pic": 90.0}
    noop = gene_floors(db_at(tmp_path, floors_json='{"schema": "gapit.floors/1"}'))
    assert noop == GeneFloors()


def test_validate_floors_unknown_gene_is_typed_error(tmp_path: Path) -> None:
    """Given a floors document naming a gene the FASTA does not carry, When
    validated at build, Then FLOORS_UNKNOWN_GENE lists the ghosts (a floor
    for a missing gene is silently dead safety config)."""
    floors = GeneFloors(genes={"geneA": 90.0, "ghost": 90.0})
    with pytest.raises(DatabaseError) as exc_info:
        validate_floors(floors, frozenset({"geneA"}), source=tmp_path / "floors.json")
    assert exc_info.value.code == "FLOORS_UNKNOWN_GENE"
    assert exc_info.value.context["genes"] == "ghost"


def test_validate_floors_passes_when_every_gene_exists(tmp_path: Path) -> None:
    """Given floors whose genes all exist in the FASTA records, When
    validated, Then nothing raises."""
    validate_floors(
        GeneFloors(genes={"geneA": 90.0}), frozenset({"geneA", "geneB"}), source=tmp_path
    )


# ------------------------------------------------------------------- gate --


def test_gate_drops_sub_floor_rows_for_floored_gene() -> None:
    """Given geneA floored at 90 with one 86% and one 95% alignment, When
    gated, Then only the 95% row survives."""
    floors = GeneFloors(genes={"geneA": 90.0})
    kept = apply_gene_floors(
        [row(identity_pct=86.0), row(qname="r2", identity_pct=95.0)], floors, default_db="db"
    )
    assert [entry.qname for entry in kept] == ["r2"]


def test_gate_keeps_row_exactly_at_the_floor() -> None:
    """Given an alignment whose identity equals the floor, When gated, Then
    it survives (>= semantics, same as --min-identity)."""
    kept = apply_gene_floors(
        [row(identity_pct=90.0)], GeneFloors(genes={"geneA": 90.0}), default_db="db"
    )
    assert len(kept) == 1


def test_gate_leaves_unfloored_genes_untouched() -> None:
    """Given only geneA floored, When gated, Then geneB's low-identity rows
    pass verbatim (floors are per-gene, not global)."""
    rows = [row(tname=GENE_B, identity_pct=50.0), row(tname=GENE_A, identity_pct=95.0)]
    kept = apply_gene_floors(rows, GeneFloors(genes={"geneA": 90.0}), default_db="db")
    assert [entry.tname for entry in kept] == [GENE_B, GENE_A]


def test_gate_default_covers_unlisted_genes() -> None:
    """Given default=90 and geneA unlisted, When gated, Then geneA's 86%
    row drops while an explicit gene-genes floor of 50 would keep it."""
    floors = GeneFloors(default=90.0)
    assert apply_gene_floors([row(identity_pct=86.0)], floors, default_db="db") == []
    relaxed = apply_gene_floors(
        [row(identity_pct=86.0)], GeneFloors(genes={"geneA": 50.0}), default_db="db"
    )
    assert len(relaxed) == 1


def test_gate_row_without_nm_counts_as_identity_100() -> None:
    """Given a row without an NM tag (identity rule: counts as 100.0), When
    gated under a 90 floor, Then it survives — same boundary as reads/2."""
    kept = apply_gene_floors([row(nm=None)], GeneFloors(genes={"geneA": 90.0}), default_db="db")
    assert len(kept) == 1


def test_gate_duplicate_records_of_one_gene_share_the_floor() -> None:
    """Given two seqids carrying the same gene (duplicate panel records),
    When gated, Then the floor keys on the decoded gene and gates both."""
    acc2 = "db~~~geneA~~~ACC2~~~RES"
    kept = apply_gene_floors(
        [row(identity_pct=86.0), row(tname=acc2, qname="r2", identity_pct=86.0)],
        GeneFloors(genes={"geneA": 90.0}),
        default_db="db",
    )
    assert kept == []


def test_gate_recomputes_breadth_from_survivors() -> None:
    """Given geneA at 100% identity over [0,50) and 86% over [50,100) with a
    90 floor, When gated then aggregated, Then breadth collapses to 50%
    (absent at min_breadth 90) — the 26ECO0084 homolog rows vanish; without
    the gate the same rows aggregate to 100% breadth (present)."""
    rows = [row(tstart=0, tend=50), row(qname="r2", tstart=50, tend=100, identity_pct=86.0)]
    floors = GeneFloors(genes={"geneA": 90.0})

    gated = aggregate_coverage(
        apply_gene_floors(rows, floors, default_db="db"), default_db="db", min_breadth=90.0
    )
    (entry,) = gated
    assert entry.breadth_pct == 50.0
    assert entry.present is False

    ungated = aggregate_coverage(rows, default_db="db", min_breadth=90.0)
    (entry,) = ungated
    assert entry.breadth_pct == 100.0
    assert entry.present is True


def test_floors_and_min_identity_are_conjunctive() -> None:
    """Given a 90 floor and the reads/2 --min-identity 95 filter, When both
    apply (the screen_reads order), Then a row must clear the max of both:
    96% survives, 93% fails the global filter, 85% fails the floor."""
    rows = [
        row(qname="r96", identity_pct=96.0),
        row(qname="r93", identity_pct=93.0),
        row(qname="r85", identity_pct=85.0),
    ]
    filtered = filter_alignments(rows, min_identity=95.0)
    gated = apply_gene_floors(filtered, GeneFloors(genes={"geneA": 90.0}), default_db="db")
    assert [entry.qname for entry in gated] == ["r96"]
