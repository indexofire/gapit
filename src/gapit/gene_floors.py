"""Per-gene identity floors for reads-mode presence (gapit.floors/1).

A gene database may carry a ``floors.json`` sidecar declaring a MINIMUM
alignment identity (%) per gene — the rightsholder-curated answer to
SPATE-family homologs that pile onto a virulence gene at ~86% identity with
~97% breadth: breadth-only presence over-calls them (26ECO0084 pic), while
the assembly/blastn path already filters on identity. Floors are DATABASE
metadata — the db-driven opt-in: a database without the sidecar screens
byte-identically to before, and the blastn contig path is untouched
(``--minid`` already gates it).

The file's presence IS the flag — no manifest field, no CLI switch:

``{"schema": "gapit.floors/1", "default": null, "genes": {"pic": 90.0}}``

``default`` applies to genes not listed (null = no floor). At build time
(``db build --floors FILE``) values are constrained to [0, 100] and every
listed gene must exist in the FASTA; malformed content raises
FLOORS_MALFORMED, a ghost gene FLOORS_UNKNOWN_GENE (the typing.json
precedent). In reads mode the gate drops sub-floor alignments BEFORE
aggregation (:func:`apply_gene_floors`), so breadth is recomputed from the
surviving rows — that collapse is the point.
"""

from collections.abc import Iterable
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gapit.db import Database
from gapit.dbcodec import decode_seqid
from gapit.errors import DatabaseError, InputError
from gapit.paf import PafRecord, alignment_identity

FLOORS_FILENAME = "floors.json"

_Pct = Annotated[float, Field(ge=0.0, le=100.0)]


class GeneFloors(BaseModel, frozen=True):
    """gapit.floors/1 — per-gene minimum alignment identity (%) for
    reads-mode presence. ``default`` covers genes not in ``genes`` (null =
    no floor); an empty document is a no-op, identical to no file."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.floors/1"] = Field(default="gapit.floors/1", alias="schema")
    default: _Pct | None = None
    genes: dict[str, _Pct] = Field(default_factory=dict)


def read_floors(path: Path) -> GeneFloors:
    """Parse one floors document; a missing/unreadable file is the typed
    INPUT_NOT_FOUND input error, malformed content (bad JSON, wrong schema
    tag, values outside [0, 100]) the FLOORS_MALFORMED DatabaseError — the
    typing.json error precedent, file in context."""
    if not path.is_file():
        raise InputError(
            f"floors file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return GeneFloors.model_validate_json(path.read_text(encoding="utf-8"))
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed floors document {path}: {' '.join(str(exc).split())}",
            code="FLOORS_MALFORMED",
            context={"file": str(path)},
        ) from exc


def gene_floors(db: Database) -> GeneFloors | None:
    """The lazy loader for a database's floors sidecar: None when the
    database ships none (today's behavior for every existing db), the parsed
    :class:`GeneFloors` otherwise. Reads mode consults this once per run —
    discovery stays cheap for the floorless majority."""
    path = db.path / FLOORS_FILENAME
    return read_floors(path) if path.is_file() else None


def validate_floors(floors: GeneFloors, genes: frozenset[str], *, source: Path) -> None:
    """``db build --floors``: every listed gene must exist in the FASTA
    records — a floor for a gene the database does not carry is silently
    dead safety config, so it is a hard FLOORS_UNKNOWN_GENE error (the
    TYPING_UNKNOWN_GENE precedent), raised before any artifact is written."""
    unknown = sorted(set(floors.genes) - genes)
    if unknown:
        raise DatabaseError(
            f"floors document {source} references unknown gene(s)"
            f" (not in the FASTA records): {', '.join(unknown)}",
            code="FLOORS_UNKNOWN_GENE",
            context={"genes": ",".join(unknown)},
        )


def apply_gene_floors(
    rows: Iterable[PafRecord], floors: GeneFloors, *, default_db: str
) -> list[PafRecord]:
    """Reads-mode presence gate: drop alignments whose target gene's
    per-alignment identity (:func:`gapit.paf.alignment_identity`) is below
    the gene's floor — ``genes.get(gene, default)``; None means no floor and
    the row survives verbatim. Runs BEFORE aggregation (breadth/depth/reads
    are then computed from the survivors), independent of the gapit.reads/2
    ``--min-identity`` filter — a row must clear both when both apply
    (max-of-both). The gene is decoded from the target seqid once per
    unique tname (duplicate records of one gene share the floor)."""
    per_target: dict[str, float | None] = {}

    def floor_of(tname: str) -> float | None:
        if tname not in per_target:
            per_target[tname] = floors.genes.get(
                decode_seqid(tname, default_db).gene, floors.default
            )
        return per_target[tname]

    return [
        row
        for row in rows
        if (floor := floor_of(row.tname)) is None or alignment_identity(row) >= floor
    ]
