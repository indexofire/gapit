"""Per-database cluster-screening profile (``cluster.json``).

The shared cluster engine defaults to the historical primary-only coverage
and its long-standing 90/90/96 thresholds.  A cluster database may opt into
a different contract with a small ``cluster.json`` sidecar; databases
without the file keep the historical behavior byte-for-byte.

This is deliberately a database-side policy, not a global engine switch:
Kaptive-compatible panels can choose secondary-inclusive coverage while a
consumer-curated cluster panel (for example ``vch_oagc``) keeps the
primary-only behavior.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gspec.errors import DatabaseError

CLUSTER_PROFILE_SCHEMA = "gapit.cluster_profile/1"
CLUSTER_PROFILE_FILENAME = "cluster.json"
CoverageMode = Literal["primary", "primary_secondary"]


class ClusterProfile(BaseModel, frozen=True):
    """The database-owned cluster-engine policy.

    ``coverage_mode`` controls whether cluster math consumes only minimap2
    primary rows (the historical default) or the union of primary and
    secondary rows (Kaptive-style coverage).  The optional threshold fields
    override the CLI defaults only while the caller leaves the corresponding
    flag at its historical default; an explicit non-default CLI value wins.
    ``locus_typeable_coverage`` / ``locus_typeable_identity`` (with
    ``fallback_evidence``) classify a detected locus as Kaptive-style
    typeable/untypeable without changing the phenotype call.  The
    ``locus_quality_*`` fields configure the native coverage/identity/gene-count
    quality verdict (annotate-only in the prototype).
    """

    model_config = ConfigDict(populate_by_name=True, extra="forbid")

    schema_name: Literal["gapit.cluster_profile/1"] = Field(
        default="gapit.cluster_profile/1", alias="schema"
    )
    coverage_mode: CoverageMode = "primary"
    min_gene_cov: float | None = Field(default=None, ge=0.0, le=100.0)
    min_gene_id: float | None = Field(default=None, ge=0.0, le=100.0)
    min_cluster_cov: float | None = Field(default=None, ge=0.0, le=100.0)
    fallback_evidence: bool = False
    locus_typeable_coverage: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_typeable_identity: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_quality_min_coverage: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_quality_min_identity: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_quality_max_missing_genes: int | None = Field(default=None, ge=0)
    locus_quality_max_extra_panel_genes: int | None = Field(default=None, ge=0)


def default_cluster_profile() -> ClusterProfile:
    return ClusterProfile()


def load_cluster_profile(db_dir: Path) -> ClusterProfile:
    """Read ``<db_dir>/cluster.json`` or return the historical default.

    A malformed profile is a typed database error, never a silent fallback:
    the file changes coverage semantics and must not be partially applied.
    """
    path = db_dir / CLUSTER_PROFILE_FILENAME
    if not path.is_file():
        return default_cluster_profile()
    try:
        return ClusterProfile.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"cannot read cluster profile {path}: {exc}",
            code="CLUSTER_PROFILE_MALFORMED",
            context={"db": db_dir.name, "file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed cluster profile {path}: {exc}",
            code="CLUSTER_PROFILE_MALFORMED",
            context={"db": db_dir.name, "file": str(path)},
        ) from exc


__all__ = [
    "CLUSTER_PROFILE_FILENAME",
    "CLUSTER_PROFILE_SCHEMA",
    "ClusterProfile",
    "CoverageMode",
    "default_cluster_profile",
    "load_cluster_profile",
]
