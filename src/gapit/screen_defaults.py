"""Per-database screening defaults (gapit.screen/1).

A database may ship a ``screen.json`` sidecar declaring the thresholds that
make sense for ITS reference set. Sequence length, allele diversity, and the
presence of a control gene all change what "a good hit" means: a 62 bp
amplicon needs a higher identity floor than an 857 bp one, and a panel with
an lpxA control can gate presence on relative depth.

The file is optional. A database without it keeps the global CLI defaults
byte-identically (the floorless-database precedent in
:mod:`gapit.gene_floors`). CLI flags always override the database defaults;
only options the user did NOT pass are filled from this document.

Schema (all keys optional)::

    {
      "schema": "gapit.screen/1",
      "reads": {
        "min_breadth": 90.0,
        "min_identity": 80.0,
        "min_mapq": 20,
        "min_depth": 5.0,
        "min_relative_depth": 0.10,
        "relative_to": "lpxA"
      },
      "contigs": {
        "minid": 80.0,
        "mincov": 80.0,
        "merge_fragments": true
      }
    }

``min_relative_depth`` is a ratio (0..1), not a percentage: 0.10 means the
gene's mean depth must be >= 10% of the ``relative_to`` control gene's mean
depth. Per-gene identity floors stay in ``floors.json``
(:mod:`gapit.gene_floors`); this sidecar owns the mode-wide defaults.
"""

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gapit import config
from gapit.db import Database
from gapit.errors import DatabaseError, InputError

SCREEN_FILENAME = "screen.json"

_Pct = Annotated[float, Field(ge=0.0, le=100.0)]
_Ratio = Annotated[float, Field(ge=0.0, le=1.0)]


class ReadsScreenDefaults(BaseModel, frozen=True):
    """Reads-mode defaults for one database; None = use the CLI default."""

    min_breadth: _Pct | None = None
    min_identity: _Pct | None = None
    min_mapq: int | None = Field(default=None, ge=0)
    min_depth: float | None = Field(default=None, ge=0.0)
    min_relative_depth: _Ratio | None = None
    relative_to: str | None = None


class ContigScreenDefaults(BaseModel, frozen=True):
    """Contig-mode defaults for one database; None = use the CLI default."""

    minid: _Pct | None = None
    mincov: _Pct | None = None
    merge_fragments: bool | None = None


class ScreenDefaults(BaseModel, frozen=True):
    """gapit.screen/1 — one database's mode-wide screening defaults."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.screen/1"] = Field(default="gapit.screen/1", alias="schema")
    reads: ReadsScreenDefaults = Field(default_factory=ReadsScreenDefaults)
    contigs: ContigScreenDefaults = Field(default_factory=ContigScreenDefaults)


def read_screen_defaults(path: Path) -> ScreenDefaults:
    """Parse one screen-defaults document; missing file = INPUT_NOT_FOUND,
    malformed content = SCREEN_DEFAULTS_MALFORMED (typing/floors precedent)."""
    if not path.is_file():
        raise InputError(
            f"screen defaults file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    try:
        return ScreenDefaults.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"could not read screen defaults {path}: {exc}",
            code="SCREEN_DEFAULTS_MALFORMED",
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed screen defaults {path}: {' '.join(str(exc).split())}",
            code="SCREEN_DEFAULTS_MALFORMED",
            context={"file": str(path)},
        ) from exc


def validate_screen_defaults(
    defaults: ScreenDefaults, genes: frozenset[str], *, source: Path
) -> None:
    """Build-time guard: a relative-depth control must exist in the database.

    A typo in ``relative_to`` would otherwise make every call absent at
    screen time (the control-missing rule), which looks like a silent
    false-negative panel. A relative gate without a control is dead config,
    so both are hard errors before any artifact is written.
    """
    if defaults.reads.min_relative_depth is None:
        return
    control = defaults.reads.relative_to or "lpxA"
    if control not in genes:
        raise DatabaseError(
            f"screen defaults {source} names unknown relative-depth control"
            f" {control!r} (not in the FASTA records)",
            code="SCREEN_DEFAULTS_UNKNOWN_GENE",
            context={"gene": control, "file": str(source)},
        )


def screen_defaults(db: Database) -> ScreenDefaults | None:
    """The lazy loader for an installed database's screen.json sidecar."""
    path = db.path / SCREEN_FILENAME
    return read_screen_defaults(path) if path.is_file() else None


def screen_defaults_for(db_name: str, datadir: Path | None) -> ScreenDefaults | None:
    """Resolve a database's screen defaults before it is opened.

    Installed databases win; a bundled database that has not been
    materialized yet is read from the wheel's copy (the zero-network first-use
    path). A malformed sidecar is never swallowed.
    """
    installed = config.datadir_path(datadir) / db_name / SCREEN_FILENAME
    if installed.is_file():
        return read_screen_defaults(installed)
    # Lazy import: bundled -> db_build_ops -> screen_defaults would otherwise
    # be a module-level cycle.
    from gapit.bundled import find_bundled

    bundled = find_bundled(db_name)
    if bundled is not None:
        bundled_path = bundled.path / SCREEN_FILENAME
        if bundled_path.is_file():
            return read_screen_defaults(bundled_path)
    return None
