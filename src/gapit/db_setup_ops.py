"""Database setup use-case: bootstrap the datadir and rebuild BLAST indices.

``perform_setup`` is the behavior behind ``gapit db setup`` (CLI) and can be
called by an agent-facing layer without duplicating the pipeline: create a
missing datadir, materialize bundled databases that are not installed yet,
then rebuild the index for every discovered database. The returned
:class:`DbSetupDocument` is the machine-readable receipt (``gapit.dbsetup/1``)
emitted by ``db setup --json``.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit import config, db
from gapit.bundled import bundled_names, materialize_bundled
from gapit.db import validate_database_name
from gapit.errors import DatabaseError
from gapit.proctools import note


class DbSetupEntry(BaseModel, frozen=True):
    """One database that ``db setup`` indexed."""

    db: str
    records: int
    dbtype: Literal["nucl", "prot"]


class DbSetupDocument(BaseModel, frozen=True):
    """gapit.dbsetup/1 — the ``gapit db setup --json`` receipt."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.dbsetup/1"] = Field(default="gapit.dbsetup/1", alias="schema")
    datadir: str
    materialized: tuple[str, ...]
    fetched: tuple[str, ...]
    databases: tuple[DbSetupEntry, ...]


def perform_setup(
    datadir: Path | None,
    *,
    debug: bool = False,
    quiet: bool = False,
    bundles: bool = True,
    fetch: bool = True,
    names: Sequence[str] | None = None,
) -> DbSetupDocument:
    """Bootstrap and index the datadir, returning the setup receipt.

    With no ``names``, every discovered database is set up: ``bundles``
    materializes wheel-shipped databases whose manifest is absent (an
    installed bundle is left untouched), then every database is reindexed in
    place. With ``names``, only those databases are touched: an installed
    name is reindexed, a missing bundled name is materialized first (unless
    ``bundles`` is false), and a missing provider name is fetched when
    ``fetch`` is true. A name that is none of those is a
    ``DATABASE_NOT_FOUND``. Progress goes to stderr through
    :func:`gapit.proctools.note` and is silenced by ``quiet``.
    """
    root = config.ensure_datadir(datadir)
    selected = list(names) if names else None
    materialized: list[str] = []
    fetched: list[str] = []

    if selected is None:
        if bundles:
            for bundled_name in bundled_names():
                if not (root / bundled_name / "gapit-manifest.json").is_file():
                    materialize_bundled(bundled_name, root, quiet=quiet)
                    materialized.append(bundled_name)
        infos = db.list_databases(root, reindex=True, debug=debug)
    else:
        # Deduplicate while preserving the caller's order; validate every
        # name before any write so a typo cannot leave a partial setup.
        selected = list(dict.fromkeys(selected))
        for name in selected:
            validate_database_name(name)
        from gapit.db_ops import perform_fetch
        from gapit.providers import REGISTRY

        installed = {database.name for database in db.discover_databases(root)}
        bundled = set(bundled_names())
        providers = set(REGISTRY)
        unknown = [
            name
            for name in selected
            if name not in installed and name not in bundled and name not in providers
        ]
        if unknown:
            available = sorted(installed | bundled | providers)
            raise DatabaseError(
                f"unknown database: {unknown[0]} (available: {', '.join(available) or '(none)'})",
                code="DATABASE_NOT_FOUND",
                context={"db": unknown[0], "datadir": str(root)},
            )
        for name in selected:
            if name in installed:
                continue
            if name in bundled and bundles:
                materialize_bundled(name, root, quiet=quiet)
                materialized.append(name)
                continue
            if name in providers:
                if not fetch:
                    raise DatabaseError(
                        f"database {name!r} is not installed; run: gapit db fetch {name}",
                        code="DATABASE_NOT_FOUND",
                        context={"db": name, "datadir": str(root)},
                    )
                for _receipt in perform_fetch(name, root, quiet=quiet, debug=debug):
                    pass
                fetched.append(name)
                continue
            raise DatabaseError(
                f"bundled database {name!r} is not installed (--no-bundles)",
                code="DATABASE_NOT_FOUND",
                context={"db": name, "datadir": str(root)},
            )
        infos = db.list_databases(root, reindex=True, debug=debug, names=frozenset(selected))
        indexed = {info.name for info in infos}
        missing = [name for name in selected if name not in indexed]
        if missing:
            raise DatabaseError(
                f"database {missing[0]!r} has no readable sequences file",
                code="DATABASE_NOT_FOUND",
                context={"db": missing[0], "datadir": str(root)},
            )

    for info in infos:
        note(quiet, f"indexed {info.name} ({info.n_sequences} sequences, {info.dbtype})")
    return DbSetupDocument(
        datadir=str(root),
        materialized=tuple(materialized),
        fetched=tuple(fetched),
        databases=tuple(
            DbSetupEntry(db=info.name, records=info.n_sequences, dbtype=info.dbtype)
            for info in infos
        ),
    )
