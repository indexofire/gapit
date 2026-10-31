"""Optional gseek bridge for BLAST and minimap2 execution.

``GAPIT_ENGINE`` selects the implementation:

- ``legacy`` (default): the existing subprocess paths, no gseek import;
- ``gseek``: gseek plans and parses the same BLAST/minimap2 work;
- ``shadow``: run both, log coarse differences, and return the legacy result.

gseek is optional and requires Python >= 3.12. Importing this module never
imports gseek; only a ``gseek``/``shadow`` call does.
"""

from __future__ import annotations

import importlib
import inspect
import logging
import os
import sys
import tempfile
import threading
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Protocol, TypeVar, cast

from gapit.db import Database
from gapit.errors import DependencyError, GapitError
from gapit.report import ScreeningParams

if TYPE_CHECKING:
    from gseek.compat import Engine
    from gseek.execution import Executor, ToolResolver
    from gseek.index import IndexStore
    from gseek.model import InputRef

    from gapit.blast import BlastRow
    from gapit.paf import PafRecord

_LOGGER = logging.getLogger("gapit.gseek")
_ENGINE_ENV = "GAPIT_ENGINE"
_VALID_MODES = frozenset({"legacy", "gseek", "shadow"})
T = TypeVar("T")


@dataclass(frozen=True)
class _GseekSession:
    """Process-wide gseek tool/index cache.

    ``gseek.align`` builds a fresh ``ToolResolver`` per call unless one is
    supplied. Its BLAST adapter probes every declared binary (including the
    optional blastp/blastx/tblastn/tblastx) by running ``-version``; with one
    call per input file that dominated the gene path. Reusing these objects
    across files also keeps the content-addressed index cache warm.
    """

    resolver: ToolResolver
    executor: Executor
    store: IndexStore
    env_key: tuple[tuple[str, str], ...]


_session_cache: _GseekSession | None = None
_SESSION_LOCK = threading.Lock()
_TARGET_REFS: dict[tuple[Path, int, int], tuple[InputRef, ...]] = {}
_TARGET_REFS_LOCK = threading.Lock()


def batch_supported() -> bool:
    """Whether the optional gseek package exposes the batch API."""
    try:
        return supports_batch(_gseek())
    except DependencyError:
        return False


def supports_batch(gseek: object) -> bool:
    """Whether gseek exposes the multi-query API with per-file annotation."""
    align_many = getattr(gseek, "align_many", None)
    if not callable(align_many):
        return False
    try:
        return "annotate_batch" in inspect.signature(align_many).parameters
    except (TypeError, ValueError):
        return False


def _supports_batch_option(gseek: object, option: str) -> bool:
    """Whether the installed align_many accepts an optional keyword."""
    align_many = getattr(gseek, "align_many", None)
    if not callable(align_many):
        return False
    try:
        return option in inspect.signature(align_many).parameters
    except (TypeError, ValueError):
        return False


def supports_api_refs(gseek: object) -> bool:
    """Whether the installed gseek accepts precomputed query/target refs."""
    align = getattr(gseek, "align", None)
    if not callable(align):
        return False
    try:
        return "target_refs" in inspect.signature(align).parameters
    except (TypeError, ValueError):
        return False


def _supports_external_blast_db() -> bool:
    """Whether the installed BLAST adapter accepts ``params['db']``.

    This is a version capability, not a private-API dependency in the normal
    path: gseek <= 0.10.2 has no such parameter and the bridge simply keeps
    its old index/subject behavior.
    """
    try:
        blast_adapter = importlib.import_module("gseek.adapters.blast")
    except ImportError:
        return False
    allowed = cast("frozenset[str]", getattr(blast_adapter, "_ALLOWED_PARAMS", frozenset[str]()))
    return "db" in allowed


def cached_target_refs(database: Database) -> tuple[InputRef, ...] | None:
    """One content ref for the database target, cached by file identity.

    Consecutive per-file gseek calls otherwise re-hash the same target for
    the index key/manifest. Database files are immutable between builds, so
    (path, mtime_ns, size) is a safe cache key; a rebuild changes it.
    """
    path = database.sequences_path
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (path, stat.st_mtime_ns, stat.st_size)
    with _TARGET_REFS_LOCK:
        cached = _TARGET_REFS.get(key)
        if cached is not None:
            return cached
    try:
        from gseek.provenance import input_ref
        from gseek.types import InputKind
    except ImportError:
        return None
    ref = input_ref(path, kind=InputKind.FASTA, compute_records=False)
    refs = (ref,)
    with _TARGET_REFS_LOCK:
        for old_key in [entry for entry in _TARGET_REFS if entry[0] == path and entry != key]:
            del _TARGET_REFS[old_key]
        _TARGET_REFS[key] = refs
    return refs


def blast_db_ready(database: Database) -> bool:
    """Whether the gene database already carries a BLAST index prefix."""
    path = database.sequences_path
    return any((path.parent / f"{path.name}{suffix}").exists() for suffix in (".nin", ".pin"))


def _session() -> _GseekSession | None:
    """The shared gseek resolver/executor/store, or None without gseek.

    The env key includes PATH and every GSEEK_* variable, so a test or caller
    that changes tool discovery starts a fresh resolver rather than reusing
    stale probes. No gseek import happens in legacy mode.
    """
    try:
        from gseek.execution import Executor, ToolResolver
        from gseek.index import IndexStore
    except ImportError:
        return None
    env_key = tuple(
        sorted(
            (key, value)
            for key, value in os.environ.items()
            if key in {"PATH", "HOME", "XDG_CACHE_HOME"} or key.startswith("GSEEK_")
        )
    )
    global _session_cache
    with _SESSION_LOCK:
        if _session_cache is None or _session_cache.env_key != env_key:
            _session_cache = _GseekSession(ToolResolver(), Executor(), IndexStore(), env_key)
        return _session_cache


class _LegacyMode:
    value = "legacy"


class _LegacyEngine:
    """No-gseek fallback so legacy installs never import the optional package."""

    mode = _LegacyMode()

    def run(
        self,
        legacy: Callable[[], T],
        gseek: Callable[[], T],
        *,
        compare: Callable[[T, T], bool] | None = None,
    ) -> T:
        return legacy()


class _GseekRecord(Protocol):
    """The subset of a gseek AlignmentRecord used by gapit."""

    query_id: str
    target_id: str
    query_length: int | None
    target_length: int | None
    query_start: int | None
    query_end: int | None
    target_start: int | None
    target_end: int | None
    identity: float | None
    aligned_length: int | None
    matches: int | None
    gap_count: int | None
    gap_openings: int | None
    evalue: float | None
    mapq: int | None
    target_description: str | None
    is_primary: bool
    tags: Mapping[str, str]

    @property
    def is_reverse(self) -> bool: ...


class _GseekResult(Protocol):
    records: tuple[_GseekRecord, ...]


class _GseekModule(Protocol):
    def align(self, **kwargs: object) -> _GseekResult: ...

    def iter_align(self, **kwargs: object) -> Iterable[_GseekRecord]: ...

    def align_many(self, **kwargs: object) -> _GseekResult: ...


def engine() -> _LegacyEngine | Engine:
    """Return the engine switch, importing gseek only for gseek/shadow modes."""
    raw = os.environ.get(_ENGINE_ENV, "legacy").strip().lower() or "legacy"
    if raw not in _VALID_MODES:
        raise GapitError(
            f"unknown engine mode {raw!r}; expected one of {sorted(_VALID_MODES)}",
            code="INVALID_ENGINE",
            context={_ENGINE_ENV: raw},
        )
    if raw == "legacy":
        return _LegacyEngine()
    try:
        from gseek.compat import Engine
    except ImportError as exc:
        raise DependencyError(
            f"{_ENGINE_ENV}={raw} requires the optional gseek package (install gapit[gseek])",
            code="MISSING_DEPENDENCY",
            context={"package": "gseek"},
        ) from exc
    return Engine.from_env(_ENGINE_ENV, logger=_LOGGER)


def engine_mode() -> str:
    """Return the selected mode as a string (backward-compatible helper)."""
    return engine().mode.value


def _gseek() -> _GseekModule:
    try:
        module = importlib.import_module("gseek")
    except ImportError as exc:
        raise DependencyError(
            f"{_ENGINE_ENV}=gseek requires the optional gseek package "
            "(Python >= 3.12; install gapit[gseek])",
            code="MISSING_DEPENDENCY",
            context={"package": "gseek"},
        ) from exc
    return cast(_GseekModule, module)


def _record_int(value: int | None, field: str) -> int:
    if value is None:
        raise GapitError(
            f"gseek {field} field is missing from an alignment record",
            code="GSEEK_RECORD_INVALID",
            context={"field": field},
        )
    return value


def _blast_row(record: _GseekRecord, *, protein: bool = False) -> BlastRow:
    """Map one gseek BLAST record back to BLAST's 15-column row contract.

    gseek stores 0-based half-open, forward-strand coordinates; BLAST emits
    1-based inclusive coordinates and reports a minus subject as sstart > send.
    gapit's dedup/swap logic depends on exactly that legacy representation.
    """
    from gapit.blast import BlastRow

    qstart = _record_int(record.query_start, "query_start") + 1
    qend = _record_int(record.query_end, "query_end")
    tstart = _record_int(record.target_start, "target_start")
    tend = _record_int(record.target_end, "target_end")
    reverse = record.is_reverse
    sstart = tend if reverse else tstart + 1
    send = tstart + 1 if reverse else tend
    return BlastRow(
        qseqid=record.query_id,
        qstart=qstart,
        qend=qend,
        qlen=_record_int(record.query_length, "query_length"),
        sseqid=record.target_id,
        sstart=sstart,
        send=send,
        slen=_record_int(record.target_length, "target_length"),
        sstrand="N/A" if protein else ("minus" if reverse else "plus"),
        evalue=record.evalue if record.evalue is not None else 0.0,
        length=_record_int(record.aligned_length, "aligned_length"),
        pident=record.identity if record.identity is not None else 0.0,
        gaps=record.gap_count or 0,
        gapopen=record.gap_openings or 0,
        stitle=record.target_description or record.target_id,
    )


_INDEX_POLICY: dict[tuple[Path, int, int], str] = {}
_INDEX_POLICY_LOCK = threading.Lock()
_MINIMAP2_NAME_LIMIT = 255


def _minimap2_index_policy(database: Database) -> str:
    """Use a direct FASTA target when minimap2 names exceed its index limit.

    minimap2's ``.mmi`` index replaces over-long target names with numeric
    sequence indices; the legacy path passes the FASTA target directly and
    keeps the names. CARD has six such names (including MexB), so an index
    changes the reported gene even though the alignment is identical.
    """
    path = database.sequences_path
    try:
        stat = path.stat()
    except OSError:
        return "auto"
    cache_key = (path, stat.st_mtime_ns, stat.st_size)
    with _INDEX_POLICY_LOCK:
        cached = _INDEX_POLICY.get(cache_key)
        if cached is not None:
            return cached
        policy = "auto"
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    name = line[1:].split(maxsplit=1)[0]
                    if len(name) > _MINIMAP2_NAME_LIMIT:
                        policy = "none"
                        break
        for old_key in [key for key in _INDEX_POLICY if key[0] == path and key != cache_key]:
            del _INDEX_POLICY[old_key]
        _INDEX_POLICY[cache_key] = policy
        return policy


def _paf_record(record: _GseekRecord) -> PafRecord:
    """Map one gseek minimap2 record to gapit's PAF boundary model."""
    from gapit.paf import PafRecord

    tags = record.tags
    nm_raw = tags.get("NM")
    try:
        nm = int(nm_raw) if nm_raw is not None else None
    except ValueError as exc:
        raise GapitError(
            f"gseek minimap2 NM tag is not an integer: {nm_raw!r}",
            code="GSEEK_RECORD_INVALID",
            context={"tag": "NM"},
        ) from exc
    return PafRecord(
        qname=record.query_id,
        qlen=_record_int(record.query_length, "query_length"),
        qstart=_record_int(record.query_start, "query_start"),
        qend=_record_int(record.query_end, "query_end"),
        strand="-" if record.is_reverse else "+",
        tname=record.target_id,
        tlen=_record_int(record.target_length, "target_length"),
        tstart=_record_int(record.target_start, "target_start"),
        tend=_record_int(record.target_end, "target_end"),
        nmatch=_record_int(record.matches, "matches"),
        alen=_record_int(record.aligned_length, "aligned_length"),
        mapq=_record_int(record.mapq, "mapq"),
        is_primary=record.is_primary,
        nm=nm,
        cs=tags.get("cs"),
    )


def run_blast_gseek(
    fasta_text: str,
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> list[BlastRow]:
    """Run one normalized query through gseek's BLAST adapter."""
    max_targets = 10000
    gseek_params: dict[str, object] = {
        "task": "blastx-fast" if dbtype == "prot" else "blastn",
        "evalue": 1e-20,
        "max_target_seqs": max_targets,
        "culling_limit": 1,
    }
    if dbtype == "prot":
        gseek_params["seg"] = False
    else:
        gseek_params["dust"] = False
        gseek_params["min_identity"] = params.minid
    if debug:
        print(
            f"gapit: gseek blast: db={database.name} dbtype={dbtype}",
            file=sys.stderr,
        )
    with tempfile.TemporaryDirectory(prefix="gapit-gseek-") as tmp:
        query_path = Path(tmp) / "query.fasta"
        query_path.write_text(fasta_text, encoding="utf-8")
        gseek = _gseek()
        supports_refs = supports_api_refs(gseek)
        target_refs = cached_target_refs(database) if supports_refs else None
        index_policy = "auto"
        if _supports_external_blast_db() and blast_db_ready(database):
            # The database already has a BLAST index (db build/fetch); use it
            # instead of making gseek materialize and index a second copy.
            gseek_params["db"] = str(database.sequences_path.absolute())
            index_policy = "none"
        session = _session()
        session_kwargs: dict[str, object] = {}
        if session is not None:
            session_kwargs = {
                "resolver": session.resolver,
                "executor": session.executor,
                "store": session.store,
            }
        ref_kwargs: dict[str, object] = {}
        if target_refs is not None:
            ref_kwargs["target_refs"] = target_refs
        from gseek.errors import ToolNotFoundError

        try:
            result = gseek.align(
                backend="blast",
                query=query_path,
                target=database.sequences_path,
                query_kind="fasta",
                params=gseek_params,
                threads=params.threads,
                index=index_policy,
                record_manifest=False,
                **session_kwargs,
                **ref_kwargs,
            )
        except ToolNotFoundError as exc:
            binary = "blastx" if dbtype == "prot" else "blastn"
            raise DependencyError(
                f"required binary not found on PATH: {binary}",
                code="MISSING_DEPENDENCY",
                context={"binary": binary, "backend": "blast"},
            ) from exc
        except Exception as exc:
            raise GapitError(
                f"gseek BLAST failed: {exc}",
                code="GSEEK_FAILED",
                context={"backend": "blast"},
            ) from exc
    return [_blast_row(record, protein=dbtype == "prot") for record in result.records]


def run_blast_many_gseek(
    queries: Sequence[Path],
    database: Database,
    params: ScreeningParams,
    *,
    dbtype: Literal["nucl", "prot"],
    debug: bool = False,
) -> list[list[BlastRow]] | None:
    """Batch several FASTA queries through gseek.align_many.

    Returns one row list per input path in input order, or None when the
    installed gseek lacks the batch API (gapit then falls back to per-file
    screening). The records carry gseek's ``batch_index`` tag after the API
    restores the original query ids, which is what lets us split them back.
    """
    gseek = _gseek()
    if not supports_batch(gseek):
        return None
    gseek_params: dict[str, object] = {
        "task": "blastx-fast" if dbtype == "prot" else "blastn",
        "evalue": 1e-20,
        "max_target_seqs": 10000,
        "culling_limit": 1,
    }
    if dbtype == "prot":
        gseek_params["seg"] = False
    else:
        gseek_params["dust"] = False
        gseek_params["min_identity"] = params.minid
    target_refs = cached_target_refs(database) if supports_api_refs(gseek) else None
    index_policy = "auto"
    if _supports_external_blast_db() and blast_db_ready(database):
        gseek_params["db"] = str(database.sequences_path.absolute())
        index_policy = "none"
    session = _session()
    session_kwargs: dict[str, object] = {}
    if session is not None:
        session_kwargs = {
            "resolver": session.resolver,
            "executor": session.executor,
            "store": session.store,
        }
    ref_kwargs: dict[str, object] = {}
    if target_refs is not None:
        ref_kwargs["target_refs"] = target_refs
    batch_kwargs: dict[str, object] = {}
    if _supports_batch_option(gseek, "compute_inputs"):
        # gapit consumes records, not gseek's input provenance; skip the full
        # query-file hash when the newer batch API allows it.
        batch_kwargs["compute_inputs"] = False
    if _supports_batch_option(gseek, "sort_records"):
        # report_from_rows() performs gapit's own stable hit ordering.
        batch_kwargs["sort_records"] = False
    if debug:
        print(
            f"gapit: gseek batch blast: db={database.name} files={len(queries)} dbtype={dbtype}",
            file=sys.stderr,
        )
    from gseek.errors import ToolNotFoundError

    try:
        result = gseek.align_many(
            backend="blast",
            queries=queries,
            target=database.sequences_path,
            params=gseek_params,
            threads=params.threads,
            index=index_policy,
            annotate_batch=True,
            **session_kwargs,
            **ref_kwargs,
            **batch_kwargs,
        )
    except ToolNotFoundError as exc:
        binary = "blastx" if dbtype == "prot" else "blastn"
        raise DependencyError(
            f"required binary not found on PATH: {binary}",
            code="MISSING_DEPENDENCY",
            context={"binary": binary, "backend": "blast"},
        ) from exc
    except Exception as exc:
        raise GapitError(
            f"gseek BLAST batch failed: {exc}",
            code="GSEEK_FAILED",
            context={"backend": "blast"},
        ) from exc
    groups: list[list[BlastRow]] = [[] for _ in queries]
    for record in result.records:
        raw_index = record.tags.get("batch_index")
        try:
            index = int(raw_index) if raw_index is not None else -1
        except ValueError as exc:
            raise GapitError(
                "gseek batch record has a non-integer batch_index",
                code="GSEEK_RECORD_INVALID",
                context={"batch_index": str(raw_index)},
            ) from exc
        if index < 0 or index >= len(groups):
            raise GapitError(
                "gseek batch record has an out-of-range batch_index",
                code="GSEEK_RECORD_INVALID",
                context={"batch_index": str(raw_index)},
            )
        groups[index].append(_blast_row(record, protein=dbtype == "prot"))
    return groups


def iter_paf_gseek(
    lanes: Sequence[tuple[Path, Path | None]],
    database: Database,
    *,
    read_type: str,
    threads: int,
    debug: bool = False,
    nm_tags: bool = False,
) -> Iterable[PafRecord]:
    """Yield gapit PAF records from gseek minimap2 for every lane.

    ``secondary="all"`` means "do not pass --secondary=no": gapit's legacy
    invocation keeps minimap2's default secondary output and filters later.
    """
    gseek = _gseek()
    from gseek.errors import ToolNotFoundError

    supports_stream = callable(getattr(gseek, "iter_align", None))
    supports_refs = supports_api_refs(gseek)
    target_refs = cached_target_refs(database) if supports_refs else None
    ref_kwargs: dict[str, object] = {}
    if target_refs is not None:
        ref_kwargs["target_refs"] = target_refs
    session = _session()
    session_kwargs: dict[str, object] = {}
    if session is not None:
        session_kwargs = {
            "resolver": session.resolver,
            "executor": session.executor,
            "store": session.store,
        }
    for r1, r2 in lanes:
        query = (r1,) if r2 is None else (r1, r2)
        gseek_params: dict[str, object] = {
            "preset": read_type,
            "secondary": "all",
            "cs": nm_tags,
        }
        if read_type == "sr":
            gseek_params.update({"k": 11, "w": 6, "B": 2})
        if debug:
            print(
                f"gapit: gseek minimap2: preset={read_type} reads={r1}",
                file=sys.stderr,
            )
        try:
            if supports_stream:
                records = gseek.iter_align(
                    backend="minimap2",
                    query=query,
                    target=database.sequences_path,
                    paired=r2 is not None,
                    params=gseek_params,
                    threads=threads,
                    index=_minimap2_index_policy(database),
                    **session_kwargs,
                    **ref_kwargs,
                )
                for record in records:
                    yield _paf_record(record)
            else:
                result = gseek.align(
                    backend="minimap2",
                    query=query,
                    target=database.sequences_path,
                    paired=r2 is not None,
                    params=gseek_params,
                    threads=threads,
                    index=_minimap2_index_policy(database),
                    record_manifest=False,
                    **session_kwargs,
                    **ref_kwargs,
                )
                for record in result.records:
                    yield _paf_record(record)
        except ToolNotFoundError as exc:
            raise DependencyError(
                "required binary not found on PATH: minimap2",
                code="MISSING_DEPENDENCY",
                context={"binary": "minimap2", "backend": "minimap2"},
            ) from exc
        except Exception as exc:
            raise GapitError(
                f"gseek minimap2 failed: {exc}",
                code="GSEEK_FAILED",
                context={"backend": "minimap2"},
            ) from exc


__all__ = [
    "batch_supported",
    "engine",
    "engine_mode",
    "iter_paf_gseek",
    "run_blast_gseek",
    "run_blast_many_gseek",
]
