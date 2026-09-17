"""Wave B0 shared provider infrastructure: the generic get_db pipeline.

Every concrete provider (Wave B1..B12) is a :class:`Provider` value — a name,
a description, a tuple of source URLs, a dbtype, and a ``transform`` that
turns downloaded files into typed ``Record``s. :func:`fetch_provider` runs
the generic abricate-get_db flow around it: download each URL into a scratch
workdir under ``db_dir``, transform, normalize sequences and function classes
(upstream load_fasta semantics), dedupe by exact normalized sequence
(first wins), sort by gene, persist ``records.jsonl``, then delegate to
:func:`gapit.dbbuild.build_database` for the ``sequences`` FASTA, the BLAST
and minimap2 indexes, and the manifest (written last, certifying the build).

Upstream's ``is_full_gene`` is deliberately NOT ported: its map result is
discarded — a no-op.
"""

import os
import re
import sys
import urllib.request
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from urllib.parse import urlparse

from gapit import __version__
from gapit.dbbuild import build_database
from gapit.errors import DatabaseError
from gapit.records import Manifest, Record, write_records

_NUCL_JUNK = re.compile(r"[^AGCT]")
_PROT_JUNK = re.compile(r"[^A-Z]")
_WHITESPACE = re.compile(r"\s+")
_CHUNK_BYTES = 1 << 20

# mgc.ac.cn (vfdb) 403s "Python-urllib" UAs specifically while serving
# browser-ish clients (Wave E diagnosis): the Mozilla compatibility token
# passes those naive UA filters and the gapit/<version> suffix still
# identifies the tool to servers that log it.
_USER_AGENT = f"Mozilla/5.0 (compatible; gapit/{__version__})"

Dbtype = Literal["nucl", "prot"]


@dataclass(frozen=True, slots=True)
class Provider:
    """One database provider: everything fetch_provider needs to know.

    ``transform`` receives the download workdir (each source URL saved under
    its basename) and yields Records; setting ``db`` to the provider name is
    the provider module's job. Sequences and function classes may arrive
    raw — fetch_provider normalizes both.
    """

    name: str
    description: str
    source_urls: tuple[str, ...]
    dbtype: Dbtype
    transform: Callable[[Path], Iterable[Record]]


def _note(quiet: bool, message: str) -> None:
    """Per-step progress on stderr when quiet is disabled (stdout stays pure).

    Mirrors gapit.dbbuild._note — reportPrivateUsage blocks importing it.
    """
    if not quiet:
        print(f"gapit: {message}", file=sys.stderr)


def _basename(url: str) -> str:
    """Final path segment of a URL, query string excluded ('' for a bare host)."""
    return urlparse(url).path.rstrip("/").rsplit("/", 1)[-1]


def _download(url: str, dest: Path) -> None:
    """Fetch ``url`` into ``dest`` atomically: urlopen a Request carrying
    the module User-Agent, stream the body into a hidden .part file beside
    ``dest`` in chunks, then os.replace. file:// URLs work (tests depend on
    it). Any URLError/OSError/ValueError becomes DatabaseError
    ``DOWNLOAD_FAILED`` with the URL in context. No shell, ever.
    """
    temp = dest.with_name(f".{dest.name}.part")
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(request) as response, temp.open("wb") as out:
            while chunk := response.read(_CHUNK_BYTES):
                out.write(chunk)
        os.replace(temp, dest)
    except (OSError, ValueError) as exc:
        raise DatabaseError(
            f"failed to download {url}: {exc}",
            code="DOWNLOAD_FAILED",
            context={"url": url},
        ) from exc
    finally:
        temp.unlink(missing_ok=True)


def _normalize_sequence(sequence: str, dbtype: Dbtype) -> str:
    """Uppercase; then nucl: non-AGCT -> 'N', prot: non-A-Z -> 'X'
    (upstream load_fasta semantics)."""
    upper = sequence.upper()
    match dbtype:
        case "nucl":
            return _NUCL_JUNK.sub("N", upper)
        case "prot":
            return _PROT_JUNK.sub("X", upper)


def _normalize_function(classes: Sequence[str]) -> tuple[str, ...]:
    """Sorted classes with each whitespace run collapsed to one '_'.

    Upstream sorts, joins with ';', then globally substitutes s/\\s+/_/g —
    ';' carries no whitespace, so per-class substitution is equivalent.
    """
    return tuple(_WHITESPACE.sub("_", drug_class) for drug_class in sorted(classes))


def _normalize_record(record: Record, dbtype: Dbtype) -> Record:
    """Apply the two per-record normalizations; db, gene, and the rest stay."""
    return record.model_copy(
        update={
            "sequence": _normalize_sequence(record.sequence, dbtype),
            "function": _normalize_function(record.function),
        }
    )


def _dedupe(records: Sequence[Record]) -> tuple[tuple[Record, ...], int]:
    """Drop later records whose (already normalized) SEQUENCE was seen; first
    wins. Returns the kept records and the dropped count (for the caller's
    stderr note). Duplicate gene NAMES are allowed and kept — upstream too.
    """
    seen: set[str] = set()
    kept: list[Record] = []
    dropped = 0
    for record in records:
        if record.sequence in seen:
            dropped += 1
        else:
            seen.add(record.sequence)
            kept.append(record)
    return tuple(kept), dropped


def fetch_provider(
    provider: Provider,
    db_dir: Path,
    *,
    fetched_at: str,
    force: bool = False,
    quiet: bool = True,
) -> Manifest:
    """Run the generic provider pipeline into ``db_dir`` (created if needed).

    Raises DatabaseError (exit 4): ``DB_ALREADY_EXISTS`` when a manifest is
    present and force is off (upstream: "Won't overwrite existing (use
    --force)"), ``DOWNLOAD_FAILED`` for a failed source download, and
    ``PROVIDER_EMPTY`` when the transform+dedupe leaves zero records.
    """
    if (db_dir / "gapit-manifest.json").is_file() and not force:
        raise DatabaseError(
            f"won't overwrite existing database {provider.name} (use --force)",
            code="DB_ALREADY_EXISTS",
            context={"db": provider.name},
        )
    db_dir.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=db_dir, prefix=".download.") as workdir_name:
        workdir = Path(workdir_name)
        for url in provider.source_urls:
            _download(url, workdir / _basename(url))
        _note(quiet, f"downloaded {len(provider.source_urls)} source file(s)")
        records = tuple(
            _normalize_record(record, provider.dbtype) for record in provider.transform(workdir)
        )
    kept, dropped = _dedupe(records)
    _note(quiet, f"read {len(records)} records from {provider.name}")
    _note(quiet, f"dropped {dropped} duplicate sequence(s), kept {len(kept)}")
    if not kept:
        raise DatabaseError(
            f"provider {provider.name} yielded no records",
            code="PROVIDER_EMPTY",
            context={"db": provider.name},
        )
    write_records(sorted(kept, key=lambda record: record.gene), db_dir / "records.jsonl")
    return build_database(
        db_dir,
        name=provider.name,
        dbtype=provider.dbtype,
        source_urls=provider.source_urls,
        fetched_at=fetched_at,
        quiet=quiet,
    )
