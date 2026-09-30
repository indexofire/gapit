"""The cluster-database provider contract and its fetch pipeline (stage 3).

The kaptive shape: GBK/GFF sources that download at fetch time and build
through the cluster pipeline (locus FASTA + features.json, ``kind:
cluster``) — no transform, no records.jsonl. Split from providers/common.py
at the 250-LOC ceiling; the gene-pipeline Provider and the shared download
seam (``download_file``/``url_basename``) stay there.
"""

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from gapit.clusterbuild import build_cluster_database, parse_cluster_input
from gapit.errors import DatabaseError
from gapit.proctools import note
from gapit.providers.common import Dbtype, download_file, url_basename
from gapit.records import Manifest


@dataclass(frozen=True, slots=True)
class ClusterProvider:
    """A GBK/GFF-backed cluster-database provider. ``license`` and ``note``
    land in the build manifest for provenance; ``snapshot`` exists only to
    mirror the Provider shape and is ALWAYS None — upstream licenses (e.g.
    Kaptive's GPL-3.0) must never ship inside the wheel. Cluster loci are
    always nucleotide, so dbtype is fixed."""

    name: str
    description: str
    source_urls: tuple[str, ...]
    license: str
    note: str
    snapshot: str | None = None
    dbtype: Dbtype = "nucl"
    kind: Literal["cluster"] = "cluster"


def fetch_cluster_provider(
    provider: ClusterProvider,
    db_dir: Path,
    *,
    fetched_at: str,
    force: bool = False,
    quiet: bool = True,
    debug: bool = False,
) -> Manifest:
    """Download the provider's GBK/GFF sources and build a kind-cluster
    database through the cluster pipeline (no typing model: phenotype calls
    stay null; users get kaptive-style locus calls). Network-only by
    contract — nothing is ever bundled. Same guards as fetch_provider
    (DB_ALREADY_EXISTS, DOWNLOAD_FAILED) plus the cluster parser's typed
    input errors."""
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
            download_file(url, workdir / url_basename(url))
        note(quiet, f"downloaded {len(provider.source_urls)} source file(s)")
        loci = parse_cluster_input(workdir / url_basename(provider.source_urls[0]))
    note(quiet, f"parsed {len(loci)} loci from {provider.name}")
    return build_cluster_database(
        db_dir,
        name=provider.name,
        loci=loci,
        typing=None,
        source_urls=provider.source_urls,
        fetched_at=fetched_at,
        content_license=provider.license,
        content_note=provider.note,
        quiet=quiet,
        debug=debug,
    )
