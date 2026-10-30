"""The cluster-database provider contract and its fetch pipeline (stage 3).

The kaptive shape: GBK/GFF sources that download at fetch time and build
through the cluster pipeline (locus FASTA + features.json, ``kind:
cluster``) — no transform, no records.jsonl. Split from providers/common.py
at the 250-LOC ceiling; the gene-pipeline Provider and the shared download
seam (``download_file``/``url_basename``) stay there.

Providers may declare MULTIPLE source URLs (every downloaded file is
parsed and the loci concatenated in ``source_urls`` order — the merge a
multi-reference provider needs) and may carry a committed
gapit.typing asset (``typing``), validated exactly like
``gapit db build --typing`` — a multi-scheme document resolves gene ids
first, then ``validate_references`` gates drift against the downloaded
content BEFORE any artifact is written.
"""

from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal

from gapit.clusterbuild import build_cluster_database, parse_cluster_input, resolve_gene_ids
from gapit.db_info_models import DbGuidance
from gapit.db_txn import staged_database
from gapit.errors import DatabaseError
from gapit.gbfeatures import FeaturesDocument, LocusSequence
from gapit.proctools import note
from gapit.providers.common import (
    Dbtype,
    download_file,
    insecure_source_urls,
    security_note,
    url_basename,
)
from gapit.records import Manifest
from gapit.typing_models import (
    read_typing_document,
    typing_schema_of,
    uses_cluster_v2,
    validate_references,
)


@dataclass(frozen=True, slots=True)
class ClusterProvider:
    """A GBK/GFF-backed cluster-database provider. ``vendor`` names the
    upstream maintainer organisation (Kaptive (klebgenomics), ...); ``name``
    is the database name users pass to ``--db``. ``license`` and ``note``
    land in the build manifest for provenance. Upstream licenses (e.g.
    Kaptive's GPL-3.0) must never ship inside the wheel — every cluster
    provider downloads on fetch. Cluster loci are always nucleotide, so
    dbtype is fixed. ``typing`` optionally points at a committed
    gapit.typing asset riding the wheel (gapit-authored, MIT) that
    installs with the database; None keeps the untyped kaptive behavior."""

    name: str
    description: str
    vendor: str
    source_urls: tuple[str, ...]
    license: str
    note: str
    # Committed extra GBK/GFF loci appended after the downloaded sources
    # (a provider-specific control locus, for example). These are wheel
    # assets, not downloads; their content is part of the provider build.
    extra_sources: tuple[Path, ...] = ()
    dbtype: Dbtype = "nucl"
    kind: Literal["cluster"] = "cluster"
    typing: Path | None = None
    # Literature references + interpretation caveats for `gapit db info`
    # (usage is synthesized from the kind).
    guidance: DbGuidance = field(default_factory=DbGuidance)


def fetch_cluster_provider(
    provider: ClusterProvider,
    db_dir: Path,
    *,
    fetched_at: str,
    force: bool = False,
    quiet: bool = True,
    debug: bool = False,
) -> Manifest:
    """Download the provider's GBK/GFF sources (every URL, loci concatenated
    in declaration order) and build a kind-cluster database through the
    cluster pipeline. A declared ``typing`` asset validates exactly like
    ``gapit db build --typing`` (multi-scheme documents resolve gene ids
    first) BEFORE any artifact is written, so a spec drifting from upstream
    content fails TYPING_UNKNOWN_GENE with no half-built database. Otherwise
    the kaptive behavior is untouched: no typing model, phenotype calls stay
    null. Same guards as fetch_provider (DB_ALREADY_EXISTS, DOWNLOAD_FAILED)
    plus the cluster parser's typed input errors."""
    if (db_dir / "gapit-manifest.json").is_file() and not force:
        raise DatabaseError(
            f"won't overwrite existing database {provider.name} (use --force)",
            code="DB_ALREADY_EXISTS",
            context={"db": provider.name},
        )
    insecure = insecure_source_urls(provider.source_urls)
    if insecure:
        security_note(f"plain HTTP source(s), no transport integrity: {', '.join(insecure)}")
    with staged_database(db_dir) as build_dir:
        if (db_dir / "gapit-manifest.json").is_file() and not force:
            raise DatabaseError(
                f"won't overwrite existing database {provider.name} (use --force)",
                code="DB_ALREADY_EXISTS",
                context={"db": provider.name},
            )
        with TemporaryDirectory(dir=build_dir, prefix=".download.") as workdir_name:
            workdir = Path(workdir_name)
            for url in provider.source_urls:
                download_file(url, workdir / url_basename(url))
            note(quiet, f"downloaded {len(provider.source_urls)} source file(s)")
            loci = _merged_loci(provider, workdir)
        note(quiet, f"parsed {len(loci)} loci from {provider.name}")
        typing_schema = ""
        if provider.typing is not None:
            document = read_typing_document(provider.typing)
            if uses_cluster_v2(document):
                loci = resolve_gene_ids(loci)
            validate_references(
                document, FeaturesDocument(loci=tuple(item.features for item in loci))
            )
            typing_schema = typing_schema_of(provider.typing)
        return build_cluster_database(
            build_dir,
            name=provider.name,
            loci=loci,
            typing=provider.typing,
            source_urls=provider.source_urls,
            fetched_at=fetched_at,
            typing_schema=typing_schema,
            content_license=provider.license,
            content_note=provider.note,
            quiet=quiet,
            debug=debug,
        )


def _merged_loci(provider: ClusterProvider, workdir: Path) -> list[LocusSequence]:
    """Parse every downloaded source into one locus list (declaration
    order); a single-URL provider parses its one file, byte-identical to the
    pre-multi-URL pipeline."""
    loci: list[LocusSequence] = []
    for url in provider.source_urls:
        loci.extend(parse_cluster_input(workdir / url_basename(url)))
    for extra in provider.extra_sources:
        loci.extend(parse_cluster_input(extra))
    return loci
