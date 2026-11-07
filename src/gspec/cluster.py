"""The cluster-screening engine core (stage 2 of the gene-cluster feature).

One batched minimap2 ``asm20`` invocation per input file (query = sample
contigs, target = the cluster db's locus FASTA — the reads-mode target
orientation, via minimap2_run), then per-gene and per-locus coverage /
identity from the cs walk (cluster_math). Loci are ranked (coverage desc,
identity desc, covered bp desc, id asc) and the rank-1 locus becomes the
file's ``best`` call when its coverage reaches ``min_cluster_cov``
(kaptive-style confidence floor; below it no call is made). Records from
multiple contigs union naturally, so loci fragmented across contigs still
screen as present. Assembly FASTA only in v1 — reads mode rejects cluster
databases; the blastn gene path is untouched (parity is sacred). The
multi-file orchestration and rendering fan-out live in screening.py beside
the gene path's.
"""

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from pydantic import BaseModel, Field, ValidationError

from gspec.cluster_math import (
    CsRun,
    Verdict,
    gene_verdict,
    intersect_coverage,
    locus_sort_key,
    parse_cs,
)
from gspec.cluster_profile import CoverageMode
from gspec.db import Database
from gspec.errors import DatabaseError, GapitError
from gspec.gbfeatures import FeaturesDocument, LocusFeatures
from gspec.gene_families import load_gene_families
from gspec.locus_map_models import LocusMap
from gspec.minimap2_run import iter_minimap2
from gspec.paf import PafRecord
from gspec.typing_models import TypingDocument, read_typing_document
from gspec.typing_results import PhenotypeDetail, SchemeCall


class ClusterParams(BaseModel, frozen=True):
    """Parameters for one cluster-screening run; ``preset`` is always asm20
    (the vpautils engine's choice) and is recorded verbatim in the params
    block of gapit.cluster/1."""

    db: str
    preset: str = "asm20"
    coverage_mode: CoverageMode = "primary"
    fallback_evidence: bool = False
    min_gene_cov: float = Field(default=90.0, ge=0.0, le=100.0)
    min_gene_id: float = Field(default=90.0, ge=0.0, le=100.0)
    min_cluster_cov: float = Field(default=96.0, ge=0.0, le=100.0)
    locus_typeable_coverage: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_typeable_identity: float | None = Field(default=None, ge=0.0, le=100.0)
    locus_quality_min_coverage: float = Field(default=80.0, ge=0.0, le=100.0)
    locus_quality_min_identity: float = Field(default=90.0, ge=0.0, le=100.0)
    locus_quality_max_missing_genes: int = Field(default=3, ge=0)
    locus_quality_max_extra_panel_genes: int = Field(default=2, ge=0)
    threads: int = Field(default=1, ge=1)


class GeneCall(BaseModel, frozen=True):
    """One gene's screen against a locus (features coordinates, 1-based)."""

    gene_id: str
    start: int
    end: int
    strand: str
    coverage_pct: float
    identity_pct: float
    verdict: Verdict
    # Table-derived typing evidence: the alignment start as a percentage of
    # the reference marker (0.0 = the 5' end is covered).  Cluster screens
    # and reads tables leave it None; only `gspec typing` populates it.
    alignment_start_pct: float | None = None


class LocusCall(BaseModel, frozen=True):
    """One locus's screen: union coverage/matches-weighted identity over all
    records' target blocks, its rank, per-gene calls, and the gene ids that
    are not fully present."""

    locus: str
    label: str
    type: str
    coverage_pct: float
    identity_pct: float
    rank: int
    genes: tuple[GeneCall, ...]
    missing: tuple[str, ...]
    kmer_hits: dict[str, float] = Field(default_factory=dict)


class BestCall(BaseModel, frozen=True):
    """The best-locus summary (rank 1, coverage >= min_cluster_cov; None
    when no locus clears the floor). ``phenotype`` and ``phenotype_detail``
    stay null until the stage-3 typing evaluator fills them (databases
    without a typing.json never annotate them)."""

    locus: str
    label: str
    type: str
    coverage_pct: float
    identity_pct: float
    genes_present: int
    genes_partial: int
    genes_absent: int
    phenotype: str | None = None
    phenotype_detail: PhenotypeDetail | None = None


class ClusterReport(BaseModel, frozen=True):
    """Result of cluster-screening one sample file. ``phenotypes`` carries
    the gapit.cluster/2 per-scheme calls (set only by the multi-scheme
    path; the /1 path annotates ``best`` instead and never sets it)."""

    file: str
    best: BestCall | None
    loci: tuple[LocusCall, ...]
    phenotypes: dict[str, SchemeCall] = Field(default_factory=dict)
    locus_maps: dict[str, LocusMap] = Field(default_factory=dict)


def load_features(database: Database) -> FeaturesDocument:
    """Read the cluster db's gapit.features/1 table; a cluster-kind dir
    without a readable one is a database error (rebuild it)."""
    path = database.path / "features.json"
    if not path.is_file():
        raise DatabaseError(
            f"cluster database {database.name} lacks features.json (rebuild it: gspec db build)",
            code="FEATURES_MISSING",
            context={"db": database.name},
        )
    try:
        return FeaturesDocument.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"cluster database {database.name} has an unreadable features.json: {exc}",
            code="FEATURES_MALFORMED",
            context={"db": database.name},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"cluster database {database.name} has a malformed features.json: {exc}",
            code="FEATURES_MALFORMED",
            context={"db": database.name},
        ) from exc


def load_typing(database: Database) -> TypingDocument | None:
    """The cluster db's typing.json as a validated document, or None when
    the database carries none (untyped screening, stage-2 output)."""
    path = database.path / "typing.json"
    return read_typing_document(path) if path.is_file() else None


def runs_by_locus(
    rows: Iterable[PafRecord], *, coverage_mode: CoverageMode = "primary"
) -> dict[str, tuple[list[tuple[CsRun, ...]], int]]:
    """Rows grouped by target locus: each record's parsed cs runs plus the
    first-seen tlen (the reads-mode denominator contract). Primary-only is
    the historical default; ``primary_secondary`` is the opt-in Kaptive-style
    union used by databases whose profile asks for it."""
    grouped: dict[str, tuple[list[tuple[CsRun, ...]], int]] = {}
    for row in rows:
        if coverage_mode == "primary" and not row.is_primary:
            continue
        if row.cs is None:
            raise GapitError(
                f"minimap2 row without a cs tag (cannot score cluster {row.tname})",
                code="CS_MISSING",
                context={"locus": row.tname},
            )
        grouped.setdefault(row.tname, ([], row.tlen))[0].append(parse_cs(row.cs, row.tstart))
    return grouped


@dataclass(frozen=True, slots=True)
class _LocusStats:
    """Rank-intermediate: the (unranked) locus call plus its covered bases."""

    call: LocusCall
    covered_bp: int


def _kmer_hit_fractions(
    query: Path, database: Database, covered_loci: set[str]
) -> dict[str, tuple[str, float]]:
    """Per-set hit fractions for a database carrying a ``gapit.kmerdb/1``
    sidecar: count the query's canonical k-mers once, intersect every
    declared set, and return ``{set_name: (of_locus, fraction)}`` filtered
    to sets anchored on covered loci (absent loci keep an empty map).
    Databases without the sidecar return ``{}`` — the alignment-only path
    is unchanged."""
    from gspec.kmerdb import count_fasta_kmers, load_kmer_set, read_kmerdb, set_hit_fraction

    document = read_kmerdb(database.path)
    if document is None:
        return {}
    relevant = {name: spec for name, spec in document.sets.items() if spec.of in covered_loci}
    if not relevant:
        return {}
    codes = count_fasta_kmers(query, document.k)
    hits: dict[str, tuple[str, float]] = {}
    for name, spec in relevant.items():
        reference = load_kmer_set(database.path, name)
        hits[name] = (spec.of, set_hit_fraction(codes, reference))
    return hits


def _locus_call(
    locus: LocusFeatures, grouped: tuple[list[tuple[CsRun, ...]], int], params: ClusterParams
) -> _LocusStats:
    """Union coverage and matches-weighted identity over every record's
    target blocks (locus denominator = first-seen tlen), per-gene calls from
    the same run set, and the not-present gene ids."""
    runs, locus_len = grouped
    covered, matched = intersect_coverage(runs, 0, locus_len)
    coverage_pct = 100.0 * covered / locus_len
    identity_pct = 100.0 * matched / covered if covered else 0.0
    genes: list[GeneCall] = []
    for gene in locus.genes:
        g_covered, g_matched = intersect_coverage(runs, gene.start - 1, gene.end)
        gene_cov = 100.0 * g_covered / (gene.end - gene.start + 1)
        gene_id_pct = 100.0 * g_matched / g_covered if g_covered else 0.0
        genes.append(
            GeneCall(
                gene_id=gene.gene_id,
                start=gene.start,
                end=gene.end,
                strand=gene.strand,
                coverage_pct=gene_cov,
                identity_pct=gene_id_pct,
                verdict=gene_verdict(
                    gene_cov,
                    gene_id_pct,
                    min_gene_cov=params.min_gene_cov,
                    min_gene_id=params.min_gene_id,
                ),
            )
        )
    return _LocusStats(
        call=LocusCall(
            locus=locus.id,
            label=locus.label,
            type=locus.type,
            coverage_pct=coverage_pct,
            identity_pct=identity_pct,
            rank=0,
            genes=tuple(genes),
            missing=tuple(gene.gene_id for gene in genes if gene.verdict != "present"),
        ),
        covered_bp=covered,
    )


def screen_cluster_file(
    query: Path,
    database: Database,
    features: FeaturesDocument,
    params: ClusterParams,
    *,
    debug: bool = False,
    detail: bool = False,
) -> ClusterReport:
    """Screen one assembly file: one minimap2 asm20 --cs invocation (query =
    the file's contigs, target = the locus FASTA), per-gene/locus union math,
    ranked loci, best call above min_cluster_cov. Rows are primary-only unless
    the database profile opts into ``primary_secondary`` coverage. With
    ``detail``, attach the rich best-locus map (alignment pieces, indels,
    quality) without changing compact renderers. When the database carries a
    ``gapit.kmerdb/1`` sidecar the query's canonical k-mers are counted once
    and every declared set's hit fraction is attached to its ``of`` locus —
    the ``key_kmers`` typing component consumes these (SNP-level calls)."""
    row_stream = iter_minimap2(
        [(query, None)],
        database,
        read_type="asm20",
        threads=params.threads,
        debug=debug,
        nm_tags=True,
    )
    if detail:
        rows = list(row_stream)
        grouped = runs_by_locus(rows, coverage_mode=params.coverage_mode)
    else:
        rows = []
        grouped = runs_by_locus(row_stream, coverage_mode=params.coverage_mode)
    stats = [
        _locus_call(locus, grouped[locus.id], params)
        for locus in features.loci
        if locus.id in grouped
    ]
    stats.sort(
        key=lambda entry: locus_sort_key(
            entry.call.locus, entry.call.coverage_pct, entry.call.identity_pct, entry.covered_bp
        )
    )
    kmer_hits = _kmer_hit_fractions(query, database, {entry.call.locus for entry in stats})
    if kmer_hits:
        stats = [
            _LocusStats(
                call=entry.call.model_copy(
                    update={
                        "kmer_hits": {
                            name: fraction[1]
                            for name, fraction in kmer_hits.items()
                            if fraction[0] == entry.call.locus
                        }
                    },
                ),
                covered_bp=entry.covered_bp,
            )
            for entry in stats
        ]
    loci = tuple(
        entry.call.model_copy(update={"rank": rank}) for rank, entry in enumerate(stats, start=1)
    )
    best: BestCall | None = None
    if loci and loci[0].coverage_pct >= params.min_cluster_cov:
        top = loci[0]
        verdicts = Counter(gene.verdict for gene in top.genes)
        best = BestCall(
            locus=top.locus,
            label=top.label,
            type=top.type,
            coverage_pct=top.coverage_pct,
            identity_pct=top.identity_pct,
            genes_present=verdicts["present"],
            genes_partial=verdicts["partial"],
            genes_absent=verdicts["absent"],
        )
    locus_maps: dict[str, LocusMap] = {}
    if detail and loci:
        from gspec.cluster_alignment import build_locus_alignment, gene_query_spans
        from gspec.locus_map import build_locus_map

        family_document = load_gene_families(database.path)
        features_by_id = {locus.id: locus for locus in features.loci}
        # Extra panel genes must come from primary alignments: a secondary
        # cross-locus hit is homology, not a second locus in the sample.
        gene_spans = {
            locus_id: gene_query_spans(features_by_id[locus_id], rows, coverage_mode="primary")
            for locus_id in grouped
            if locus_id in features_by_id
        }
        locus_maps = {
            call.locus: build_locus_map(
                call,
                build_locus_alignment(
                    features_by_id[call.locus], rows, coverage_mode=params.coverage_mode
                ),
                features,
                loci,
                family_document,
                gene_spans,
                min_coverage=params.locus_quality_min_coverage,
                min_identity=params.locus_quality_min_identity,
                max_missing_genes=params.locus_quality_max_missing_genes,
                max_extra_panel_genes=params.locus_quality_max_extra_panel_genes,
            )
            for call in loci
        }
    return ClusterReport(file=str(query), best=best, loci=loci, locus_maps=locus_maps)
