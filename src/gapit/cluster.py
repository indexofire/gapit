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
from dataclasses import dataclass
from pathlib import Path

from pydantic import BaseModel, Field, ValidationError

from gapit.cluster_math import (
    CsRun,
    Verdict,
    gene_verdict,
    intersect_coverage,
    locus_sort_key,
    parse_cs,
)
from gapit.db import Database
from gapit.errors import DatabaseError, GapitError
from gapit.gbfeatures import FeaturesDocument, LocusFeatures
from gapit.minimap2_run import run_minimap2
from gapit.paf import PafRecord
from gapit.typing_models import PhenotypeDetail, TypingDocument, read_typing_document


class ClusterParams(BaseModel, frozen=True):
    """Parameters for one cluster-screening run; ``preset`` is always asm20
    (the vpautils engine's choice) and is recorded verbatim in the params
    block of gapit.cluster/1."""

    db: str
    preset: str = "asm20"
    min_gene_cov: float = Field(default=90.0, ge=0.0, le=100.0)
    min_gene_id: float = Field(default=90.0, ge=0.0, le=100.0)
    min_cluster_cov: float = Field(default=96.0, ge=0.0, le=100.0)
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
    """Result of cluster-screening one sample file."""

    file: str
    best: BestCall | None
    loci: tuple[LocusCall, ...]


def load_features(database: Database) -> FeaturesDocument:
    """Read the cluster db's gapit.features/1 table; a cluster-kind dir
    without a readable one is a database error (rebuild it)."""
    path = database.path / "features.json"
    if not path.is_file():
        raise DatabaseError(
            f"cluster database {database.name} lacks features.json (rebuild it: gapit db build)",
            code="FEATURES_MISSING",
            context={"db": database.name},
        )
    try:
        return FeaturesDocument.model_validate_json(path.read_text(encoding="utf-8"))
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


def _runs_by_locus(rows: list[PafRecord]) -> dict[str, tuple[list[tuple[CsRun, ...]], int]]:
    """Primary rows grouped by target locus: each record's parsed cs runs
    plus the first-seen tlen (the reads-mode denominator contract)."""
    grouped: dict[str, tuple[list[tuple[CsRun, ...]], int]] = {}
    for row in rows:
        if not row.is_primary:
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
) -> ClusterReport:
    """Screen one assembly file: one minimap2 asm20 --cs invocation (query =
    the file's contigs, target = the locus FASTA), primary rows only, per-
    gene/locus union math, ranked loci, best call above min_cluster_cov."""
    rows = run_minimap2(
        [(query, None)],
        database,
        read_type="asm20",
        threads=params.threads,
        debug=debug,
        nm_tags=True,
    )
    grouped = _runs_by_locus(rows)
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
    return ClusterReport(file=str(query), best=best, loci=loci)
