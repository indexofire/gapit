"""Panel gene-family and locus-containment sidecar (gspec.gene_families/1).

For panel accounting we do not claim phylogenetic orthology: the unit is a
homology family built from a one-time all-vs-all gene search.  A locus's
family set drives containment (``A`` is a subset of ``B``) and diagnostic
gene checks, while the raw family assignments also give the visualization a
stable colour/key across loci.

Pure functions only.  The expensive BLAST step lives in
``scripts/build_gene_families.py``; the cluster engine only reads the small
sidecar.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from gspec.errors import DatabaseError
from gspec.gbfeatures import FeaturesDocument

GENE_FAMILY_SCHEMA = "gspec.gene_families/1"
GENE_FAMILY_FILENAME = "gene_families.json"


class GeneFamily(BaseModel, frozen=True):
    """One homology family and its member gene ids."""

    family_id: str
    genes: tuple[str, ...]


class LocusFamilySet(BaseModel, frozen=True):
    """The family membership of one locus."""

    locus: str
    families: tuple[str, ...]


class ContainmentEdge(BaseModel, frozen=True):
    """``subset`` is a strict family-set subset of ``superset``; the
    diagnostic families are the ones that distinguish the two."""

    subset: str
    superset: str
    diagnostic: tuple[str, ...]


class GeneFamilyDocument(BaseModel, frozen=True):
    """gspec.gene_families/1 — families, locus sets, containment DAG."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: str = Field(default=GENE_FAMILY_SCHEMA, alias="schema")
    families: tuple[GeneFamily, ...]
    loci: tuple[LocusFamilySet, ...]
    containment: tuple[ContainmentEdge, ...]


class HomologyHit(BaseModel, frozen=True):
    """One filtered all-vs-all gene hit (BLAST tabular input)."""

    query: str
    subject: str
    identity_pct: float
    alignment_len: int
    query_len: int
    subject_len: int
    evalue: float


def filter_homology_hits(
    hits: Iterable[HomologyHit],
    *,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    max_evalue: float = 1e-5,
) -> tuple[HomologyHit, ...]:
    """Keep reciprocal-candidate hits at the family threshold.

    Coverage is measured against the shorter gene so a fragment can still
    join its family; identity and E-value stay strict.  Reciprocal-best-hit
    selection is deliberately left to the caller (``cluster_families`` uses
    the filtered graph's connected components, which is stable for
    multi-copy families).
    """
    kept: list[HomologyHit] = []
    for hit in hits:
        if hit.evalue > max_evalue or hit.identity_pct < min_identity:
            continue
        shorter = min(hit.query_len, hit.subject_len)
        if shorter <= 0 or 100.0 * hit.alignment_len / shorter < min_coverage:
            continue
        kept.append(hit)
    return tuple(kept)


class _UnionFind:
    """Deterministic union-find over gene ids."""

    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            item, self.parent[item] = self.parent[item], root
        return root

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            # Lexicographic root keeps the output independent of edge order.
            first, second = sorted((left_root, right_root))
            self.parent[second] = first


def _family_id(index: int) -> str:
    return f"GF{index:06d}"


def cluster_families(gene_ids: Iterable[str], hits: Iterable[HomologyHit]) -> dict[str, str]:
    """Gene id -> deterministic family id from connected homology components."""
    union = _UnionFind()
    for gene_id in sorted(set(gene_ids)):
        union.add(gene_id)
    for hit in hits:
        if hit.query != hit.subject:
            union.union(hit.query, hit.subject)
    components: dict[str, list[str]] = {}
    for gene_id in sorted(union.parent):
        components.setdefault(union.find(gene_id), []).append(gene_id)
    ordered = sorted(components.values(), key=lambda genes: genes[0])
    assignments: dict[str, str] = {}
    for index, genes in enumerate(ordered, start=1):
        family_id = _family_id(index)
        for gene_id in genes:
            assignments[gene_id] = family_id
    return assignments


def build_family_document(
    features: FeaturesDocument, assignments: Mapping[str, str]
) -> GeneFamilyDocument:
    """Assemble families, locus sets, and strict containment edges.

    A gene missing from ``assignments`` is a hard error: a partially applied
    family table would silently distort containment and extra-gene counts.
    """
    by_family: dict[str, list[str]] = {}
    for locus in features.loci:
        for gene in locus.genes:
            family = assignments.get(gene.gene_id)
            if family is None:
                raise DatabaseError(
                    f"gene family assignment missing for {gene.gene_id!r}",
                    code="GENE_FAMILIES_INCOMPLETE",
                    context={"gene": gene.gene_id, "locus": locus.id},
                )
            by_family.setdefault(family, []).append(gene.gene_id)
    families = tuple(
        GeneFamily(family_id=family, genes=tuple(sorted(genes)))
        for family, genes in sorted(by_family.items())
    )
    loci = tuple(
        LocusFamilySet(
            locus=locus.id,
            families=tuple(sorted({assignments[gene.gene_id] for gene in locus.genes})),
        )
        for locus in features.loci
    )
    edges: list[ContainmentEdge] = []
    for subset in loci:
        subset_set = set(subset.families)
        for superset in loci:
            if subset.locus == superset.locus:
                continue
            superset_set = set(superset.families)
            if subset_set < superset_set:
                edges.append(
                    ContainmentEdge(
                        subset=subset.locus,
                        superset=superset.locus,
                        diagnostic=tuple(sorted(superset_set - subset_set)),
                    )
                )
    return GeneFamilyDocument(families=families, loci=loci, containment=tuple(edges))


def write_gene_families(document: GeneFamilyDocument, path: Path) -> None:
    """Deterministic JSON bytes: indent 2, by alias, trailing newline."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document.model_dump_json(indent=2, by_alias=True) + "\n", encoding="utf-8")


def read_gene_families(path: Path) -> GeneFamilyDocument:
    """Read and validate a gene-family sidecar (typed DB error on failure)."""
    try:
        return GeneFamilyDocument.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"cannot read gene families {path}: {exc}",
            code="GENE_FAMILIES_MALFORMED",
            context={"file": str(path)},
        ) from exc
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed gene families {path}: {exc}",
            code="GENE_FAMILIES_MALFORMED",
            context={"file": str(path)},
        ) from exc


def load_gene_families(db_dir: Path) -> GeneFamilyDocument | None:
    """The optional sidecar next to a cluster database, or None."""
    path = db_dir / GENE_FAMILY_FILENAME
    return read_gene_families(path) if path.is_file() else None


def locus_family_sets(document: GeneFamilyDocument) -> dict[str, frozenset[str]]:
    """Locus -> family set, convenient for containment checks."""
    return {entry.locus: frozenset(entry.families) for entry in document.loci}


def gene_family_map(document: GeneFamilyDocument) -> dict[str, str]:
    """Gene id -> family id (the inverse of the family member lists)."""
    return {gene_id: family.family_id for family in document.families for gene_id in family.genes}


def containment_for(document: GeneFamilyDocument, locus: str) -> tuple[ContainmentEdge, ...]:
    """Edges where ``locus`` is either side (subset or superset)."""
    return tuple(
        edge for edge in document.containment if edge.subset == locus or edge.superset == locus
    )


__all__ = [
    "GENE_FAMILY_FILENAME",
    "GENE_FAMILY_SCHEMA",
    "ContainmentEdge",
    "GeneFamily",
    "GeneFamilyDocument",
    "HomologyHit",
    "LocusFamilySet",
    "build_family_document",
    "cluster_families",
    "containment_for",
    "filter_homology_hits",
    "gene_family_map",
    "load_gene_families",
    "locus_family_sets",
    "read_gene_families",
    "write_gene_families",
]
