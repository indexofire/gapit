"""gapit.kmerdb/1 — the k-mer substrate for cluster databases.

A cluster database may carry an optional ``kmer.json`` sidecar declaring
k-mer sets over its loci, materialized at build time as sorted arrays of
canonical k-mer codes (little-endian uint64) under ``kmers/``:

    {"schema": "gapit.kmerdb/1",
     "k": 51,
     "sets": {"06A.diff": {"of": "06A", "against": ["06B", "06F", "06G"], "n": 1182}}}

``of`` names the anchor locus; ``against`` the sibling loci whose union is
subtracted. Differential sets (SNP-spanning k-mers unique to one locus
against its family) are the k-mer granularity analogue of ``key_genes`` and
resolve reference pairs that locus-level coverage/identity cannot (a single
SNP moves the identity component by ~0.001 but flips the differential set
from 1.0 to 0.0).

The canonical code is the lexicographically smaller of a k-mer and its
reverse complement, packed two bits per base (A=0, C=1, G=2, T=3) — the
same convention as the gmlst cgMLST prefilter. Counting is a single
streaming pass; query codes are counted once per file and intersected with
each declared set (sets are small: hundreds to a few thousand members).
"""

from __future__ import annotations

import struct
from collections.abc import Iterator
from pathlib import Path

from pydantic import BaseModel, Field, field_validator

from gapit.errors import DatabaseError

_BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def iter_canonical_kmer_codes(sequence: str, k: int) -> Iterator[int]:
    """Yield the canonical (min of forward/reverse-complement) two-bit
    packed code for every non-ambiguous k-mer window. Runs with a rolling
    forward code and a rolling reverse code; any N (or other IUPAC letter)
    resets both windows so windows spanning it are skipped."""
    if k <= 0 or k > 31:
        raise ValueError(f"k must be in [1, 31]: got {k}")
    mask = (1 << (2 * k)) - 1
    top = 2 * (k - 1)
    fwd = 0
    rev = 0
    seen = 0
    for char in sequence:
        code = _BASE_INDEX.get(char.upper())
        if code is None:
            fwd = rev = seen = 0
            continue
        fwd = ((fwd << 2) | code) & mask
        rev = (rev >> 2) | (3 - code) << top
        seen += 1
        if seen >= k:
            yield fwd if fwd <= rev else rev


def count_sequence_kmers(sequence: str, k: int) -> set[int]:
    """The set of distinct canonical k-mer codes of one sequence."""
    return set(iter_canonical_kmer_codes(sequence, k))


def count_fasta_kmers(path: Path, k: int) -> set[int]:
    """Distinct canonical k-mer codes over every record of a FASTA file
    (plain or gzip/bzip2 via the shared text opener). Sequence lines are
    upper-cased; the whole file is one set."""
    from gapit.fasta import open_text

    codes: set[int] = set()
    chunks: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if chunks:
                    codes.update(iter_canonical_kmer_codes("".join(chunks), k))
                    chunks.clear()
            else:
                chunks.append(line.strip())
    if chunks:
        codes.update(iter_canonical_kmer_codes("".join(chunks), k))
    return codes


class KmerSetSpec(BaseModel, frozen=True):
    """One declared set: the ``of`` locus's k-mers minus the union of the
    ``against`` loci. ``n`` is filled at materialization time."""

    of: str = Field(min_length=1)
    against: tuple[str, ...] = ()
    n: int = Field(default=0, ge=0)


class KmerDbDocument(BaseModel, frozen=True):
    """The ``kmer.json`` sidecar: one k size for every set in the document
    (one k per database keeps counting single-pass)."""

    schema_name: str = Field(default="gapit.kmerdb/1", alias="schema")
    k: int = Field(ge=1, le=31)
    sets: dict[str, KmerSetSpec] = Field(default_factory=dict)

    @field_validator("schema_name")
    @classmethod
    def _schema(cls, value: str) -> str:
        if value != "gapit.kmerdb/1":
            raise ValueError(f"unsupported kmerdb schema: {value}")
        return value


def read_kmerdb(db_dir: Path) -> KmerDbDocument | None:
    """Load ``<db_dir>/kmer.json``; None when the database carries no
    k-mer substrate. A malformed file is a typed database error (the
    typing rules may reference its sets, so a silent fallback would
    under-call rules)."""
    path = db_dir / "kmer.json"
    if not path.is_file():
        return None
    try:
        return KmerDbDocument.model_validate_json(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise DatabaseError(
            f"malformed kmer sidecar: {path}: {exc}",
            code="KMERDB_MALFORMED",
            context={"file": str(path)},
        ) from exc


def _pack(codes: set[int]) -> bytes:
    return struct.pack(f"<{len(codes)}Q", *sorted(codes))


def materialize_kmer_sets(db_dir: Path, document: KmerDbDocument, sequences: Path) -> None:
    """Compute every declared set from the locus FASTA and write it as a
    sorted little-endian uint64 array under ``kmers/``, updating each
    spec's ``n`` in the stored ``kmer.json``. Validates that every ``of``
    and ``against`` locus exists (a stale set name is a build error, not a
    silent empty set)."""
    from gapit.fasta import iter_fasta

    loci: dict[str, str] = {}
    for record in iter_fasta(sequences):
        loci[record.id] = record.sequence.upper()
    missing = sorted(
        {name for spec in document.sets.values() for name in (spec.of, *spec.against)} - set(loci)
    )
    if missing:
        raise DatabaseError(
            f"kmer sets reference unknown loci: {missing}",
            code="KMERDB_UNKNOWN_LOCUS",
            context={"loci": ", ".join(missing)},
        )
    kmers_dir = db_dir / "kmers"
    kmers_dir.mkdir(parents=True, exist_ok=True)
    cache: dict[str, set[int]] = {}
    sized: dict[str, KmerSetSpec] = {}
    for name, spec in document.sets.items():
        anchor = cache.setdefault(spec.of, count_sequence_kmers(loci[spec.of], document.k))
        codes = set(anchor)
        for other in spec.against:
            codes -= cache.setdefault(other, count_sequence_kmers(loci[other], document.k))
        (kmers_dir / f"{name}.u64").write_bytes(_pack(codes))
        sized[name] = KmerSetSpec(of=spec.of, against=spec.against, n=len(codes))
    updated = KmerDbDocument(k=document.k, sets=sized)
    (db_dir / "kmer.json").write_text(
        updated.model_dump_json(indent=1, by_alias=True), encoding="utf-8"
    )


def load_kmer_set(db_dir: Path, name: str) -> frozenset[int]:
    """Read one materialized set back. Missing files are a typed error
    (the db validated set names at build/screen time)."""
    path = db_dir / "kmers" / f"{name}.u64"
    if not path.is_file():
        raise DatabaseError(
            f"kmer set not materialized: {name}",
            code="KMERDB_SET_MISSING",
            context={"set": name},
        )
    raw = path.read_bytes()
    return frozenset(struct.unpack(f"<{len(raw) // 8}Q", raw))


def set_hit_fraction(query: set[int], reference: frozenset[int]) -> float:
    """|query ∩ reference| / |reference| — the fraction of the declared
    set the query carries. An empty declared set scores 0.0 (a curator
    artifact: the locus shares every k-mer with its siblings)."""
    if not reference:
        return 0.0
    return sum(1 for code in reference if code in query) / len(reference)


__all__ = [
    "KmerDbDocument",
    "KmerSetSpec",
    "count_fasta_kmers",
    "count_sequence_kmers",
    "iter_canonical_kmer_codes",
    "load_kmer_set",
    "materialize_kmer_sets",
    "read_kmerdb",
    "set_hit_fraction",
]
