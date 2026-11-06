"""Tests for the k-mer containment screening path (screening_kmer.py).

Covers: FASTQ k-mer counting, per-locus containment → LocusCall conversion,
typing evaluation on k-mer evidence, and an end-to-end simulated-reads
serotype call — the full cluster+FASTQ pipeline that bypasses minimap2."""

import json
import random
import tempfile
from pathlib import Path

from gapit.cluster import LocusCall
from gapit.gbfeatures import FeaturesDocument, LocusFeatures
from gapit.kmerdb import (
    KmerDbDocument,
    KmerSetSpec,
    count_fastq_kmers,
    count_sequence_kmers,
    materialize_kmer_sets,
)
from gapit.screening_kmer import best_phenotype, locus_calls_from_kmers
from gapit.typing_models import TypingDocument, read_typing_document

_K = 31


def _make_loci(seed: int = 42) -> tuple[str, str]:
    """Two 1200 bp sequences sharing the first 720 bp."""
    rng = random.Random(seed)
    shared = "".join(rng.choice("ACGT") for _ in range(720))
    tail_a = "".join(rng.choice("ACGT") for _ in range(480))
    tail_b = "".join(rng.choice("ACGT") for _ in range(480))
    return shared + tail_a, shared + tail_b


def _mini_features() -> FeaturesDocument:
    return FeaturesDocument(
        loci=(
            LocusFeatures(id="typeA", label="typeA", type="typeA", genes=()),
            LocusFeatures(id="typeB", label="typeB", type="typeB", genes=()),
        )
    )


def _materialize(tmp_path: Path, locus_a: str, locus_b: str) -> KmerDbDocument:
    sequences = tmp_path / "sequences"
    sequences.write_text(f">typeA\n{locus_a}\n>typeB\n{locus_b}\n")
    doc = KmerDbDocument(
        k=_K,
        sets={
            "typeA.diff": KmerSetSpec(of="typeA", against=("typeB",)),
            "typeB.diff": KmerSetSpec(of="typeB", against=("typeA",)),
        },
    )
    materialize_kmer_sets(tmp_path, doc, sequences)
    return doc


class TestLocusCallsFromKmers:
    def test_perfect_query_hits_own_locus(self, tmp_path: Path):
        locus_a, _ = _make_loci()
        doc = _materialize(tmp_path, locus_a, _make_loci()[1])
        query = count_sequence_kmers(locus_a, _K)
        calls = locus_calls_from_kmers(query, tmp_path, _mini_features(), doc)
        assert len(calls) >= 1
        assert calls[0].locus == "typeA"
        assert calls[0].coverage_pct > 99.0
        assert calls[0].identity_pct > 99.0
        assert calls[0].kmer_hits.get("typeA.diff") == 1.0

    def test_no_match_yields_empty(self, tmp_path: Path):
        rng = random.Random(99)
        unrelated = "".join(rng.choice("ACGT") for _ in range(1000))
        a, b = _make_loci()
        doc = _materialize(tmp_path, a, b)
        query = count_sequence_kmers(unrelated, _K)
        calls = locus_calls_from_kmers(query, tmp_path, _mini_features(), doc)
        assert len(calls) == 0

    def test_rank_ordering(self, tmp_path: Path):
        locus_a, locus_b = _make_loci()
        doc = _materialize(tmp_path, locus_a, locus_b)
        query = count_sequence_kmers(locus_a, _K)
        calls = locus_calls_from_kmers(query, tmp_path, _mini_features(), doc)
        for i, call in enumerate(calls):
            assert call.rank == i + 1


class TestBestPhenotype:
    def _load_typing(self) -> TypingDocument:
        doc: dict[str, object] = {
            "schema": "gapit.typing/3",
            "schemes": [
                {
                    "name": "test",
                    "cutoff": 0.70,
                    "ambiguity_margin": 0.05,
                    "fallback": "untypeable",
                    "rules": [
                        {
                            "model": "cluster_match",
                            "phenotype": "Type A",
                            "coverage": {"weight": 0.5, "floor": 50.0},
                            "identity": {"weight": 0.2, "floor": 80.0},
                            "key_genes": {"weight": 0.05, "genes": []},
                            "key_kmers": {"weight": 0.25, "k": 31, "set": "typeA.diff"},
                        },
                        {
                            "model": "cluster_match",
                            "phenotype": "Type B",
                            "coverage": {"weight": 0.5, "floor": 50.0},
                            "identity": {"weight": 0.2, "floor": 80.0},
                            "key_genes": {"weight": 0.05, "genes": []},
                            "key_kmers": {"weight": 0.25, "k": 31, "set": "typeB.diff"},
                        },
                    ],
                }
            ],
        }
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
            json.dump(doc, f)
            f.flush()
            return read_typing_document(Path(f.name))

    def test_perfect_match_calls_phenotype(self):
        typing = self._load_typing()
        locus = LocusCall(
            locus="typeA",
            label="typeA",
            type="typeA",
            coverage_pct=99.5,
            identity_pct=99.0,
            rank=1,
            genes=(),
            missing=(),
            kmer_hits={"typeA.diff": 1.0},
        )
        assert best_phenotype([locus], typing) == "Type A"


class TestFastqEndToEnd:
    def test_simulated_reads_serotype_call(self, tmp_path: Path):
        locus_a, locus_b = _make_loci(seed=7)
        reads = [locus_a[i : i + 250] for i in range(0, len(locus_a) - 250, 12)]
        fastq = tmp_path / "reads.fastq"
        with fastq.open("w") as fh:
            for i, seq in enumerate(reads):
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
        codes = count_fastq_kmers(fastq, k=_K)
        assert len(codes) > 500
        full_a = count_sequence_kmers(locus_a, _K)
        full_b = count_sequence_kmers(locus_b, _K)
        conta = len(codes & full_a) / len(full_a)
        contb = len(codes & full_b) / len(full_b)
        assert conta > 0.95
        assert conta > contb
