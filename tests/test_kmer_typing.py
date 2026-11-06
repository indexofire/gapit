"""Tests for the gapit.kmerdb/1 k-mer substrate and key_kmers typing.

The canonical counter, differential set materialization, and the SNP-level
discrimination that locus coverage/identity cannot deliver (a single SNP
moves the identity component by ~0.001 but flips the differential set from
1.0 to 0.0)."""

import json
from pathlib import Path

import pytest

from gapit.kmerdb import (
    KmerDbDocument,
    count_sequence_kmers,
    iter_canonical_kmer_codes,
    load_kmer_set,
    materialize_kmer_sets,
    read_kmerdb,
    set_hit_fraction,
)


class TestCanonicalCounter:
    def test_single_kmer_forward(self):
        assert list(iter_canonical_kmer_codes("ACGT", 4)) == [0b00011011]

    def test_reverse_complement_symmetry(self):
        forward = "ACGTTGCA"
        reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        assert list(iter_canonical_kmer_codes(forward, 8)) == list(
            iter_canonical_kmer_codes(reverse, 8)
        )

    def test_n_resets_window(self):
        seq = "ACGTNACGT"
        assert len(list(iter_canonical_kmer_codes(seq, 5))) == 0

    def test_k_beyond_31_rejected(self):
        with pytest.raises(ValueError):
            list(iter_canonical_kmer_codes("A" * 32, 32))


class TestDifferentialSet:
    def test_single_snp_separation(self, tmp_path: Path):
        """Two 1200 bp sequences differing by one SNP: each one's
        differential k-mer set fully covers itself and zero-covers the
        other."""
        import random
        rng = random.Random(42)
        base = "".join(rng.choice("ACGT") for _ in range(1200))
        pos = 600
        snip = base[:pos] + ("A" if base[pos] != "A" else "C") + base[pos + 1 :]
        assert len(base) == len(snip)
        k = 31
        codes_a = count_sequence_kmers(base, k)
        codes_b = count_sequence_kmers(snip, k)
        diff_a = codes_a - codes_b
        diff_b = codes_b - codes_a
        assert len(diff_a) > 0
        assert len(diff_b) > 0
        assert set_hit_fraction(codes_a, frozenset(diff_a)) == 1.0
        assert set_hit_fraction(codes_b, frozenset(diff_a)) == 0.0
        assert set_hit_fraction(codes_b, frozenset(diff_b)) == 1.0
        assert set_hit_fraction(codes_a, frozenset(diff_b)) == 0.0


class TestMaterializeSets:
    def test_roundtrip(self, tmp_path: Path):
        """Materialize two loci with a shared block + unique tails, read
        the .u64 back, verify set sizes and membership."""
        shared = "ACGTTGCAAG" * 60
        locus_a = shared + "TTTTAAAACCCCGGGG" * 20
        locus_b = shared + "GGGGCCCCAAAATTTT" * 20
        sequences = tmp_path / "sequences"
        sequences.write_text(f">locA\n{locus_a}\n>locB\n{locus_b}\n")
        doc = KmerDbDocument(
            k=31,
            sets={
                "locA.diff": {"of": "locA", "against": ["locB"]},
                "locB.diff": {"of": "locB", "against": ["locA"]},
            },
        )
        materialize_kmer_sets(tmp_path, doc, sequences)
        stored = read_kmerdb(tmp_path)
        assert stored is not None
        assert stored.k == 31
        assert stored.sets["locA.diff"].n > 0
        assert stored.sets["locB.diff"].n > 0
        set_a = load_kmer_set(tmp_path, "locA.diff")
        codes_a = count_sequence_kmers(locus_a, 31)
        assert set_hit_fraction(codes_a, set_a) == 1.0
        codes_b = count_sequence_kmers(locus_b, 31)
        assert set_hit_fraction(codes_b, set_a) == 0.0

    def test_unknown_locus_fails(self, tmp_path: Path):
        sequences = tmp_path / "sequences"
        sequences.write_text(">locA\nACGT\n")
        doc = KmerDbDocument(k=31, sets={"x": {"of": "ghost", "against": []}})
        from gapit.errors import DatabaseError

        with pytest.raises(DatabaseError):
            materialize_kmer_sets(tmp_path, doc, sequences)
