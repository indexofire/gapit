"""Tests for the gapit.kmerdb/1 k-mer substrate and key_kmers typing.

The canonical counter, differential set materialization, and the SNP-level
discrimination that locus coverage/identity cannot deliver (a single SNP
moves the identity component by ~0.001 but flips the differential set from
1.0 to 0.0)."""

from pathlib import Path

import pytest

from gapit.kmerdb import (
    KmerDbDocument,
    KmerSetSpec,
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
                "locA.diff": KmerSetSpec(of="locA", against=("locB",)),
                "locB.diff": KmerSetSpec(of="locB", against=("locA",)),
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
        doc = KmerDbDocument(k=31, sets={"x": KmerSetSpec(of="ghost")})
        from gapit.errors import DatabaseError

        with pytest.raises(DatabaseError):
            materialize_kmer_sets(tmp_path, doc, sequences)


class TestFastqKmerCounting:
    def test_fastq_kmers_match_fasta(self, tmp_path: Path):
        """The same sequence counted as FASTA and as simulated FASTQ
        produces the same canonical k-mer set."""
        import random
        rng = random.Random(7)
        seq = "".join(rng.choice("ACGT") for _ in range(3000))
        fa = tmp_path / "seq.fna"
        fa.write_text(f">seq\n{seq}\n")
        fq = tmp_path / "seq.fastq"
        reads = [seq[i : i + 250] for i in range(0, len(seq) - 250, 125)]
        with fq.open("w") as fh:
            for j, read in enumerate(reads):
                fh.write(f"@r{j}\n{read}\n+\n{'I' * len(read)}\n")
        from gapit.kmerdb import count_fastq_kmers, count_fasta_kmers

        fa_codes = count_fasta_kmers(fa, 31)
        fq_codes = count_fastq_kmers(fq, k=31)
        # Reads tile with overlap; edge k-mers may differ slightly
        assert fq_codes <= fa_codes  # no spurious k-mers
        assert len(fq_codes & fa_codes) / len(fa_codes) > 0.95  # >95% recall

    def test_paired_fastq(self, tmp_path: Path):
        """Two FASTQ files concatenated count the union."""
        seq = "ACGT" * 100
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text(f"@a\n{seq[:200]}\n+\n{'I' * 200}\n")
        r2.write_text(f"@b\n{seq[200:]}\n+\n{'I' * len(seq[200:])}\n")
        from gapit.kmerdb import count_fastq_kmers, count_sequence_kmers

        union = count_fastq_kmers(r1, r2, k=31)
        full = count_sequence_kmers(seq, 31)
        assert union == full


class TestContainmentIdentity:
    def test_perfect_match(self):
        from gapit.kmerdb import containment_to_identity

        assert containment_to_identity(1.0, 31) == 100.0

    def test_zero(self):
        from gapit.kmerdb import containment_to_identity

        assert containment_to_identity(0.0, 31) == 0.0

    def test_snp_level_containment(self):
        """A single SNP in a 1200bp sequence at k=31 reduces containment
        to ~0.97 but the identity estimate stays above 98% — enough for
        family-level floors while key_kmers delivers the SNP call."""
        import random
        from gapit.kmerdb import containment_to_identity, count_sequence_kmers

        rng = random.Random(3)
        base = "".join(rng.choice("ACGT") for _ in range(1200))
        snp = base[:600] + ("A" if base[600] != "A" else "C") + base[601:]
        codes_a = count_sequence_kmers(base, 31)
        codes_b = count_sequence_kmers(snp, 31)
        containment = len(codes_a & codes_b) / len(codes_a)
        ident = containment_to_identity(containment, 31)
        assert containment > 0.95
        assert ident > 98.0
