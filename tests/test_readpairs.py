"""Unit tests for positional-FASTQ sample auto-pairing (gapit.readpairs).

The wildcard workflow: ``gapit screen -d db *.fastq.gz`` hands gapit a flat
glob of reads files; pair_samples infers per-sample lanes from filename
conventions (rightsholder request: sample_1/2.fastq.gz, R1/R2 styles).
"""

from collections.abc import Callable
from pathlib import Path

from gapit.readpairs import SampleLanes, pair_samples


def collect_warnings() -> tuple[list[str], Callable[[str], None]]:
    """A warn sink plus the list it fills."""

    messages: list[str] = []

    def warn(message: str) -> None:
        messages.append(message)

    return messages, warn


def test_pair_samples_pairs_r1_r2_markers() -> None:
    """Given sample_R1/sample_R2 files, When paired, Then one sample with one
    paired lane."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("sample_R1.fq"), Path("sample_R2.fq")], warn=warn)
    assert samples == [
        SampleLanes(sample="sample", lanes=((Path("sample_R1.fq"), Path("sample_R2.fq")),))
    ]
    assert messages == []


def test_pair_samples_pairs_underscore_1_2_with_gz_extensions() -> None:
    """Given the rightsholder's sample_1/sample_2 .fq.gz shape, When paired,
    Then the extensions strip and one sample 's1' remains."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("s1_1.fq.gz"), Path("s1_2.fq.gz")], warn=warn)
    assert samples == [SampleLanes(sample="s1", lanes=((Path("s1_1.fq.gz"), Path("s1_2.fq.gz")),))]
    assert messages == []


def test_pair_samples_pairs_bcl2fastq_r1_001_markers() -> None:
    """Given bcl2fastq _R1_001/_R2_001 suffixes, When paired, Then the marker
    (not the sample) strips."""
    messages, warn = collect_warnings()
    samples = pair_samples(
        [Path("run_S1_R1_001.fastq.gz"), Path("run_S1_R2_001.fastq.gz")], warn=warn
    )
    assert samples == [
        SampleLanes(
            sample="run_S1",
            lanes=((Path("run_S1_R1_001.fastq.gz"), Path("run_S1_R2_001.fastq.gz")),),
        )
    ]
    assert messages == []


def test_pair_samples_keeps_lane_number_in_sample_key() -> None:
    """Given _L001_R1_001 lanes, When paired, Then L001 stays part of the
    sample key (one sample per lane directory, per the bcl2fastq layout)."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("iso_L001_R1_001.fq"), Path("iso_L001_R2_001.fq")], warn=warn)
    assert samples == [
        SampleLanes(
            sample="iso_L001",
            lanes=((Path("iso_L001_R1_001.fq"), Path("iso_L001_R2_001.fq")),),
        )
    ]
    assert messages == []


def test_pair_samples_pairs_dot_1_dot_2_markers() -> None:
    """Given sample.1/sample.2 markers, When paired, Then one sample."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("a.1.fastq"), Path("a.2.fastq")], warn=warn)
    assert samples == [SampleLanes(sample="a", lanes=((Path("a.1.fastq"), Path("a.2.fastq")),))]
    assert messages == []


def test_pair_samples_markers_and_extensions_are_case_insensitive() -> None:
    """Given lowercase r1 / mixed-case R2 and uppercase extensions, When
    paired, Then both the marker and the extension still strip."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("x_r1.FASTQ.GZ"), Path("x_R2.Fq")], warn=warn)
    assert samples == [SampleLanes(sample="x", lanes=((Path("x_r1.FASTQ.GZ"), Path("x_R2.Fq")),))]
    assert messages == []


def test_pair_samples_missing_mate_becomes_single_end_with_warning() -> None:
    """Given an R1 whose R2 is absent from the glob, When paired, Then a
    single-end lane plus one warning naming the file."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("solo_R1.fq.gz")], warn=warn)
    assert samples == [SampleLanes(sample="solo", lanes=((Path("solo_R1.fq.gz"), None),))]
    assert messages == ["no mate found for solo_R1.fq.gz — screening single-end"]


def test_pair_samples_lone_r2_becomes_single_end_with_warning() -> None:
    """Given an R2 with no R1, When paired, Then it still screens as its own
    single-end lane with a warning (warning, never an error)."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("z_R2.fq")], warn=warn)
    assert samples == [SampleLanes(sample="z", lanes=((Path("z_R2.fq"), None),))]
    assert messages == ["no mate found for z_R2.fq — screening single-end"]


def test_pair_samples_duplicate_names_merge_as_extra_lanes() -> None:
    """Given identical basenames in two directories (multi-lane sequencing),
    When paired, Then one sample with two lanes sorted by name then path."""
    messages, warn = collect_warnings()
    samples = pair_samples(
        [
            Path("lane2/iso_R1.fq"),
            Path("lane1/iso_R1.fq"),
            Path("lane2/iso_R2.fq"),
            Path("lane1/iso_R2.fq"),
        ],
        warn=warn,
    )
    assert samples == [
        SampleLanes(
            sample="iso",
            lanes=(
                (Path("lane1/iso_R1.fq"), Path("lane1/iso_R2.fq")),
                (Path("lane2/iso_R1.fq"), Path("lane2/iso_R2.fq")),
            ),
        )
    ]
    assert messages == []


def test_pair_samples_extra_r1_beyond_r2_count_warns() -> None:
    """Given two R1 files but one R2 in a sample (duplicate basenames across
    lane directories), When paired, Then one paired lane plus one single-end
    lane with a warning for the unpaired R1."""
    messages, warn = collect_warnings()
    samples = pair_samples(
        [
            Path("lane1/m_R1.fq"),
            Path("lane2/m_R1.fq"),
            Path("lane1/m_R2.fq"),
        ],
        warn=warn,
    )
    assert samples == [
        SampleLanes(
            sample="m",
            lanes=(
                (Path("lane1/m_R1.fq"), Path("lane1/m_R2.fq")),
                (Path("lane2/m_R1.fq"), None),
            ),
        )
    ]
    assert messages == ["no mate found for lane2/m_R1.fq — screening single-end"]


def test_pair_samples_marker_less_files_are_own_single_end_samples() -> None:
    """Given files with no detectable mate marker, When paired, Then each is
    its own single-end sample keyed by its stripped stem, with warnings."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("tank.fastq.gz"), Path("pond.fq")], warn=warn)
    assert samples == [
        SampleLanes(sample="pond", lanes=((Path("pond.fq"), None),)),
        SampleLanes(sample="tank", lanes=((Path("tank.fastq.gz"), None),)),
    ]
    assert messages == [
        "no mate found for pond.fq — screening single-end",
        "no mate found for tank.fastq.gz — screening single-end",
    ]


def test_pair_samples_sorts_samples_and_lanes_deterministically() -> None:
    """Given samples in non-sorted input order, When paired, Then samples sort
    lexicographically and each sample's lanes sort by original name."""
    messages, warn = collect_warnings()
    samples = pair_samples(
        [
            Path("beta_2.fq"),
            Path("alpha_2.fq"),
            Path("beta_1.fq"),
            Path("alpha_1.fq"),
        ],
        warn=warn,
    )
    assert [sample.sample for sample in samples] == ["alpha", "beta"]
    assert samples[0].lanes == ((Path("alpha_1.fq"), Path("alpha_2.fq")),)
    assert samples[1].lanes == ((Path("beta_1.fq"), Path("beta_2.fq")),)
    assert messages == []


def test_pair_samples_bz2_extensions_strip() -> None:
    """Given .fq.bz2 inputs in the pairing layer, When paired, Then the
    double extension strips the same as .gz (stem parsing only; whether the
    reads engine can read the file is a separate contract)."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("b_1.fq.bz2"), Path("b_2.fq.bz2")], warn=warn)
    assert samples == [SampleLanes(sample="b", lanes=((Path("b_1.fq.bz2"), Path("b_2.fq.bz2")),))]
    assert messages == []


def test_pair_samples_nested_markers_take_the_longest_suffix() -> None:
    """Given _R1_001 competing with the shorter _R1/_1 markers, When paired,
    Then the longest marker wins so the sample key keeps nothing of it."""
    messages, warn = collect_warnings()
    samples = pair_samples([Path("s_R1_001.fq"), Path("s_R2_001.fq")], warn=warn)
    assert samples == [SampleLanes(sample="s", lanes=((Path("s_R1_001.fq"), Path("s_R2_001.fq")),))]
    assert messages == []
