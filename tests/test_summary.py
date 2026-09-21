"""Unit tests for the summary core: table parsing + matrix aggregation (SPEC.md §6).

Behaviors verified against real abricate 1.4.0 on the synthetic fixtures in
tests/data/summary (labels are feature_a/feature_b — no biological data).
"""

from collections.abc import Callable
from pathlib import Path

import pytest

from gapit.errors import InputError
from gapit.summary import SummaryMatrix, SummaryParams, build_summary

FIXTURES = Path(__file__).parent / "data" / "summary"


class Recorder:
    """Collects warn() messages for assertion."""

    def __init__(self) -> None:
        self.messages: list[str] = []

    def __call__(self, message: str) -> None:
        self.messages.append(message)


def summarize(
    names: list[str | Path],
    identity: bool = False,
    nopath: bool = False,
    warn: Callable[[str], None] | None = None,
) -> tuple[SummaryMatrix, list[str]]:
    """build_summary over the given paths (relative names resolve via cwd),
    returning (matrix, warn messages)."""
    recorder = Recorder()
    matrix = build_summary(
        [Path(name) for name in names],
        SummaryParams(identity=identity, nopath=nopath),
        warn=warn if warn is not None else recorder,
    )
    return matrix, recorder.messages


def write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def summarize_text(
    tmp_path: Path, name: str, text: str, *, nopath: bool = False, identity: bool = False
) -> tuple[SummaryMatrix, list[str]]:
    path = write(tmp_path, name, text)
    return summarize([path], identity=identity, nopath=nopath)


def test_dutch_single_file_rows_keyed_by_file_column(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given exactly one report, When summarized, Then rows are keyed by that
    report's FILE column values (abricate "dutch mode", issue #32)."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["multi_sample.tsv"])
    assert [row.file for row in matrix.rows] == [
        "aa_assembly.fa",
        "mm_assembly.fa",
        "zz_assembly.fa",
    ]


def test_multi_file_rows_keyed_by_input_filename(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given several reports, When summarized, Then rows are keyed by the input
    filenames (as given), sorted lexicographically — not argument order."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["sample_b.tsv", "sample_a.tsv", "empty.tsv"])
    assert [row.file for row in matrix.rows] == [
        "empty.tsv",
        "sample_a.tsv",
        "sample_b.tsv",
    ]


def test_gene_columns_are_sorted_union() -> None:
    """Given reports mentioning different genes, When summarized, Then the gene
    universe is the sorted union across all files."""
    matrix, _ = summarize([FIXTURES / "sample_a.tsv", FIXTURES / "sample_b.tsv"])
    assert matrix.genes == ("feature_a", "feature_b")


def test_duplicate_gene_hits_are_semicolon_joined_in_row_order() -> None:
    """Given two feature_a hits in one file, When summarized, Then the cell is
    both %COVERAGE strings ';'-joined in file order and NUM_FOUND counts
    DISTINCT genes (2, not 3)."""
    matrix, _ = summarize([FIXTURES / "sample_a.tsv"])
    (row,) = matrix.rows
    assert row.num_found == 2
    assert row.cells["feature_a"] == ("99.50", "52.00")
    assert row.cells["feature_b"] == ("76.00",)


def test_cell_values_preserve_original_number_strings(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given cells, When summarized, Then values stay the original report
    strings verbatim (no re-formatting, no float round-trip)."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["sample_b.tsv"])
    (sample_b_row, sample_c_row) = matrix.rows
    assert sample_b_row.file == "sampleB.fa"
    assert sample_b_row.cells["feature_b"] == ("100.00",)
    assert sample_c_row.cells["feature_a"] == ("90.00",)


def test_zero_hit_file_appears_with_num_found_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given a header-only report among several inputs, When summarized, Then
    it still appears with NUM_FOUND 0 (non-dutch mode)."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["sample_a.tsv", "empty.tsv"])
    empty = matrix.rows[0]
    assert empty.file == "empty.tsv"
    assert empty.num_found == 0
    assert empty.cells == {}


def test_dutch_header_only_report_has_no_rows() -> None:
    """Given a single header-only report, When summarized, Then zero rows
    (dutch mode invents no keys)."""
    matrix, _ = summarize([FIXTURES / "empty.tsv"])
    assert matrix.rows == ()
    assert matrix.genes == ()


def test_duplicate_input_paths_warn_and_skip(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given the same path twice, When summarized, Then a warning names the
    file and the second occurrence contributes nothing (key = path as given)."""
    monkeypatch.chdir(FIXTURES)
    matrix, warnings = summarize(["sample_a.tsv", "sample_a.tsv"])
    assert warnings == ["Skipping duplicate file: sample_a.tsv"]
    assert [row.file for row in matrix.rows] == ["sample_a.tsv"]


def test_identity_param_selects_identity_column(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given --identity, When summarized, Then cells carry the %IDENTITY
    strings instead of %COVERAGE."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["sample_a.tsv", "sample_b.tsv"], identity=True)
    sample_a = next(row for row in matrix.rows if row.file == "sample_a.tsv")
    assert sample_a.cells["feature_a"] == ("98.75", "91.00")
    assert sample_a.cells["feature_b"] == ("95.10",)


def test_nopath_basenames_dutch_keys(tmp_path: Path) -> None:
    """Given dutch mode and --nopath, When summarized, Then keys are the
    basename of the report's FILE column values."""
    text = (
        (FIXTURES / "multi_sample.tsv")
        .read_text(encoding="utf-8")
        .replace("aa_assembly.fa", "/deep/dir/aa_assembly.fa")
    )
    matrix, _ = summarize_text(tmp_path, "deep.tsv", text, nopath=True)
    assert [row.file for row in matrix.rows] == [
        "aa_assembly.fa",
        "mm_assembly.fa",
        "zz_assembly.fa",
    ]


def test_nopath_multi_file_sorts_by_full_key_not_label(tmp_path: Path) -> None:
    """Given --nopath with keys 1dir/zeta.tsv and 2dir/mid.tsv, When
    summarized, Then rows keep as-given-key sort order even though the labels
    are basenames (verified upstream quirk: labels can look unsorted)."""
    header = (FIXTURES / "empty.tsv").read_text(encoding="utf-8")
    zeta = write(tmp_path, "1dir/zeta.tsv", header)
    mid = write(tmp_path, "2dir/mid.tsv", header)
    matrix, _ = summarize([zeta, mid], nopath=True)
    assert [row.file for row in matrix.rows] == ["zeta.tsv", "mid.tsv"]


def test_hash_prefixed_lines_after_first_are_skipped(tmp_path: Path) -> None:
    """Given a concatenated report (two #FILE headers), When summarized, Then
    only the first line builds the column map and later '#' lines are skipped
    as data."""
    text = (FIXTURES / "sample_a.tsv").read_text(encoding="utf-8") + (
        FIXTURES / "sample_b.tsv"
    ).read_text(encoding="utf-8")
    matrix, _ = summarize_text(tmp_path, "cat.tsv", text)
    assert [row.file for row in matrix.rows] == ["sampleA.fa", "sampleB.fa", "sampleC.fa"]


def test_first_row_without_hash_is_header_and_data(tmp_path: Path) -> None:
    """Given a noheader report whose first data row lacks '#', When
    summarized, Then that row is BOTH the column map and a data row
    (upstream quirk: zip(@hdr, @col) on the header line itself)."""
    text = "GENE\t%COVERAGE\t%IDENTITY\tFILE\nfeature_a\t88.00\t87.00\tk.fa\n"
    matrix, _ = summarize_text(tmp_path, "nh.tsv", text)
    assert matrix.genes == ("GENE", "feature_a")
    assert [(row.file, row.num_found) for row in matrix.rows] == [("GENE", 1), ("feature_a", 1)]
    assert matrix.rows[0].cells["GENE"] == ("%COVERAGE",)
    assert matrix.rows[1].cells["feature_a"] == ("88.00",)


def test_csv_separator_autodetected_per_file(monkeypatch: pytest.MonkeyPatch) -> None:
    """Given comma-separated reports (no --csv flag), When summarized, Then
    each file's separator is auto-detected [gapit-extension] and parsing
    matches the tab equivalents."""
    monkeypatch.chdir(FIXTURES)
    matrix, _ = summarize(["sample_a.csv", "sample_b.csv"])
    sample_a = next(row for row in matrix.rows if row.file == "sample_a.csv")
    assert sample_a.cells["feature_a"] == ("99.50", "52.00")


def test_missing_file_raises_input_error(tmp_path: Path) -> None:
    """Given a nonexistent report path, When summarized, Then InputError
    INPUT_NOT_FOUND (exit 5) with the file in context."""
    missing = tmp_path / "nope.tsv"
    with pytest.raises(InputError) as excinfo:
        summarize([missing])
    assert excinfo.value.exit_code == 5
    assert excinfo.value.code == "INPUT_NOT_FOUND"
    assert excinfo.value.context["file"] == str(missing)


def test_short_data_row_raises_input_error(tmp_path: Path) -> None:
    """Given a report with a data row shorter than the mapped columns, When
    summarized, Then InputError SUMMARY_MALFORMED with file and line context
    (gapit refuses upstream's silent undef)."""
    text = (FIXTURES / "empty.tsv").read_text(encoding="utf-8") + "sampleA.fa\tcontig1\t1\n"
    path = write(tmp_path, "short.tsv", text)
    with pytest.raises(InputError) as excinfo:
        summarize([path])
    assert excinfo.value.code == "SUMMARY_MALFORMED"
    assert excinfo.value.context["file"] == str(path)
    assert excinfo.value.context["line"] == "2"


def test_header_without_gene_column_raises_input_error(tmp_path: Path) -> None:
    """Given a first row lacking GENE/%COVERAGE/%IDENTITY columns, When
    summarized, Then InputError SUMMARY_MALFORMED naming the header."""
    path = write(tmp_path, "badhdr.tsv", "FILE\tNOTE\nk.fa\tnone\n")
    with pytest.raises(InputError) as excinfo:
        summarize([path])
    assert excinfo.value.code == "SUMMARY_MALFORMED"
    assert "GENE" in str(excinfo.value)


def test_empty_zero_byte_file_is_valid_with_zero_rows(tmp_path: Path) -> None:
    """Given a 0-byte report alongside another, When summarized, Then it still
    appears with NUM_FOUND 0 (upstream: zero rows, no error)."""
    zero = write(tmp_path, "zero.tsv", "")
    filled = write(tmp_path, "aaa.tsv", (FIXTURES / "sample_a.tsv").read_text(encoding="utf-8"))
    matrix, _ = summarize([filled, zero])
    assert [row.file for row in matrix.rows] == [str(filled), str(zero)]
    assert matrix.rows[1].num_found == 0


def test_dutch_empty_zero_byte_file_has_no_rows(tmp_path: Path) -> None:
    """Given a single 0-byte report, When summarized, Then no rows and no
    genes (verified upstream: header-only output)."""
    matrix, _ = summarize([write(tmp_path, "zero.tsv", "")])
    assert matrix.rows == ()
    assert matrix.genes == ()


def test_num_found_equals_distinct_gene_count() -> None:
    """Given one file with feature_a twice and feature_b once, When
    summarized, Then NUM_FOUND is the distinct-gene count 2."""
    matrix, _ = summarize([FIXTURES / "sample_a.tsv"])
    (row,) = matrix.rows
    assert row.num_found == len(row.cells) == 2


def test_nopath_dutch_same_basename_dirs_merge_into_one_row(tmp_path: Path) -> None:
    """Given one dutch report whose FILE column contains dirA/x.fa and
    dirB/x.fa, When summarized with --nopath, Then both keys merge into ONE
    row labelled x.fa (setdefault collision semantics, intentional): cells
    carry the union of genes in row order and NUM_FOUND counts the union's
    distinct genes."""
    header = (FIXTURES / "empty.tsv").read_text(encoding="utf-8")
    text = header + (
        "dirA/x.fa\tc1\t1\t80\t+\tfeature_a\t1-80\t================\t0\t99.50\t98.75\tdb\tA\tp\tR\n"
        "dirA/x.fa\tc1\t1\t80\t+\tfeature_c\t1-80\t================\t0\t99.00\t97.00\tdb\tA\tp\tR\n"
        "dirB/x.fa\tc1\t1\t80\t+\tfeature_b\t1-80\t================\t0\t76.00\t95.10\tdb\tA\tp\tR\n"
        "dirB/x.fa\tc1\t1\t80\t+\tfeature_c\t1-80\t================\t0\t88.00\t93.00\tdb\tA\tp\tR\n"
    )
    matrix, _ = summarize_text(tmp_path, "collide.tsv", text, nopath=True)
    (row,) = matrix.rows
    assert row.file == "x.fa"
    assert row.num_found == 3
    assert row.cells == {
        "feature_a": ("99.50",),
        "feature_c": ("99.00", "88.00"),
        "feature_b": ("76.00",),
    }
