"""CLI tests for the positional-FASTQ wildcard workflow (gapit screen *.gz).

`gapit screen -d db *.fastq.gz` routes to the reads engine without --r1/--r2
and auto-pairs samples from filenames (gapit.readpairs); one files[] entry
per sample, keyed by the sample name. The --r1/--r2 path is untouched (its
own tests in test_cli_reads.py / test_cli_reads2.py stay the parity proof).
"""

import glob
import json
import shutil
from pathlib import Path

import pytest
from pydantic import TypeAdapter
from typer.testing import CliRunner, Result

from gapit.cli import app
from gapit.cmd_screen_reads_args import reads_positional
from gapit.formats.reads_json import Reads2Document, ReadsDocument

READS_DB_DIR = Path(__file__).parent / "data" / "reads_db"
READS = Path(__file__).parent / "data" / "reads"
CONTIGS = Path(__file__).parent / "data" / "contigs"

runner = CliRunner()
reads_adapter = TypeAdapter(ReadsDocument)
reads2_adapter = TypeAdapter(Reads2Document)

# The pairing table for the filename conventions lives in test_readpairs.py;
# these lock the CLI surface on top of it.


def envelope(stderr: str) -> dict[str, str]:
    return json.loads([line for line in stderr.splitlines() if line.strip()][-1])


@pytest.fixture()
def datadir(tmp_path: Path) -> Path:
    target = tmp_path / "datadir"
    shutil.copytree(READS_DB_DIR, target)
    return target


@pytest.fixture()
def wildcard_dir(tmp_path: Path) -> Path:
    """A directory holding exactly the rightsholder's glob shape: two samples
    in two marker styles, all .fq.gz."""
    target = tmp_path / "glob"
    target.mkdir()
    for name in ("s1_1.fq.gz", "s1_2.fq.gz", "s2_R1.fq.gz", "s2_R2.fq.gz"):
        shutil.copy(READS / name, target / name)
    return target


def screen(datadir: Path, *extra: str) -> Result:
    return runner.invoke(app, ["screen", "--db", "tinyreads", "--datadir", str(datadir), *extra])


def test_reads_positional_classification(tmp_path: Path) -> None:
    """Given the extension/sniff classifier, When probed, Then reads and
    contig extensions answer outright and ambiguous names get one sniff."""
    assert reads_positional(Path("a.fastq.gz")) is True
    assert reads_positional(Path("a.FQ")) is True
    assert reads_positional(Path("a.fq.GZ")) is True
    assert reads_positional(Path("contig.fa")) is False
    assert reads_positional(Path("x.gbk.bz2")) is False
    fastq_txt = tmp_path / "run.txt"
    fastq_txt.write_text("@r001\nACGT\n+\nIIII\n", encoding="utf-8")
    assert reads_positional(fastq_txt) is True
    fasta_txt = tmp_path / "asm.txt"
    fasta_txt.write_text(">c1\nACGT\n", encoding="utf-8")
    assert reads_positional(fasta_txt) is False
    garbage_txt = tmp_path / "junk.txt"
    garbage_txt.write_text("not sequencing data\n", encoding="utf-8")
    assert reads_positional(garbage_txt) is False
    assert reads_positional(tmp_path / "missing.txt") is False


def test_wildcard_glob_screens_two_samples(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the user's exact shape (a shell glob of gz FASTQ handing four
    files to screen), When run with --format json (the single-document
    opt-in), Then one gapit.reads/1 document with one files[] entry per
    sample: s1 (the _1/_2 pair, tetX present, 12 reads) and s2 (the R1/R2
    pair, sub-threshold on its own, 6 reads) — no coverage leaks between
    samples."""
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--format", "json", *sorted(glob.glob("*.gz")))
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/1"
    assert document.params.read_type == "sr"
    assert [entry.reads for entry in document.files] == [["s1"], ["s2"]]
    (s1_gene,) = document.files[0].genes
    assert (s1_gene.gene, s1_gene.present, s1_gene.reads_mapped) == ("tetX", True, 12)
    (s2_gene,) = document.files[1].genes
    assert (s2_gene.gene, s2_gene.present, s2_gene.reads_mapped) == ("tetX", False, 6)
    assert s2_gene.breadth_pct < s1_gene.breadth_pct
    assert "Screening sample s1 reads:" in result.stderr
    assert "Detected 1 present genes in sample s1" in result.stderr
    assert "Detected 0 present genes in sample s2" in result.stderr


def test_wildcard_default_format_is_streaming_tsv(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the user's exact shape with NO --format, When run, Then stdout
    is the streaming table: the #SAMPLE header, then one row per PRESENT
    gene per sample (the default hides absent calls — s2's low-breadth
    tetX joins only under -A); --format csv answers with the comma
    spelling of the same table."""
    monkeypatch.chdir(wildcard_dir)
    files = sorted(glob.glob("*.gz"))
    result = screen(datadir, *files)
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    assert lines[0] == (
        "#SAMPLE\tGENE\tBREADTH%\tDEPTH\tREADS\tPRESENT\tDATABASE\tACCESSION\tPRODUCT"
    )
    assert [line.split("\t")[0] for line in lines[1:]] == ["s1"]
    assert lines[1].split("\t")[:2] == ["s1", "tetX"]
    assert lines[1].split("\t")[4:6] == ["12", "yes"]
    all_result = screen(datadir, "-A", *files)
    assert all_result.exit_code == 0
    all_lines = all_result.stdout.splitlines()
    assert [line.split("\t")[0] for line in all_lines[1:]] == ["s1", "s2"]
    assert all_lines[2].split("\t")[4:6] == ["6", "no"]
    csv_result = screen(datadir, "--format", "csv", *files)
    assert csv_result.exit_code == 0
    assert csv_result.stdout.startswith("#SAMPLE,GENE,BREADTH%,")
    assert csv_result.stdout.count("\n") == result.stdout.count("\n")


def test_wildcard_md_keys_sections_by_sample(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the wildcard run with --format md, When rendered, Then each
    sample gets its own `## <sample>` section and the frontmatter is STATIC
    (thresholds knowable before sample 1; no run totals — those live in the
    JSON document), so the path streams one chunk per completed sample."""
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--format", "md", *sorted(glob.glob("*.gz")))
    assert result.exit_code == 0
    assert result.stdout.startswith("---\nschema: gapit.reads/1\n")
    assert "## `s1`" in result.stdout
    assert "## `s2`" in result.stdout
    assert "files:" not in result.stdout and "genes_found:" not in result.stdout


def test_wildcard_unpaired_file_warns_and_screens_single_end(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a glob that caught only one mate, When run, Then the sample
    screens single-end with a stderr WARNING, and --quiet silences the
    warning while stdout stays the same document (created_at aside)."""
    monkeypatch.chdir(READS)
    args = ["--format", "json", "tetx_R1.fq"]
    plain = screen(datadir, *args)
    assert plain.exit_code == 0
    assert "WARNING: no mate found for tetx_R1.fq — screening single-end" in plain.stderr
    document = reads_adapter.validate_json(plain.stdout)
    assert document.files[0].reads == ["tetx"]
    quiet = screen(datadir, "--quiet", *args)
    assert quiet.exit_code == 0
    assert "no mate found" not in quiet.stderr
    normalized = plain.stdout.replace(json.loads(plain.stdout)["created_at"], "TS")
    assert normalized == quiet.stdout.replace(json.loads(quiet.stdout)["created_at"], "TS")


def test_wildcard_marker_less_file_is_its_own_sample(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a marker-less FASTQ in the glob, When run, Then it screens as
    its own single-end sample keyed by its stem."""
    monkeypatch.chdir(READS)
    result = screen(datadir, "--format", "json", "tetx_full.fq.gz")
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.files[0].reads == ["tetx_full"]
    (gene,) = document.files[0].genes
    assert gene.present is True


def test_mixed_assembly_and_reads_positionals_exit_2(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given a FASTA contig mixed into a FASTQ glob, When run, Then usage
    error exit 2 naming the reads files as the offenders."""
    monkeypatch.chdir(READS)
    result = screen(datadir, str(CONTIGS / "full.fa"), "s1_1.fq.gz", "s1_2.fq.gz")
    assert result.exit_code == 2
    error = envelope(result.stderr)
    assert error["code"] == "USAGE_ERROR"
    assert error["message"] == (
        "mixed assembly and reads inputs; screen them separately: s1_1.fq.gz, s1_2.fq.gz"
    )


def test_wildcard_rejects_blastn_only_flags(datadir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Given positional FASTQ with the blastn-only flags, When run, Then the
    reads-mode guards fire exactly as on the --r1/--r2 path — except --jobs,
    which is LEGAL here (it parallelizes the wildcard's samples; its own
    tests live in test_cli_reads_glob_jobs.py)."""
    monkeypatch.chdir(READS)
    for extra, message in [
        (["--merge-fragments", "s1_1.fq.gz"], "--merge-fragments is not available in reads mode"),
        (
            ["--minid", "90", "s1_1.fq.gz"],
            (
                "--minid/--mincov apply to blastn only; use --min-identity/--min-breadth"
                " with positional reads"
            ),
        ),
    ]:
        result = screen(datadir, *extra)
        assert result.exit_code == 2, extra
        error = envelope(result.stderr)
        assert error["code"] == "USAGE_ERROR"
        assert error["message"] == message


def test_wildcard_accepts_reads2_filters(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the wildcard run with --min-identity, When run, Then the
    reads-mode-only guard no longer fires (positional FASTQ IS reads mode)
    and the document is gapit.reads/2."""
    monkeypatch.chdir(wildcard_dir)
    result = screen(datadir, "--format", "json", "--min-identity", "95", *sorted(glob.glob("*.gz")))
    assert result.exit_code == 0
    document = reads2_adapter.validate_json(result.stdout)
    assert document.schema_name == "gapit.reads/2"
    assert document.params.min_identity == 95.0


def test_explicit_aligner_blastn_keeps_contig_pipeline(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given an explicit --aligner blastn over FASTQ positionals, When run,
    Then the blastn contig pipeline keeps them (explicit engine beats the
    auto-detect; today's normalization path unchanged)."""
    from gapit.db import make_blast_db

    make_blast_db(datadir / "tinyreads" / "sequences", "tinyreads")
    monkeypatch.chdir(READS)
    result = screen(datadir, "--aligner", "blastn", "--nopath", "tetx_full.fq")
    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    # The 100 bp reads screen as sub-coverage contig fragments (no rows
    # survive --mincov 80); the point is the engine: abricate TSV shape, not
    # a reads JSON document.
    assert lines[0].startswith("#FILE\t")
    assert "gapit.reads" not in result.stdout


def test_explicit_aligner_minimap2_over_fastq_still_exits_2(
    datadir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given --aligner minimap2 over a positional FASTQ, When run, Then the
    frozen assemblies-path guard still fires (the auto mode is a default-
    engine behavior only)."""
    monkeypatch.chdir(READS)
    result = screen(datadir, "--aligner", "minimap2", "tetx_full.fq")
    assert result.exit_code == 2
    assert envelope(result.stderr)["message"] == "minimap2 engine requires FASTA assemblies"


def test_wildcard_sniffs_ambiguous_extension_into_reads_mode(datadir: Path, tmp_path: Path) -> None:
    """Given a FASTQ file whose extension does not say so, When screened
    positionally, Then the content sniff routes it into reads mode."""
    sniffed = tmp_path / "illumina_run.data"
    shutil.copy(READS / "tetx_full.fq", sniffed)
    result = screen(datadir, "--format", "json", str(sniffed))
    assert result.exit_code == 0
    document = reads_adapter.validate_json(result.stdout)
    assert document.files[0].reads == ["illumina_run.data"]


def test_r1_path_output_is_byte_identical_to_positional_equivalent(
    datadir: Path, wildcard_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Given the same s1 pair via --r1/--r2 and via the positional wildcard,
    When both run, Then each document carries the same gene call for the
    sample (reads[] differs by contract: lane files vs sample key)."""
    monkeypatch.chdir(wildcard_dir)
    positional = screen(datadir, "--format", "json", "s1_1.fq.gz", "s1_2.fq.gz")
    assert positional.exit_code == 0
    explicit = runner.invoke(
        app,
        [
            "screen",
            "--r1",
            "s1_1.fq.gz",
            "--r2",
            "s1_2.fq.gz",
            "--db",
            "tinyreads",
            "--datadir",
            str(datadir),
            "--format",
            "json",
        ],
    )
    assert explicit.exit_code == 0
    positional_doc = reads_adapter.validate_json(positional.stdout)
    explicit_doc = reads_adapter.validate_json(explicit.stdout)
    assert explicit_doc.files[0].reads == ["s1_1.fq.gz", "s1_2.fq.gz"]
    assert positional_doc.files[0].reads == ["s1"]
    assert positional_doc.files[0].genes == explicit_doc.files[0].genes
