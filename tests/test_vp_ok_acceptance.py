"""vp_ok acceptance suite (Phase 3): the literature-locked RIMD 2210633
call plus structural checks over unlabelled vpasero genomes.

The database comes from the wheel's bundled cluster snapshot (no provider
anymore — materialization is the only install path). Everything here is
skipif-guarded on the local presence of the real genomes (private
upstream; try the working repos, then the developer scratch copies) — the
offline suite never touches this file's fixtures. Only RIMD 2210633
(O3:K6, literature-locked) gets a pinned call; the vpasero panel carries
no ground truth in this repo, so those genomes assert structural validity
(vocabulary + confidence fields) and leave the observed calls to the
calibration report.
"""

import json
import re
from pathlib import Path

import pytest
from click.testing import CliRunner

from gapit.bundled import find_bundled
from gapit.cli import app
from gapit.fasta import iter_fasta

runner = CliRunner()

_VPAUTILS = Path.home() / "repos/github/vpautils"
_VPASERO_DIR = _VPAUTILS.parent / "vpasero/test/genomes"
_SCRATCH = Path("/tmp/opencode/vp_provider_test")


def _first(*candidates: Path) -> Path | None:
    return next((candidate for candidate in candidates if candidate.is_file()), None)


RIMD = _first(
    _VPAUTILS / "tests/data/GCF_000196095.1.fna",
    _SCRATCH / "GCF_000196095.1.fna",
)
VPASERO = sorted(
    path
    for path in (
        _first(_VPASERO_DIR / "06VPA0006.fna", _SCRATCH / "06VPA0006.fna"),
        _first(_VPASERO_DIR / "13VPA0055.fna", _SCRATCH / "13VPA0055.fna"),
        _first(_VPASERO_DIR / "17VPA0078.fna", _SCRATCH / "17VPA0078.fna"),
    )
    if path is not None
)
REAL_DATA = pytest.mark.skipif(
    RIMD is None,
    reason="local RIMD genome absent (private upstream)",
)


@pytest.fixture(scope="module")
def real_db(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The vp_ok database built ONCE by materializing the wheel's bundled
    cluster snapshot (one synthetic-carrier screen, zero network); returns
    its datadir."""
    datadir = tmp_path_factory.mktemp("vp-accept-dd")
    bundled = find_bundled("vp_ok")
    assert bundled is not None
    carrier = datadir.parent / "carrier.fna"
    wanted = {"OL3", "KL6"}
    with carrier.open("w", encoding="utf-8") as handle:
        for record in iter_fasta(bundled.sequences_path):
            if record.id in wanted:
                handle.write(f">{record.id} sampled\n{record.sequence}\n")
    result = runner.invoke(
        app,
        ["screen", str(carrier), "--db", "vp_ok", "--datadir", str(datadir), "--quiet"],
    )
    assert result.exit_code == 0, result.stderr
    return datadir


def phenotypes_of(datadir: Path, sample: Path) -> dict[str, dict[str, object]]:
    """The per-scheme phenotype calls of one screened sample (cluster/2)."""
    result = runner.invoke(
        app,
        [
            "screen",
            str(sample),
            "--db",
            "vp_ok",
            "--datadir",
            str(datadir),
            "--format",
            "json",
            "--quiet",
        ],
    )
    assert result.exit_code == 0, result.stderr
    report = json.loads(result.stdout)
    assert report["schema"] == "gapit.cluster/2"
    calls: dict[str, dict[str, object]] = report["files"][0]["phenotypes"]
    return calls


@REAL_DATA
class TestRimdAcceptance:
    def test_genome_identity_is_rimd_2210633(self) -> None:
        """Given the RIMD candidate file, When its first FASTA header is
        read, Then it names RIMD 2210633 — the acceptance ground truth is
        the literature-locked isolate, not a lookalike."""
        assert RIMD is not None
        header = RIMD.read_text(encoding="utf-8").splitlines()[0]
        assert "RIMD 2210633" in header

    def test_rimd_calls_o3_k6_with_sensible_confidence(self, real_db: Path) -> None:
        """Given RIMD 2210633 screened against the bundled vp_ok database,
        When the cluster/2 report renders, Then o_group calls O3, k_group
        calls K6, serotype composes O3:K6 — every scheme at score >= 0.9
        with high confidence (the locked acceptance)."""
        assert RIMD is not None
        phenotypes = phenotypes_of(real_db, RIMD)
        for name, expected in (("o_group", "O3"), ("k_group", "K6"), ("serotype", "O3:K6")):
            call = phenotypes[name]
            assert call["phenotype"] == expected, name
            score = call["score"]
            assert isinstance(score, (int, float)) and score >= 0.9, name
            assert call["confidence"] == "high", name


@REAL_DATA
@pytest.mark.skipif(len(VPASERO) < 2, reason="local vpasero genomes absent")
class TestVpaseroStructural:
    def test_unlabelled_genomes_render_valid_vocabulary(self, real_db: Path) -> None:
        """Given vpasero panel genomes (no ground truth in this repo), When
        screened, Then every scheme renders a phenotype inside the declared
        vocabulary (O-group O<n>/OUT, K-group K<n>/KUT, serotype
        "<o>:<k>") with a confidence field present — labels stay unpinned
        (observed calls live in the calibration report)."""
        for genome in VPASERO[:3]:
            phenotypes = phenotypes_of(real_db, genome)
            assert re.fullmatch(r"O\d+|OUT", str(phenotypes["o_group"]["phenotype"])), genome
            assert re.fullmatch(r"K\d+|KUT", str(phenotypes["k_group"]["phenotype"])), genome
            assert re.fullmatch(
                r"(O\d+|OUT):(K\d+|KUT)", str(phenotypes["serotype"]["phenotype"])
            ), genome
            for call in phenotypes.values():
                assert call["confidence"] in ("high", "low", "ambiguous"), genome
