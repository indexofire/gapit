"""Unit tests for the nme_sero provider (N. meningitidis serogroup/ctrA/sodC)."""

import json
from pathlib import Path

from gapit.gene_floors import read_floors
from gapit.providers.nme_sero import PROVIDER, finalize, transform
from gapit.screen_defaults import read_screen_defaults


def _workdir(tmp_path: Path) -> Path:
    work = tmp_path / "work"
    work.mkdir()
    (work / "seroALLELES.fa").write_text(
        ">A\nACGTACGTACGT\n"
        ">B\nACGTACGTACGA\n"
        ">C\nACGTACGTACGG\n"
        ">W\nACGTACGTACGC\n"
        ">X\nACGTACGTACGT\n"
        ">Y\nACGTACGTACGA\n",
        encoding="utf-8",
    )
    (work / "ctrA.fa").write_text(">ctrA\nACGTACGTACGT\n", encoding="utf-8")
    return work


def test_transform_yields_only_serogroup_ctrA_and_sodC(tmp_path: Path) -> None:
    records = list(transform(_workdir(tmp_path)))
    genes = {record.gene for record in records}
    assert genes == {
        "SERO_A",
        "SERO_B",
        "SERO_C",
        "SERO_W",
        "SERO_X",
        "SERO_Y",
        "ctrA",
        "sodC",
    }
    assert len(records) == 8


def test_finalize_writes_typing_screen_and_floors_without_sidecars(tmp_path: Path) -> None:
    build = tmp_path / "build"
    build.mkdir()
    finalize(_workdir(tmp_path), build)

    assert (build / "typing.json").is_file()
    assert (build / "floors.json").is_file()
    assert (build / "screen.json").is_file()
    assert not (build / "meningotype").exists()

    screen = read_screen_defaults(build / "screen.json")
    assert screen is not None
    assert screen.reads.relative_to == "sodC"
    assert screen.reads.min_relative_depth == 0.10
    floors = read_floors(build / "floors.json")
    assert floors.genes["sodC"] == 80.0
    typing = json.loads((build / "typing.json").read_text(encoding="utf-8"))
    assert typing["schemes"][0]["name"] == "men_serogroup"
    assert typing["schemes"][0]["control_gene"] == "sodC"
    wy_rule = next(rule for rule in typing["schemes"][0]["rules"] if rule["phenotype"] == "W/Y")
    assert any("EX7E" in note for note in wy_rule["notes"])
    assert any("NEXT STEP" in note for note in wy_rule["notes"])
    assert any("menwy" in note for note in wy_rule["notes"])


def test_provider_is_serogroup_only_fetch_on_demand() -> None:
    assert PROVIDER.name == "nme_sero"
    assert PROVIDER.dbtype == "nucl"
    assert "GPL-3.0" in (PROVIDER.license or "")
    assert any("seroALLELES.fa" in url for url in PROVIDER.source_urls)
    assert any("ctrA.fa" in url for url in PROVIDER.source_urls)
    assert all("PorB" not in url and "BAST" not in url for url in PROVIDER.source_urls)
    assert all("mlst" not in url.lower() for url in PROVIDER.source_urls)
