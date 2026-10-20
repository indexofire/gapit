"""Acceptance for the vc_serogroup bundled database (handed over from gside)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gapit.typing_models import TypingDocument

_DB = Path(__file__).resolve().parents[1] / "src" / "gapit" / "data" / "dbs" / "vc_serogroup"


def test_bundled_files_present():
    for name in ("bundled.json", "records.jsonl", "sequences", "typing.json", "floors.json"):
        assert (_DB / name).is_file(), name


def test_records_and_sequences_agree():
    import json

    records = [json.loads(line) for line in (_DB / "records.jsonl").read_text().splitlines()]
    genes = {r["gene"] for r in records}
    assert genes == {"wbeN", "wbfR", "ctxA"}
    seq_genes = {
        line.split("gene=")[1].split("|")[0]
        for line in (_DB / "sequences").read_text().splitlines()
        if line.startswith(">gapit|")
    }
    assert seq_genes == genes


def test_typing_document_valid_and_exclusive():
    doc = TypingDocument.model_validate_json((_DB / "typing.json").read_text())
    scheme = doc.schemes[0]
    assert scheme.name == "serogroup"
    phenotypes = {r.phenotype: r for r in scheme.rules}
    assert phenotypes["O1"].requires == ["wbeN"]
    assert phenotypes["O139"].requires == ["wbfR"]
    assert "wbfR" in (phenotypes["O1"].excludes or [])
    assert "wbeN" in (phenotypes["O139"].excludes or [])
