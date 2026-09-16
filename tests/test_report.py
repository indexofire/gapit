"""Unit tests for the ScreeningParams and Report models."""

import pytest
from pydantic import ValidationError

from gapit.blast import BlastRow
from gapit.hits import Hit, process_rows
from gapit.report import Report, ScreeningParams


def test_default_params() -> None:
    """Given no overrides, When constructed, Then SPEC.md §1 defaults apply."""
    params = ScreeningParams(db="ncbi")
    assert params.minid == 80.0
    assert params.mincov == 80.0
    assert params.threads == 1


def test_minid_bounds() -> None:
    """Given minid outside (0, 100], When constructed, Then ValidationError."""
    with pytest.raises(ValidationError):
        ScreeningParams(db="ncbi", minid=0.0)
    with pytest.raises(ValidationError):
        ScreeningParams(db="ncbi", minid=100.5)


def test_mincov_bounds() -> None:
    """Given mincov outside [0, 100], When constructed, Then ValidationError."""
    with pytest.raises(ValidationError):
        ScreeningParams(db="ncbi", mincov=-0.1)
    with pytest.raises(ValidationError):
        ScreeningParams(db="ncbi", mincov=100.01)


def test_threads_lower_bound() -> None:
    """Given threads < 1, When constructed, Then ValidationError."""
    with pytest.raises(ValidationError):
        ScreeningParams(db="ncbi", threads=0)


def test_report_wraps_hits_as_tuple() -> None:
    """Given a hit list, When wrapped in a Report, Then hits become a tuple."""
    row = BlastRow(
        qseqid="c",
        qstart=1,
        qend=100,
        qlen=100,
        sseqid="db~~~g~~~a~~~r",
        sstart=1,
        send=100,
        slen=100,
        sstrand="plus",
        evalue=1e-40,
        length=100,
        pident=100.0,
        gaps=0,
        gapopen=0,
        stitle="db~~~g~~~a~~~r product",
    )
    (hit,) = process_rows([row], mincov=80.0, default_db="db")
    report = Report(file="x.fa", hits=(hit,))
    assert isinstance(report.hits, tuple)
    assert isinstance(report.hits[0], Hit)


def test_report_is_frozen() -> None:
    """Given a Report, When a field is assigned, Then pydantic rejects it."""
    report = Report(file="x.fa", hits=())
    attribute = "file"
    with pytest.raises(ValidationError):
        setattr(report, attribute, "y.fa")
