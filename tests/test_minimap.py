"""Unit tests for the minimap / COVERAGE_MAP arithmetic (SPEC.md §4).

Every expected string below was computed BY HAND from the Perl algorithm:
width = 15 - (broken?1:0); scale = length/width; xi = int(x/scale);
yi = int(y/scale); '=' where xi <= i <= yi; '/' after box width//2 when broken.
"""

from gapit.minimap import minimap


def test_full_length_hit_on_short_gene_all_equals() -> None:
    """Given L=60, x=1, y=60: scale=60/15=4.0 -> xi=int(0.25)=0, yi=int(15.0)=15
    -> every box 0..14 satisfies 0 <= i <= 15."""
    assert minimap(1, 60, 60, 0) == "==============="


def test_tiny_gene_full_hit_all_equals() -> None:
    """Given L=100, x=1, y=100: scale=100/15~6.667 -> xi=int(0.15)=0,
    yi=int(15.0)=15 -> all 15 boxes '='."""
    assert minimap(1, 100, 100, 0) == "==============="


def test_partial_hit_at_gene_end_marks_last_box_only() -> None:
    """Given L=15000, x=14001, y=15000: scale=1000.0 -> xi=14, yi=15 ->
    only box 14 (and the phantom 15) are '='."""
    assert minimap(14001, 15000, 15000, 0) == "..............="


def test_late_hit_marks_trailing_boxes() -> None:
    """Given L=15000, x=10001, y=15000: xi=10, yi=15 -> boxes 10..14 '='."""
    assert minimap(10001, 15000, 15000, 0) == "..........====="


def test_yi_truncation_can_exceed_last_box() -> None:
    """Given full hit x=1, y=15000, L=15000: yi=int(15.0)=15 > width-1=14 —
    the 1-based quirk; output is still 15 chars, all '='."""
    result = minimap(1, 15000, 15000, 0)
    assert result == "==============="
    assert len(result) == 15


def test_full_length_hit_on_long_gene_keeps_box_zero_empty() -> None:
    """Given x=1001, y=15000, L=15000: xi=int(1.001)=1 -> box 0 stays '.' even
    though coverage is ~93% (documented quirk — do not fix)."""
    assert minimap(1001, 15000, 15000, 0) == ".=============="


def test_broken_map_full_length() -> None:
    """Given L=1000, x=1, y=1000, broken=1: width=14, scale=1000/14~71.43 ->
    xi=0, yi>=13 -> 14 '=' with '/' inserted after box 7; total 15 chars."""
    result = minimap(1, 1000, 1000, 1)
    assert result == "========/======"
    assert len(result) == 15


def test_broken_map_partial_hit() -> None:
    """Given L=15000, x=1, y=5000, broken=1: width=14, scale~1071.43 ->
    xi=0, yi=int(4.667)=4 -> boxes 0..4 '=', slash after box 7."""
    assert minimap(1, 5000, 15000, 1) == "=====.../......"
