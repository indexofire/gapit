"""COVERAGE_MAP minimap — exact abricate arithmetic (SPEC.md §4)."""


def minimap(x: int, y: int, length: int, broken: int) -> str:
    """15-char coverage bar: '=' over [x, y] scaled onto the map width.

    Replicates the Perl int() truncation quirks verbatim — do not "fix" the
    1-based coordinate artifacts (yi may exceed width-1). ``broken`` is the
    gap-open count; a broken map is 14 boxes with '/' after box width//2.
    """
    width = 15 - (1 if broken > 0 else 0)
    scale = length / width
    xi = int(x / scale)
    yi = int(y / scale)
    chars: list[str] = []
    for i in range(width):
        chars.append("=" if xi <= i <= yi else ".")
        if broken > 0 and i == width // 2:
            chars.append("/")
    return "".join(chars)
