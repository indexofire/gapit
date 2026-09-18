r"""ecoli_vf provider (phac-nml E. coli virulence factors) — transform only.

Upstream ``get_ecoli_vf`` (abricate-get_db 1.4.0) parses
``repaired_ecoli_vfs_shortnames.ffn`` in three steps per record::

    >VFG000748(gi:2865308) (espF) EspF [EspF (VF0182)] [Escherichia coli ...]

    1. id  ``^(\w+)(?:\((.*?)\))?$``  -> base id + optional paren accession
       (upstream ``die``s on a non-match; we SKIP the record instead —
       same deviation family as vfdb/ecoh, one bad header must not kill
       a fetch)
    2. accession = ``$2 || $1``        -> paren content; an absent OR EMPTY
       capture falls back to the base id (perl ``||`` is falsy-based)
    3. description: repeatedly strip trailing bracket groups
       (``s/\s\[.*?\]$//g``), then ``^(?:\((.*?)\)\s+)?(.*)$`` -> a
       leading paren group RENAMES the gene (overrides the base id), the
       remainder becomes the product (``$DESC || $ID`` fallback).
"""

import re
from collections.abc import Iterable
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

_NAME = "ecoli_vf"
_SOURCE_FILE = "repaired_ecoli_vfs_shortnames.ffn"

# Upstream regexes verbatim (the /x id regex ignores whitespace; joined here).
_ID = re.compile(r"^(\w+)(?:\((.*?)\))?$")
_TRAILING_BRACKETS = re.compile(r"\s\[.*?\]$")
_DESC = re.compile(r"^(?:\((.*?)\)\s+)?(.*)$")


def transform(workdir: Path) -> Iterable[Record]:
    """Yield Records from ``workdir/repaired_ecoli_vfs_shortnames.ffn``.

    Sequence and the leading/trailing-stripped description arrive via
    gapit.fasta; gene/accession/product derive per the docstring steps.
    """
    for fasta in iter_fasta(workdir / _SOURCE_FILE):
        id_match = _ID.match(fasta.id)
        if id_match is None:
            continue
        base = id_match.group(1)
        accession = id_match.group(2) or base  # perl $2 || $1: '' is falsy too
        description = fasta.description
        while (stripped := _TRAILING_BRACKETS.sub("", description)) != description:
            description = stripped
        desc_match = _DESC.match(description)
        assert desc_match is not None  # (.*) matches any string — cannot fail
        gene = desc_match.group(1) or base
        product = desc_match.group(2)
        yield Record(
            db=_NAME,
            gene=gene,
            sequence=fasta.sequence,
            accession=accession,
            function=("virulence",),
            product=product or gene,  # save_fasta: -desc => ($DESC || $ID)
            source_id=base,
        )


PROVIDER = Provider(
    name=_NAME,
    description="E. coli virulence factors (phac-nml)",
    source_urls=("https://github.com/phac-nml/ecoli_vf/raw/master/data/" + _SOURCE_FILE,),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
