"""Victors provider (virulence factors, nucleotide) — transform only.

Upstream ``get_victors`` (abricate-get_db 1.4.0) cross-references two
phidias.us downloads, both served under ``.php`` URLs:

- ``gen_downloads.php`` — nucleotide CDS content (.ffn), ids like
  ``gi|115534241:2616-3152`` carrying the source GI and coordinates;
- ``gen_downloads_protein.php`` — protein content (.faa), headers like
  ``>gi|115534244|ref|YP_783826.1| hypothetical protein pCJ01p4
  [Campylobacter jejuni]``.

Protein headers are keyed by gi; each .ffn record looks its gi up in that
map for accession and product, falling back to ``gi|<gi>:<start>-<stop>``
and ``hypothetical protein``. The gene stays the ORIGINAL .ffn id —
upstream never renames it here.

Perl quirks (phase 7 notepad): upstream's ``.`` in both regexes matches
the literal ``|`` (escaped here); a protein header failing the regex is
skipped (``next unless``); an .ffn id failing its regex would reuse the
PREVIOUS match's stale capture variables — such records are skipped, the
same policy as the vfdb provider. Upstream's product line
(``$s->{DESC} =~ ... || 'hypothetical protein'``) is a void-context match
whose result is discarded; the notepad-locked digest semantics (map
product or ``hypothetical protein``) are implemented instead.
"""

import re
from collections.abc import Iterable
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

_NAME = "victors"
_FFN_NAME = "gen_downloads.php"  # nucleotide CDS content
_FAA_NAME = "gen_downloads_protein.php"  # protein content
_HYPOTHETICAL = "hypothetical protein"
_FUNCTION = ("virulence",)  # locked gapit/v1 func vocabulary (Wave F2c)

# Upstream m"^>gi.(\d+).ref.([^|]+). ([^[]+)" — its dots matched the pipes.
# \s+ eats the id/product separator; the product runs to end of line.
_FAA_HEADER = re.compile(r">gi\|(\d+)\|ref\|([^|]+)\|\s+(.+)")
# Upstream m/gi.(\d+):(\d+)-(\d+)/ searched against the .ffn record id.
_FFN_ID = re.compile(r"gi\|(\d+):(\d+)-(\d+)")


def _read_protein_map(path: Path) -> dict[str, tuple[str, str]]:
    """gi -> (accession, product) from .faa header lines.

    Headers failing the pattern are skipped (upstream ``next unless``);
    the product stops at the first `` [`` strain bracket (upstream
    ``([^[]+)``). A later duplicate gi overwrites an earlier one, exactly
    like upstream's hash assignment.
    """
    gi_map: dict[str, tuple[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            header = _FAA_HEADER.match(line)
            if header is None:
                continue
            gi, accession, product = (
                header.group(1),
                header.group(2),
                header.group(3).split(" [", 1)[0],
            )
            gi_map[gi] = (accession, product)
    return gi_map


def transform(workdir: Path) -> Iterable[Record]:
    """Yield Records from the two phidias.us downloads in ``workdir``.

    Each .ffn record keeps its ORIGINAL id as the gene; the gi carved out
    of that id looks accession/product up in the protein map, falling back
    to ``gi|<gi>:<start>-<stop>`` / ``hypothetical protein`` for an unknown
    gi. Records whose id carries no gi coordinates are skipped (upstream
    would build them from the previous match's stale captures).
    """
    gi_map = _read_protein_map(workdir / _FAA_NAME)
    for fasta in iter_fasta(workdir / _FFN_NAME):
        coords = _FFN_ID.search(fasta.id)
        if coords is None:
            continue
        gi = coords.group(1)
        fallback = f"gi|{gi}:{coords.group(2)}-{coords.group(3)}"
        accession, product = gi_map.get(gi, (fallback, _HYPOTHETICAL))
        yield Record(
            db=_NAME,
            gene=fasta.id,
            sequence=fasta.sequence,
            accession=accession,
            function=_FUNCTION,
            product=product,
            source_id=gi,
        )


PROVIDER = Provider(
    name=_NAME,
    description="Victors virulence factors",
    source_urls=(
        "http://phidias.us/victors/downloads/gen_downloads.php",
        "http://phidias.us/victors/downloads/gen_downloads_protein.php",
    ),
    dbtype="nucl",
    transform=transform,
)
