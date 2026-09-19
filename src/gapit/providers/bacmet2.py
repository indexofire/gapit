"""BacMet2 provider (experimentally confirmed, PROTEIN) — transform only.

Upstream ``get_bacmet2`` (abricate-get_db 1.4.0) downloads
``BacMet2_EXP_database.fasta`` — a protein file, hence the set's one
``dbtype="prot"`` provider (screening uses blastx). Headers look like::

    >BAC0098|ctpC|sp|P0A502|CTPC_MYCTU Probable manganese/zinc-exporting

The id token splits on ``|``: ``bac_id|gene|db|uniprot`` plus an optional
5th field (the UniProt entry name) that upstream ignores. Upstream
recombines ``gene-bac_id`` as the gene name and ``db:uniprot`` as the
accession; product falls back to the gene when the description is empty
(save_fasta's ``DESC || ID``). Function is the locked ``biocide`` constant:
BacMet covers biocides + metals but the 4-field id carries no class, so
``biocide`` is the user-approved approximation (Wave F2c).
"""

from collections.abc import Iterable
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

_NAME = "bacmet2"
_FUNCTION = ("biocide",)  # locked gapit/v1 func vocabulary (Wave F2c)


def transform(workdir: Path) -> Iterable[Record]:
    """Yield Records from ``workdir/BacMet2_EXP_database.fasta``.

    A record whose id has fewer than 4 pipe-separated fields is skipped
    (upstream perl would emit undef-indexed garbage). Sequence arrives raw
    — fetch_provider X-normalizes it for the prot dbtype.
    """
    for fasta in iter_fasta(workdir / "BacMet2_EXP_database.fasta"):
        parts = fasta.id.split("|")
        if len(parts) < 4:
            continue
        gene = f"{parts[1]}-{parts[0]}"
        yield Record(
            db=_NAME,
            gene=gene,
            sequence=fasta.sequence,
            accession=f"{parts[2]}:{parts[3]}",
            function=_FUNCTION,
            product=fasta.description or gene,
            source_id=fasta.id,
        )


PROVIDER = Provider(
    name=_NAME,
    description="BacMet2 experimentally confirmed biocide/resistance genes (protein)",
    source_urls=("http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta",),
    dbtype="prot",
    transform=transform,
    snapshot=None,
)
