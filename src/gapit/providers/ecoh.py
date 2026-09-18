"""ecoh database provider — E. coli O and H antigen genes (srst2 EcOH).

Transform-only Wave B10 module: ``PROVIDER`` wires the pinned metadata into
the frozen B0 :class:`gapit.providers.common.Provider` contract and
:func:`transform` parses the downloaded ``EcOH.fasta`` into typed Records.
Functional categories are allele-prefix derived — ``fliC*`` → H-antigen,
``wzx``/``wzy``/``wzt``/``wzm`` → O-antigen — fallback ``antigen``.
"""

from collections.abc import Iterator
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "ecoh"
SOURCE_FILE = "EcOH.fasta"


def transform(workdir: Path) -> Iterator[Record]:
    """Parse ``workdir/EcOH.fasta`` into Records (db = ecoh).

    srst2 header convention: ``>cluster__gene__allele__seq_id acc;word;word``
    (abricate-get_db ``get_ecoh``). The id splits on ``__`` into exactly 4
    parts and the gene is part 3 — the allele, e.g. ``fliC-H1``. The
    description splits on ``;``: the first piece is the accession, the
    space-joined rest is the product. Ids without exactly 4 ``__``-parts are
    skipped: upstream perl would read an undefined gene field there, and a
    typed Record cannot carry one. ``function`` is allele-prefix derived
    per the module-docstring rule (fallback ``antigen``).
    """
    for fasta in iter_fasta(workdir / SOURCE_FILE):
        parts = fasta.id.split("__")
        if len(parts) != 4:
            continue
        gene = parts[2]
        if gene.startswith("fliC"):
            function: tuple[str, ...] = ("H-antigen",)
        elif gene.startswith(("wzx", "wzy", "wzt", "wzm")):
            function = ("O-antigen",)
        else:
            function = ("antigen",)
        description = fasta.description.split(";")
        yield Record(
            db=NAME,
            gene=gene,
            accession=description[0],
            function=function,
            product=" ".join(description[1:]),
            sequence=fasta.sequence,
            source_id=fasta.id,
        )


PROVIDER = Provider(
    name=NAME,
    description="E. coli O and H antigens (srst2 EcOH)",
    source_urls=("https://raw.githubusercontent.com/katholt/srst2/master/data/EcOH.fasta",),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
