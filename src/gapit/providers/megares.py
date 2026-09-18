"""megares database provider — MEGARes v3 antimicrobial resistance genes.

Transform-only Wave B9 module: ``PROVIDER`` wires the pinned metadata into
the frozen B0 :class:`gapit.providers.common.Provider` contract and
:func:`transform` unpacks the downloaded ``megares_v3.00.zip`` and parses
every ``megares_drugs_*.fasta`` inside it into typed Records.
"""

import zipfile
from collections.abc import Iterator
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "megares"
ARCHIVE = "megares_v3.00.zip"
_DATABASE_GLOB = "megares_drugs_*.fasta"


def transform(workdir: Path) -> Iterator[Record]:
    """Extract ``workdir/megares_v3.00.zip`` and parse every nested
    ``megares_drugs_*.fasta`` into Records (db = megares).

    MEGARes v3 header convention (abricate-get_db ``get_megares``):
    ``>id|type|class|mech|group|note`` where the note field exists only on
    records requiring SNP confirmation — a non-empty note skips the record
    (SPEC §8). Keepers map gene=group, accession=source_id=id (the MEG_
    number), product=colon-joined type/class/mech/group; function = the
    ``class`` field (x[2], e.g. Tetracyclines) as a 1-tuple — upstream
    sets no ABX, the class IS the functional category. Archives often nest
    the fasta one directory deep — upstream
    ``unzip -j`` flattens, here the glob is recursive (sorted, so multi-file
    archives parse deterministically).

    Documented deviation: upstream perl splits into six list variables and
    keeps records even when fields are missing (undef gene, empty product
    pieces); gapit skips ids with fewer than 5 pipe-fields or an empty among
    the 5 leading ones. ``split('|', 5)`` mirrors perl's implicit
    list-assignment limit: any pipes past the group fold into the note, and
    a folded note is non-empty, so over-long ids skip exactly like upstream.
    """
    with zipfile.ZipFile(workdir / ARCHIVE) as archive:
        archive.extractall(workdir)
    for fasta_path in sorted(workdir.rglob(_DATABASE_GLOB)):
        for fasta in iter_fasta(fasta_path):
            parts = fasta.id.split("|", 5)
            if len(parts) < 5 or "" in parts[:5]:
                continue
            if len(parts) == 6 and parts[5]:
                continue
            yield Record(
                db=NAME,
                gene=parts[4],
                accession=parts[0],
                function=(parts[2],),
                product=":".join(parts[1:5]),
                sequence=fasta.sequence,
                source_id=parts[0],
            )


PROVIDER = Provider(
    name=NAME,
    description="MEGARes antimicrobial resistance genes",
    source_urls=("https://www.meglab.org/downloads/megares_v3.00.zip",),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
