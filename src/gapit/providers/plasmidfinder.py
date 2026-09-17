"""plasmidfinder database provider — CGE PlasmidFinder replicons (transform only).

Upstream ``get_plasmidfinder`` (abricate-get_db 1.4.0) downloads the bitbucket
HEAD.zip — an archive with an arbitrary top-level directory — unzips every
``*.fsa`` and, per record, splits the accession suffix off the id::

    >IncFII_1_NC_004631.1  ->  gene IncFII_1, accession NC_004631.1

Perl semantics preserved exactly (verified against the parity env): the
product is the ORIGINAL full id (DESC is assigned before the regex munging),
trailing underscore runs are stripped from the captured id, and a non-matching
id — or one that strips to empty — keeps the ORIGINAL id as gene with an empty
accession. Upstream additionally warns on the empty-id case; the transform has
no diagnostic channel, so the record is kept silently (``fetch_provider`` owns
stderr). The ``_1`` copy number STAYS on the gene: only the accession suffix
is removed.
"""

import re
import zipfile
from collections.abc import Iterator
from pathlib import Path

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "plasmidfinder"
ARCHIVE = "HEAD.zip"

# Upstream regex verbatim: group 2 is the accession ([A-Z]+ prefix or the
# literal NC_, digits, optional .version); group 1 is everything before the
# last viable underscore — greedy, so the copy number stays in group 1.
_ID = re.compile(r"^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)$")
_FUNCTION = ("replicon",)  # locked gapit/v1 func vocabulary (Wave F2c)


def transform(workdir: Path) -> Iterator[Record]:
    """Yield Records from every ``*.fsa`` inside ``workdir/HEAD.zip``.

    Members are extracted into the workdir (bitbucket archives nest under an
    arbitrary top-level directory) and matched with a recursive glob, in
    sorted order for determinism.
    """
    with zipfile.ZipFile(workdir / ARCHIVE) as archive:
        archive.extractall(workdir)
    for fsa in sorted(workdir.glob("**/*.fsa")):
        for fasta in iter_fasta(fsa):
            id_match = _ID.match(fasta.id)
            captured = id_match.group(1).rstrip("_") if id_match is not None else ""
            yield Record(
                db=NAME,
                gene=captured or fasta.id,
                accession=id_match.group(2) if id_match is not None else "",
                function=_FUNCTION,
                product=fasta.id,
                sequence=fasta.sequence,
                source_id=fasta.id,
            )


PROVIDER = Provider(
    name=NAME,
    description="CGE PlasmidFinder replicons",
    source_urls=("https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip",),
    dbtype="nucl",
    transform=transform,
)
