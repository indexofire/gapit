"""resfinder database provider — CGE ResFinder acquired resistance genes.

Upstream ``get_resfinder`` (abricate-get_db 1.4.0) git-clones the resfinder_db
repo; gapit downloads the same tree as the bitbucket ``HEAD.zip`` archive (no
git dependency) whose members sit under an arbitrary top-level directory.
``phenotypes.txt`` maps each gene to its function categories; every ``*.fsa``
carries ids like ``demoA_1_FAKE0001`` — gene prefix, copy number, accession.

Issue #62 repair: resfinder .fsa files can glue a record header onto the end
of the previous sequence line (a letter immediately followed by ``>``); a
newline is inserted there in the raw text before parsing, exactly like
upstream's ``sed s/([A-Z])>/\\1\\n>/gi`` (the /i makes it both cases).

Perl semantics kept: Class-cell pieces are filtered of unknown/notes/none
markers BEFORE stripping (upstream greps, then trims), the LAST phenotypes
row for a gene wins (plain hash assign), and the product is always the gene
prefix (upstream's ``$anno{$id}{DESC}`` is never assigned). Deviation: a
record whose id has no ``_<digits>_<accession>`` structure is SKIPPED —
upstream would emit it with undef fields (same policy as vfdb).
"""

import re
import zipfile
from collections.abc import Iterator
from pathlib import Path
from tempfile import TemporaryDirectory

from gapit.errors import DatabaseError
from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "resfinder"
ARCHIVE = "HEAD.zip"

_PHENOTYPES = "phenotypes.txt"
_GENE = re.compile(r"^(.*?)_\w+$")  # col0 minus its final _<word> chunk
_CLASS_SPLIT = re.compile(r",\s*")
_CLASS_FILTER = re.compile(r"(unknown|notes|^none)", re.IGNORECASE)
# abricate issue #62: a letter directly followed by '>' is a glued header.
_INLINE_HEADER = re.compile(r"([A-Za-z])>")
_ID = re.compile(r"^(.*?)_(\d+)_(\S+)$")  # base _ copy _ accession


def _read_classes(root: Path) -> dict[str, tuple[str, ...]]:
    """phenotypes.txt under the extracted archive: {gene prefix: classes}.

    Comment lines start with '#'; rows split on tabs; the gene prefix is
    col0 minus its final ``_<word>`` chunk; classes are the col2 pieces
    (split on comma + whitespace) that do not match unknown/notes/none —
    filtered BEFORE stripping, like upstream — each stripped. A later row
    for the same gene overwrites the earlier one (plain hash assign), and a
    missing Class cell means no classes (upstream splits undef into an
    empty list).
    """
    found = sorted(root.glob(f"**/{_PHENOTYPES}"))
    if not found:
        raise DatabaseError(
            f"{ARCHIVE} contains no {_PHENOTYPES}",
            code="PROVIDER_INVALID",
            context={"db": NAME},
        )
    classes: dict[str, tuple[str, ...]] = {}
    for line in found[0].read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        row = line.split("\t")
        gene = _GENE.match(row[0])
        if gene is None:
            continue
        cell = row[2] if len(row) > 2 else ""
        classes[gene.group(1)] = tuple(
            piece.strip()
            for piece in (_CLASS_SPLIT.split(cell) if cell else [])
            if not _CLASS_FILTER.search(piece)
        )
    return classes


def transform(workdir: Path) -> Iterator[Record]:
    """Yield Records from ``workdir/HEAD.zip`` (the bitbucket archive).

    The archive is extracted into a scratch directory (members sit under an
    arbitrary top-level directory); each ``*.fsa`` is repaired for glued
    headers into a scratch copy before FASTA parsing, so the download is
    never modified. Files are visited in sorted order for determinism.
    """
    with (
        TemporaryDirectory(prefix=".resfinder-extract.") as extract_name,
        TemporaryDirectory(prefix=".resfinder-repair.") as repair_name,
    ):
        extract_dir = Path(extract_name)
        with zipfile.ZipFile(workdir / ARCHIVE) as archive:
            archive.extractall(extract_dir)
        classes = _read_classes(extract_dir)
        for index, fsa in enumerate(sorted(extract_dir.glob("**/*.fsa"))):
            repaired = _INLINE_HEADER.sub(r"\1\n>", fsa.read_text(encoding="utf-8"))
            temp = Path(repair_name) / f"{index:06d}.fsa"
            temp.write_text(repaired, encoding="utf-8")
            for fasta in iter_fasta(temp):
                id_match = _ID.match(fasta.id)
                if id_match is None:
                    continue
                base = id_match.group(1)
                yield Record(
                    db=NAME,
                    gene=f"{base}_{id_match.group(2)}",
                    sequence=fasta.sequence,
                    accession=id_match.group(3),
                    function=classes.get(base, ()),
                    product=base,
                    source_id=fasta.id,
                )


PROVIDER = Provider(
    name=NAME,
    description="CGE ResFinder acquired resistance genes",
    source_urls=("https://bitbucket.org/genomicepidemiology/resfinder_db/get/HEAD.zip",),
    dbtype="nucl",
    transform=transform,
)
