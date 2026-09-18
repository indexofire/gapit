"""VFDB provider (set A, nucleotide) — transform only.

Upstream ``get_vfdb`` (abricate-get_db 1.4.0) decompresses
``VFDB_setA_nt.fas.gz`` and, per record, pulls the accession from the
``<gene>(<db>|<acc>.<version>)`` id suffix and renames the gene to the
leading paren group of the description::

    >VFG000676(gb|AAD32411) (lef) anthrax toxin lethal factor precursor ...

Perl quirk: a record failing either regex keeps the PREVIOUS record's $1/$2
(stale fields). We skip such records instead — see the phase 7 notepad.

Encoding quirk (Wave E): the real ``VFDB_setA_nt.fas.gz`` is not UTF-8-clean
(a latin-1 0xA0 nbsp crashes ``gapit.fasta``'s strict decoder), so the
decompressed bytes are decoded as latin-1 before parsing — see the comment
in :func:`transform`.
"""

import gzip
import re
from collections.abc import Iterable
from pathlib import Path
from tempfile import NamedTemporaryFile

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

_NAME = "vfdb"

# Upstream regexes verbatim (abricate issue #64 comment by @VGalata):
# accession = group 2; group 3 (".2" version suffix) is consumed, not used.
_ID = re.compile(r"^(\w+)\(\w+\|(\w+)(\.\d+)?\)$")
_DESC = re.compile(r"^\((.*?)\)")


def transform(workdir: Path) -> Iterable[Record]:
    """Yield Records from ``workdir/VFDB_setA_nt.fas.gz``.

    The gunzipped text is decoded as latin-1 into a scratch UTF-8 file for
    :func:`gapit.fasta.iter_fasta` (it takes a path); the scratch file is
    removed on completion or early generator close (argannot precedent).

    A record failing EITHER regex is skipped; product keeps the full
    description (upstream leaves DESC untouched) and source_id keeps the
    original id.
    """
    # Byte-preserving decode: upstream perl is byte-oriented here, and
    # latin-1 maps every byte 1:1 — descriptions keep the character and
    # sequence junk flows on to fetch_provider's N-normalization. Our
    # JSONL/FASTA outputs are UTF-8 by contract, so writing the scratch
    # file as UTF-8 re-encodes the character deterministically. Never
    # errors="replace": that would silently rewrite the data.
    text = gzip.decompress((workdir / "VFDB_setA_nt.fas.gz").read_bytes()).decode("latin-1")
    with NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        prefix=".vfdb.",
        suffix=".fas",
        dir=workdir,
        delete=False,
    ) as handle:
        handle.write(text)
        scratch = Path(handle.name)
    try:
        for fasta in iter_fasta(scratch):
            id_match = _ID.match(fasta.id)
            desc_match = _DESC.match(fasta.description)
            if id_match is None or desc_match is None:
                continue
            yield Record(
                db=_NAME,
                gene=desc_match.group(1),
                sequence=fasta.sequence,
                accession=id_match.group(2),
                function=("virulence",),
                product=fasta.description,
                source_id=fasta.id,
            )
    finally:
        scratch.unlink(missing_ok=True)


PROVIDER = Provider(
    name=_NAME,
    description="VFDB virulence factors (set A, nucleotide)",
    source_urls=("http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz",),
    dbtype="nucl",
    transform=transform,
    snapshot="vfdb.tar.gz",
)
