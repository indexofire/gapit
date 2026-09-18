"""ARG-ANNOT provider (acquired resistance genes) — transform only.

Upstream ``get_argannot`` (abricate-get_db 1.4.0) downloads
``ARG-ANNOT_NT_V6_July2019.txt`` — a nucleotide FASTA whose raw text is
littered with stray backslashes (upstream calls this "fix syntax errors in
the FASTA file"), so every backslash is stripped BEFORE parsing. Headers
look like::

    >(AGly)aac2-Ie:NC_011896:3039059-3039607:549

The id token splits on ``:``: gene is ``x[0]``, accession is ``x[1]:x[2]``
(later fields, e.g. the trailing length, are ignored — upstream never reads
``x[3]``). Product falls back to the gene when the header carries no
description (save_fasta's ``DESC || ID``); real ARG-ANNOT headers carry
none. Records with fewer than 3 colon-separated fields are SKIPPED
(documented deviation): upstream perl would splice an undefined ``x[2]``
into the accession; a typed Record refuses to carry that.
"""

from collections.abc import Iterator
from pathlib import Path
from tempfile import NamedTemporaryFile

from gapit.fasta import iter_fasta
from gapit.providers.common import Provider
from gapit.records import Record

NAME = "argannot"
SOURCE_FILE = "ARG-ANNOT_NT_V6_July2019.txt"


def transform(workdir: Path) -> Iterator[Record]:
    """Yield Records from ``workdir/ARG-ANNOT_NT_V6_July2019.txt``.

    The raw text is repaired first — every backslash stripped, upstream's
    global ``s/\\\\//g`` — and written to a scratch file inside the workdir
    so :func:`gapit.fasta.iter_fasta` can parse it (it takes a path). The
    scratch file is removed on completion or early generator close.
    Sequence arrives raw — fetch_provider N-normalizes it for the nucl
    dbtype.
    """
    repaired = (workdir / SOURCE_FILE).read_text(encoding="utf-8").replace("\\", "")
    with NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="\n",
        prefix=".argannot.",
        suffix=".txt",
        dir=workdir,
        delete=False,
    ) as handle:
        handle.write(repaired)
        scratch = Path(handle.name)
    try:
        for fasta in iter_fasta(scratch):
            parts = fasta.id.split(":")
            if len(parts) < 3:
                continue
            gene = parts[0]
            yield Record(
                db=NAME,
                gene=gene,
                sequence=fasta.sequence,
                accession=":".join(parts[1:3]),
                function=(),
                product=fasta.description or gene,
                source_id=fasta.id,
            )
    finally:
        scratch.unlink(missing_ok=True)


PROVIDER = Provider(
    name=NAME,
    description="ARG-ANNOT acquired resistance genes",
    # Upstream 301-migrated; the deep link 404s on both domains (verified 2026-09-17).
    # Wayback CDX: 2020-06-26 + 2026-01-14 captures share digest 5ZDVKAPZ4ZGIUDSYFVTFUYYY4CS5MUXZ
    # (2024 differs: suspect partial). The 2026-01-14 memento intermittently serves a 9KB
    # Wayback outage interstitial (LoadShardBlock datanode, observed 2026-09-18); the
    # 2020-06-26 memento probed good (text/plain, 2,139,519 bytes, same stable digest
    # family). Snapshot pinned, not latest; veto to upstream if it returns.
    # The id_ suffix requests the original WARC bytes with no Wayback rewrite/redirect
    # layer — the right form for machine fetching: verified 2026-09-18 (200, text/plain,
    # 2139519B, FASTA head ">(AGly)aac:AJ628983:1985-2539:555") while the plain
    # memento form flapped (interstitial twice, 45s apart).
    source_urls=(
        "https://web.archive.org/web/20200626214628id_/"
        "https://www.mediterranee-infection.com/wp-content/uploads/2019/09/"
        "ARG-ANNOT_NT_V6_July2019.txt",
    ),
    dbtype="nucl",
    transform=transform,
    snapshot=None,
)
