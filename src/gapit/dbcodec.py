"""gapit/v1 database header codec: tagged, percent-encoded sequence ids.

Native header format (Phase 7): ``gapit|db=<v>|gene=<v>|acc=<v>|func=<v>`` with
fixed key order and lowercase-hex percent escapes inside values. Encode is
total for any input; decode is strict (``DatabaseError HEADER_MALFORMED``).
Seqids without the ``gapit|`` prefix delegate to the legacy ``~~~`` rules in
``gapit.db`` (SPEC.md §4 step 5).

The fourth slot is named ``func`` because gapit answers presence/absence — the
values it carries are functional categories (AMR classes, virulence, O-antigen,
replicon, biocide), not only resistance (Wave F1 rename; outputs keep the
frozen ``RESISTANCE``/``resistance`` names).
"""

import re
from collections.abc import Sequence

from gapit.db import DbHeader, parse_db_header
from gapit.errors import DatabaseError

PREFIX = "gapit|"

# Values-only escape alphabet (probe-verified byte-identical through
# makeblastdb/blastn/minimap2); every other char travels verbatim.
_ENCODE_TABLE = str.maketrans(
    {
        "%": "%25",
        "|": "%7C",
        "=": "%3D",
        ";": "%3B",
        " ": "%20",
        "\t": "%09",
        "\n": "%0A",
        "\r": "%0D",
    }
)

_ESCAPE_RE = re.compile("%([0-9a-fA-F]{2})")
_BAD_ESCAPE_RE = re.compile(r"%(?![0-9a-fA-F]{2})")


def is_gapit_header(seqid: str) -> bool:
    """True iff the seqid carries the native ``gapit|`` tagged prefix."""
    return seqid.startswith(PREFIX)


def encode_seqid(db: str, gene: str, accession: str, function: Sequence[str]) -> str:
    """Render one record as a native gapit/v1 seqid (total: never raises).

    The function list is joined with ``;`` first, then encoded as one value,
    so it decodes back to the abricate display string ``a;b``.
    """
    return "|".join(
        (
            "gapit",
            f"db={db.translate(_ENCODE_TABLE)}",
            f"gene={gene.translate(_ENCODE_TABLE)}",
            f"acc={accession.translate(_ENCODE_TABLE)}",
            f"func={';'.join(function).translate(_ENCODE_TABLE)}",
        )
    )


def decode_seqid(seqid: str, default_db: str) -> DbHeader:
    """Parse a seqid into a ``DbHeader``.

    ``gapit|``-prefixed seqids use the strict native format: bad escapes,
    missing/empty gene, a missing ``db``/``acc`` KEY, and duplicate keys raise
    ``HEADER_MALFORMED``; unknown keys are skipped (forward compat); an empty
    ``db`` value falls back to ``default_db`` (symmetric with legacy) and an
    empty ``acc`` value is allowed (unpublished sequences); the ``func`` key is
    intentionally optional and defaults to ``""``. Anything else delegates to
    ``gapit.db.parse_db_header``.
    """
    if not is_gapit_header(seqid):
        return parse_db_header(seqid, default_db)
    fields: dict[str, str] = {}
    for segment in seqid.split("|")[1:]:
        key, sep, raw_value = segment.partition("=")
        if not sep:
            raise _malformed(seqid, "segment_without_key")
        if not key:
            raise _malformed(seqid, "empty_key")
        if key in fields:
            raise _malformed(seqid, f"duplicate_key:{key}")
        fields[key] = _decode_value(raw_value, key=key, seqid=seqid)
    gene = fields.get("gene", "")
    if not gene:
        raise _malformed(seqid, "missing_gene")
    if "db" not in fields:
        raise _malformed(seqid, "missing_db")
    if "acc" not in fields:
        raise _malformed(seqid, "missing_acc")
    return DbHeader(
        database=fields["db"] or default_db,
        gene=gene,
        accession=fields["acc"],
        # func was encoded as one ;-joined value; the decoded string already is
        # the abricate display form "a;b".
        function=fields.get("func", ""),
    )


def _decode_value(raw: str, *, key: str, seqid: str) -> str:
    """Percent-decode one raw value; every ``%`` must precede two hex digits."""
    if _BAD_ESCAPE_RE.search(raw):
        raise _malformed(seqid, f"invalid_percent_escape:{key}")
    return _ESCAPE_RE.sub(lambda match: chr(int(match.group(1), 16)), raw)


def _malformed(seqid: str, reason: str) -> DatabaseError:
    """Build the strict-decode failure (exit 4) with seqid + reason context."""
    return DatabaseError(
        f"malformed gapit/v1 sequence header: {reason}",
        code="HEADER_MALFORMED",
        context={"seqid": seqid, "reason": reason},
    )
