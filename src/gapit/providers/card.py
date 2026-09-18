"""CARD provider (Wave B2): transform card.mcmaster.ca data into Records.

Upstream quirk: the source URL ``https://card.mcmaster.ca/latest/data`` has no
extension but serves a tar.bz2, so the download lands at ``<workdir>/data``.
Only the root-level ``card.json`` is extracted (upstream:
``tar xf card.tar.bz2 card.json``), matched by root-normalized member name
because the re-tarred archive prefixes members with ``./``. Every dict-valued
top-level entry whose
``model_type`` is ``"protein homolog model"`` becomes one Record; a homolog
model carrying ``model_param.snp`` is a hard error (upstream ``err()``); every
other model type — and any non-dict entry, like the ``_version`` string — is
    skipped. Sequence/function normalization is NOT done here: that is
    fetch_provider's job (B0 contract).

``Any`` is confined to the parsed-JSON values flowing through these helpers:
CARD node shapes are heterogeneous, and checking them into pydantic models
would dwarf the module. Every value crossing out to a Record is a str/int.
"""

import json
import re
import tarfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from gapit.errors import DatabaseError
from gapit.providers.common import Dbtype, Provider
from gapit.records import Record

NAME = "card"
DESCRIPTION = "CARD protein homolog resistance models"
SOURCE_URLS = ("https://card.mcmaster.ca/latest/data",)
DBTYPE: Dbtype = "nucl"

_MEMBER = "card.json"  # the one member we extract from the tarball root
_CHAR_SPACE = re.compile(r"\s")  # per character, like perl s/\s/_/g (drug classes)
_RUN_SPACE = re.compile(r"\s+")  # per run, like perl s/\s+/_/g (model names)


def _root_name(member_name: str) -> str:
    """Member path with a leading './' and a single leading '/' stripped:
    card.mcmaster.ca re-tarred its archive with './'-prefixed members
    (Wave E), so exact-name lookups miss."""
    return member_name.removeprefix("./").removeprefix("/")


def _models(tarball: Path) -> Iterator[dict[str, Any]]:
    """Extract card.json from the tarball; yield its dict-valued entries.

    The member is matched by root-normalized name (first hit in archive
    order). Non-dict top-level values (the ``_version`` metadata string) are
    skipped — upstream's ``next unless ref($g) eq 'HASH'``.
    """
    with tarfile.open(tarball, "r:bz2") as tar:
        member = next(
            (m for m in tar.getmembers() if _root_name(m.name) == _MEMBER),
            None,
        )
        if member is None:
            raise DatabaseError(
                f"no {_MEMBER} member in {tarball}",
                code="PROVIDER_INVALID",
                context={"archive": str(tarball), "expected": _MEMBER},
            )
        extracted = tar.extractfile(member)
        if extracted is None:
            raise DatabaseError(
                f"{_MEMBER} in {tarball} is not a regular file",
                code="PROVIDER_INVALID",
                context={"db": NAME},
            )
        card: Any = json.loads(extracted.read())
    for model in card.values():
        if isinstance(model, dict):
            yield model


def _drug_classes(model: dict[str, Any]) -> tuple[str, ...]:
    """Drug Class category names: ' antibiotic' stripped once, then each
    whitespace character -> '_' (upstream's two s/// on $abx), in category
    order — sorting is fetch_provider's job."""
    categories: Any = model.get("ARO_category", {})
    classes: list[str] = []
    for category in categories.values():
        if category.get("category_aro_class_name") == "Drug Class":
            name: str = category.get("category_aro_name", "")
            classes.append(_CHAR_SPACE.sub("_", name.replace(" antibiotic", "", 1)))
    return tuple(classes)


def _first_dna(model: dict[str, Any]) -> Any:
    """The dna_sequence dict of the lexicographically-first key under
    model_sequences.sequence (perl: ``my ($key) = sort keys %$dna``).

    Its siblings carry accession/strand/fmin/fmax; the sequence itself is the
    ``sequence`` key of that same dict.
    """
    sequences: Any = model.get("model_sequences", {}).get("sequence", {})
    first = sorted(sequences)[0]
    return sequences[first].get("dna_sequence", {})


def _model_record(model: dict[str, Any]) -> Record:
    """One 'protein homolog model' -> one Record (db=card)."""
    model_name: str = model.get("model_name", "")
    model_param: Any = model.get("model_param") or {}
    if "snp" in model_param:
        raise DatabaseError(
            f"{model_name} has model_param.snp",
            code="PROVIDER_INVALID",
            context={"model": model_name},
        )
    dna: Any = _first_dna(model)
    strand: Any = dna.get("strand", "+")
    fmin: Any = dna.get("fmin", 0)
    fmax: Any = dna.get("fmax", 0)
    start, stop = (fmax, fmin) if strand == "-" else (fmin, fmax)
    product: Any = model.get("ARO_description") or model.get("ARO_accession", "")
    return Record(
        db=NAME,
        gene=_RUN_SPACE.sub("_", model_name),
        sequence=dna.get("sequence", ""),
        accession=f"{dna.get('accession', '')}:{start}-{stop}",
        function=_drug_classes(model),
        product=product,
        source_id=model.get("ARO_accession", ""),
    )


def transform(workdir: Path) -> Iterator[Record]:
    """Provider transform: the tarball at ``<workdir>/data`` -> Records.

    Non-homolog model types are skipped BEFORE the model_param.snp check
    (upstream's statement order), so a variant model carrying an snp is
    skipped, not an error.
    """
    for model in _models(workdir / "data"):
        if model.get("model_type") == "protein homolog model":
            yield _model_record(model)


PROVIDER = Provider(
    name=NAME,
    description=DESCRIPTION,
    source_urls=SOURCE_URLS,
    dbtype=DBTYPE,
    transform=transform,
    snapshot="card.tar.gz",
)
