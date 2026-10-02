"""gapit.typing — the declarative scoring schema of a typed database.

A typing document is a DATABASE-SIDE artifact (``typing.json`` inside a db
directory, installed via ``gapit db build --typing``): it declares how hits
against the database map to phenotype calls. This module owns the document
container and its validation (readers, reference checks, the cluster
single-scheme guard); the RULE models live in gapit.typing_rules, the
evaluation RESULT models in gapit.typing_results, the arithmetic in
gapit.typing_engine (cluster path) and gapit.typing_gene (gene path).

Two document versions parse:

- ``gapit.typing/2`` (current) — named multi-scheme: ``schemes`` is a list
  of :class:`TypingScheme` (``name``/``rules``/``cutoff``/
  ``ambiguity_margin``/``fallback``); the gene path evaluates EVERY scheme
  and keys the report's phenotype calls by name. Scheme names are
  ``[a-z0-9_]+`` and unique within the document.
- ``gapit.typing/1`` (legacy, cluster dbs) — the flat rules/cutoff/margin/
  fallback shape. It loads as ONE anonymous scheme named ``default``, so
  existing typed cluster databases keep their byte-identical behavior.

Rule model kinds (discriminated by ``model``, defined in gapit.typing_rules):

- ``weighted_genes`` — gene-presence score: positive/negative gene weight
  maps, an identity floor a gene must clear to count, and an any-of gene set
  the locus must match at least once.
- ``cluster_match`` — locus-level component score: weighted coverage,
  identity, and key-genes components, each with its own floor. Component
  weights are documented as summing near 1 but are NOT forced to (the total
  may exceed 1; the cutoff compares the raw score).
- ``learned_linear`` — a trained linear model over named features with a
  bias and free-text training provenance.

Validation failures raise typed DatabaseErrors (TYPING_MALFORMED) carrying
the file path; a missing file is the INPUT_NOT_FOUND input error (the --tsv
precedent).
"""

from pathlib import Path
from string import Formatter
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    TypeAdapter,
    ValidationError,
    model_validator,
)

from gapit.errors import DatabaseError, InputError
from gapit.gbfeatures import FeaturesDocument
from gapit.typing_rules import TypingRule, check_rule_references

# The scheme name a legacy gapit.typing/1 document degrades to (one
# anonymous scheme; report/1 phenotype keys on it for v1-typed gene dbs).
V1_SCHEME_NAME = "default"

_SCHEME_NAME_PATTERN = r"^[a-z0-9_]+$"


class TypingScheme(BaseModel, frozen=True):
    """One named scheme: rules plus the decision thresholds the evaluator
    applies to them (identical semantics for both db kinds).

    Stage-2 optional fields: ``control_gene`` (prs/ipaH-style gate — when
    set and that gene is absent by verdict semantics the whole scheme
    outputs its fallback with the note ``control gene absent``),
    ``unique_group`` + ``mixed_phenotype`` (when more than one gene of a
    group is present the scheme outputs ``mixed_phenotype`` with the pair
    in ``ambiguous`` and a mixed-infection note; both fields must be set
    together), and ``compose`` (``"{o_group}:{k_group}"`` — a rule-less
    scheme rendered from sibling schemes' calls after they evaluate).
    """

    name: str = Field(pattern=_SCHEME_NAME_PATTERN)
    rules: tuple[TypingRule, ...] = ()
    cutoff: float = Field(ge=0.0, le=1.0)
    ambiguity_margin: float = Field(ge=0.0)
    fallback: str
    control_gene: str | None = None
    unique_group: dict[str, list[str]] = Field(default_factory=dict)
    mixed_phenotype: str | None = None
    compose: str | None = None

    @model_validator(mode="after")
    def _rules_xor_compose(self) -> "TypingScheme":
        if self.compose is not None:
            if self.rules:
                raise ValueError(f"scheme {self.name!r}: a compose scheme declares no rules")
            if self.control_gene is not None or self.unique_group:
                raise ValueError(
                    f"scheme {self.name!r}: a compose scheme takes no control_gene/unique_group"
                )
        elif not self.rules:
            raise ValueError(f"scheme {self.name!r}: declare rules or compose")
        if bool(self.unique_group) != (self.mixed_phenotype is not None):
            raise ValueError(
                f"scheme {self.name!r}: unique_group and mixed_phenotype must be set together"
            )
        return self


class TypingDocumentV1(BaseModel, frozen=True):
    """gapit.typing/1 — the legacy flat document (one anonymous scheme);
    parsed then degraded via :meth:`as_document`."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.typing/1"] = Field(default="gapit.typing/1", alias="schema")
    rules: tuple[TypingRule, ...] = Field(min_length=1)
    cutoff: float = Field(ge=0.0, le=1.0)
    ambiguity_margin: float = Field(ge=0.0)
    fallback: str

    def as_document(self) -> "TypingDocument":
        """The same spec as one gapit.typing/2 default scheme."""
        return TypingDocument.model_validate(
            {
                "schemes": [
                    {
                        "name": V1_SCHEME_NAME,
                        "rules": [rule.model_dump() for rule in self.rules],
                        "cutoff": self.cutoff,
                        "ambiguity_margin": self.ambiguity_margin,
                        "fallback": self.fallback,
                    }
                ]
            }
        )


class TypingDocument(BaseModel, frozen=True):
    """gapit.typing/2 — the typing spec document of a typed database: named
    schemes, each independently evaluated. A v1 document is equivalent to a
    single ``default`` scheme."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.typing/2"] = Field(default="gapit.typing/2", alias="schema")
    schemes: tuple[TypingScheme, ...] = Field(min_length=1)

    @model_validator(mode="after")
    def _unique_scheme_names(self) -> "TypingDocument":
        names = [scheme.name for scheme in self.schemes]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ValueError(f"duplicate scheme names: {duplicates}")
        return self

    @model_validator(mode="after")
    def _compose_references_siblings(self) -> "TypingDocument":
        by_name = {scheme.name: scheme for scheme in self.schemes}
        for scheme in self.schemes:
            if scheme.compose is None:
                continue
            for placeholder in template_placeholders(scheme.compose):
                sibling = by_name.get(placeholder)
                if sibling is None or sibling is scheme or sibling.compose is not None:
                    raise ValueError(
                        f"compose scheme {scheme.name!r} placeholder {placeholder!r}"
                        " must name a non-composing sibling scheme"
                    )
        return self


def template_placeholders(template: str) -> tuple[str, ...]:
    """The placeholder names of a ``{name}`` compose template in
    first-occurrence order (deduplicated); any other brace syntax — empty
    placeholders, format specs/conversions, unclosed braces — is a
    ValueError (surfaced as TYPING_MALFORMED by the readers)."""
    names: list[str] = []
    for _literal, field, spec, conversion in Formatter().parse(template):
        if field is None:
            continue
        if not field or spec or conversion:
            raise ValueError(f"compose template placeholders must be plain {{name}}: {template!r}")
        if field not in names:
            names.append(field)
    return tuple(names)


def scheme_genes(scheme: TypingScheme) -> tuple[str, ...]:
    """A scheme's scheme-level gene references: the control gene plus every
    unique-group member (rules' references are walked separately)."""
    controlled = () if scheme.control_gene is None else (scheme.control_gene,)
    return (
        *controlled,
        *(member for members in scheme.unique_group.values() for member in members),
    )


_document_adapter: TypeAdapter[TypingDocument] = TypeAdapter(TypingDocument)
_v1_document_adapter: TypeAdapter[TypingDocumentV1] = TypeAdapter(TypingDocumentV1)
_json_object_adapter: TypeAdapter[dict[str, object]] = TypeAdapter(dict[str, object])


def validate_references(document: TypingDocument, features: FeaturesDocument) -> None:
    """Every gene/locus id a rule references must exist in the database's
    feature table (TYPING_UNKNOWN_GENE, one code for both kinds) — enforced
    at ``db build --typing`` time and again at screen time. Checks every
    scheme's rules and its scheme-level genes (control gene, unique
    groups)."""
    gene_ids = {gene.gene_id for locus in features.loci for gene in locus.genes}
    locus_ids = {locus.id for locus in features.loci}
    for scheme in document.schemes:
        for gene_id in scheme_genes(scheme):
            if gene_id not in gene_ids:
                raise DatabaseError(
                    f"typing scheme {scheme.name!r} references unknown gene {gene_id!r}"
                    " (not in the database features)",
                    code="TYPING_UNKNOWN_GENE",
                    context={"gene": gene_id, "scheme": scheme.name},
                )
        for rule in scheme.rules:
            check_rule_references(rule, gene_ids, locus_ids)


def single_scheme(document: TypingDocument) -> TypingScheme:
    """gapit.cluster/1 carries ONE phenotype slot per file: a cluster db's
    typing document must hold exactly one scheme (a v1 document parses as
    one default scheme; single-scheme v2 documents are fine)."""
    if len(document.schemes) != 1:
        names = ", ".join(scheme.name for scheme in document.schemes)
        raise DatabaseError(
            f"a cluster database supports exactly one typing scheme, got {len(document.schemes)}"
            f" ({names}); merge the schemes or screen with a gene database",
            code="TYPING_MALFORMED",
            context={"schemes": names},
        )
    return document.schemes[0]


def _reason(exc: ValidationError) -> str:
    """The validation failure rendered as one line (str(exc) is multiline)."""
    return " ".join(str(exc).split())


def _read_text(path: Path) -> str:
    """The typing file's text; unreadable bytes are TYPING_MALFORMED."""
    try:
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise DatabaseError(
            f"could not read typing document {path}: {exc}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        ) from exc


def _parse_json(text: str, path: Path) -> dict[str, object]:
    try:
        return _json_object_adapter.validate_json(text)
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed typing document {path}: {_reason(exc)}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        ) from exc


def read_typing_document(path: Path) -> TypingDocument:
    """Read and validate a gapit.typing/1 or /2 file. A /1 document loads
    as one anonymous ``default`` scheme (byte-equivalent evaluation); a
    missing file is the typed INPUT_NOT_FOUND input error, malformed
    content the TYPING_MALFORMED DatabaseError — both carrying the file in
    context."""
    if not path.is_file():
        raise InputError(
            f"typing file not found or unreadable: {path}",
            code="INPUT_NOT_FOUND",
            context={"file": str(path)},
        )
    data = _parse_json(_read_text(path), path)
    try:
        if data.get("schema") == "gapit.typing/1":
            return _v1_document_adapter.validate_python(data).as_document()
        return _document_adapter.validate_python(data)
    except ValidationError as exc:
        raise DatabaseError(
            f"malformed typing document {path}: {_reason(exc)}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        ) from exc


def typing_schema_of(path: Path) -> str:
    """The file's literal schema tag (``gapit.typing/1`` or ``/2``) — the
    value a build records in its manifest's ``typing_schema``. Assumes
    :func:`read_typing_document` accepted the file."""
    tag = _parse_json(_read_text(path), path).get("schema")
    if not isinstance(tag, str) or tag not in ("gapit.typing/1", "gapit.typing/2"):
        raise DatabaseError(
            f"malformed typing document {path}: unknown schema tag {tag!r}",
            code="TYPING_MALFORMED",
            context={"file": str(path)},
        )
    return tag
