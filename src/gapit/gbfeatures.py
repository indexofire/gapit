"""GenBank FEATURES/ORIGIN parsing for cluster databases (no Biopython).

One GenBank record = one locus: the ``LOCUS`` name is the locus id, ``CDS``
features become genes (coordinates 1-based inclusive, strand from
``complement(...)``), and the ``ORIGIN`` block becomes the locus sequence
(uppercased — the seqconvert contract). Two real-world styles are supported:

- Bakta/modern (the rightsholder's VP databases): the source feature carries
  at most a ``K locus`` note; the label falls back to that note's value (the
  locus id when absent) and the type ALWAYS falls back to the locus id.
- kaptive-style: ``/note="K locus: KL1"`` + ``/note="K type: K1"`` on the
  ``source`` feature; note labels auto-detect (kaptive ``find_label``
  semantics, simplified) — any note whose key part (before ``:``) contains
  ``locus`` sets the label, any key containing ``type`` sets the type.

Deliberate limits, enforced with typed errors: compound CDS locations
(``join()``/``order()``/anything but ``[complement(][<>]?int..[<]?int[)]``)
are REJECTED — unsupported compound CDS; ``<``/``>`` partial-end prefixes
are stripped (the gene is kept, the imprecision is not modelled); records
without ORIGIN sequence and loci without CDS are rejected. CDS naming
follows the kaptive convention: ``/gene`` -> ``/locus_tag`` -> ``/note``
(the real databases name some orf CDS only in a note) -> a positional
``<locus>_NN`` id (kaptive's own numbering for unnamed CDS). Also owns the
locus/gene models shared with the GFF parser, the ``gapit.features/1``
document, and its writer.
"""

import re
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from gapit.errors import InputError
from gapit.fasta import open_text


class GeneFeature(BaseModel, frozen=True):
    """One CDS: 1-based inclusive coordinates on the locus sequence."""

    gene_id: str
    start: int = Field(ge=1)
    end: int = Field(ge=1)
    strand: Literal["+", "-"]
    product: str = ""


class LocusFeatures(BaseModel, frozen=True):
    """One locus and its genes — the feature table of a cluster database."""

    id: str
    label: str
    type: str
    genes: tuple[GeneFeature, ...]


class LocusSequence(BaseModel, frozen=True):
    """A parsed locus: its feature table plus the ORIGIN/FASTA sequence."""

    features: LocusFeatures
    sequence: str


class FeaturesDocument(BaseModel, frozen=True):
    """gapit.features/1 — the feature table of a cluster database
    (``features.json``, written next to ``sequences``; no sequence data)."""

    model_config = ConfigDict(populate_by_name=True)

    schema_name: Literal["gapit.features/1"] = Field(default="gapit.features/1", alias="schema")
    loci: tuple[LocusFeatures, ...]


def write_features_document(loci: tuple[LocusFeatures, ...], path: Path) -> None:
    """Write gapit.features/1 as indented JSON (by alias) + trailing newline."""
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(FeaturesDocument(loci=loci).model_dump_json(indent=2, by_alias=True))
        handle.write("\n")


def _fail(code: str, message: str, **context: str) -> InputError:
    """Uniform typed failure; context always carries the file."""
    return InputError(message, code=code, context=context)


_SIMPLE_LOCATION = re.compile(r"^([<>]?)(\d+)\.\.([<>]?)(\d+)$")


def _parse_location(location: str, locus: str, path: Path) -> tuple[int, int, Literal["+", "-"]]:
    """CDS location -> (start, end, strand). Only ``[complement()]a..b`` is
    supported; anything compound (join/order/multi-interval) is rejected."""
    strand: Literal["+", "-"] = "+"
    inner = location.strip()
    if inner.startswith("complement(") and inner.endswith(")"):
        strand = "-"
        inner = inner[len("complement(") : -1]
    match = _SIMPLE_LOCATION.fullmatch(inner)
    if match is None:
        raise _fail(
            "UNSUPPORTED_LOCATION",
            f"{path}: locus {locus!r} has a compound CDS location (unsupported): {location!r}",
            locus=locus,
            location=location,
        )
    return int(match.group(2)), int(match.group(4)), strand


def _first(values: list[str]) -> str:
    """First occurrence of a repeated qualifier (parse order; repeats are
    rare for the single-valued keys we consume)."""
    return values[0] if values else ""


def _note_value(note: str, marker: str) -> str | None:
    """kaptive find_label, simplified: a note whose key part (before ``:``)
    contains the marker yields the stripped value part; else None."""
    key, sep, value = note.partition(":")
    if not sep or marker not in key.lower():
        return None
    return value.strip() or None


class _Record:
    """Accumulator for one GenBank record between ``LOCUS`` and ``//``."""

    def __init__(self, locus_id: str) -> None:
        self.id = locus_id
        self.notes: list[str] = []
        self.cds: list[tuple[str, dict[str, list[str]]]] = []
        self.sequence: list[str] = []
        # the feature between its header line and the next feature/section
        self.key = ""
        self.location = ""
        self.quals: dict[str, list[str]] = {}

    def close_feature(self) -> None:
        """Route the finished feature into notes (source) or cds (CDS)."""
        if self.key == "source":
            self.notes.extend(self.quals.get("note", ()))
        elif self.key == "CDS":
            self.cds.append((self.location, dict(self.quals)))
        self.key, self.location, self.quals = "", "", {}

    def start_feature(self, key: str, location: str) -> None:
        self.close_feature()
        self.key, self.location = key, location

    def qualifier(self, content: str) -> None:
        """``/key=value`` (quoted or bare) into the qualifier map (repeats
        append in order; a leading quote is stripped here, the closing
        partner by whichever continuation line carries it)."""
        key, _, value = content[1:].partition("=")
        value = value.strip()
        if value.startswith('"'):
            value = value[1:]
            if value.endswith('"'):
                value = value[:-1]  # complete single-line quoted value
        self.quals.setdefault(key, []).append(value)

    def continuation(self, content: str) -> None:
        """A column-22 line after the header: a wrapped location while no
        qualifier was seen, else a wrapped qualifier value."""
        piece = content.strip()
        if not self.quals:
            self.location += piece
            return
        last_key = next(reversed(self.quals))
        merged = self.quals[last_key][-1] + piece
        if merged.endswith('"'):
            merged = merged[:-1]
        self.quals[last_key][-1] = merged

    def _note_back(self, marker: str) -> str:
        """First source-note value whose key contains the marker, else the
        locus id (the Bakta fallback)."""
        for note in self.notes:
            value = _note_value(note, marker)
            if value is not None:
                return value
        return self.id

    def features(self, path: Path) -> LocusSequence:
        """Finalize into a LocusSequence, enforcing every documented rule."""
        sequence = "".join(self.sequence).upper()
        if not sequence:
            raise _fail(
                "INVALID_GENBANK",
                f"{path}: locus {self.id!r} has no ORIGIN sequence",
                locus=self.id,
            )
        if not self.cds:
            raise _fail(
                "LOCUS_WITHOUT_GENES",
                f"{path}: locus {self.id!r} carries no CDS features",
                locus=self.id,
            )
        label = self._note_back("locus")
        kind = self._note_back("type")
        return LocusSequence(
            features=LocusFeatures(
                id=self.id, label=label, type=kind, genes=self._genes(sequence, path)
            ),
            sequence=sequence,
        )

    def _genes(self, sequence: str, path: Path) -> tuple[GeneFeature, ...]:
        return tuple(
            self._gene(location, quals, sequence, path, ordinal)
            for ordinal, (location, quals) in enumerate(self.cds, start=1)
        )

    def _gene(
        self,
        location: str,
        quals: dict[str, list[str]],
        sequence: str,
        path: Path,
        ordinal: int,
    ) -> GeneFeature:
        start, end, strand = _parse_location(location, self.id, path)
        # /gene -> /locus_tag -> /note (kaptive names some orf CDS only in a
        # note) -> kaptive's positional <locus>_NN id for unnamed CDS.
        gene_id = (
            _first(quals.get("gene", []))
            or _first(quals.get("locus_tag", []))
            or _first(quals.get("note", []))
            or f"{self.id}_{ordinal:02d}"
        )
        if not 1 <= start <= end <= len(sequence):
            raise _fail(
                "GENE_OUT_OF_BOUNDS",
                f"{path}: locus {self.id!r} gene {gene_id!r} falls outside the locus sequence",
                locus=self.id,
                gene=gene_id,
            )
        return GeneFeature(
            gene_id=gene_id,
            start=start,
            end=end,
            strand=strand,
            product=_first(quals.get("product", [])),
        )


def parse_genbank_features(path: Path) -> list[LocusSequence]:
    """Parse a (plain/.gz/.bz2) GenBank file into one LocusSequence per
    record; every documented malformation raises a typed InputError."""
    loci: list[LocusSequence] = []
    record: _Record | None = None
    in_origin = False
    try:
        with open_text(path) as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if in_origin:
                    if line.startswith("//"):
                        assert record is not None
                        loci.append(record.features(path))
                        record, in_origin = None, False
                    else:
                        assert record is not None
                        # seqconvert ORIGIN rule: drop the 10-column coordinate
                        # prefix, remove remaining whitespace, keep other chars.
                        record.sequence.append(re.sub(r"\s+", "", line[10:]))
                    continue
                if line.startswith("ORIGIN"):
                    if record is not None:
                        record.close_feature()
                    in_origin = True
                    continue
                if line.startswith("//"):
                    if record is not None:
                        record.close_feature()
                        loci.append(record.features(path))
                        record = None
                    continue
                if line.startswith("LOCUS") and line[5:6] in (" ", "\t"):
                    if record is not None:
                        raise _fail(
                            "INVALID_GENBANK",
                            f"{path}: record {record.id!r} is not terminated by //",
                            locus=record.id,
                        )
                    record = _Record(line.split()[1])
                    continue
                if record is None or not line.startswith(" "):
                    continue  # pre-FEATURES keywords and column-1 header lines
                if line[:5] == "     " and line[5] != " ":  # feature header, column 6
                    assert record is not None
                    record.start_feature(line[5:21].strip(), line[21:].strip())
                    continue
                assert record is not None
                content = line[21:]
                if content.startswith("/"):
                    record.qualifier(content)
                else:
                    record.continuation(content)
            if record is not None:
                raise _fail(
                    "INVALID_GENBANK",
                    f"{path}: last record {record.id!r} is not terminated by //",
                    locus=record.id,
                )
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        raise _fail("INVALID_GENBANK", f"{path}: could not read input: {exc}") from exc
    if not loci:
        raise _fail("INVALID_GENBANK", f"{path}: no GenBank records found")
    return loci
