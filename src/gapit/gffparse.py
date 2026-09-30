"""GFF3 cluster parsing: CDS rows plus embedded ``##FASTA`` or sidecar FASTA.

Each FASTA sequence is one locus (id = the header's first token, label/type
fall back to the id — GFF carries no locus labels) and the ``type=CDS`` rows
anchored on that sequence are its genes. GFF coordinates are 1-based
inclusive and pass through unchanged. Attribute lookup order is ``gene``,
``Name``, ``locus_tag``; ``product`` is optional. Percent-escapes in
attribute values are NOT decoded (none of the target databases use them).

The FASTA comes from the embedded ``##FASTA`` directive or from a sidecar
``<same-stem>.fa/.fna/.fasta`` (each also tried with ``.gz``); when neither
exists a typed error lists every searched path. ``##gff-version 3`` must be
the first line. CDS rows on sequences the FASTA does not carry, loci without
CDS, CDS without a name attribute, and genes outside their locus bounds are
rejected with typed errors (codes shared with gapit.gbfeatures).
"""

from pathlib import Path
from typing import IO

from gapit.errors import InputError
from gapit.fasta import iter_fasta, open_text
from gapit.gbfeatures import GeneFeature, LocusFeatures, LocusSequence

_SIDECAR_SUFFIXES = (".fa", ".fna", ".fasta")
_GFF_SUFFIXES = (".gff3", ".gff")


def _fail(code: str, message: str, **context: str) -> InputError:
    return InputError(message, code=code, context=context)


def _attributes(field: str) -> dict[str, str]:
    """``k=v;k=v`` pairs; pieces without ``=`` and empty keys are dropped."""
    pairs: dict[str, str] = {}
    for part in field.split(";"):
        key, sep, value = part.partition("=")
        if sep and key.strip():
            pairs[key.strip()] = value.strip()
    return pairs


def _gene(row: list[str], path: Path, line_no: int) -> tuple[str, GeneFeature]:
    """One CDS data row -> (seqid, gene); typed errors carry the line."""
    seqid, _source, _type, start_text, end_text, _score, strand, _phase, attrs_text = row
    try:
        start, end = int(start_text), int(end_text)
    except ValueError as exc:
        raise _fail(
            "INVALID_GFF",
            f"{path}: non-integer CDS coordinates at line {line_no}: {start_text}..{end_text}",
            line=str(line_no),
        ) from exc
    if strand not in ("+", "-"):
        raise _fail(
            "INVALID_GFF",
            f"{path}: CDS at line {line_no} lacks a +/- strand",
            line=str(line_no),
        )
    attrs = _attributes(attrs_text)
    gene_id = attrs.get("gene") or attrs.get("Name") or attrs.get("locus_tag") or ""
    if not gene_id:
        raise _fail(
            "GENE_WITHOUT_NAME",
            f"{path}: CDS at line {line_no} has no gene/Name/locus_tag attribute",
            line=str(line_no),
        )
    gene = GeneFeature(
        gene_id=gene_id, start=start, end=end, strand=strand, product=attrs.get("product", "")
    )
    return seqid, gene


def _embedded_fasta(handle: IO[str], path: Path) -> list[tuple[str, str]]:
    """(id, sequence) pairs from the lines after ``##FASTA``."""
    records: list[tuple[str, str]] = []
    seqid: str | None = None
    chunks: list[str] = []
    for raw in handle:
        line = raw.rstrip("\n")
        if line.startswith(">"):
            if seqid is not None:
                records.append((seqid, "".join(chunks)))
            parts = line[1:].split(None, 1)
            seqid, chunks = (parts[0] if parts else ""), []
        elif line.strip() and seqid is None:
            raise _fail("INVALID_GFF", f"{path}: FASTA content before the first '>' after ##FASTA")
        elif line.strip():
            chunks.append(line.strip())
    if seqid is not None:
        records.append((seqid, "".join(chunks)))
    return records


def _sidecar_candidates(path: Path) -> list[Path]:
    """``<stem>.fa/.fna/.fasta`` (each plain and ``.gz``) for a GFF path; the
    stem drops the .gz wrapper first, then the .gff3/.gff suffix."""
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[: -len(".gz")]
    stem = name
    for suffix in _GFF_SUFFIXES:
        if name.lower().endswith(suffix):
            stem = name[: -len(suffix)]
            break
    return [
        path.parent / f"{stem}{suffix}{gz}" for suffix in _SIDECAR_SUFFIXES for gz in ("", ".gz")
    ]


def _sidecar_records(path: Path) -> list[tuple[str, str]] | None:
    """First existing sidecar's records, or None (caller lists the paths)."""
    for candidate in _sidecar_candidates(path):
        if candidate.is_file():
            return [(record.id, record.sequence) for record in iter_fasta(candidate)]
    return None


def _build_loci(
    records: list[tuple[str, str]], genes: dict[str, list[GeneFeature]], path: Path
) -> list[LocusSequence]:
    """FASTA records + CDS rows -> loci; enforces the shared typed rules."""
    loci: list[LocusSequence] = []
    seen: set[str] = set()
    for seqid, sequence in records:
        if seqid in seen:
            raise _fail("INVALID_GFF", f"{path}: duplicate FASTA sequence {seqid!r}", seqid=seqid)
        seen.add(seqid)
        upper = sequence.upper()
        locus_genes = genes.pop(seqid, None)
        if not locus_genes:
            raise _fail(
                "LOCUS_WITHOUT_GENES",
                f"{path}: locus {seqid!r} carries no CDS features",
                locus=seqid,
            )
        for gene in locus_genes:
            if not 1 <= gene.start <= gene.end <= len(upper):
                raise _fail(
                    "GENE_OUT_OF_BOUNDS",
                    f"{path}: locus {seqid!r} gene {gene.gene_id!r} is outside the locus",
                    locus=seqid,
                    gene=gene.gene_id,
                )
        loci.append(
            LocusSequence(
                features=LocusFeatures(id=seqid, label=seqid, type=seqid, genes=tuple(locus_genes)),
                sequence=upper,
            )
        )
    if genes:
        orphan = next(iter(genes))
        raise _fail(
            "INVALID_GFF",
            f"{path}: CDS rows for sequence {orphan!r} missing from the FASTA",
            seqid=orphan,
        )
    if not loci:
        raise _fail("INVALID_GFF", f"{path}: no FASTA sequences found")
    return loci


def parse_gff3_features(path: Path) -> list[LocusSequence]:
    """Parse a (plain/.gz/.bz2) GFF3 file into one LocusSequence per FASTA
    sequence; every documented malformation raises a typed InputError."""
    genes: dict[str, list[GeneFeature]] = {}
    embedded: list[tuple[str, str]] | None = None
    try:
        with open_text(path) as handle:
            first = handle.readline().rstrip("\r\n")
            if first != "##gff-version 3":
                raise _fail(
                    "INVALID_GFF",
                    f"{path}: first line must be '##gff-version 3', found: {first!r}",
                )
            for line_no, raw in enumerate(handle, start=2):
                line = raw.rstrip("\n")
                if line.rstrip("\r") == "##FASTA":
                    embedded = _embedded_fasta(handle, path)
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                row = line.split("\t")
                if len(row) != 9:
                    raise _fail(
                        "INVALID_GFF",
                        f"{path}: expected 9 tab-separated columns at line {line_no}",
                        line=str(line_no),
                    )
                if row[2] != "CDS":
                    continue
                seqid, gene = _gene(row, path, line_no)
                genes.setdefault(seqid, []).append(gene)
    except (OSError, UnicodeDecodeError, EOFError) as exc:
        raise _fail("INVALID_GFF", f"{path}: could not read input: {exc}") from exc
    records = embedded if embedded is not None else _sidecar_records(path)
    if records is None:
        searched = ";".join(str(candidate) for candidate in _sidecar_candidates(path))
        raise _fail(
            "GFF_WITHOUT_FASTA",
            f"{path}: no ##FASTA section and no sidecar FASTA (searched: {searched})",
            searched=searched,
        )
    return _build_loci(records, genes, path)
