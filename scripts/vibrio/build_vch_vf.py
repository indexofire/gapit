#!/usr/bin/env python3
"""Build the bundled V. cholerae CTX/tcpA virulence marker database.

``vch_vf`` is a gene database: CTX prophage core genes, the RS2 regulatory
module, allele-specific ctxB markers, and tcpA allele markers.  The source
sequences are public INSDC records, cut from their CDS FASTA products:

- CTXphi NC_015209.1 (classical ctxB; complete phage annotation);
- V. cholerae N16961 NC_002505.1 (El Tor ctxB/rstR/tcpA);
- V. cholerae O395 NC_012582.1 (classical tcpA);
- V. cholerae O139 FC1225 NZ_NGQY01000008.1 (O139-associated tcpA variant).

Usage::

    python scripts/vibrio/build_vch_vf.py \
        --ctxphi NC_015209.1_cds.fna --n16961 NC_002505.1_cds.fna \
        --o395 NC_012582.1_cds.fna --o139 NZ_NGQY01000008.1_cds.fna \
        --dest src/gapit/data/dbs/vch_vf
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import UTC, datetime
from pathlib import Path

from gapit.dbcodec import encode_seqid

#: gene name -> (source label, header regex, product)
MARKERS: tuple[tuple[str, str, str, str], ...] = (
    (
        "ctxA",
        "ctxphi",
        r"\[gene=ctxA\]",
        "cholera enterotoxin catalytic subunit CtxA",
    ),
    (
        "ctxB_classical",
        "ctxphi",
        r"\[gene=ctxB\]",
        "cholera enterotoxin binding subunit CtxB (classical allele)",
    ),
    (
        "ctxB_eltor",
        "n16961",
        r"\[gene=ctxB\]",
        "cholera enterotoxin binding subunit CtxB (El Tor allele)",
    ),
    ("zot", "ctxphi", r"\[gene=zot\]", "zonula occludens toxin Zot"),
    ("ace", "ctxphi", r"\[gene=ace\]", "accessory cholera enterotoxin Ace"),
    ("cep", "ctxphi", r"\[gene=cep\]", "core-encoded pilin Cep"),
    (
        "orfU",
        "ctxphi",
        r"\[locus_tag=AV950_gp09\]",
        "minor coat protein pIII (OrfU)",
    ),
    (
        "rstR",
        "ctxphi",
        r"\[gene=rstR\]",
        "CTXphi transcriptional repressor RstR",
    ),
    ("rstA", "ctxphi", r"\[gene=rstA\]", "CTXphi replication initiation factor RstA"),
    ("rstB", "ctxphi", r"\[gene=rstB\]", "CTXphi RstB"),
    (
        "tcpA_classical",
        "o395",
        r"\[locus_tag=VC395_RS04335\]",
        "toxin-coregulated pilus major pilin TcpA (classical allele)",
    ),
    (
        "tcpA_eltor",
        "n16961",
        r"\[gene=tcpA\]",
        "toxin-coregulated pilus major pilin TcpA (El Tor allele)",
    ),
    (
        "tcpA_o139",
        "o139",
        r"\[gene=tcpA\]",
        "toxin-coregulated pilus major pilin TcpA (O139-associated variant)",
    ),
)


def _records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header = line[1:]
            sequence = []
        else:
            sequence.append(line.strip())
    if header:
        records.append((header, "".join(sequence)))
    return records


def _pick(records: list[tuple[str, str]], pattern: str, marker: str) -> tuple[str, str]:
    match = re.compile(pattern)
    for header, sequence in records:
        if match.search(header):
            return header, sequence
    raise SystemExit(f"no CDS record for {marker}: /{pattern}/")


def _accession(header: str) -> str:
    location = re.search(r"\[location=([^\]]+)\]", header)
    protein = re.search(r"\[protein_id=([^\]]+)\]", header)
    return protein.group(1) if protein else (location.group(1) if location else "")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("ctxphi", "n16961", "o395", "o139"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--dest", type=Path, required=True)
    args = parser.parse_args(argv)

    sources = {
        "ctxphi": _records(args.ctxphi),
        "n16961": _records(args.n16961),
        "o395": _records(args.o395),
        "o139": _records(args.o139),
    }
    dest: Path = args.dest
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "sequences").open("w", encoding="utf-8", newline="\n") as handle:
        for gene, source, pattern, product in MARKERS:
            header, sequence = _pick(sources[source], pattern, gene)
            seqid = encode_seqid(
                "vch_vf",
                gene,
                _accession(header),
                ("virulence",),
            )
            handle.write(f">{seqid} {product}\n")
            for offset in range(0, len(sequence), 60):
                handle.write(sequence[offset : offset + 60] + "\n")
    (dest / "bundled.json").write_text(
        json.dumps(
            {
                "name": "vch_vf",
                "description": ("V. cholerae CTX prophage core/RS2 genes + ctxB/tcpA allele panel"),
                "vendor": "gspec-curated (INSDC)",
                "dbtype": "nucl",
                "kind": "gene",
                "snapshotted": datetime.now(UTC).strftime("%Y-%m-%d"),
                "license": "Public domain (INSDC sequences)",
                "note": (
                    "CTX core (ctxA, ctxB, zot, ace, cep, orfU), RS2"
                    " (rstR, rstA, rstB), allele-specific ctxB classical/El"
                    " Tor, and tcpA classical/El Tor/O139 markers cut from"
                    " NC_015209.1 (CTXphi), N16961 (NC_002505.1), O395"
                    " (NC_012582.1), and O139 FC1225"
                    " (NZ_NGQY01000008.1). No typing scheme."
                ),
                "guidance": {
                    "references": [
                        {
                            "citation": (
                                "Waldor MK, Mekalanos JJ. Lysogenic conversion"
                                " by a filamentous phage encoding cholera toxin."
                                " Science. 1996;272(5270):1910-1914."
                            ),
                            "doi": "10.1126/science.272.5270.1910",
                            "note": "CTXphi biology and ctxAB carriage.",
                        },
                        {
                            "citation": (
                                "Taylor RK, Miller VL, Furlong DB, Mekalanos JJ."
                                " Use of phoA gene fusions to identify a pilus"
                                " colonization factor coordinately regulated with"
                                " cholera toxin. Proc Natl Acad Sci USA."
                                " 1987;84(9):2833-2837."
                            ),
                            "doi": "10.1073/pnas.84.9.2833",
                            "note": "tcpA/TCP as the CTXphi receptor and colonization factor.",
                        },
                        {
                            "citation": (
                                "Li X, Han Y, Zhao W, et al. Diversity and"
                                " Complexity of CTXPhi and Pre-CTXPhi Families in"
                                " Vibrio cholerae from Seventh Pandemic."
                                " Microorganisms. 2024;12(10):1935."
                            ),
                            "doi": "10.3390/microorganisms12101935",
                            "note": "ctxB and tcpA allele/type framework used by this panel.",
                        },
                    ],
                    "caveats": (
                        "Marker presence is genetic potential, not toxin"
                        " expression: CTXphi can be partial, defective, or"
                        " lost, and tcpA presence does not prove a functional"
                        " TCP.",
                        "ctxB_classical and ctxB_eltor are allele markers;"
                        " O139 strains use the El Tor ctxB allele, so there is"
                        " deliberately no separate ctxB_o139 record.",
                        "tcpA_o139 is the O139-associated variant represented"
                        " by strain FC1225 (WP_050926100.1); it is ~99.7%"
                        " identical to the canonical El Tor tcpA and can"
                        " cross-hit. Read identity/coverage, not just presence.",
                        "cep is the core-encoded pilin; orfU is carried as the"
                        " NC_015209.1 minor coat protein pIII annotation. Both"
                        " are CTXphi structural/assembly genes, not toxin genes.",
                        "rstR/rstA/rstB are RS2 phage regulatory/replication"
                        " genes; the RstR sequence is shared by the reference"
                        " CTXphi and N16961 records, so it is a single marker."
                        " RS1 (rstC) and CTXphi integration-site context are"
                        " not modelled here.",
                        "This panel does not assign an O serogroup, biotype, or"
                        " pathotype; combine it with the O-antigen and typing"
                        " panels.",
                    ),
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (dest / "screen.json").write_text(
        json.dumps(
            {
                "schema": "gapit.screen/1",
                "contigs": {"minid": 80.0, "mincov": 50.0},
                "reads": {"min_breadth": 50.0, "min_identity": 80.0},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {dest}: {len(MARKERS)} markers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
