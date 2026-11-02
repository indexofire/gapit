# Vibrio database builders

Developer scripts behind the bundled `vpa_pai` and `vch_vf` databases.

## `vpa_pai` (cluster)

`build_vpa_pai.py` reconstructs VPaI-1..VPaI-10 as ten cluster loci from public
INSDC records:

- RIMD 2210633 (`NC_004603.1`, `NC_004605.1`) for VPaI-1..VPaI-7 and the T6SS
  VPaI-10;
- AQ3810 (`GCA_000154045.1`, scaffolds `DS180230` and `DS180088`) for VPaI-8
  and VPaI-9.

The island nomenclature and boundaries come from Boyd et al. 2008 (BMC
Microbiol 8:110, PMID 18590559) and Hurley et al. 2006 (BMC Genomics 7:104).
VPaI-10 is the RIMD T6SS region VP1387-VP1420; see the `bundled.json` caveat.

```bash
python scripts/vibrio/build_vpa_pai.py \
  --rimd GCF_000196095.1_genomic.gbff \
  --aq3810 GCA_000154045.1_genomic.gbff \
  --dest src/gapit/data/dbs/vpa_pai
```

## `vch_vf` (gene)

`build_vch_vf.py` cuts 13 screen-only markers from public CDS FASTA products:

- CTXphi `NC_015209.1`: ctxA, classical ctxB, zot, ace, cep, orfU, rstR,
  rstA, rstB;
- N16961 `NC_002505.1`: El Tor ctxB and tcpA;
- O395 `NC_012582.1`: classical tcpA;
- O139 FC1225 `NZ_NGQY01000008.1`: O139-associated tcpA variant.

```bash
python scripts/vibrio/build_vch_vf.py \
  --ctxphi NC_015209.1_cds.fna --n16961 NC_002505.1_cds.fna \
  --o395 NC_012582.1_cds.fna --o139 NZ_NGQY01000008.1_cds.fna \
  --dest src/gapit/data/dbs/vch_vf
```

Both builders emit only public-domain INSDC content; no GPL or non-commercial
database files are copied.
