# Shigella serotype database (branch prototype)

This branch starts the gapit-native Shigella molecular serotype panel.

## Current slice

- `src/gapit/data/dbs/shi_sero/`
  - `sequences`: 95 ShigaTyper marker sequences + one split `ipaH_c` control
    record (96 FASTA records; gapit/v1 headers)
  - `markers.json`: marker → category / serotype-group metadata
  - `serotype.json`: first `gapit.serotype/1` rule set
  - `bundled.json`: bundle metadata + provenance
  - `screen.json`: contig/reads marker thresholds
- `scripts/shigella/build_shigella_sero.py`: deterministic builder from
  ShigaTyper 2.0.5 `ShigellaRef5.fasta`.

The bundle now carries `serotype.json` (`gapit.serotype/1`) and a first
conservative marker-set resolver: `gapit typing` reads the screen table and
returns one `serotype` call. It covers the direct Sd/Sb wzx markers, S.
sonnei forms I/II, Sf6, the S. flexneri POAC combinations, and Sd1 rfp+/rfp-.
S. boydii 6 vs 10 uses the `wbaM` 5' alignment start (IS629 interruption vs
intact type 10). The `ipaH_c` control is split into two non-overlapping
presence fragments (1-640 and 641-780) so divergent ipaH alleles remain
detectable; reads tables carry no alignment geometry and leave the boydii
6/10 pair unresolved.

## Licensing

- ShigaTyper marker FASTA: **Public Domain** (package metadata).
- ShigaPass is **GPL3** and its database content is deliberately **not copied**
  into this MIT bundle. ShigaPass is used only as a method/reference for future
  clean-room curation.
- If the panel grows beyond the ShigaTyper public-domain set, new sequences must
  be re-extracted from public-domain INSDC/GenBank/RefSeq records and curated
  here, with provenance recorded in `bundled.json`.

## Next slices

1. ~~Add `gapit.serotype/1`/profile schema and the first ShigaTyper marker
   resolver.~~ Done in the first slice.
2. Expand/validate the ShigaTyper decision-tree port: 144-genome S. boydii
   ShigaTyper/ShigaPass validation plus 669 S. flexneri/S. sonnei/S.
   dysenteriae assemblies against ShigeiFinder (see `VALIDATION.md`).
3. Add a clean-room ShigaPass-style profile layer if needed:
   RFB 150-mers + fliC + CRISPR + POAC; MLST optional/fallback only.
4. Evaluate gseek/KMA as an optional marker backend for reads/low-quality
   assemblies; the assembly path can stay on BLASTn.
5. Carry the remaining known gaps: 2av/2a refinement, mixed Sf/Ss sonnei
   ambiguity, `cadA`-negative form II, and reads-level `wbaM` junction.

## Non-goals

- Copying GPL ShigaPass database files into the wheel.
- Calling read-level `wbaM` junction/variant accuracy without the required
  read-level engine.
- Pretending marker presence alone is a complete serotype call.
