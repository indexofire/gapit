# Summarizing reports

`gapit summary` collapses one or more abricate-format report tables (including gapit screen
output) into a gene presence/absence matrix: one row per file, one column per gene, cells
carrying the hit metric. It replaces `abricate --summary`.

Input reports are TSV or CSV with the standard 15-column header (`#FILE  SEQUENCE  ...`),
exactly what `gapit screen` writes. See [./screen.md](./screen.md) for producing them and
[./outputs.md](./outputs.md) for the output schemas.

## Options

Transcribed from `gapit summary --help` (gapit 0.1.0):

| Flag | Type | Default | Description |
|---|---|---|---|
| `FILE...` | path(s) | required | Abricate-format report file(s) to summarize. At least one. |
| `--identity` | flag | off | Cells show %IDENTITY instead of %COVERAGE. |
| `--nopath` | flag | off | Basename row keys (FILE values / input filenames). |
| `--quiet` | flag | off | Silence stderr diagnostics. |
| `--format` | tsv\|csv\|json\|md | `tsv` | Output format. |

## Dutch mode vs multi-file mode

Row keys depend on the input count, matching abricate's behavior:

- **Exactly one report (dutch mode):** rows are keyed by that report's `FILE` column, so you get
  one row per assembly inside the report.
- **More than one report:** rows are keyed by the input filenames as given, one row per file.
  A file with no hits still appears, with `NUM_FOUND 0`.

```console
$ gapit summary tests/data/summary/multi_sample.tsv
#FILE	NUM_FOUND	feature_a	feature_b
aa_assembly.fa	1	99.00	.
mm_assembly.fa	1	.	50.00
zz_assembly.fa	1	91.00	.
```

The same three assemblies summarized as separate reports instead would key by filename:

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv tests/data/summary/empty.tsv
#FILE	NUM_FOUND	feature_a	feature_b
tests/data/summary/empty.tsv	0	.	.
tests/data/summary/sample_a.tsv	2	99.50;52.00	76.00
tests/data/summary/sample_b.tsv	2	90.00	100.00
```

Note the zero-hit file: it appears with `NUM_FOUND 0` and `.` in every gene column. Multiple
hits on the same gene in one file join with `;` in report order (see `99.50;52.00`).

## Cells and `--identity`

By default cells hold each hit's %COVERAGE. With `--identity` they hold %IDENTITY instead; gapit
notes the switch on stderr:

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv --identity --nopath
Using %IDENTITY for the summary table instead of %COVERAGE
#FILE	NUM_FOUND	feature_a	feature_b
sample_a.tsv	2	98.75;91.00	95.10
sample_b.tsv	2	97.00	99.99
```

Cells keep the original report strings verbatim; no reformatting or averaging happens.
`NUM_FOUND` counts distinct genes. The gene universe is the union of all `GENE` values across
all inputs, sorted lexicographically; a gene absent from a row shows `.`.

## Input parsing

- **Separator auto-detection, per file.** Each input is read as TSV or CSV by sniffing its own
  first line, so you can mix tab and comma reports in one call without flags:

```console
$ gapit summary tests/data/summary/sample_a.csv tests/data/summary/sample_b.csv --nopath
#FILE	NUM_FOUND	feature_a	feature_b
sample_a.csv	2	99.50;52.00	76.00
sample_b.csv	2	90.00	100.00
```

- **Duplicate inputs are skipped.** A path listed twice (compared as given, before any
  basenaming) warns on stderr and is processed once:

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_a.tsv
WARNING: Skipping duplicate file: tests/data/summary/sample_a.tsv
#FILE	NUM_FOUND	feature_a	feature_b
tests/data/summary/sample_a.tsv	2	99.50;52.00	76.00
```

- **Malformed input is a typed error, not silence.** A missing file exits 5 with an
  `INPUT_NOT_FOUND` envelope; a row too short for the header map exits 5 with
  `SUMMARY_MALFORMED`.
- **Header handling matches abricate's parser.** The first row of the first file becomes the
  column-name map. Lines whose first column starts with `#` are then skipped as headers, so a
  conventional report contributes data rows only. Two quirks fall out of that rule: a noheader
  report's first row doubles as header map and data row, and a file whose key starts with `#`
  is skipped entirely (the check runs before `--nopath` basenaming).

## Output formats

`--format tsv` (default) and `--format csv` print the matrix above. `--format json` and
`--format md` emit the `gapit.summary/1` document, which keeps every cell string and adds
machine-readable params:

```console
$ gapit summary tests/data/summary/sample_a.tsv tests/data/summary/sample_b.tsv --format json --quiet
{
  "schema": "gapit.summary/1",
  "tool": {
    "name": "gapit",
    "version": "0.1.0"
  },
  "created_at": "2026-09-19T01:12:48Z",
  "params": {
    "metric": "%COVERAGE",
    "nopath": false
  },
  "genes": [
    "feature_a",
    "feature_b"
  ],
  "rows": [
    {
      "file": "tests/data/summary/sample_a.tsv",
      "num_found": 2,
      "cells": {
        "feature_a": [
          "99.50",
          "52.00"
        ],
        "feature_b": [
          "76.00"
        ]
      }
    },
    {
      "file": "tests/data/summary/sample_b.tsv",
      "num_found": 2,
      "cells": {
        "feature_a": [
          "90.00"
        ],
        "feature_b": [
          "100.00"
        ]
      }
    }
  ]
}
```

Introspect the schema with `gapit schema summary`; field contracts live in
[./outputs.md](./outputs.md).

## Notes

- **Sort order is by key, not label.** Rows sort by the row key exactly as given, before
  `--nopath` basenaming. With files spread over directories the basenamed labels can therefore
  appear unsorted:

```console
$ gapit summary /tmp/gapit-demo/sortdemo/1dir/zeta.tsv /tmp/gapit-demo/sortdemo/2dir/mid.tsv --nopath
#FILE	NUM_FOUND	feature_a	feature_b
zeta.tsv	2	91.00;99.00	50.00
mid.tsv	2	99.50;52.00	76.00
```

  The order comes from `1dir/zeta.tsv` < `2dir/mid.tsv`, so `zeta.tsv` prints first even though
  `mid.tsv` sorts earlier alphabetically. Without `--nopath` the full keys print in obvious
  order:

```console
$ gapit summary /tmp/gapit-demo/sortdemo/1dir/zeta.tsv /tmp/gapit-demo/sortdemo/2dir/mid.tsv
#FILE	NUM_FOUND	feature_a	feature_b
/tmp/gapit-demo/sortdemo/1dir/zeta.tsv	2	91.00;99.00	50.00
/tmp/gapit-demo/sortdemo/2dir/mid.tsv	2	99.50;52.00	76.00
```

- **Dutch mode keys come from the report.** With one input, the row labels are the `FILE` values
  found inside the report, not the report's own filename.
- **`--quiet` suppresses the `--identity` note and duplicate warnings**, never stdout.
