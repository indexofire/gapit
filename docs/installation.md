# Installing gapit

Two install paths: the PyPI wheel (you provide the external binaries) or a git clone with
[pixi](https://pixi.sh) (binaries managed for you).

## From PyPI

```bash
pip install gapit
```

The wheel bundles the Python package and the `gapit`/`gapit-mcp` console scripts — no
database content (every provider downloads from upstream at fetch time; several upstream
licenses forbid redistribution). It does **not** bundle BLAST+ or minimap2 either — install
them first ([External binaries](#external-binaries), e.g. `conda create -n gapit-env -c
bioconda blast minimap2`).

## From source (pixi)

### Prerequisites

- [git](https://git-scm.com) and [pixi](https://pixi.sh). On macOS/Linux:
  `curl -fsSL https://pixi.sh/install.sh | bash`
- Python 3.11+ if you install the package outside pixi (pip/pyproject). With pixi this is
  moot: the environment ships its own Python (the dev env pins 3.14).

### Install

```bash
git clone https://github.com/indexofire/gapit.git
cd gapit
pixi install
```

Run the CLI through pixi, or put the env's bin directory on PATH:

```bash
pixi run gapit --version
# or, equivalent:
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
gapit --version
```

## Verify

```console
$ gapit --version
gapit 0.5.2
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.5.2"}
```

## External binaries

gapit shells out to BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`.
Both come from the pixi environment (conda-forge and bioconda), so there is nothing to install
by hand. If you run gapit outside pixi, make sure these binaries are on PATH yourself. Input
normalization (plain/gz/bz2 FASTA, FASTQ, GenBank, EMBL) is native — no `any2fasta` needed.

## Database bootstrap

Screening needs at least one database in the datadir (`$GAPIT_DATADIR`, then
`~/.local/share/gapit/db`; override per call with `--datadir`). Nothing is bundled:
every provider downloads from upstream when fetched (network required), because several
upstream licenses — CARD's McMaster terms, VFDB's CC BY-NC, Kaptive's GPL-3.0 — forbid
redistribution inside an MIT-licensed package. See [Databases](./databases.md) for the
full database table.

```bash
gapit db fetch all         # installs the default set (card, vfdb) into the default datadir
```

The same install into a scratch datadir, so you can see what a fetch prints. Progress
lines go to stderr, one JSON receipt per database to stdout:

```console
$ gapit db fetch all --datadir /tmp/gapit-docs/dd
gapit: downloaded 1 source file(s)
gapit: read 6059 records from card
gapit: generated /tmp/gapit-docs/dd/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/card"}
gapit: downloaded 1 source file(s)
gapit: read 4769 records from vfdb
gapit: generated /tmp/gapit-docs/dd/vfdb/sequences
gapit: self-check passed for vfdb
gapit: BLAST index built (nucl)
{"db":"vfdb","records":4769,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/vfdb"}
```

Confirm the database catalog and what's installed any time with `gapit db list` (or
`gapit db list --json`).

## Shell completions

```bash
gapit --install-completion   # bash, zsh, or fish; installs for the current shell
gapit --show-completion      # print the completion script to copy or customize
```

## Upgrading

```bash
git pull
pixi install
```

Database content does not update with the code: re-run `gapit db fetch <name> --force`
to rebuild an installed database from upstream.
