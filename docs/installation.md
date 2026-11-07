# Installing gspec

Three install paths: the PyPI wheel (you provide the external binaries), the repository's
[`environment.yml`](https://github.com/indexofire/gspec/blob/main/environment.yml) (conda
installs the binaries for you), or a git clone with [pixi](https://pixi.sh) (binaries
managed for you).

## From PyPI

```bash
pip install gspec
```

The wheel bundles the Python package and the `gspec`/`gspec-mcp` console scripts, plus
seven audited, permissively licensed, gspec-curated database bundles that materialize into
the datadir on first use (see
[Bundled databases](./databases.md#bundled-databases-install-time-ready)); every other
provider downloads from upstream at fetch time, because several upstream licenses forbid
redistribution. It does **not** bundle BLAST+ or minimap2 either — install
them first ([External binaries](#external-binaries), e.g. `conda create -n gspec-env -c
conda-forge -c bioconda blast minimap2`).

## From conda (environment.yml)

If you want conda to manage the external binaries too, use the repository's
[`environment.yml`](https://github.com/indexofire/gspec/blob/main/environment.yml). It pins
Python 3.12, installs BLAST+ and minimap2 from conda-forge/bioconda, and installs gspec
from PyPI in the same environment:

```bash
curl -fsSLO https://raw.githubusercontent.com/indexofire/gspec/main/environment.yml
conda env create -f environment.yml
conda activate gspec-env
```

gspec is not on bioconda yet (recipe submission is deferred), so it arrives through the
file's `pip:` subsection; BLAST+ and minimap2 have no PyPI distribution, which is why they
must be conda dependencies.

## From source (pixi)

### Prerequisites

- [git](https://git-scm.com) and [pixi](https://pixi.sh). On macOS/Linux:
  `curl -fsSL https://pixi.sh/install.sh | bash`
- Python 3.12 if you install the package outside pixi (pip/pyproject). With pixi this is
  moot: the environment ships its own Python (the dev env pins 3.12).

### Install

```bash
git clone https://github.com/indexofire/gspec.git
cd gspec
pixi install
```

Run the CLI through pixi, or put the env's bin directory on PATH:

```bash
pixi run gspec --version
# or, equivalent:
export PATH="$PWD/.pixi/envs/default/bin:$PATH"
gspec --version
```

## Verify

```console
$ gspec --version
gspec 0.10.0
$ gspec --version --json
{"schema":"gspec.version/1","name":"gspec","version":"0.10.0"}
```

## Optional gseek engine

On Python >= 3.12, `pip install "gspec[gseek]"` enables
`GSPEC_ENGINE=gseek|shadow`; the pixi `default` environment (which enables the `py312` feature) includes
gseek; `py312` is the 3.12 baseline feature, not a separate environment. See [gseek engine](gseek-engine.md).

## External binaries

gspec shells out to BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`.
Both come from the pixi environment or `environment.yml` (conda-forge and bioconda), so
there is nothing to install by hand on those paths. If you installed gspec with plain
`pip`, make sure these binaries are on PATH yourself. Input normalization (plain/gz/bz2
FASTA, FASTQ, GenBank, EMBL) is native — no `any2fasta` needed.

## Database bootstrap

Screening needs at least one database in the datadir (`$GSPEC_DATADIR`, then
`~/.local/share/gspec/db`; override per call with `--datadir`). On a fresh install,
`gspec db path` shows the effective location, `gspec db list` works before any database
exists (it bootstraps the datadir and lists the catalog), and `gspec db setup` creates the
datadir and materializes the bundled databases.
Four audited, permissively licensed databases ship inside the wheel and need no network:
the first screen against a bundled name materializes it automatically, and `gspec db setup`
materializes all four at once. Every other provider downloads from upstream when fetched
(network required), because several
upstream licenses — CARD's McMaster terms, VFDB's CC BY-NC, Kaptive's GPL-3.0 — forbid
redistribution inside an MIT-licensed package. See [Databases](./databases.md) for the
full database table and the [bundled section](./databases.md#bundled-databases-install-time-ready).

To keep databases inside an activated conda environment instead of your home directory,
set `GSPEC_DATADIR="${CONDA_PREFIX}/db/gspec"` (the home default survives environment
rebuilds and is shared by every env).

```bash
gspec db fetch -d all         # installs the default set (card, vfdb) into the default datadir
```

The same install into a scratch datadir, so you can see what a fetch prints. Progress
lines go to stderr, one JSON receipt per database to stdout:

```console
$ gspec db fetch -d all --datadir /tmp/gspec-docs/dd
gspec: downloaded 1 source file(s)
gspec: read 6059 records from card
gspec: generated /tmp/gspec-docs/dd/card/sequences
gspec: self-check passed for card
gspec: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/gspec-docs/dd/card"}
gspec: downloaded 1 source file(s)
gspec: read 4769 records from vfdb
gspec: generated /tmp/gspec-docs/dd/vfdb/sequences
gspec: self-check passed for vfdb
gspec: BLAST index built (nucl)
{"db":"vfdb","records":4769,"dbtype":"nucl","destination":"/tmp/gspec-docs/dd/vfdb"}
```

Confirm the database catalog and what's installed any time with `gspec db list` (or
`gspec db list --json`).

## Shell completions

The CLI is built with click: completions activate through the `_GSPEC_COMPLETE`
environment variable — there is no `--install-completion` flag. Pick your shell:

Bash — in `~/.bashrc`:

```bash
eval "$(_GSPEC_COMPLETE=bash_source gspec)"
```

Zsh — in `~/.zshrc` (after `compinit`), or drop the emitted `#compdef` script on your
`fpath` as `~/.zfunc/_gspec`:

```zsh
eval "$(_GSPEC_COMPLETE=zsh_source gspec)"
```

Fish:

```fish
_GSPEC_COMPLETE=fish_source gspec > ~/.config/fish/completions/gspec.fish
```

The emitted script (on bash it ends with `complete -o nosort -F _gspec_completion gspec`)
calls back into gspec at completion time — `_GSPEC_COMPLETE=bash_complete gspec` plus the
words under the cursor — so completions always follow the installed version. click also
ships a powershell emitter (`_GSPEC_COMPLETE=powershell_source gspec`).

## Upgrading

From a git clone with pixi:

```bash
git pull
pixi install
```

With `environment.yml`, upgrade gspec from PyPI inside the activated environment and refresh
the conda-managed binaries as needed:

```bash
conda activate gspec-env
pip install --upgrade gspec
conda update blast minimap2
```

Database content does not update with the code: re-run `gspec db fetch -d <name> --force`
to rebuild an installed database from upstream.
