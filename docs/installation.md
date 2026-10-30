# Installing gapit

Three install paths: the PyPI wheel (you provide the external binaries), the repository's
[`environment.yml`](https://github.com/indexofire/gapit/blob/main/environment.yml) (conda
installs the binaries for you), or a git clone with [pixi](https://pixi.sh) (binaries
managed for you).

## From PyPI

```bash
pip install gapit
```

The wheel bundles the Python package and the `gapit`/`gapit-mcp` console scripts, plus
seven audited, permissively licensed database bundles that materialize into the datadir on
first use (see
[Bundled databases](./databases.md#bundled-databases-install-time-ready)); every other
provider downloads from upstream at fetch time, because several upstream licenses forbid
redistribution. It does **not** bundle BLAST+ or minimap2 either — install
them first ([External binaries](#external-binaries), e.g. `conda create -n gapit-env -c
conda-forge -c bioconda blast minimap2`).

## From conda (environment.yml)

If you want conda to manage the external binaries too, use the repository's
[`environment.yml`](https://github.com/indexofire/gapit/blob/main/environment.yml). It pins
Python 3.12, installs BLAST+ and minimap2 from conda-forge/bioconda, and installs gapit
from PyPI in the same environment:

```bash
curl -fsSLO https://raw.githubusercontent.com/indexofire/gapit/main/environment.yml
conda env create -f environment.yml
conda activate gapit-env
```

gapit is not on bioconda yet (recipe submission is deferred), so it arrives through the
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
gapit 0.6.6
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.6.6"}
```

## Optional gseek engine

On Python >= 3.12, `pip install "gapit[gseek]"` enables
`GAPIT_ENGINE=gseek|shadow`; the pixi `default` environment (which enables the `py312` feature) includes
gseek; `py312` is the 3.12 baseline feature, not a separate environment. See [gseek engine](gseek-engine.md).

## External binaries

gapit shells out to BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`) and `minimap2`.
Both come from the pixi environment or `environment.yml` (conda-forge and bioconda), so
there is nothing to install by hand on those paths. If you installed gapit with plain
`pip`, make sure these binaries are on PATH yourself. Input normalization (plain/gz/bz2
FASTA, FASTQ, GenBank, EMBL) is native — no `any2fasta` needed.

## Database bootstrap

Screening needs at least one database in the datadir (`$GAPIT_DATADIR`, then
`~/.local/share/gapit/db`; override per call with `--datadir`). On a fresh install,
`gapit db path` shows the effective location, `gapit db list` works before any database
exists (it bootstraps the datadir and lists the catalog), and `gapit db setup` creates the
datadir and materializes the bundled databases.
Seven audited, permissively licensed databases ship inside the wheel and need no network:
the first screen against a bundled name materializes it automatically, and `gapit db setup`
materializes all seven at once. Every other provider downloads from upstream when fetched
(network required), because several
upstream licenses — CARD's McMaster terms, VFDB's CC BY-NC, Kaptive's GPL-3.0 — forbid
redistribution inside an MIT-licensed package. See [Databases](./databases.md) for the
full database table and the [bundled section](./databases.md#bundled-databases-install-time-ready).

To keep databases inside an activated conda environment instead of your home directory,
set `GAPIT_DATADIR="${CONDA_PREFIX}/db/gapit"` (the home default survives environment
rebuilds and is shared by every env).

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

The CLI is built with click: completions activate through the `_GAPIT_COMPLETE`
environment variable — there is no `--install-completion` flag. Pick your shell:

Bash — in `~/.bashrc`:

```bash
eval "$(_GAPIT_COMPLETE=bash_source gapit)"
```

Zsh — in `~/.zshrc` (after `compinit`), or drop the emitted `#compdef` script on your
`fpath` as `~/.zfunc/_gapit`:

```zsh
eval "$(_GAPIT_COMPLETE=zsh_source gapit)"
```

Fish:

```fish
_GAPIT_COMPLETE=fish_source gapit > ~/.config/fish/completions/gapit.fish
```

The emitted script (on bash it ends with `complete -o nosort -F _gapit_completion gapit`)
calls back into gapit at completion time — `_GAPIT_COMPLETE=bash_complete gapit` plus the
words under the cursor — so completions always follow the installed version. click also
ships a powershell emitter (`_GAPIT_COMPLETE=powershell_source gapit`).

## Upgrading

From a git clone with pixi:

```bash
git pull
pixi install
```

With `environment.yml`, upgrade gapit from PyPI inside the activated environment and refresh
the conda-managed binaries as needed:

```bash
conda activate gapit-env
pip install --upgrade gapit
conda update blast minimap2
```

Database content does not update with the code: re-run `gapit db fetch <name> --force`
to rebuild an installed database from upstream.
