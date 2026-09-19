# Installing gapit

gapit is installed from a git clone: [pixi](https://pixi.sh) creates the environment, so
you don't manage Python or external tools yourself.

## Prerequisites

- [git](https://git-scm.com) and [pixi](https://pixi.sh). On macOS/Linux:
  `curl -fsSL https://pixi.sh/install.sh | bash`
- Python 3.11+ if you install the package outside pixi (pip/pyproject). With pixi this is
  moot: the environment ships its own Python (the dev env pins 3.14).

## Install

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
gapit 0.1.0
$ gapit --version --json
{"schema":"gapit.version/1","name":"gapit","version":"0.1.0"}
```

## External binaries

gapit shells out to BLAST+ (`blastn`, `blastx`, `makeblastdb`, `blastdbcmd`),
`any2fasta`, and `minimap2`. All three come from the pixi environment (conda-forge and
bioconda), so there is nothing to install by hand. If you run gapit outside pixi, make
sure these binaries are on PATH yourself.

## Database bootstrap

Screening needs at least one database in the datadir (`$GAPIT_DATADIR`, then
`~/.local/share/gapit/db`; override per call with `--datadir`). `card` and `vfdb` ship
inside the package and install with no network; the other providers download from
upstream when fetched. See [Databases](./databases.md) for the full provider table.

```bash
gapit db fetch            # installs the default set (card, vfdb) into the default datadir
```

The same install into a scratch datadir, so you can see what a fetch prints. Progress
lines go to stderr, one JSON receipt per database to stdout:

```console
$ gapit db fetch --datadir /tmp/gapit-docs/dd
gapit: installed card from bundled snapshot card.tar.gz
gapit: generated /tmp/gapit-docs/dd/card/sequences
gapit: self-check passed for card
gapit: BLAST index built (nucl)
{"db":"card","records":6059,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/card"}
gapit: installed vfdb from bundled snapshot vfdb.tar.gz
gapit: generated /tmp/gapit-docs/dd/vfdb/sequences
gapit: self-check passed for vfdb
gapit: BLAST index built (nucl)
{"db":"vfdb","records":4769,"dbtype":"nucl","destination":"/tmp/gapit-docs/dd/vfdb"}
```

Confirm what's installed any time with `gapit list` (or `gapit list --json`).

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
to rebuild an installed database from upstream (or its bundled snapshot).
