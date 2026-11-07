# Data directory and database layout

gspec keeps **tool environments** and **database state** separate. This page explains
where databases live, what a database directory contains, why the default is not inside
a conda environment, and how sharing and containers behave.

If you only need the command: `gspec db path` prints the effective datadir, and
`gspec db list` lists the catalog after bootstrapping that directory.

## Quick decision table

| Situation | Recommended datadir |
|---|---|
| One user, several conda/envs | Default user-level datadir (shared cache); refresh after gspec upgrades |
| Project needs reproducible DB versions | Explicit `--datadir` or `GSPEC_DATADIR` per project |
| Persistent machine / workstation | `$XDG_DATA_HOME/gspec/db`, else `~/.local/share/gspec/db` |
| CI | Cache or mount a workspace datadir; pin DB versions |
| Docker / Podman | Named volume or bind mount + `GSPEC_DATADIR`; never rely on container `$HOME` |
| Multi-user HPC | Per-user datadir; shared read-only references need an explicit permission model |
| Strict safety boundary | Do not share a writable datadir across users or untrusted jobs |

## What is a gspec database?

A database is a **directory under the datadir**, one subdirectory per database name:

```text
<datadir>/<name>/
├── sequences                  # required: FASTA used by BLAST/minimap2
├── gapit-manifest.json        # provenance + certification (gapit.manifest/1)
├── records.jsonl              # gene truth store (gspec-built gene databases)
├── features.json              # locus/gene table (cluster databases)
├── typing.json                # optional gapit.typing/1 or /2 rules
├── floors.json                # optional per-gene identity floors
├── sequences.asm20.mmi        # cluster DBs: prebuilt minimap2 asm20 index
└── sequences.n* / sequences.p* / sequences.nhr ...   # BLAST index
```

`sequences` is the screening input. Its headers carry the gene metadata: native gspec
databases use the tagged `gapit/v1` form (`gapit|db=...|gene=...|acc=...|func=...` with
the product after a space); legacy abricate databases use `~~~`-separated headers.
gspec can read both.

The BLAST index (`.nin/.nsq/.nhr` for nucleotide, `.pin/.psq/.phr` for protein) is a
generated artifact, not source data. A missing or stale index is rebuilt by
`gspec db setup` or by the normal build pipeline. Cluster databases additionally carry a
`sequences.asm20.mmi` minimap2 target index, which `db setup` also refreshes; legacy
cluster screening uses it when fresh.

`gapit-manifest.json` is written **last**. It certifies the artifacts (`sha256`,
`n_records`, `dbtype`, `kind`, `source`, `typing_schema`, licenses, tool versions).
Discovery and listing use it, but a database directory remains readable without one for
legacy abricate compatibility.

### Gene databases

A gene database holds one FASTA record per known gene/marker. Contig screening uses
`blastn` (or `blastx` for protein databases); reads screening uses `minimap2`. A gene
database may carry an optional `typing.json` so `gspec typing` can designate phenotypes
from a screen table.

### Cluster databases

A cluster database (`kind: cluster`) holds one FASTA record per **locus** plus
`features.json`, which maps loci to genes and coordinates. `gspec screen` runs the
minimap2 cluster engine and scores genes/loci from the `cs` walk. A multi-scheme
`typing.json` makes the cluster output `gapit.cluster/2`; otherwise the database renders
`gapit.cluster/1`. `vpa_ok` is the bundled typed-cluster example.

### Where databases come from

| Source | Meaning |
|---|---|
| `db fetch` | Provider catalog: downloads from upstream at fetch time |
| Bundled | Content shipped inside the wheel, materialized on first use |
| `db build` | Custom FASTA, GBK, or GFF/GBK cluster database |
| `db install` | SHA256-verified local bytes, atomically installed |
| abricate datadir | Legacy `~~~` headers, no manifest required |

Bundled databases are still user data after materialization: they are copied and indexed
under the datadir, not executed from the wheel.

## Where the datadir lives

Resolution order, highest first:

1. `--datadir PATH` / `-D PATH`
2. `$GSPEC_DATADIR`
3. `$XDG_DATA_HOME/gspec/db` (only when `XDG_DATA_HOME` is absolute)
4. `~/.local/share/gspec/db`

`gspec db path` prints the effective path **without creating it**. Add `--json` for a
`gapit.datadir/1` document with `datadir`, `exists`, and `source` (`cli`, `env`, `xdg`,
or `default`) — useful in Docker entrypoints and agent workflows.

Two rename accommodations keep pre-0.8 installs working: `$GAPIT_DATADIR` is still read
when `$GSPEC_DATADIR` is unset, and when only the pre-0.8 `~/.local/share/gapit/db`
exists it keeps being used as the home default. Run
`mv ~/.local/share/gapit ~/.local/share/gspec` to adopt the new path.

Commands differ on whether they create the directory:

| Command | Missing datadir |
|---|---|
| `gspec db list` | Creates the empty root, then lists the catalog |
| `gspec db setup` | Creates the root, materializes bundled databases, indexes all |
| `gspec db fetch` / `db build` | Creates the root as part of the write path |
| `screen` / `typing` against a bundled name | Creates the root and materializes that bundle |
| `screen` / `typing` against a fetched name | Requires it (`DATADIR_NOT_FOUND`, exit 4) |
| `db search` / `db outdated` | Requires it (`DATADIR_NOT_FOUND`, exit 4) |

`db list` and `db setup` are the explicit first-run bootstrap paths. `db path` is pure
introspection.

### Bootstrapping with `gspec db setup`

`gspec db setup` creates a missing datadir, materializes bundled databases that are
not installed yet, and rebuilds the BLAST index for every database it finds. An
already-installed bundle (manifest present) is left untouched; every discovered
database is reindexed in place, so the command is a repair path as well as a first-run
path.

Select databases explicitly:

```bash
gspec db setup                 # same as --all (kept for compatibility)
gspec db setup -a              # all installed + bundled databases
gspec db setup -d nme_sero   # only nme_sero
gspec db setup -d eco_dec -d cje_cps
```

`--all/-a` never fetches providers: it materializes missing bundled databases and
rebuilds the index for everything installed or bundled. `--db/-d NAME` is explicit
per database: an installed name is reindexed, a missing bundled name is materialized,
and a missing provider name is fetched with the normal `db fetch` pipeline before it
is indexed. A name that is none of those is a `DATABASE_NOT_FOUND` error before
anything is written.

| Option | Effect |
|---|---|
| `--all`, `-a` | Set up every installed/bundled database (never fetches providers) |
| `--db`, `-d NAME` | Set up NAME (repeatable); missing providers are fetched, missing bundles materialized |
| `--datadir`, `-D PATH` | Use PATH instead of the resolved datadir |
| `--no-bundles`, `-B` | Skip bundled materialization (named providers may still be fetched) |
| `--json`, `-J` | Print one `gapit.dbsetup/1` receipt on stdout and suppress progress |
| `--quiet`, `-q` | Silence progress diagnostics on stderr |
| `--debug`, `-v` | Echo external command lines (`makeblastdb`) to stderr |

The JSON receipt names the bundles that were materialized, the providers that were
fetched, and every database that was indexed:

```json
{"schema":"gapit.dbsetup/1","datadir":"/path/to/db","materialized":["eco_dec"],"fetched":[],"databases":[{"db":"eco_dec","records":16,"dbtype":"nucl"}]}
```

Introspect it with `gspec schema dbsetup`.

## Why not `$CONDA_PREFIX` by default?

Conda environments are **toolchain environments**, not user-data stores:

- an env can be read-only or shared;
- it can be deleted and recreated frequently;
- runtime writes into `$CONDA_PREFIX` are generally discouraged;
- every env would duplicate large reference databases;
- a wheel installed in an env is not automatically the database version you want at
  runtime (bundles materialize once into the shared datadir).

For **immutable reference data shipped by a conda package**, `$CONDA_PREFIX/share/...`
can be a reasonable location. gspec's databases are fetched, built, refreshed, and
customized at runtime, so they are user state and belong under the user data directory.

If you do want environment-local data, opt in explicitly:

```bash
export GSPEC_DATADIR="$CONDA_PREFIX/db/gspec"
gspec db list
```

The caveat is the opposite of the default: deleting the env deletes the databases; each
env gets its own copy; and nested shells must restore or unset the variable.

## Sharing a datadir across environments

Sharing is a feature, not an accident: one download, one index build, and identical
reference content across envs. It is **not** an isolation or security boundary.

### What is safe

- Concurrent **reads** (`screen`, `typing`, `db list`) of immutable database files.
- One writer at a time per database.
- Different custom database names; they coexist.

### What needs care

**Version skew.** Environment A may ship a newer bundled snapshot than the datadir
already contains. `db list` shows installed state, but `db outdated` does not compare a
wheel's bundled snapshot with the installed copy. After upgrading gspec in one env,
either refresh explicitly (`gspec db fetch -d <name> --force`) or use a per-project datadir.

**Name collisions.** Two envs building the same custom name do not silently overwrite
without `--force`, but `--force` affects every env sharing the datadir. Use unique names
or explicit project datadirs.

**Writer concurrency.** Concurrent `db fetch`/`db build` writers now serialize on a
hidden per-database `flock`, and builds happen in a hidden staging directory before a
rename swap. Readers normally see the old database or the new one; there is still a
short rename window where the canonical path is being replaced. Do not schedule
concurrent writes to the same database against a shared datadir if a reader cannot
tolerate that window. `gspec db setup` rebuilds BLAST indices in place, so do not run it
concurrently with `db fetch`/`db build` for the same database.

**Permissions.** New datadirs and transaction staging directories are created `0700`.
A shared multi-user datadir needs an explicit group/ACL policy; the default is per-user.

**Trust.** Some providers are HTTP-only and downloads are not checksum-pinned. A shared
datadir broadcasts the impact of one compromised fetch. For high-assurance use, pin
database versions/checksums and prefer a read-only shared reference directory.

## Docker

Container `$HOME` and `$CONDA_PREFIX` writes are normally **ephemeral**: they live in the
container writable layer and disappear with the container. A read-only root filesystem
makes the default path fail entirely.

Use an explicit mount and environment variable:

```bash
docker volume create gspec-data

docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -v gspec-data:/data/gspec \
  -e GSPEC_DATADIR=/data/gspec \
  gspec:0.8.0 gspec db list
```

For CI, mount or cache a workspace directory and pin versions:

```bash
docker run --rm \
  -v "$PWD/.gspec-db:/data/gspec" \
  -e GSPEC_DATADIR=/data/gspec \
  gspec:0.8.0 gspec screen sample.fna --db eco_dec
```

You can bake reference databases into an image for **read-only** screening:

```dockerfile
RUN gspec db setup --datadir /opt/gspec/db
ENV GSPEC_DATADIR=/opt/gspec/db
```

Runtime `db fetch`/`db build` writes then land in the container layer (and vanish) or
fail with a read-only rootfs. If the container must write custom databases, mount a
writable volume; if it also needs baked references, seed the volume from the image or use
a project datadir instead.

## Security and reproducibility checklist

- Keep the datadir under the user's home or an explicit project/volume path.
- Use `gspec db path --json` in entrypoints rather than guessing paths.
- Pin database versions for published analyses; record the manifest.
- Refresh consciously after a gspec upgrade.
- Do not share a writable datadir across users without ACLs.
- Treat provider downloads as a trust boundary; verify upstream content where possible.
- Back up `records.jsonl`/custom inputs, not the generated BLAST index.

## Related commands

```text
gspec db path [--json]      # resolved datadir, no creation
gspec db list [--json]      # catalog + install state; bootstraps root
gspec db setup              # bootstrap + materialize + index
gspec db fetch -d NAME         # provider download/build
gspec db build NAME INPUT   # custom gene or cluster database
gspec db outdated           # installed DB age (strict; requires datadir)
```
