# Data directory and database layout

gapit keeps **tool environments** and **database state** separate. This page explains
where databases live, what a database directory contains, why the default is not inside
a conda environment, and how sharing and containers behave.

If you only need the command: `gapit db path` prints the effective datadir, and
`gapit db list` lists the catalog after bootstrapping that directory.

## Quick decision table

| Situation | Recommended datadir |
|---|---|
| One user, several conda/envs | Default user-level datadir (shared cache); refresh after gapit upgrades |
| Project needs reproducible DB versions | Explicit `--datadir` or `GAPIT_DATADIR` per project |
| Persistent machine / workstation | `$XDG_DATA_HOME/gapit/db`, else `~/.local/share/gapit/db` |
| CI | Cache or mount a workspace datadir; pin DB versions |
| Docker / Podman | Named volume or bind mount + `GAPIT_DATADIR`; never rely on container `$HOME` |
| Multi-user HPC | Per-user datadir; shared read-only references need an explicit permission model |
| Strict safety boundary | Do not share a writable datadir across users or untrusted jobs |

## What is a gapit database?

A database is a **directory under the datadir**, one subdirectory per database name:

```text
<datadir>/<name>/
├── sequences                  # required: FASTA used by BLAST/minimap2
├── gapit-manifest.json        # provenance + certification (gapit.manifest/1)
├── records.jsonl              # gene truth store (gapit-built gene databases)
├── features.json              # locus/gene table (cluster databases)
├── typing.json                # optional gapit.typing/1 or /2 rules
├── floors.json                # optional per-gene identity floors
└── sequences.n* / sequences.p* / sequences.nhr ...   # BLAST index
```

`sequences` is the screening input. Its headers carry the gene metadata: native gapit
databases use the tagged `gapit/v1` form (`gapit|db=...|gene=...|acc=...|func=...` with
the product after a space); legacy abricate databases use `~~~`-separated headers.
gapit can read both.

The BLAST index (`.nin/.nsq/.nhr` for nucleotide, `.pin/.psq/.phr` for protein) is a
generated artifact, not source data. A missing or stale index is rebuilt by
`gapit setupdb` or by the normal build pipeline.

`gapit-manifest.json` is written **last**. It certifies the artifacts (`sha256`,
`n_records`, `dbtype`, `kind`, `source`, `typing_schema`, licenses, tool versions).
Discovery and listing use it, but a database directory remains readable without one for
legacy abricate compatibility.

### Gene databases

A gene database holds one FASTA record per known gene/marker. Contig screening uses
`blastn` (or `blastx` for protein databases); reads screening uses `minimap2`. A gene
database may carry an optional `typing.json` so `gapit typing` can designate phenotypes
from a screen table.

### Cluster databases

A cluster database (`kind: cluster`) holds one FASTA record per **locus** plus
`features.json`, which maps loci to genes and coordinates. `gapit screen` runs the
minimap2 cluster engine and scores genes/loci from the `cs` walk. A multi-scheme
`typing.json` makes the cluster output `gapit.cluster/2`; otherwise the database renders
`gapit.cluster/1`. `vp_ok` is the bundled typed-cluster example.

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
2. `$GAPIT_DATADIR`
3. `$XDG_DATA_HOME/gapit/db` (only when `XDG_DATA_HOME` is absolute)
4. `~/.local/share/gapit/db`

`gapit db path` prints the effective path **without creating it**. Add `--json` for a
`gapit.datadir/1` document with `datadir`, `exists`, and `source` (`cli`, `env`, `xdg`,
or `default`) — useful in Docker entrypoints and agent workflows.

Commands differ on whether they create the directory:

| Command | Missing datadir |
|---|---|
| `gapit db list` | Creates the empty root, then lists the catalog |
| `gapit setupdb` | Creates the root, materializes bundled databases, indexes all |
| `gapit db fetch` / `db build` | Creates the root as part of the write path |
| `screen` / `typing` against a bundled name | Creates the root and materializes that bundle |
| `screen` / `typing` against a fetched name | Requires it (`DATADIR_NOT_FOUND`, exit 4) |
| `db search` / `db outdated` | Requires it (`DATADIR_NOT_FOUND`, exit 4) |

`db list` and `setupdb` are the explicit first-run bootstrap paths. `db path` is pure
introspection.

## Why not `$CONDA_PREFIX` by default?

Conda environments are **toolchain environments**, not user-data stores:

- an env can be read-only or shared;
- it can be deleted and recreated frequently;
- runtime writes into `$CONDA_PREFIX` are generally discouraged;
- every env would duplicate large reference databases;
- a wheel installed in an env is not automatically the database version you want at
  runtime (bundles materialize once into the shared datadir).

For **immutable reference data shipped by a conda package**, `$CONDA_PREFIX/share/...`
can be a reasonable location. gapit's databases are fetched, built, refreshed, and
customized at runtime, so they are user state and belong under the user data directory.

If you do want environment-local data, opt in explicitly:

```bash
export GAPIT_DATADIR="$CONDA_PREFIX/db/gapit"
gapit db list
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
wheel's bundled snapshot with the installed copy. After upgrading gapit in one env,
either refresh explicitly (`gapit db fetch <name> --force`) or use a per-project datadir.

**Name collisions.** Two envs building the same custom name do not silently overwrite
without `--force`, but `--force` affects every env sharing the datadir. Use unique names
or explicit project datadirs.

**Writer concurrency.** Concurrent `db fetch`/`db build` writers now serialize on a
hidden per-database `flock`, and builds happen in a hidden staging directory before a
rename swap. Readers normally see the old database or the new one; there is still a
short rename window where the canonical path is being replaced. Do not schedule
concurrent writes to the same database against a shared datadir if a reader cannot
tolerate that window. `gapit setupdb` rebuilds BLAST indices in place, so do not run it
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
docker volume create gapit-data

docker run --rm -it \
  --user "$(id -u):$(id -g)" \
  -v gapit-data:/data/gapit \
  -e GAPIT_DATADIR=/data/gapit \
  gapit:0.6.5 gapit db list
```

For CI, mount or cache a workspace directory and pin versions:

```bash
docker run --rm \
  -v "$PWD/.gapit-db:/data/gapit" \
  -e GAPIT_DATADIR=/data/gapit \
  gapit:0.6.5 gapit screen sample.fna --db ecoli_dec
```

You can bake reference databases into an image for **read-only** screening:

```dockerfile
RUN gapit setupdb --datadir /opt/gapit/db
ENV GAPIT_DATADIR=/opt/gapit/db
```

Runtime `db fetch`/`db build` writes then land in the container layer (and vanish) or
fail with a read-only rootfs. If the container must write custom databases, mount a
writable volume; if it also needs baked references, seed the volume from the image or use
a project datadir instead.

## Security and reproducibility checklist

- Keep the datadir under the user's home or an explicit project/volume path.
- Use `gapit db path --json` in entrypoints rather than guessing paths.
- Pin database versions for published analyses; record the manifest.
- Refresh consciously after a gapit upgrade.
- Do not share a writable datadir across users without ACLs.
- Treat provider downloads as a trust boundary; verify upstream content where possible.
- Back up `records.jsonl`/custom inputs, not the generated BLAST index.

## Related commands

```text
gapit db path [--json]      # resolved datadir, no creation
gapit db list [--json]      # catalog + install state; bootstraps root
gapit setupdb               # bootstrap + materialize + index
gapit db fetch NAME         # provider download/build
gapit db build NAME INPUT   # custom gene or cluster database
gapit db outdated           # installed DB age (strict; requires datadir)
```
