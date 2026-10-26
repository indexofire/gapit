# Security policy

## Supported versions

Security fixes are applied to the current `0.6.x` line. Older releases are
best-effort; users should upgrade with `pip install -U gapit`.

## Reporting a vulnerability

Please do **not** open a public issue for a suspected vulnerability. Use
GitHub's private security advisory flow:

<https://github.com/indexofire/gapit/security/advisories/new>

Include the gapit version (`gapit --version --json`), a minimal reproducer,
and whether the issue affects CLI, MCP, database fetching, or a bundled asset.
We will acknowledge and coordinate disclosure through the advisory.

## Trust boundaries

- **Provider downloads.** `gapit db fetch` downloads content from upstream
  providers. Some upstreams are HTTP-only and not checksum-pinned; treat fetched
  database content as third-party data. Pin versions/checksums for regulated or
  reproducible work, and prefer a read-only shared reference directory.
- **User-provided databases.** `gapit db build`/`db install` trust the local
  paths and bytes supplied by the user. A database can influence screening and
  typing results; review custom content before distributing it.
- **Shared datadir.** A writable datadir shared across users or untrusted jobs is
  not a security boundary. Writers serialize with a per-database lock, but
  readers/users share the same database state and trust assumptions. Use
  per-project `--datadir` or `GAPIT_DATADIR` where isolation matters.
- **External binaries.** gapit invokes BLAST+ and minimap2 with argv lists (no
  shell). Their integrity and PATH resolution are the host's responsibility.
- **MCP.** `gapit mcp` is a local stdio server. It does not authenticate clients;
  do not expose it as a network service without an authenticating proxy.

## Out of scope

- Incorrect biological interpretation caused by upstream database content or
  user-curated rules (report data issues to the upstream/curator).
- Denial of service from deliberately unbounded local inputs where no documented
  size guard exists; report it if a realistic guard is missing.
