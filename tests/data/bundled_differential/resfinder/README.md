# resfinder differential fixture

`HEAD.zip` is a byte-exact copy of the upstream ResFinder database archive
(https://bitbucket.org/genomicepidemiology/resfinder_db/get/HEAD.zip, downloaded
2026-10-02, sha256 526186759ae4e182fbb8b4a12e05dcaca9f5acf920db7a498a8fec9d78109109)
— the exact bytes the wheel's bundled `resfinder` snapshot was produced from.

License: Apache-2.0 (resfinder_db). Redistribution permitted; kept here as a test
fixture only (never packaged in the wheel — the wheel carries the built snapshot
under `src/gapit/data/dbs/resfinder/`).

`tests/test_bundled.py::test_resfinder_snapshot_is_byte_identical_to_a_fetch_build`
monkeypatches the provider URL to this file and proves the bundled-built and
fetch-built `sequences` are byte-identical. If you regenerate the bundled snapshot
from a newer upstream, refresh this archive in the same change or that test fails.
