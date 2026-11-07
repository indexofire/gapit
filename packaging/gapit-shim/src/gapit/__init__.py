"""Compatibility package for the ``gapit`` import name.

The software is now distributed as ``gspec`` and its implementation lives in
the ``gspec`` package. This stub keeps the pre-0.8 import surface working:

    import gapit
    from gapit.errors import GapitError      # same class object as gspec
    from gapit.cli import app
    import gapit.typing_engine

``gapit.__path__`` points at the ``gspec`` package directory, so
``gapit.<module>`` loads from the same source tree. A small set of modules
whose *identity* matters is aliased eagerly, so ``except
gapit.errors.GapitError`` still catches what ``gspec`` raises.

Nothing else lives in this distribution: the commands (``gspec``, ``gapit``,
``gspec-mcp``, ``gapit-mcp``) and the implementation are installed by
``gspec``, which this package depends on.
"""

from __future__ import annotations

import sys

import gspec

__version__ = gspec.__version__

#: The distribution that now carries the implementation.
RENAMED_TO = "gspec"

#: Last release published under the old project name.
DEPRECATED_AFTER = "0.7.3"

# Submodules resolve out of the gspec package tree.
__path__ = list(gspec.__path__)

# Modules whose identity must be shared: exception types, dataclasses, and
# anything else a caller may compare or catch across the two names.
_SHARED = (
    "errors",
    "config",
    "db",
    "records",
    "hits",
    "typing_models",
    "typing_rules",
    "serotype_models",
    "cluster_profile",
    "locus_map_models",
)

for _name in _SHARED:
    try:
        _module = __import__(f"gspec.{_name}", fromlist=["_"])
    except ImportError:  # optional or moved module: leave it to the path alias
        continue
    sys.modules[f"{__name__}.{_name}"] = _module
    setattr(sys.modules[__name__], _name, _module)

del _name, _module, sys
