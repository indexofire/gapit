"""Canonical schema-id handling for the ``gspec.*`` document names.

0.9.0 renamed every document id from ``gapit.<doc>/<n>`` to
``gspec.<doc>/<n>``. Documents written by older releases are still valid
input, so every model canonicalises the id on validation: a ``gapit.*``
value is accepted and re-emitted as ``gspec.*``. Nothing else about the
documents changes.
"""

from __future__ import annotations

from pydantic import BeforeValidator

#: Prefix used by 0.8 and earlier.
LEGACY_PREFIX = "gapit."

#: Prefix emitted by every current document.
CURRENT_PREFIX = "gspec."


def canonical_schema_id(value: object) -> object:
    """Map a legacy ``gapit.*`` schema id to its ``gspec.*`` spelling.

    Unknown shapes pass through untouched so pydantic's ``Literal`` check
    still produces the usual validation error.
    """
    if isinstance(value, str) and value.startswith(LEGACY_PREFIX):
        return CURRENT_PREFIX + value[len(LEGACY_PREFIX) :]
    return value


#: Attach to a ``Literal["gspec.<doc>/<n>"]`` annotation to accept the
#: pre-0.9 ``gapit.<doc>/<n>`` spelling as input.
Canonical = BeforeValidator(canonical_schema_id)

__all__ = ["CURRENT_PREFIX", "LEGACY_PREFIX", "Canonical", "canonical_schema_id"]
