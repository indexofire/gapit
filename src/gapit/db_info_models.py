"""Shared models for `gapit db info` provenance (a leaf module, no imports).

Two consumers need these shapes: :mod:`gapit.bundled` (the wheel's
``bundled.json`` declares the literature and caveats that ship with a
database) and :mod:`gapit.db_info_ops` (the ``gapit.dbinfo/1`` document).
Keeping them here means neither has to import the other.
"""

from pydantic import BaseModel, Field


class DbReference(BaseModel, frozen=True):
    """One citable source behind a database's content or scheme.

    ``citation`` is the short human form (authors, year, title, venue);
    ``doi`` and ``url`` are optional machine-followable pointers. Empty
    strings stay empty rather than None so JSON consumers get a stable
    shape.
    """

    citation: str = Field(min_length=1)
    doi: str = ""
    url: str = ""
    note: str = ""


class DbGuidance(BaseModel, frozen=True):
    """The database's own guidance: how it may be used and where it misleads.

    ``caveats`` are the interpretation limits a caller must not skip
    (unresolvable subtypes, documented assay artefacts, threshold effects);
    ``usage`` is the how-to (the commands that constitute the intended
    workflow). Both serialize as plain string lists.
    """

    references: tuple[DbReference, ...] = ()
    caveats: tuple[str, ...] = ()
    usage: tuple[str, ...] = ()
