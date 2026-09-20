"""The schema registry: every document name ``gapit schema`` can introspect.

Lives beside the document models but one level up: it needs the summary
document (formats/summary.py) and the error envelope (errors.py), and
formats/summary.py already imports from formats/json.py — so this module, not
formats/json.py, is the cycle-free single home. Insertion order is part of the
public surface (``gapit schema`` unknown-name message iterates it).
"""

from pydantic import BaseModel

from gapit.errors import ErrorEnvelope
from gapit.formats.json import (
    ListDocument,
    Reads2Document,
    ReadsDocument,
    ReportDocument,
    VersionDocument,
)
from gapit.formats.summary import SummaryDocument

SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "report": ReportDocument,
    "reads": ReadsDocument,
    "reads2": Reads2Document,
    "summary": SummaryDocument,
    "list": ListDocument,
    "error": ErrorEnvelope,
    "version": VersionDocument,
}
