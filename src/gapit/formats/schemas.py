"""The schema registry: every document name ``gapit schema`` can introspect.

Lives beside the document models but one level up: it needs the summary
document (formats/summary.py) and the error envelope (errors.py), and
formats/summary.py already imports from formats/json.py — so this module, not
formats/json.py, is the cycle-free single home. Insertion order is part of the
public surface (``gapit schema`` unknown-name message iterates it).
"""

from pydantic import BaseModel

from gapit.errors import ErrorEnvelope
from gapit.formats.cluster import ClusterDocument
from gapit.formats.json import ReportDocument, VersionDocument
from gapit.formats.reads_json import Reads2Document, ReadsDocument
from gapit.formats.summary import SummaryDocument
from gapit.formats.typing_result import TypingResultDocument
from gapit.gbfeatures import FeaturesDocument
from gapit.typing_models import TypingDocument

SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "report": ReportDocument,
    "reads": ReadsDocument,
    "reads2": Reads2Document,
    "cluster": ClusterDocument,
    "summary": SummaryDocument,
    "error": ErrorEnvelope,
    "version": VersionDocument,
    "features": FeaturesDocument,
    "typing": TypingDocument,
    "typing_result": TypingResultDocument,
}
