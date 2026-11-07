"""The schema registry: every document name ``gspec schema`` can introspect.

Lives beside the document models but one level up: it needs the summary
document (formats/summary.py) and the error envelope (errors.py), and
formats/summary.py already imports from formats/json.py — so this module, not
formats/json.py, is the cycle-free single home. Insertion order is part of the
public surface (``gspec schema`` unknown-name message iterates it).
"""

from pydantic import BaseModel

from gspec.db_info_ops import DbInfoDocument
from gspec.db_ops import DatadirDocument, DbListDocument
from gspec.db_query_ops import DbOutdatedDocument
from gspec.db_setup_ops import DbSetupDocument
from gspec.errors import ErrorEnvelope
from gspec.formats.cluster import ClusterDocument
from gspec.formats.cluster2 import Cluster2Document
from gspec.formats.json import ReportDocument, VersionDocument
from gspec.formats.locus_map import LocusMapDocument
from gspec.formats.reads_json import Reads2Document, ReadsDocument
from gspec.formats.summary import SummaryDocument
from gspec.formats.typing_result import TypingResultDocument
from gspec.gbfeatures import FeaturesDocument
from gspec.gene_floors import GeneFloors
from gspec.screen_defaults import ScreenDefaults
from gspec.serotype_models import SerotypeDocument
from gspec.typing_models import TypingDocument

SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "report": ReportDocument,
    "reads": ReadsDocument,
    "reads2": Reads2Document,
    "cluster": ClusterDocument,
    "cluster2": Cluster2Document,
    "locus_map": LocusMapDocument,
    "summary": SummaryDocument,
    "error": ErrorEnvelope,
    "version": VersionDocument,
    "dblist": DbListDocument,
    "dbsetup": DbSetupDocument,
    "dboutdated": DbOutdatedDocument,
    "datadir": DatadirDocument,
    "dbinfo": DbInfoDocument,
    "features": FeaturesDocument,
    "typing": TypingDocument,
    "typing_result": TypingResultDocument,
    "floors": GeneFloors,
    "screen": ScreenDefaults,
    "serotype": SerotypeDocument,
}
