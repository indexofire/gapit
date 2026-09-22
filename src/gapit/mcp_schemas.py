"""MCP tools/list declarations: names, descriptions, inputSchemas (protocol: gapit.mcp).

Pure data — the roster mirrors the handler map in :mod:`gapit.mcp_tools`.
Descriptions are wire-visible wording; treat any edit like a schema change
(AGENTS.md §5).
"""

from gapit.cmd_db_outdated import DEFAULT_STALE_DAYS
from gapit.cmd_db_search import DEFAULT_LIMIT, SearchField
from gapit.formats.schemas import SCHEMA_MODELS
from gapit.reads import ReadTypeEnum

_FILES: dict[str, object] = {"type": "array", "items": {"type": "string"}}
_FLAG: dict[str, object] = {"type": "boolean", "default": False}
_STR: dict[str, object] = {"type": "string"}
_DAYS: dict[str, object] = {"type": "integer", "minimum": 0, "default": DEFAULT_STALE_DAYS}
_SEARCH_FIELDS: list[str] = [field.value for field in SearchField]
_READ_TYPES: list[str] = [preset.value for preset in ReadTypeEnum]


def _tool_entry(
    name: str,
    description: str,
    properties: dict[str, object] | None = None,
    required: list[str] | None = None,
) -> dict[str, object]:
    schema = {"type": "object", "properties": properties or {}, "required": required or []}
    return {"name": name, "description": description, "inputSchema": schema}


TOOLS: list[dict[str, object]] = [
    _tool_entry(
        "screen",
        "Screen contig files for AMR/virulence genes (json = gapit.report/1;"
        " aligner minimap2 = fast assembly survey emitting gapit.reads/1).",
        dict(
            files=_FILES,
            db={"type": "string", "default": "ncbi"},
            minid={"type": "number"},
            mincov={"type": "number"},
            format={"type": "string", "enum": ["json", "tsv", "md"], "default": "json"},
            aligner={"type": "string", "enum": ["blastn", "minimap2"], "default": "blastn"},
            min_breadth={"type": "number", "minimum": 0, "maximum": 100, "default": 90},
            min_identity={"type": "number", "minimum": 0, "maximum": 100, "default": 0},
            min_mapq={"type": "integer", "minimum": 0, "default": 0},
            datadir=_STR,
        ),
        ["files"],
    ),
    _tool_entry(
        "screen_reads",
        "Screen FASTQ reads for genes via minimap2 (json = gapit.reads/1;"
        " min_identity/min_mapq > 0 emits gapit.reads/2). Returns the rendered"
        " document.",
        dict(
            r1=_FILES,
            r2=_FILES,
            read_type={"type": "string", "enum": _READ_TYPES, "default": "sr"},
            min_breadth={"type": "number", "minimum": 0, "maximum": 100, "default": 90},
            min_identity={"type": "number", "minimum": 0, "maximum": 100, "default": 0},
            min_mapq={"type": "integer", "minimum": 0, "default": 0},
            format={"type": "string", "enum": ["json", "md"], "default": "json"},
            db={"type": "string", "default": "ncbi"},
            datadir=_STR,
        ),
        ["r1"],
    ),
    _tool_entry(
        "summary",
        "Summarize report tables into a gapit.summary/1 matrix.",
        dict(files=_FILES, identity={"type": "boolean"}, nopath={"type": "boolean"}),
        ["files"],
    ),
    _tool_entry(
        "schema",
        "Print the JSON Schema of a gapit output document.",
        dict(name={"type": "string", "enum": sorted(SCHEMA_MODELS)}),
        ["name"],
    ),
    _tool_entry("db_list", "List database providers and their installed state (gapit.dblist/1)."),
    _tool_entry(
        "db_fetch",
        "Fetch provider database(s) into the datadir (name omitted: card+vfdb from bundled"
        " snapshots; network installs can take minutes). One JSON receipt line per database"
        " (db, records, dbtype, destination).",
        dict(name=_STR, datadir=_STR, force=_FLAG),
    ),
    _tool_entry(
        "db_build",
        "Build a custom database from a LOCAL FASTA filesystem path (plain, abricate ~~~, or"
        " gapit| headers; .gz/.bz2 accepted). Returns a JSON receipt (db, records, dbtype,"
        " destination).",
        dict(
            name=_STR,
            fasta=_STR,
            tsv=_STR,
            dbtype={"type": "string", "enum": ["nucl", "prot"]},
            description=_STR,
            datadir=_STR,
            force=_FLAG,
        ),
        ["name", "fasta"],
    ),
    _tool_entry(
        "db_search",
        "Search installed databases for a term (case-insensitive substring, or exact"
        " full-field equality). Hit rows are TSV: DB, GENE, ACCESSION, FUNCTION, PRODUCT, LENGTH.",
        dict(
            term=_STR,
            db=_STR,
            field={"type": "string", "enum": _SEARCH_FIELDS, "default": "any"},
            exact=_FLAG,
            limit={"type": "integer", "minimum": 0, "default": DEFAULT_LIMIT},
            datadir=_STR,
        ),
        ["term"],
    ),
    _tool_entry(
        "db_outdated",
        "Report installed database ages against the staleness threshold (days) and newer"
        " bundled snapshots. Rows are TSV with columns NAME, FETCHED_AT, AGE_DAYS, STATUS.",
        dict(days=_DAYS, datadir=_STR),
    ),
]
