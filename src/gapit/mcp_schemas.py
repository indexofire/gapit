"""MCP tools/list declarations: names, descriptions, inputSchemas (protocol: gapit.mcp).

Pure data — the roster mirrors the handler map in :mod:`gapit.mcp_tools`.
Descriptions are wire-visible wording; treat any edit like a schema change
(AGENTS.md §5).
"""

from gapit.db_query_ops import DEFAULT_LIMIT, DEFAULT_STALE_DAYS, SearchField
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
        "Screen contig files for known genes (json = gapit.report/1;"
        " aligner minimap2 = fast assembly survey emitting gapit.reads/1)."
        " A cluster-kind database dispatches to the cluster engine"
        " (minimap2 asm20; json = gapit.cluster/1 with a best-locus call,"
        " per-locus gene verdicts, the phenotype call + score breakdown"
        " when the db carries a single-scheme typing.json, or"
        " gapit.cluster/2 with a phenotypes map of every scheme's call"
        " when the typing.json holds more than one scheme, and the"
        " min_gene_cov/min_gene_id/min_cluster_cov thresholds).",
        dict(
            files=_FILES,
            db=_STR,
            minid={"type": "number"},
            mincov={"type": "number"},
            format={"type": "string", "enum": ["json", "tsv", "csv", "md"], "default": "json"},
            aligner={"type": "string", "enum": ["blastn", "minimap2"], "default": "blastn"},
            merge_fragments=_FLAG,
            min_breadth={"type": "number", "minimum": 0, "maximum": 100, "default": 90},
            min_identity={"type": "number", "minimum": 0, "maximum": 100, "default": 0},
            min_mapq={"type": "integer", "minimum": 0, "default": 0},
            all_genes=_FLAG,
            min_gene_cov={"type": "number", "minimum": 0, "maximum": 100, "default": 90},
            min_gene_id={"type": "number", "minimum": 0, "maximum": 100, "default": 90},
            min_cluster_cov={"type": "number", "minimum": 0, "maximum": 100, "default": 96},
            datadir=_STR,
        ),
        ["files", "db"],
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
            all_genes=_FLAG,
            format={"type": "string", "enum": ["json", "tsv", "csv", "md"], "default": "json"},
            db=_STR,
            datadir=_STR,
        ),
        ["r1", "db"],
    ),
    _tool_entry(
        "summary",
        "Summarize report tables into a gapit.summary/1 matrix.",
        dict(files=_FILES, identity={"type": "boolean"}, nopath={"type": "boolean"}),
        ["files"],
    ),
    _tool_entry(
        "typing",
        "Designate phenotypes from a gapit/abricate screen result table (TSV or CSV"
        " text). Returns gapit.typing_result/1 JSON by default; its notes carry"
        " any follow-up guidance (for meningococcal W/Y, the synG EX7E check).",
        dict(
            table=_STR,
            datadir=_STR,
            format={"type": "string", "enum": ["json", "tsv", "md"], "default": "json"},
        ),
        ["table"],
    ),
    _tool_entry(
        "schema",
        "Print the JSON Schema of a gapit output document.",
        dict(name={"type": "string", "enum": sorted(SCHEMA_MODELS)}),
        ["name"],
    ),
    _tool_entry(
        "db_list",
        "List known databases with NAME, upstream PROVIDER, and install state (gapit.dblist/1).",
        dict(datadir=_STR),
    ),
    _tool_entry(
        "db_path",
        "Resolve the effective datadir and its source without creating it (gapit.datadir/1).",
        dict(datadir=_STR),
    ),
    _tool_entry(
        "db_setup",
        "Bootstrap the datadir, materialize bundled databases, and rebuild BLAST"
        " indices (gapit.dbsetup/1). names selects databases (missing providers are"
        " fetched, missing bundles materialized); no_bundles skips bundled"
        " materialization and indexes only what is already installed.",
        dict(datadir=_STR, no_bundles=_FLAG, names=_FILES),
    ),
    _tool_entry(
        "db_fetch",
        "Fetch database(s) by NAME into the datadir (name 'all': the card+vfdb"
        " default set; every database downloads from its upstream provider, so installs"
        " need network and can take minutes). One JSON receipt line per database"
        " (db, records, dbtype, destination).",
        dict(name=_STR, datadir=_STR, force=_FLAG),
        ["name"],
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
            kind={"type": "string", "enum": ["gene", "cluster"]},
            typing=_STR,
            floors=_STR,
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
        "Report installed database ages against the staleness threshold (days)."
        " Rows are TSV with columns NAME, FETCHED_AT, AGE_DAYS, STATUS.",
        dict(days=_DAYS, datadir=_STR),
    ),
]
