"""vch_oagc provider: Vibrio cholerae O-antigen gene clusters (OAGC).

The 190 OL reference loci come from Murase et al. 2022 (Microb Genom
000860; NCBI BioProject PRJDB8600), fetched as one NCBI efetch GenBank
batch. A committed control locus (V. cholerae N16961 ompW) travels with
the wheel; the committed gapit.typing/2 asset defines the O-group and
species schemes. OLU loci are not in this build yet.
"""

from pathlib import Path

from gapit.db_info_models import DbGuidance, DbReference
from gapit.providers.cluster_common import ClusterProvider

_ASSETS = Path(__file__).parent / "assets" / "vch_oagc"

# The 190 O-antigen locus records used by the Murase et al. reference set.
# One comma-separated efetch request is the stable public download path.
_ACCESSIONS = (
    "LC594800",
    "LC594801",
    "LC594802",
    "LC594803",
    "LC594804",
    "LC594805",
    "LC594806",
    "LC594807",
    "LC594808",
    "LC594809",
    "LC594810",
    "LC594811",
    "LC594812",
    "LC594813",
    "LC594814",
    "LC594815",
    "LC594816",
    "LC594817",
    "LC594818",
    "LC594819",
    "LC594820",
    "LC594821",
    "LC594822",
    "LC594823",
    "LC594824",
    "LC594825",
    "LC594826",
    "LC594827",
    "LC594828",
    "LC594829",
    "LC594830",
    "LC594831",
    "LC594832",
    "LC594833",
    "LC594834",
    "LC594835",
    "LC594836",
    "LC594837",
    "LC594838",
    "LC594839",
    "LC594840",
    "LC594841",
    "LC594842",
    "LC594843",
    "LC594844",
    "LC594845",
    "LC594846",
    "LC594847",
    "LC594848",
    "LC594849",
    "LC594850",
    "LC594851",
    "LC594852",
    "LC594853",
    "LC594855",
    "LC594856",
    "LC594857",
    "LC594858",
    "LC594859",
    "LC594860",
    "LC594861",
    "LC594862",
    "LC594863",
    "LC594864",
    "LC594865",
    "LC594866",
    "LC594867",
    "LC594868",
    "LC594869",
    "LC594870",
    "LC594871",
    "LC594872",
    "LC594873",
    "LC594874",
    "LC594875",
    "LC594876",
    "LC594877",
    "LC594878",
    "LC594879",
    "LC594880",
    "LC594881",
    "LC594882",
    "LC594883",
    "LC594884",
    "LC594885",
    "LC594886",
    "LC594887",
    "LC594888",
    "LC594889",
    "LC594890",
    "LC594891",
    "LC594892",
    "LC594893",
    "LC594894",
    "LC594895",
    "LC594897",
    "LC594898",
    "LC594899",
    "LC594900",
    "LC594901",
    "LC594902",
    "LC594903",
    "LC594904",
    "LC594905",
    "LC594906",
    "LC594907",
    "LC594908",
    "LC594909",
    "LC594910",
    "LC594911",
    "LC594912",
    "LC594913",
    "LC594914",
    "LC594915",
    "LC594916",
    "LC594917",
    "LC594918",
    "LC594919",
    "LC594920",
    "LC594921",
    "LC594922",
    "LC594923",
    "LC594924",
    "LC594925",
    "LC594926",
    "LC594927",
    "LC594928",
    "LC594929",
    "LC594930",
    "LC594931",
    "LC594932",
    "LC594933",
    "LC594934",
    "LC594935",
    "LC594936",
    "LC594937",
    "LC594938",
    "LC594939",
    "LC594940",
    "LC594941",
    "LC594942",
    "LC594943",
    "LC594944",
    "LC594945",
    "LC594946",
    "LC594947",
    "LC594948",
    "LC594949",
    "LC594950",
    "LC594951",
    "LC594952",
    "LC594953",
    "LC594954",
    "LC594955",
    "LC594956",
    "LC594957",
    "LC594958",
    "LC594959",
    "LC594960",
    "LC594961",
    "LC594962",
    "LC594963",
    "LC594964",
    "LC594965",
    "LC594966",
    "LC594967",
    "LC594968",
    "LC594969",
    "LC594970",
    "LC594971",
    "LC594972",
    "LC594973",
    "LC594974",
    "LC594975",
    "LC594976",
    "LC594977",
    "LC594978",
    "LC594979",
    "LC594980",
    "LC594981",
    "LC594982",
    "LC594983",
    "LC594984",
    "LC594985",
    "LC594986",
    "LC594987",
    "LC594988",
    "LC594989",
    "LC594990",
    "LC594991",
)
_SOURCE_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=nuccore&id=" + ",".join(_ACCESSIONS) + "&rettype=gb&retmode=text"
)

PROVIDER = ClusterProvider(
    name="vch_oagc",
    description=(
        "Vibrio cholerae O-antigen gene clusters (OAGC; 190 OL loci) + ompW species control"
    ),
    vendor="gapit(Hangzhou CDC)",
    source_urls=(_SOURCE_URL,),
    license="CC BY 4.0 (paper) / INSDC public domain (sequences)",
    note=(
        "190 O-antigen biosynthetic locus records from Murase et al. 2022"
        " (Microb Genom 000860; NCBI BioProject PRJDB8600), fetched from NCBI,"
        " plus a committed V. cholerae N16961 ompW control locus; build uses"
        " a committed gapit.typing/2 O-group + species asset."
    ),
    typing=_ASSETS / "typing.json",
    extra_sources=(_ASSETS / "control.gbk",),
    guidance=DbGuidance(
        references=(
            DbReference(
                citation=(
                    "Murase K, Arakawa E, Izumiya H, Iguchi A, Takemura T,"
                    " Kikuchi T, Nakagawa I, Thomson NR, Ohnishi M, Morita M."
                    " Genomic dissection of the Vibrio cholerae O-serogroup global"
                    " reference strains: reassessing our view of diversity and"
                    " plasticity between two chromosomes."
                    " Microb Genom 2022;8(8):mgen000860."
                ),
                doi="10.1099/mgen.0.000860",
                note=(
                    "Source of the OL reference loci (NCBI BioProject PRJDB8600;"
                    " LC594800-LC594991 subset used by this panel)."
                ),
            ),
            DbReference(
                citation=(
                    "Nandi B, Nandy RK, Mukhopadhyay S, Nair GB, Shimada T,"
                    " Ghose AC. Rapid method for species-specific identification"
                    " of Vibrio cholerae using primers targeted to the gene of"
                    " outer membrane protein OmpW. J Clin Microbiol"
                    " 2000;38(11):4145-4151."
                ),
                doi="10.1128/JCM.38.11.4145-4151.2000",
                note=(
                    "ompW species-control marker; control locus is the N16961"
                    " ompW CDS (NC_002506.1:819609-820262)."
                ),
            ),
        ),
        caveats=(
            "The OL reference set contains many clusters that have not been"
            " serologically verified. An OL/O call is a locus match, not a"
            " confirmed O serogroup; confirm with agglutination or a validated"
            " serotyping scheme.",
            "OLU (untypeable/unassigned) loci are not yet included. An OUT call"
            " means no reference locus matched, not that the isolate is"
            " non-O1/non-O139; additional OLU references will be added from"
            " public databases.",
            "OL1-Inaba and OL1-Ogawa both map to O1 in this build; the panel"
            " does not resolve Inaba vs Ogawa (wbeT phase variation) or other"
            " subserotypes.",
            "O-antigen loci can recombine and assemblies can be partial. A weak"
            " or partial locus match should be confirmed; coverage/identity and"
            " the cluster evidence columns are part of the call.",
            "ctxA/ctxB are deliberately absent from this database; they belong"
            " to the separate vch_vf panel.",
            "ompW is a species-control marker, not an absolute species or"
            " serogroup determinant; interpret it with the isolate context.",
        ),
    ),
)
