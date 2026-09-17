"""Database providers: the Provider contract and per-provider modules (Wave B).

Each concrete provider (Wave B1..B12) is a module exporting a ``Provider``
value consumed by :func:`gapit.providers.common.fetch_provider`.
``REGISTRY`` maps every provider name to its ``Provider`` — the CLI's
``gapit db fetch`` / ``gapit db list`` lookup table (Wave C-a). Imports are
static and alphabetized; no dynamic import magic.
"""

from gapit.providers import (
    argannot,
    bacmet2,
    card,
    ecoh,
    ecoli_vf,
    megares,
    ncbi,
    plasmidfinder,
    resfinder,
    upec_expec_vf,
    vfdb,
    victors,
)
from gapit.providers.common import Provider

REGISTRY: dict[str, Provider] = {
    argannot.PROVIDER.name: argannot.PROVIDER,
    bacmet2.PROVIDER.name: bacmet2.PROVIDER,
    card.PROVIDER.name: card.PROVIDER,
    ecoh.PROVIDER.name: ecoh.PROVIDER,
    ecoli_vf.PROVIDER.name: ecoli_vf.PROVIDER,
    megares.PROVIDER.name: megares.PROVIDER,
    ncbi.PROVIDER.name: ncbi.PROVIDER,
    plasmidfinder.PROVIDER.name: plasmidfinder.PROVIDER,
    resfinder.PROVIDER.name: resfinder.PROVIDER,
    upec_expec_vf.PROVIDER.name: upec_expec_vf.PROVIDER,
    vfdb.PROVIDER.name: vfdb.PROVIDER,
    victors.PROVIDER.name: victors.PROVIDER,
}
