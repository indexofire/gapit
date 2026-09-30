"""Database providers: the Provider contracts and per-provider modules (Wave B).

Each concrete provider is a module exporting a ``Provider`` (gene pipeline,
providers.common) or ``ClusterProvider`` (kaptive-style GBK providers,
providers.cluster_common) value consumed by ``fetch_provider`` /
``fetch_cluster_provider``. ``REGISTRY`` maps every provider name to its
value — the CLI's ``gapit db fetch`` / ``gapit db list`` lookup table (Wave
C-a). Imports are static and alphabetized; no dynamic import magic.
"""

from gapit.providers import (
    ab_k,
    ab_o,
    argannot,
    bacmet2,
    card,
    ecoh,
    ecoli_kps,
    ecoli_vf,
    kosc_k,
    kosc_o,
    kpsc_k,
    kpsc_o,
    megares,
    ncbi,
    plasmidfinder,
    resfinder,
    upec_expec_vf,
    vfdb,
    victors,
)
from gapit.providers.cluster_common import ClusterProvider
from gapit.providers.common import Provider

AnyProvider = Provider | ClusterProvider

REGISTRY: dict[str, AnyProvider] = {
    ab_k.PROVIDER.name: ab_k.PROVIDER,
    ab_o.PROVIDER.name: ab_o.PROVIDER,
    argannot.PROVIDER.name: argannot.PROVIDER,
    bacmet2.PROVIDER.name: bacmet2.PROVIDER,
    card.PROVIDER.name: card.PROVIDER,
    ecoh.PROVIDER.name: ecoh.PROVIDER,
    ecoli_kps.PROVIDER.name: ecoli_kps.PROVIDER,
    ecoli_vf.PROVIDER.name: ecoli_vf.PROVIDER,
    kosc_k.PROVIDER.name: kosc_k.PROVIDER,
    kosc_o.PROVIDER.name: kosc_o.PROVIDER,
    kpsc_k.PROVIDER.name: kpsc_k.PROVIDER,
    kpsc_o.PROVIDER.name: kpsc_o.PROVIDER,
    megares.PROVIDER.name: megares.PROVIDER,
    ncbi.PROVIDER.name: ncbi.PROVIDER,
    plasmidfinder.PROVIDER.name: plasmidfinder.PROVIDER,
    resfinder.PROVIDER.name: resfinder.PROVIDER,
    upec_expec_vf.PROVIDER.name: upec_expec_vf.PROVIDER,
    vfdb.PROVIDER.name: vfdb.PROVIDER,
    victors.PROVIDER.name: victors.PROVIDER,
}

__all__ = ["REGISTRY", "AnyProvider", "ClusterProvider", "Provider"]
