"""Database providers: the Provider contracts and per-provider modules (Wave B).

Each concrete provider is a module exporting a ``Provider`` (gene pipeline,
providers.common) or ``ClusterProvider`` (kaptive-style GBK providers,
providers.cluster_common) value consumed by ``fetch_provider`` /
``fetch_cluster_provider``. ``REGISTRY`` maps every provider name to its
value — the CLI's ``gapit db fetch`` / ``gapit db list`` lookup table (Wave
C-a). Imports are static and alphabetized; no dynamic import magic.
"""

from gapit.providers import (
    argannot,
    bacmet2,
    card,
    ecoh,
    ecoli_vf,
    kaptive_ak,
    kaptive_k,
    kaptive_o,
    kaptive_oc,
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
    argannot.PROVIDER.name: argannot.PROVIDER,
    bacmet2.PROVIDER.name: bacmet2.PROVIDER,
    card.PROVIDER.name: card.PROVIDER,
    ecoh.PROVIDER.name: ecoh.PROVIDER,
    ecoli_vf.PROVIDER.name: ecoli_vf.PROVIDER,
    kaptive_ak.PROVIDER.name: kaptive_ak.PROVIDER,
    kaptive_k.PROVIDER.name: kaptive_k.PROVIDER,
    kaptive_o.PROVIDER.name: kaptive_o.PROVIDER,
    kaptive_oc.PROVIDER.name: kaptive_oc.PROVIDER,
    megares.PROVIDER.name: megares.PROVIDER,
    ncbi.PROVIDER.name: ncbi.PROVIDER,
    plasmidfinder.PROVIDER.name: plasmidfinder.PROVIDER,
    resfinder.PROVIDER.name: resfinder.PROVIDER,
    upec_expec_vf.PROVIDER.name: upec_expec_vf.PROVIDER,
    vfdb.PROVIDER.name: vfdb.PROVIDER,
    victors.PROVIDER.name: victors.PROVIDER,
}

__all__ = ["REGISTRY", "AnyProvider", "ClusterProvider", "Provider"]
