import numpy as np
import pandas as pd

from scpagwas.config import PagwasConfig
from scpagwas.io.gwas import clean_gwas


def _raw():
    return pd.DataFrame({
        "chrom": ["1", "1", "chr6", "X", "1"],
        "pos": [1000, 30000000, 30000000, 500, 1000],
        "rsid": ["rs1", "rs2", "rs3", "rs4", "rs1"],
        "beta": [0.2, 0.5, 0.1, 0.3, 0.2],
        "se": [0.1, 0.1, 0.1, 0.1, 0.1],
        "maf": [0.2, 0.1, 0.3, 0.4, 0.7],
    })


def test_clean_gwas_dedup_and_filters():
    out = clean_gwas(_raw(), PagwasConfig())
    assert out["rsid"].is_unique
    assert "rs4" not in set(out["rsid"])
    assert "rs3" not in set(out["rsid"])
    assert out["chrom"].dtype.kind in "iu"
    assert (out["maf"] <= 0.5).all()


def test_clean_gwas_maf_flip():
    out = clean_gwas(_raw(), PagwasConfig(mhc_filter=False, sex_filter=False))
    assert 0.3 in set(np.round(out["maf"], 6))


def test_clean_gwas_drops_unknown_contigs():
    raw = _raw()
    raw.loc[len(raw)] = ["GL000191.1", 500, "rs9", 0.1, 0.1, 0.2]
    out = clean_gwas(raw, PagwasConfig(sex_filter=False))
    assert "rs9" not in set(out["rsid"])
