import numpy as np

from scpagwas.config import PagwasConfig
from scpagwas.utils import winsorize, storey_qvalue, marchenko_pastur_edge


def test_config_defaults():
    cfg = PagwasConfig()
    assert cfg.min_pathway_size == 5
    assert cfg.n_topgenes == 1000
    assert cfg.gwas_z_filter is None


def test_winsorize_clips_and_keeps_middle():
    x = np.array([-100.0, 0.0, 1.0, 2.0, 100.0])
    out = winsorize(x, 0.2, 0.8)
    assert out.min() > -100.0 and out.max() < 100.0
    assert out[2] == 1.0


def test_winsorize_nan_to_number():
    out = winsorize(np.array([np.nan, 1.0, 2.0, 3.0]))
    assert not np.isnan(out).any()


def test_storey_qvalue_bounded_and_pi0():
    rng = np.random.default_rng(0)
    p = np.concatenate([rng.uniform(0, 0.05, 50), rng.uniform(0, 1, 950)])
    q, pi0 = storey_qvalue(p)
    assert q.shape == p.shape
    assert 0.0 <= q.min() and q.max() <= 1.0
    assert 0.0 < pi0 <= 1.0


def test_mp_edge_grows_with_genes():
    assert marchenko_pastur_edge(1000, 200) < marchenko_pastur_edge(1000, 800)
