import numpy as np
import pandas as pd

from scpagwas.io.ld import LDIndex
from scpagwas.link import link_pathway, make_ld_matrix


def test_make_ld_matrix_symmetric():
    ld = pd.DataFrame({"SNP_A": ["a", "b"], "SNP_B": ["b", "c"], "R": [0.5, -0.2]})
    M = make_ld_matrix(["a", "b", "c"], ld)
    assert M.shape == (3, 3)
    assert np.allclose(M, M.T)
    assert np.allclose(np.diag(M), 1.0)
    assert M[0, 1] == 0.5 and M[1, 2] == -0.2


def test_ld_index_missing_chrom_returns_empty():
    idx = LDIndex.from_frames({})
    assert list(idx.get(1).columns) == ["SNP_A", "SNP_B", "R"]


def test_ld_index_restrict_keeps_only_internal_pairs():
    idx = LDIndex.from_frames({1: pd.DataFrame({
        "SNP_A": ["a", "b", "c"], "SNP_B": ["b", "c", "d"], "R": [0.1, 0.2, 0.3]})})
    restricted = idx.restrict({"a", "b"})
    out = restricted.get(1)
    assert len(out) == 1
    assert out.iloc[0]["SNP_A"] == "a" and out.iloc[0]["SNP_B"] == "b"


def test_link_pathway_builds_y_and_design():
    gwas_by_chrom = {1: pd.DataFrame({
        "chrom": [1, 1], "pos": [100, 200], "rsid": ["rs1", "rs2"],
        "beta": [0.2, 0.3], "se": [0.1, 0.1], "maf": [0.2, 0.2]})}
    snp_gene = pd.DataFrame({"rsid": ["rs1", "rs2"], "label": ["G1", "G2"],
                             "pos": [100, 200]})
    anno = pd.DataFrame({"chrom": [1, 1], "start": [100, 200], "end": [150, 250],
                         "label": ["G1", "G2"]})
    ld = LDIndex.from_frames({1: pd.DataFrame(
        {"SNP_A": ["rs1"], "SNP_B": ["rs2"], "R": [0.3]})})
    pl = link_pathway("P", {1: anno}, gwas_by_chrom, snp_gene, ld)
    assert pl is not None
    assert np.allclose(pl.y, [0.04, 0.09])
    assert pl.snps.shape[0] == 2
    assert np.allclose(pl.ldsq[0, 1], 0.09)


def test_link_pathway_none_when_no_snps():
    pl = link_pathway("P", {1: pd.DataFrame({"chrom": [1], "start": [1],
                                            "end": [2], "label": ["X"]})},
                      {}, pd.DataFrame({"rsid": [], "label": [], "pos": []}),
                      LDIndex.from_frames({}))
    assert pl is None
