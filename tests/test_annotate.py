import pandas as pd

from scpagwas.annotate import build_pathway_blocks, map_snps_to_genes


def _gwas():
    return pd.DataFrame({
        "chrom": [1, 1, 2],
        "pos": [1500, 9000, 100],
        "rsid": ["rs1", "rs2", "rs3"],
        "beta": [0.1, 0.2, 0.3],
        "se": [0.1, 0.1, 0.1],
        "maf": [0.2, 0.2, 0.2],
    })


def _anno():
    return pd.DataFrame({
        "chrom": [1, 1],
        "start": [2000, 12000],
        "end": [3000, 13000],
        "label": ["GENEA", "GENEB"],
    })


def test_map_window():
    sg = map_snps_to_genes(_gwas(), _anno(), marg=10000)
    assert set(sg.columns) == {"rsid", "label", "pos"}
    assert set(sg.loc[sg["rsid"] == "rs1", "label"]) == {"GENEA"}
    assert "rs3" not in set(sg["rsid"])  # chromosome 2 has no annotation


def test_map_snp_in_multiple_windows():
    sg = map_snps_to_genes(_gwas(), _anno(), marg=10000)
    assert set(sg.loc[sg["rsid"] == "rs2", "label"]) == {"GENEA", "GENEB"}


def test_build_blocks_min():
    sg = map_snps_to_genes(_gwas(), _anno(), marg=10000)
    blocks = build_pathway_blocks(
        {"P": ["GENEA", "GENEB"], "Q": ["GENEA"]}, _anno(), sg, min_blocks=2
    )
    assert "P" in blocks and "Q" not in blocks
    assert 1 in blocks["P"]
