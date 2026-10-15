import numpy as np
import pandas as pd

from scpagwas.genes import gene_pcc


def test_gene_pcc_matches_numpy():
    genes, cells = pd.Index(["g1", "g2"]), pd.Index(["c1", "c2", "c3"])
    data = np.array([[1.0, 2.0, 3.0], [3.0, 1.0, 2.0]])
    gpas = pd.Series([1.0, 2.0, 3.0], index=cells)
    pcc = gene_pcc(gpas, data, genes, cells)
    assert np.isclose(pcc.loc["g1"], np.corrcoef(data[0], gpas)[0, 1])
    assert np.isclose(pcc.loc["g2"], np.corrcoef(data[1], gpas)[0, 1])
