import pandas as pd

from scpagwas.io.annotation import load_block_annotation


def _anno():
    return pd.DataFrame({
        "chrom": ["chr1", "1", "2"],
        "start": [100, 200, 300],
        "end": [150, 250, 350],
        "label": ["A", "B", "C"],
        "extra": [1, 2, 3],
    })


def test_load_from_tsv(tmp_path):
    p = tmp_path / "anno.tsv"
    _anno().to_csv(p, sep="\t", index=False)
    out = load_block_annotation(p)
    assert list(out.columns) == ["chrom", "start", "end", "label"]
    assert out["chrom"].dtype.kind in "iu"
    assert set(out["chrom"]) == {1, 2}


def test_load_from_parquet(tmp_path):
    p = tmp_path / "anno.parquet"
    _anno().to_parquet(p)
    out = load_block_annotation(p)
    assert len(out) == 3


def test_load_handles_sex_and_mito_chroms(tmp_path):
    df = pd.DataFrame({
        "chrom": ["chrX", "chrY", "chrM", "chr1"],
        "start": [1, 2, 3, 4],
        "end": [4, 5, 6, 7],
        "label": ["A", "B", "C", "D"],
    })
    p = tmp_path / "anno.tsv"
    df.to_csv(p, sep="\t", index=False)
    out = load_block_annotation(p)
    assert dict(zip(out["label"], out["chrom"])) == {"A": 23, "B": 24, "C": 26, "D": 1}


def test_load_reads_gzipped_tsv(tmp_path):
    p = tmp_path / "anno.tsv.gz"
    _anno().to_csv(p, sep="\t", index=False, compression="gzip")
    out = load_block_annotation(p)
    assert len(out) == 3
