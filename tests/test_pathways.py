from scpagwas.io.pathways import read_gmt, write_gmt, example_pathways


def test_gmt_roundtrip(tmp_path):
    pw = {"p1": ["A", "B"], "p2": ["C"]}
    f = tmp_path / "x.gmt"
    write_gmt(pw, f)
    assert read_gmt(f) == pw


def test_example_pathways_nonempty():
    pw = example_pathways()
    assert isinstance(pw, dict) and len(pw) >= 3
