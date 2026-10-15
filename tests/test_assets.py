from scpagwas.io import assets


def test_builtin_pathways_kegg():
    pw = assets.load_builtin_pathways("kegg")
    assert len(pw) > 100
    assert all(isinstance(v, list) for v in pw.values())


def test_builtin_block_annotation():
    anno = assets.load_builtin_block_annotation("hg38")
    assert list(anno.columns) == ["chrom", "start", "end", "label"]
    assert len(anno) > 50000
    assert 23 in set(anno["chrom"])  # chrX mapped to 23


def test_example_gwas_path_exists():
    assert assets.example_gwas_path().exists()


def test_builtin_ld_one_chromosome():
    ld = assets.load_builtin_ld(chromosomes=[22])
    assert ld.chromosomes == [22]
    assert list(ld.get(22).columns) == ["SNP_A", "SNP_B", "R"]
