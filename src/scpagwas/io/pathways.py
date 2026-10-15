"""Pathway gene-set input/output (GMT)."""

from pathlib import Path

_EXAMPLE = {
    "EXAMPLE_TCELL_A": [f"GENEA{i}" for i in range(1, 21)],
    "EXAMPLE_TCELL_B": [f"GENEB{i}" for i in range(1, 21)],
    "EXAMPLE_MYELOID": [f"GENEC{i}" for i in range(1, 21)],
}


def read_gmt(path):
    """Read a GMT file into ``{pathway: [gene, ...]}``."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        out[parts[0]] = [g for g in parts[2:] if g]
    return out


def write_gmt(pathways, path):
    """Write ``{pathway: [gene, ...]}`` to a GMT file."""
    lines = [f"{n}\tna\t" + "\t".join(g) for n, g in pathways.items()]
    Path(path).write_text("\n".join(lines) + "\n")


def example_pathways():
    """A tiny synthetic pathway set for demos and tests."""
    return {k: list(v) for k, v in _EXAMPLE.items()}
