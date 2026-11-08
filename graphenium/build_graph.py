#!/usr/bin/env python3
"""Build Graphenium graph JSON for the MUB-6 embedded Datalog layer.

The MUB-6 facts are plain relational tables (`ortho(i,j)` and `mub(i,j)`),
while Graphenium's Datalog EDB is derived from a graph.  This bridge creates a
small synthetic graph whose nodes are the vector indices and whose edges encode
the three relations used by `graphenium/mub_cliques.dl`:

    orthogonality        -> `calls`    (mirrored in both directions)
    mutual unbiasedness   -> `imports`  (mirrored in both directions)
    canonical order i < j -> `contains` (one direction only)

Graphenium's embedded EDB exposes `calls`, `imports`, and `contains` directly,
which keeps the Datalog evaluation fast while still running the clique logic
entirely inside `gm`.  All paths written into the graph are relative to the
repository root.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "graphenium-out"

VECTOR_COUNTS = {"f6": 48, "s6": 90}


def read_pairs(path: Path) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) != 2:
            raise SystemExit(f"{path}:{line_number}: expected two integers, got {raw!r}")
        a, b = map(int, parts)
        pairs.append((min(a, b), max(a, b)))
    return sorted(set(pairs))


def build(dataset: str) -> Path:
    expected = VECTOR_COUNTS[dataset]
    ortho = read_pairs(DATA / f"{dataset}_ortho.facts")
    mub = read_pairs(DATA / f"{dataset}_mub.facts")

    max_index = max((i for pair in ortho + mub for i in pair), default=-1)
    if max_index >= expected:
        raise SystemExit(
            f"{dataset}: vector index {max_index} exceeds expected range 0..{expected - 1}"
        )

    nodes = [
        {
            "id": f"v{i}",
            "label": f"v{i}",
            "file_type": "code",
            "source_file": f"data/{dataset}_mub.facts",
            "source_location": "",
            "community": 0,
        }
        for i in range(expected)
    ]

    links: list[dict[str, object]] = []

    def add_edge(a: int, b: int, relation: str, source_file: str) -> None:
        links.append(
            {
                "source": f"v{a}",
                "target": f"v{b}",
                "relation": relation,
                "confidence": "EXTRACTED",
                "confidence_score": 1.0,
                "weight": 1.0,
                "source_file": source_file,
            }
        )

    ortho_file = f"data/{dataset}_ortho.facts"
    mub_file = f"data/{dataset}_mub.facts"
    for a, b in ortho:
        add_edge(a, b, "calls", ortho_file)
        add_edge(b, a, "calls", ortho_file)
    for a, b in mub:
        add_edge(a, b, "imports", mub_file)
        add_edge(b, a, "imports", mub_file)
    for a in range(expected):
        for b in range(a + 1, expected):
            add_edge(a, b, "contains", "graphenium/order")

    graph = {
        "nodes": nodes,
        "links": links,
        "hyperedges": [],
        "metadata": {
            "ast_only": True,
            "schema_version": "0.2.0",
            "project_root": ".",
            "extraction_modes": ["mub-facts"],
            "languages": ["mub-6"],
        },
    }

    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / f"{dataset}_graph.json"
    target.write_text(json.dumps(graph, indent=2) + "\n", encoding="utf-8")
    print(
        f"{dataset}: {len(nodes)} vectors, {len(ortho)} orthogonal pairs, "
        f"{len(mub)} MUB pairs, {len(links)} graph edges -> {target.relative_to(ROOT)}"
    )
    return target


def main(argv: list[str]) -> int:
    datasets = argv[1:] or ["f6", "s6"]
    if datasets == ["all"]:
        datasets = ["f6", "s6"]
    for dataset in datasets:
        if dataset not in VECTOR_COUNTS:
            raise SystemExit(f"unknown dataset {dataset!r}; expected f6, s6, or all")
        build(dataset)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
