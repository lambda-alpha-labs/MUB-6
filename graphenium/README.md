# Graphenium embedded Datalog layer

This directory contains the Graphenium-native version of the MUB-6 clique
program. It replaces the optional Soufflé reference path: no external Datalog
engine is required.

## Files

| File | Purpose |
|---|---|
| `build_graph.py` | Converts `data/{f6,s6}_{ortho,mub}.facts` into small Graphenium graph JSON files under `graphenium-out/`. |
| `mub_cliques.dl` | Graphenium Datalog rules for all 6-cliques and complete bipartite `K_{6,6}` subgraphs. |
| `verify_embedded.sh` | Builds the graphs, runs both datasets through `gm`, and asserts 16 F6 bases with zero MUB pairs and an empty S6 obstruction. |

## Relation encoding

Graphenium's Datalog EDB is derived from a graph. The bridge maps the MUB-6
relations onto existing EDB predicates:

| Mathematical relation | Graphenium relation | Direction |
|---|---|---|
| orthogonality | `calls` | both directions |
| mutual unbiasedness | `imports` | both directions |
| canonical order `i < j` | `contains` | one direction |

Using the existing EDB predicates keeps the graph indexed by relation, which
makes the 6-clique join fast enough for ordinary CI. All `source_file` values
in the generated graph are relative paths such as `data/f6_ortho.facts`.

## Run

```sh
# From the repository root:
make datalog

# Or directly:
./graphenium/verify_embedded.sh
```

The script requires `gm` on `PATH`, or `GM=/path/to/gm`. The top-level
`make verify` target runs this layer along with the Python generator, the C
exhaustive search, and the Lean 4 build.
