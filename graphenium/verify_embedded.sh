#!/usr/bin/env bash
# Verify the MUB-6 graph obstruction with Graphenium's embedded Datalog engine.
#
# Usage:
#   ./graphenium/verify_embedded.sh
#   GM=/path/to/gm ./graphenium/verify_embedded.sh
#
# Requires `gm` on PATH (or GM set), plus Python 3 for the fact -> graph bridge.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
gm="${GM:-gm}"
out="$root/graphenium-out"
program="$root/graphenium/mub_cliques.dl"
python="${PYTHON:-python3}"

if ! command -v "$gm" >/dev/null 2>&1; then
  echo "error: Graphenium CLI '$gm' not found; install it or set GM=/path/to/gm" >&2
  exit 1
fi

mkdir -p "$out"

echo "[MUB-6] building Graphenium graphs from exact facts..."
"$python" "$root/graphenium/build_graph.py" all

datalog="$(cat "$program")"

run_datalog() {
  local graph="$1"
  local log="$2"
  "$gm" query "" --datalog "$datalog" --graph "$graph" | tee "$log"
}

echo "[MUB-6] running embedded Datalog on F6..."
f6_log="$out/f6_datalog.txt"
run_datalog "$out/f6_graph.json" "$f6_log" >/dev/null
f6_cliques="$(grep -cE '^\| [0-9]+ \|' "$f6_log" || true)"
if ! grep -q 'Results for `clique6`' "$f6_log"; then
  echo "error: F6 Datalog run did not report any 6-cliques" >&2
  exit 1
fi
if [ "$f6_cliques" -ne 16 ]; then
  echo "error: expected 16 F6 bases, got $f6_cliques" >&2
  exit 1
fi
if grep -q 'Results for `mub_clique_pair`' "$f6_log"; then
  echo "error: unexpected F6 MUB pair between bases" >&2
  exit 1
fi
echo "[MUB-6] F6: $f6_cliques bases, 0 MUB pairs between bases"

echo "[MUB-6] running embedded Datalog on S6..."
s6_log="$out/s6_datalog.txt"
run_datalog "$out/s6_graph.json" "$s6_log" >/dev/null
if ! grep -q '(no results)' "$s6_log"; then
  echo "error: expected no S6 basis and no S6 MUB pair" >&2
  exit 1
fi
echo "[MUB-6] S6: 0 bases, 0 orthogonal pairs, 0 MUB pairs"

# Keep generated artifacts portable: no developer-machine absolute paths.
if grep -RInE '/Users/|/home/|/root/|[A-Za-z]:[\\/]' "$out" \
    --include='*.json' --include='*.txt' --include='*.log' >/dev/null 2>&1; then
  echo "error: generated artifacts contain absolute paths" >&2
  exit 1
fi
echo "[MUB-6] path sanitation check passed"

echo "[MUB-6] embedded Datalog verification passed."
