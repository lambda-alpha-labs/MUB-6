#!/usr/bin/env python3
"""
Generate exact facts and Lean certificates for the MUB-6 showcase.

Two data sets are produced:

F6 (the Fourier matrix):
  48 biunimodular vectors, i.e. vectors unbiased with respect to the
  standard basis and F6.  The vectors are generated from Bjorck's
  algebraic construction.  The alphabet is

      A = { d^m * zeta^k : -2 <= m <= 2, 0 <= k < 12 },
      d = (1 - sqrt(3))/2 + i*sqrt(sqrt(3)/2),  zeta = exp(2*pi*i/12).

  The 48 code words below are the complete set of length-6 codes whose
  DFT is again unimodular; they were obtained by exhaustive enumeration
  of A^6 / global phase (with the first coordinate fixed to 1).  Every
  code gives one vector; the Fourier transform permutes the set.

S6 (Tao's spectral matrix):
  90 vectors from the 18 seeds supplied by the project.  The arithmetic
  here is exact integer arithmetic in the order

      Z[sqrt(5), i*sqrt(3), i*sqrt(15)],

  scaled by 8.  The U9 table is the corrected table (the two entries
  involving omega^2 * beta and omega^2 * beta-bar had the i*sqrt(15)
  sign flipped in the original hand-written table).

The script writes:
  data/f6_ortho.facts, data/f6_mub.facts
  data/s6_ortho.facts, data/s6_mub.facts
  lean4/MUB6/Certificates.lean
  data/summary.json

It also asserts every structural property used by the Lean layer.  The
F6 graph is computed in double precision and checked against a generous
margin; the S6 graph is computed in exact integer arithmetic.
"""

from __future__ import annotations

import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
LEAN = ROOT / "lean4" / "MUB6"

# ---------------------------------------------------------------------------
# Exact arithmetic in Z[sqrt(5), i*sqrt(3), i*sqrt(15)]
# ---------------------------------------------------------------------------

Q = tuple[int, int, int, int]
ZERO: Q = (0, 0, 0, 0)
ONE: Q = (8, 0, 0, 0)
OMEGA: Q = (-4, 0, 4, 0)
OMEGA2: Q = (-4, 0, -4, 0)
UNIT_NORM: tuple[int, int] = (64, 0)
VECTOR_NORM: tuple[int, int] = (384, 0)
UNBIASED_NORM: tuple[int, int] = (24576, 0)


def qadd(x: Q, y: Q) -> Q:
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3])


def qconj(x: Q) -> Q:
    return (x[0], x[1], -x[2], -x[3])


def qmul(x: Q, y: Q) -> Q:
    a1, b1, c1, d1 = x
    a2, b2, c2, d2 = y
    return (
        a1 * a2 + 5 * b1 * b2 - 3 * c1 * c2 - 15 * d1 * d2,
        a1 * b2 + b1 * a2 - 3 * c1 * d2 - 3 * d1 * c2,
        a1 * c2 + 5 * b1 * d2 + c1 * a2 + 5 * d1 * b2,
        a1 * d2 + b1 * c2 + c1 * b2 + d1 * a2,
    )


def qnorm(x: Q) -> tuple[int, int]:
    a, b, c, d = x
    return ((a * a + 5 * b * b) + 3 * (c * c + 5 * d * d), 2 * a * b + 6 * c * d)


# Corrected 8 * U9.  U9 = {1, omega, omega^2} * {1, beta, beta-bar},
# omega = (-1 + i sqrt(3))/2, beta = (1 + i sqrt(15))/4.
U9: list[Q] = [
    (8, 0, 0, 0),      # 1
    (2, 0, 0, 2),      # beta
    (2, 0, 0, -2),     # beta-bar
    (-4, 0, 4, 0),     # omega
    (-1, -3, 1, -1),   # omega * beta
    (-1, 3, 1, 1),     # omega * beta-bar
    (-4, 0, -4, 0),    # omega^2
    (-1, 3, -1, -1),   # omega^2 * beta       (corrected)
    (-1, -3, -1, 1),   # omega^2 * beta-bar   (corrected)
]

# Provided S6 seeds.  Each seed has five U9 indices.  Its five cyclic
# shifts, each prefixed by index 0 (the scalar 1), give five vectors.
S6_SEEDS: list[list[int]] = [
    [0, 1, 3, 3, 4], [0, 2, 3, 3, 5], [0, 4, 3, 3, 1], [0, 5, 3, 3, 2],
    [0, 6, 1, 7, 6], [0, 6, 2, 8, 6], [0, 6, 4, 4, 6], [0, 6, 5, 5, 6],
    [0, 6, 7, 1, 6], [0, 6, 8, 2, 6], [0, 7, 3, 3, 7], [0, 7, 4, 4, 7],
    [0, 8, 3, 3, 8], [0, 8, 5, 5, 8], [1, 1, 4, 3, 4], [1, 6, 1, 7, 7],
    [2, 2, 5, 3, 5], [2, 6, 2, 8, 8],
]

# The row permutation of the standard spectral matrix that matches the
# supplied seeds.  The standard (thesis/paper) matrix is S; we use P*S
# with P = (0,1,2,5,4,3).  Any of the 60 equivalent working row orders
# would give the same combinatorial obstruction.
S6_ROW_ORDER = (0, 1, 2, 5, 4, 3)


def build_s6_vectors() -> list[list[Q]]:
    vectors: list[list[Q]] = []
    for seed in S6_SEEDS:
        for shift in range(5):
            code = [0] + [seed[(shift + k) % 5] for k in range(5)]
            vectors.append([U9[i] for i in code])
    return vectors


def build_s6_matrix() -> list[list[Q]]:
    one, w, w2 = U9[0], U9[3], U9[6]
    standard = [
        [one, one, one, one, one, one],
        [one, one, w, w, w2, w2],
        [one, w, one, w2, w2, w],
        [one, w, w2, one, w, w2],
        [one, w2, w2, w, one, w],
        [one, w2, w, w2, w, one],
    ]
    return [standard[i] for i in S6_ROW_ORDER]


def dot_exact(u: Sequence[Q], v: Sequence[Q]) -> Q:
    acc = ZERO
    for x, y in zip(u, v):
        acc = qadd(acc, qmul(x, qconj(y)))
    return acc


# ---------------------------------------------------------------------------
# F6: Bjorck biunimodular construction
# ---------------------------------------------------------------------------

# Complete set of 48 length-6 codes over the 60-element alphabet
#   A[(m+2)*12 + k] = d^m * zeta^k,  -2 <= m <= 2, 0 <= k < 12.
# The first coordinate is the scalar 1 (index 24).  These codes were
# found by exhaustive search and are listed here so the script is
# deterministic and does not depend on a search order.
F6_CODES: list[list[int]] = [
    [24, 13, 8, 3, 22, 35], [24, 13, 14, 27, 46, 41],
    [24, 15, 0, 9, 18, 33], [24, 15, 18, 33, 42, 39],
    [24, 17, 4, 3, 14, 31], [24, 17, 22, 27, 38, 37],
    [24, 19, 8, 9, 22, 29], [24, 19, 14, 33, 46, 47],
    [24, 21, 0, 3, 18, 27], [24, 21, 18, 27, 42, 45],
    [24, 23, 4, 9, 14, 25], [24, 23, 22, 33, 38, 43],
    [24, 25, 14, 9, 4, 23], [24, 25, 24, 33, 28, 33],
    [24, 25, 28, 33, 28, 25], [24, 25, 38, 57, 52, 47],
    [24, 27, 18, 3, 0, 21], [24, 27, 28, 27, 24, 31],
    [24, 27, 32, 27, 24, 35], [24, 27, 42, 51, 48, 45],
    [24, 29, 22, 9, 8, 19], [24, 29, 24, 33, 32, 33],
    [24, 29, 32, 33, 32, 29], [24, 29, 46, 57, 56, 43],
    [24, 31, 14, 3, 4, 17], [24, 31, 24, 27, 28, 27],
    [24, 31, 28, 27, 28, 31], [24, 31, 38, 51, 52, 41],
    [24, 33, 18, 9, 0, 15], [24, 33, 28, 33, 24, 25],
    [24, 33, 32, 33, 24, 29], [24, 33, 42, 57, 48, 39],
    [24, 35, 22, 3, 8, 13], [24, 35, 24, 27, 32, 27],
    [24, 35, 32, 27, 32, 35], [24, 35, 46, 51, 56, 37],
    [24, 37, 38, 27, 22, 17], [24, 37, 56, 51, 46, 35],
    [24, 39, 42, 33, 18, 15], [24, 39, 48, 57, 42, 33],
    [24, 41, 46, 27, 14, 13], [24, 41, 52, 51, 38, 31],
    [24, 43, 38, 33, 22, 23], [24, 43, 56, 57, 46, 29],
    [24, 45, 42, 27, 18, 21], [24, 45, 48, 51, 42, 27],
    [24, 47, 46, 33, 14, 19], [24, 47, 52, 57, 38, 25],
]


def f6_alphabet() -> list[complex]:
    zeta = cmath.exp(2j * cmath.pi / 12)
    d = (1 - math.sqrt(3)) / 2 + 1j * math.sqrt(math.sqrt(3) / 2)
    alphabet: list[complex] = []
    for m in range(-2, 3):
        for k in range(12):
            alphabet.append((d ** m) * (zeta ** k))
    return alphabet


def build_f6_vectors() -> list[list[complex]]:
    A = f6_alphabet()
    return [[A[i] for i in code] for code in F6_CODES]


def f6_matrix() -> list[list[complex]]:
    q = cmath.exp(2j * cmath.pi / 6)
    return [[q ** (i * j) for j in range(6)] for i in range(6)]


def dot_float(u: Sequence[complex], v: Sequence[complex]) -> complex:
    return sum(x * y.conjugate() for x, y in zip(u, v))


# ---------------------------------------------------------------------------
# Graphs, cliques, certificates
# ---------------------------------------------------------------------------


def cliques6(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    found: set[tuple[int, ...]] = set()

    def extend(cur: list[int], cand: set[int]) -> None:
        if len(cur) == 6:
            found.add(tuple(sorted(cur)))
            return
        if len(cand) < 6 - len(cur):
            return
        ordered = sorted(cand)
        for pos, v in enumerate(ordered):
            nxt = set(ordered[pos + 1:])
            nxt &= {x for x in range(n) if adj[v][x]}
            extend(cur + [v], nxt)

    for i in range(n):
        extend([i], {x for x in range(n) if adj[i][x] and x > i})
    return sorted(found)


def write_facts(path: Path, pairs: Iterable[tuple[int, int]]) -> int:
    rows = list(pairs)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for i, j in rows:
            handle.write(f"{i}\t{j}\n")
    return len(rows)


def mask_words(n: int, pred: Callable[[int, int], bool]) -> list[int]:
    words = [0] * ((n * n + 63) // 64)
    for i in range(n):
        for j in range(n):
            if pred(i, j):
                bit = i * n + j
                words[bit // 64] |= 1 << (bit % 64)
    return words


def lean_array(words: Sequence[int]) -> str:
    return "#[" + ", ".join(str(w) for w in words) + "]"


def lean_list_of_lists(items: Sequence[Sequence[int]]) -> str:
    return "[" + ",\n   ".join("[" + ", ".join(str(x) for x in row) + "]" for row in items) + "]"


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    LEAN.mkdir(parents=True, exist_ok=True)

    # --- S6 exact data -----------------------------------------------------
    s6_vectors = build_s6_vectors()
    s6_matrix = build_s6_matrix()
    assert len(s6_vectors) == 90
    assert all(len(v) == 6 for v in s6_vectors)
    assert all(qnorm(entry) == UNIT_NORM for v in s6_vectors for entry in v)
    assert all(
        qnorm(dot_exact(v, col)) == UNBIASED_NORM
        for v in s6_vectors
        for col in ([s6_matrix[r][c] for r in range(6)] for c in range(6))
    )

    s6_ortho: list[list[bool]] = [[False] * 90 for _ in range(90)]
    s6_mub: list[list[bool]] = [[False] * 90 for _ in range(90)]
    s6_ortho_pairs: list[tuple[int, int]] = []
    s6_mub_pairs: list[tuple[int, int]] = []
    for i in range(90):
        for j in range(i + 1, 90):
            n = qnorm(dot_exact(s6_vectors[i], s6_vectors[j]))
            if n == (0, 0):
                s6_ortho[i][j] = s6_ortho[j][i] = True
                s6_ortho_pairs.append((i, j))
            if n == UNBIASED_NORM:
                s6_mub[i][j] = s6_mub[j][i] = True
                s6_mub_pairs.append((i, j))

    assert s6_ortho_pairs == [], "S6 must have no orthogonal pair (Brierley/Weigert)"
    assert len(s6_mub_pairs) == 1125

    # --- F6 numeric data ---------------------------------------------------
    f6_vectors = build_f6_vectors()
    f6_mat = f6_matrix()
    assert len(f6_vectors) == 48
    assert all(len(v) == 6 for v in f6_vectors)
    assert all(abs(abs(entry) - 1.0) < 1e-12 for v in f6_vectors for entry in v)
    assert all(
        abs(abs(dot_float(v, col)) ** 2 - 6.0) < 1e-9
        for v in f6_vectors
        for col in ([f6_mat[r][c] for r in range(6)] for c in range(6))
    )

    f6_ortho: list[list[bool]] = [[False] * 48 for _ in range(48)]
    f6_mub: list[list[bool]] = [[False] * 48 for _ in range(48)]
    f6_ortho_pairs: list[tuple[int, int]] = []
    f6_mub_pairs: list[tuple[int, int]] = []
    f6_margin = float("inf")
    f6_edge_error = 0.0
    for i in range(48):
        for j in range(i + 1, 48):
            x = abs(dot_float(f6_vectors[i], f6_vectors[j])) ** 2
            if abs(x) < 1e-9:
                f6_edge_error = max(f6_edge_error, abs(x))
                f6_ortho[i][j] = f6_ortho[j][i] = True
                f6_ortho_pairs.append((i, j))
            elif abs(x - 6.0) < 1e-9:
                f6_edge_error = max(f6_edge_error, abs(x - 6.0))
                f6_mub[i][j] = f6_mub[j][i] = True
                f6_mub_pairs.append((i, j))
            else:
                # Every non-edge must stay well away from both thresholds.
                f6_margin = min(f6_margin, abs(x), abs(x - 6.0))

    # The numerical certification margin for the F6 graph: non-edges are
    # at least 1e-6 away from both 0 and 6, and edges are within 1e-9.
    assert f6_margin > 1e-6, f"F6 non-edge margin too small: {f6_margin}"
    assert f6_edge_error < 1e-9, f"F6 edge error too large: {f6_edge_error}"

    f6_bases = cliques6(f6_ortho)
    assert len(f6_bases) == 16, f"expected 16 F6 bases, got {len(f6_bases)}"
    assert all(
        not all(f6_mub[i][j] for i in b1 for j in b2)
        for pos, b1 in enumerate(f6_bases)
        for b2 in f6_bases[pos + 1:]
    )

    (DATA / "f6_bases.csv").write_text(
        "\n".join(",".join(str(v) for v in base) for base in f6_bases) + "\n",
        encoding="utf-8",
    )
    (DATA / "s6_bases.csv").write_text("", encoding="utf-8")

    # --- Facts -------------------------------------------------------------
    def both(pairs: Sequence[tuple[int, int]]) -> list[tuple[int, int]]:
        out: list[tuple[int, int]] = []
        for i, j in pairs:
            out.append((i, j))
            out.append((j, i))
        return out

    counts = {
        "f6_ortho": write_facts(DATA / "f6_ortho.facts", both(f6_ortho_pairs)),
        "f6_mub": write_facts(DATA / "f6_mub.facts", both(f6_mub_pairs)),
        "s6_ortho": write_facts(DATA / "s6_ortho.facts", both(s6_ortho_pairs)),
        "s6_mub": write_facts(DATA / "s6_mub.facts", both(s6_mub_pairs)),
    }

    # --- Lean certificates -------------------------------------------------
    f6_ortho_mask = mask_words(48, lambda i, j: f6_ortho[i][j])
    f6_mub_mask = mask_words(48, lambda i, j: f6_mub[i][j])
    s6_mub_mask = mask_words(90, lambda i, j: s6_mub[i][j])

    certificates = f"""-- AUTO-GENERATED by python/generate_facts.py.  Do not edit by hand.
--
-- F6: 48 Bjorck biunimodular vectors, 16 orthonormal sextuples, and the
--     orthogonality / unbiasedness masks used by MasterProof.lean.
-- S6: 90 seeds-derived vectors.  The exact arithmetic produces no
--     orthogonal pair, so the S6 orthogonality mask is identically zero.
--     The S6 MUB mask is included for completeness of the certificate.

namespace MUB6

/-- The 48 length-6 codes over the 60-element Bjorck alphabet. -/
def f6Codes : List (List Nat) :=
  {lean_list_of_lists(F6_CODES)}

/-- All 16 orthonormal sextuples (bases) in the F6 orthogonality graph. -/
def f6Bases : List (List Nat) :=
  {lean_list_of_lists(f6_bases)}

/-- Bit (i*48+j) is set iff vectors i and j are orthogonal. -/
def f6OrthoMask : Array UInt64 := {lean_array(f6_ortho_mask)}

/-- Bit (i*48+j) is set iff vectors i and j are mutually unbiased. -/
def f6MubMask : Array UInt64 := {lean_array(f6_mub_mask)}

/-- Bit (i*90+j) is set iff S6 vectors i and j are orthogonal.
    The exact computation gives the zero mask. -/
def s6OrthoMask : Array UInt64 := Array.replicate 127 0

/-- Bit (i*90+j) is set iff S6 vectors i and j are mutually unbiased. -/
def s6MubMask : Array UInt64 := {lean_array(s6_mub_mask)}

end MUB6
"""
    (LEAN / "Certificates.lean").write_text(certificates, encoding="utf-8")

    summary = {
        "f6": {
            "vectors": len(f6_vectors),
            "orthogonal_pairs": len(f6_ortho_pairs),
            "mub_pairs": len(f6_mub_pairs),
            "bases": len(f6_bases),
            "no_mub_pair_between_bases": True,
            "non_edge_margin": f6_margin,
            "edge_error": f6_edge_error,
        },
        "s6": {
            "vectors": len(s6_vectors),
            "orthogonal_pairs": len(s6_ortho_pairs),
            "mub_pairs": len(s6_mub_pairs),
            "bases": 0,
            "reason": "no orthogonal pair exists among the 90 candidates",
        },
        "facts": counts,
        "sha256": {
            name: sha256_of(DATA / name)
            for name in ("f6_ortho.facts", "f6_mub.facts", "s6_ortho.facts", "s6_mub.facts")
        },
        "lean_certificate_sha256": sha256_of(LEAN / "Certificates.lean"),
    }
    (DATA / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
