# MUB-6 — exact graph obstructions in dimension six

[![Verify](https://github.com/lambda-alpha-labs/MUB-6/actions/workflows/verify.yml/badge.svg)](https://github.com/lambda-alpha-labs/MUB-6/actions/workflows/verify.yml)

This repository is the Graphenium showcase for two finite, exact
obstructions in the dimension-six MUB problem:

- the Fourier pair \(\{I_6,F_6\}\) has **16** orthonormal bases that are
  unbiased to both members, and **no two of those 16 bases are mutually
  unbiased**;
- Tao's spectral matrix \(S_6\) has **90** vectors unbiased to both
  \(I_6\) and \(S_6\), and **no two of those 90 vectors are orthogonal**
  — so the pair cannot even be extended to a third basis.

The pipeline is split into three layers:

1. `python/generate_facts.py` — exact/numerically certified generation of
   the two vector sets, the orthogonality and unbiasedness graphs, the
   clique enumeration, and the Lean certificates.
2. `graphenium/` — Graphenium's embedded Datalog layer. `build_graph.py`
   bridges the exact facts into a small Graphenium graph, and
   `mub_cliques.dl` finds all 6-cliques and all complete bipartite
   \(K_{6,6}\) subgraphs. No Soufflé installation is required.
3. `lean4/MUB6/` — Lean 4 checks of the algebra and the finite graph
   certificates.

The important point is that the repository does **not** claim to prove
the full dimension-six MUB classification from scratch. The completeness
of the 48/90 vector sets is the published input; the code verifies the
exact coordinates, the graph construction, the clique enumeration, and
the combinatorial obstruction. See [Trust boundary](#trust-boundary).

This showcase is also documented as a
[Graphenium worked case study](https://github.com/lambda-alpha-labs/Graphenium/blob/main/worked/mub6-dimension-six/README.md).

---

## Results

| Data set | Vectors | Orthogonal pairs | MUB pairs | 6-cliques (bases) | MUB pairs between bases |
|---|---:|---:|---:|---:|---:|
| \(F_6\) | 48 | 300 | 144 | 16 | **0** |
| Tao \(S_6\) | 90 | **0** | 1125 | **0** | **0** |

For \(F_6\), every vector is exactly unbiased to \(I_6\) and \(F_6\).
For \(S_6\), every one of the 90 vectors is exactly unbiased to \(I_6\)
and the chosen \(S_6\) representative, and the exact orthogonality graph
is empty.

The \(F_6\) graph is computed in double precision with a certified
separation margin:

- every non-edge is at least `0.4307...` away from both thresholds
  `0` and `6`;
- every edge is within `9.8e-15` of its threshold.

The \(S_6\) graph is computed in exact integer arithmetic (no floating
point at all).

---

## Verification

```sh
make verify
```

runs the Python generator, the exhaustive C search, Graphenium's embedded
Datalog layer, and the Mathlib-free Lean 4 core build in sequence. The
GitHub Actions workflow in `.github/workflows/verify.yml` runs the same
command on a fresh Ubuntu runner with Python, a C compiler, Lean via
`elan`, and the published Graphenium CLI.

---

## Trust boundary

### Inputs assumed from the literature

The following are **assumed**, not re-proved here:

- The 48 vectors unbiased to \((I_6,F_6)\) are exactly the Björck
  biunimodular vectors used below. This is the classical Grassl/Björck
  computation; see Bengtsson–Bruzda–Ericsson–Larsson–Tadej–Życzkowski,
  *J. Math. Phys.* **48**, 052106 (2007), Sec. VI.
- The 90 vectors unbiased to \((I_6,S_6)\) are exactly the supplied
  seed set. This is the Brierley–Weigert computation; see
  *Phys. Rev. A* **79**, 052316 (2009) and S. Brierley's thesis,
  University of York (2009), Table 3.1.
- A third basis extending \((I_6,S_6)\) would have to consist of six
  mutually orthogonal vectors from the 90-vector set. This is the
  standard MUB completeness reduction.

### Verified here

- The \(S_6\) algebra is exact: `Algebra.lean` checks, by pure kernel
  `decide`, that all 90 vectors have unit entries and that all
  \(90\times 6 = 540\) inner products with \(S_6\) have squared modulus
  exactly \(6\) (after the common scaling).
- `Algebra.lean` also checks, directly from the exact vectors, that no
  two of the 90 \(S_6\) vectors are orthogonal. This is an 8100-fold
  kernel computation, not a mask lookup.
- `MasterProof.lean` checks the finite graph certificates:
  - the \(S_6\) orthogonality mask cross-check is identically zero;
  - the 16 \(F_6\) bases are six-element cliques with distinct vertices;
  - no two of the 16 \(F_6\) bases are completely joined in the MUB
    relation.
- The Python driver independently recomputes the same graphs and
  cliques, and asserts all counts and margins.
- Graphenium's embedded Datalog layer runs entirely inside `gm query
  --datalog`. `graphenium/build_graph.py` maps the exact `ortho` / `mub`
  facts onto Graphenium's graph EDB, and `graphenium/mub_cliques.dl`
  recovers the 16 F6 bases with no MUB pair, and no S6 basis or
  orthogonal pair. No external Datalog binary is needed.

### Not verified here

- The classification of **all** vectors unbiased to a given pair is not
  formalized. That is the hard analytic/algebraic input.
- The \(F_6\) coordinates are not re-formalized in Lean; Lean checks the
  graph certificate produced from them. The \(S_6\) coordinates are
  formalized exactly.
- The Python driver remains the independent executable reference for the
  fact generation and clique counts.

No `native_decide`, no `sorry`, and no Mathlib are used in
`Algebra.lean`, `Certificates.lean`, or `MasterProof.lean`.
`Eisenstein.lean` uses Mathlib for `nlinarith`, `interval_cases`, and
`norm_num`; it is intentionally not part of the Mathlib-free core
build.

The compiled core theorems were audited with `#print axioms`: they
depend only on `propext` (and `Quot.sound` for the derived `Fin`
theorem), not on `sorryAx` and not on `Lean.ofReduceBool` from
`native_decide`.  The Mathlib `Eisenstein.lean` theorems likewise
depend only on `propext`, `Classical.choice`, and `Quot.sound`.

---

## Important corrections to the supplied data

Three corrections were necessary while implementing the exact pipeline.

### 1. The \(U_9\) table had two sign errors

The supplied table contained

```lean
(-1,  3, -1,  1), -- claimed 8 * omega^2 * beta
(-1, -3, -1, -1)  -- claimed 8 * omega^2 * beta-bar
```

Both have the wrong \(i\sqrt{15}\) sign. The correct values are

```lean
(-1,  3, -1, -1), -- 8 * omega^2 * beta
(-1, -3, -1,  1)  -- 8 * omega^2 * beta-bar
```

The corrected table is the one used in `Algebra.lean` and in
`python/generate_facts.py`. With the uncorrected table the seed vectors
do not have unit-modulus entries.

### 2. The \(S_6\) representative must be row-permuted

The standard thesis/paper representative is (with
\(\omega=e^{2\pi i/3}\))

\[
S=\begin{pmatrix}
1&1&1&1&1&1\\
1&1&\omega&\omega&\omega^2&\omega^2\\
1&\omega&1&\omega^2&\omega^2&\omega\\
1&\omega&\omega^2&1&\omega&\omega^2\\
1&\omega^2&\omega^2&\omega&1&\omega\\
1&\omega^2&\omega&\omega^2&\omega&1
\end{pmatrix}.
\]

The supplied 90 seeds are not unbiased to this representative. They are
unbiased to the row permutation

\[
P=(0,1,2,5,4,3)
\]

of \(S\). The matrix \(P S\) is Hadamard-equivalent to \(S\) and is the
representative used in `Algebra.lean`. There are 60 equivalent working
row orders; any one gives the same obstruction.

### 3. The \(F_6\) vectors are not in the \(r_\pm\) alphabet

The \(r_\pm\) alphabet in the original brief produces only 12 of the
48 vectors (the equilateral family), not the full set. The complete set
is the Bjorck biunimodular set. Its alphabet is

\[
\mathcal A=\{d^m\zeta^k:-2\le m\le 2,\ 0\le k<12\},
\qquad
d=\frac{1-\sqrt3}{2}+i\sqrt{\frac{\sqrt3}{2}},
\quad \zeta=e^{2\pi i/12}.
\]

The 48 length-6 codes over this alphabet are listed in
`f6Codes` in `Certificates.lean`. Each code gives one vector; the
Fourier transform permutes the set. This reproduces the known counts
\(N_v=48\), \(N_t=16\) from Brierley's Table 3.1. The exhaustive search
that produced the codes is included as `tools/f6_bjorck_search.c`; it
reports exactly 48 solutions.

---

## Repository layout

```
Makefile                   # one-command verification: make verify
.github/workflows/
  verify.yml               # fresh Ubuntu CI: Python, C, Graphenium, Lean

python/
  generate_facts.py        # exact/numerically certified generator + driver
tools/
  f6_bjorck_search.c       # exhaustive search for the 48 Bjorck codes
graphenium/
  build_graph.py           # facts -> Graphenium graph JSON bridge
  mub_cliques.dl           # Graphenium embedded 6-clique / K_{6,6} program
  verify_embedded.sh       # runs both datasets through `gm`
  README.md                # relation encoding and usage notes
data/
  f6_ortho.facts           # 600 directed orthogonal pairs
  f6_mub.facts             # 288 directed unbiased pairs
  s6_ortho.facts           # empty (no orthogonal pair)
  s6_mub.facts             # 2250 directed unbiased pairs
  f6_bases.csv             # the 16 bases
  s6_bases.csv             # empty (no orthogonal pair, hence no basis)
  summary.json             # counts, margins, hashes
lean4/
  lakefile.toml
  MUB6.lean                # imports the Mathlib-free core
  MUB6/
    Algebra.lean           # exact Q(sqrt(-3),sqrt(5)) arithmetic + 90 vectors
    Certificates.lean      # auto-generated graph certificates
    MasterProof.lean       # decide checks of the graph obstructions
    Eisenstein.lean        # Eisenstein norm obstruction (Mathlib)
```

---

## Running it

### Quick verification (all layers)

```sh
make verify
```

This runs the Python generator, the exhaustive C search, Graphenium's
embedded Datalog layer, and the Mathlib-free Lean core build in sequence.
Prerequisites are Python 3, a C compiler, Lean 4 with `lake` (for example
via `elan`), and the Graphenium CLI `gm`.  Set `GM=/path/to/gm` to override
the default `gm` lookup.  The CI workflow uses the same command on a fresh
Ubuntu runner.

### Python driver

```sh
python3 python/generate_facts.py
```

This rewrites `data/*.facts`, `data/summary.json`, and
`lean4/MUB6/Certificates.lean`, and asserts every structural property
used by Lean. It uses only the Python standard library.

### Lean core build

```sh
cd lean4
lake build
```

This builds `Algebra.lean`, `Certificates.lean`, and `MasterProof.lean`
without Mathlib. The `decide` proofs are pure kernel computations.

### Eisenstein proof

`Eisenstein.lean` imports Mathlib.  It was checked against Mathlib
`v4.34.0-rc2` (the same Lean toolchain as the core build).  The default
`lakefile.toml` is deliberately Mathlib-free, so to build this file add

```toml
[[require]]
name = "mathlib"
scope = "leanprover-community"
rev = "v4.34.0-rc2"
```

to a Mathlib-enabled Lake project and run

```sh
lake env lean MUB6/Eisenstein.lean
```

The theorem proved there is

```lean
theorem eisNorm_ne_six (a b : ℤ) : a ^ 2 - a * b + b ^ 2 ≠ 6
```

together with the finite residue check
`∀ a b : ZMod 4, a^2 - a*b + b^2 ≠ 2`.

### Graphenium embedded Datalog layer

```sh
./graphenium/verify_embedded.sh
```

This builds two small Graphenium graphs from the exact facts, then runs
`graphenium/mub_cliques.dl` through `gm query --datalog`. It asserts 16
F6 bases with zero MUB pairs between them and an empty S6 obstruction.
No Soufflé installation is required.  See `graphenium/README.md` for the
relation encoding and the `calls` / `imports` / `contains` mapping.

---

## What the Lean theorems say

`Algebra.lean`:

```lean
theorem s6_length : s6Vectors.length = 90
theorem s6_units : s6Vectors.all isUnitVector = true
theorem s6_all_unbiased :
    s6Vectors.all (fun v => isUnbiasedTo v s6) = true

-- 8100 exact inner products; pure kernel `decide`, about 50 s.
theorem s6_no_ortho_exact :
    (List.range 90).all (fun i =>
      (List.range 90).all (fun j =>
        i == j || !(isOrthoPair (List.getD s6Vectors i [])
          (List.getD s6Vectors j [])))) = true
```

`MasterProof.lean`:

```lean
theorem s6_no_orthogonal_pair_exact (i j : Nat)
    (hi : i < 90) (hj : j < 90) (hij : i ≠ j) :
    isOrthoPair (List.getD s6Vectors i []) (List.getD s6Vectors j []) = false

theorem s6_no_orthogonal_pair :
    ∀ i j : Fin 90, i ≠ j → bit s6OrthoMask (i.val * 90 + j.val) = false

theorem f6_bases_valid :
    bases16.all
      (fun b => b.length == 6 && b.eraseDups.length == 6 && isClique b) = true

theorem f6_no_mub_pair :
    bases16.all
      (fun b1 => bases16.all
        (fun b2 => b1 == b2 || !(isMubPair b1 b2))) = true
```

The \(S_6\) theorem is the strong form: not only are there no 6-cliques,
there is not even an orthogonal pair. A third basis would require six
mutually orthogonal vectors, hence an orthogonal pair, so the pair
\((I_6,S_6)\) cannot be extended. The 8100-fold exact computation is why
the Mathlib-free core build takes about a minute; this is a deliberate
trade for a kernel-checked result rather than a `native_decide` shortcut.

For \(F_6\), the Datalog/Python layer enumerates all 6-cliques; Lean
checks that the injected list contains 16 valid bases and that no two
are MUB. A quartet containing the standard/Fourier pair would require
two of those bases to be MUB.

---

## References

- I. Bengtsson, W. Bruzda, Å. Ericsson, J.-Å. Larsson, W. Tadej,
  K. Życzkowski, *Mutually unbiased bases and Hadamard matrices of
  order six*, J. Math. Phys. **48**, 052106 (2007).
- S. Brierley, S. Weigert, *Constructing mutually unbiased bases in
  dimension six*, Phys. Rev. A **79**, 052316 (2009).
- S. Brierley, *Mutually Unbiased Bases in Low Dimensions*, PhD thesis,
  University of York (2009).
- T. Tao, *Fuglede's conjecture is false in 5 and higher dimensions*,
  Math. Res. Lett. **11**, 251 (2004).
