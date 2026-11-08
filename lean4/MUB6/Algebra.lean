/-
  Exact algebra for Tao's spectral matrix S6.

  Every entry of S6 and of the 90 bi-unbiased vectors lies in
  Q(sqrt(-3), sqrt(5)).  We work with the integer basis

      1, sqrt(5), i*sqrt(3), i*sqrt(15)

  and store 8*z as a quadruple (a,b,c,d):

      8*z = a + b*sqrt(5) + c*(i*sqrt(3)) + d*(i*sqrt(15)).

  The factor 8 clears the denominators in omega and beta.  The
  multiplication and norm formulas are the ones supplied by the
  project, with the U9 table corrected at the two omega^2 entries.

  This file is Mathlib-free: it only uses core Lean 4 definitions and
  the kernel `decide` tactic.  The algebraic checks below are therefore
  pure kernel computations.
-/

namespace MUB6

abbrev Q3Q5Int := Int × Int × Int × Int

def qAdd (x y : Q3Q5Int) : Q3Q5Int :=
  (x.1 + y.1, x.2.1 + y.2.1, x.2.2.1 + y.2.2.1, x.2.2.2 + y.2.2.2)

def qConj (x : Q3Q5Int) : Q3Q5Int :=
  (x.1, x.2.1, -x.2.2.1, -x.2.2.2)

/-- Multiplication in the integer lattice spanned by
`1, sqrt(5), i*sqrt(3), i*sqrt(15)`. -/
def qMul (x y : Q3Q5Int) : Q3Q5Int :=
  let (a1, b1, c1, d1) := (x.1, x.2.1, x.2.2.1, x.2.2.2)
  let (a2, b2, c2, d2) := (y.1, y.2.1, y.2.2.1, y.2.2.2)
  ( a1*a2 + 5*b1*b2 - 3*c1*c2 - 15*d1*d2,
    a1*b2 + b1*a2 - 3*c1*d2 - 3*d1*c2,
    a1*c2 + 5*b1*d2 + c1*a2 + 5*d1*b2,
    a1*d2 + b1*c2 + c1*b2 + d1*a2 )

/-- Squared modulus, as an element `p + q*sqrt(5)` of `Z[sqrt(5)]`. -/
def qNormSq (x : Q3Q5Int) : Int × Int :=
  let (a, b, c, d) := (x.1, x.2.1, x.2.2.1, x.2.2.2)
  ( (a*a + 5*b*b) + 3*(c*c + 5*d*d), 2*a*b + 6*c*d )

/-- `8 * U9`, where
`U9 = {1, omega, omega^2} * {1, beta, beta-bar}`,
`omega = (-1 + i*sqrt(3))/2`, `beta = (1 + i*sqrt(15))/4`.

The last two entries are the corrected ones: the `i*sqrt(15)` sign in
the supplied table was flipped for `omega^2 * beta` and
`omega^2 * beta-bar`. -/
def u9 : List Q3Q5Int :=
  [(8, 0, 0, 0),      -- 0: 1
   (2, 0, 0, 2),      -- 1: beta
   (2, 0, 0, -2),     -- 2: beta-bar
   (-4, 0, 4, 0),     -- 3: omega
   (-1, -3, 1, -1),   -- 4: omega * beta
   (-1, 3, 1, 1),     -- 5: omega * beta-bar
   (-4, 0, -4, 0),    -- 6: omega^2
   (-1, 3, -1, -1),   -- 7: omega^2 * beta
   (-1, -3, -1, 1)]   -- 8: omega^2 * beta-bar

def u9At (i : Nat) : Q3Q5Int := List.getD u9 i (0, 0, 0, 0)

/-- The 18 supplied seeds.  Each seed has five U9 indices. -/
def seeds : List (List Nat) :=
  [[0, 1, 3, 3, 4], [0, 2, 3, 3, 5], [0, 4, 3, 3, 1], [0, 5, 3, 3, 2],
   [0, 6, 1, 7, 6], [0, 6, 2, 8, 6], [0, 6, 4, 4, 6], [0, 6, 5, 5, 6],
   [0, 6, 7, 1, 6], [0, 6, 8, 2, 6], [0, 7, 3, 3, 7], [0, 7, 4, 4, 7],
   [0, 8, 3, 3, 8], [0, 8, 5, 5, 8], [1, 1, 4, 3, 4], [1, 6, 1, 7, 7],
   [2, 2, 5, 3, 5], [2, 6, 2, 8, 8]]

/-- The five cyclic shifts of a length-5 seed. -/
def shifts5 : List Nat → List (List Nat)
  | [a, b, c, d, e] =>
      [[a, b, c, d, e], [b, c, d, e, a], [c, d, e, a, b],
       [d, e, a, b, c], [e, a, b, c, d]]
  | _ => []

/-- The 90 length-6 U9 codes.  The leading 0 is the scalar 1. -/
def s6Codes : List (List Nat) :=
  seeds.flatMap (fun s => (shifts5 s).map (fun t => 0 :: t))

/-- The 90 exact S6 vectors. -/
def s6Vectors : List (List Q3Q5Int) :=
  s6Codes.map (fun c => c.map (fun i => u9At i))

/-- Exact Hermitian inner product (conjugate-linear in the second slot). -/
def dot (u v : List Q3Q5Int) : Q3Q5Int :=
  (u.zipWith (fun x y => qMul x (qConj y)) v).foldl qAdd (0, 0, 0, 0)

/-- The six columns of a 6x6 matrix represented as a list of rows. -/
def cols6 (M : List (List Q3Q5Int)) : List (List Q3Q5Int) :=
  (List.range 6).map (fun j => M.map (fun row => List.getD row j (0, 0, 0, 0)))

/-- `v` is unbiased to the columns of `M`.

For the scaled representation, `v` and the columns of `M` have entries
`8 * unit`, so the squared inner product is `64^2 * 6 = 24576`. -/
def isUnbiasedTo (v : List Q3Q5Int) (M : List (List Q3Q5Int)) : Bool :=
  (cols6 M).all (fun col => qNormSq (dot v col) == (24576, 0))

/-- A length-6 vector all of whose entries have squared modulus 64. -/
def isUnitVector (v : List Q3Q5Int) : Bool :=
  v.length == 6 && v.all (fun x => qNormSq x == (64, 0))

def one : Q3Q5Int := u9At 0
def w : Q3Q5Int := u9At 3
def w2 : Q3Q5Int := u9At 6

/-- Tao's spectral matrix in the representative matching the seeds.

The standard thesis/paper matrix is row-permuted by
`(0,1,2,5,4,3)`; the resulting matrix is Hadamard-equivalent to the
standard one and is the representative for which the supplied 90 seeds
are exactly the bi-unbiased vectors. -/
def s6 : List (List Q3Q5Int) :=
  [ [one, one, one, one, one, one],
    [one, one, w, w, w2, w2],
    [one, w, one, w2, w2, w],
    [one, w2, w, w2, w, one],
    [one, w2, w2, w, one, w],
    [one, w, w2, one, w, w2] ]

theorem s6_length : s6Vectors.length = 90 := by
  decide

theorem s6_units : s6Vectors.all isUnitVector = true := by
  decide

/-- Every one of the 90 vectors is exactly unbiased to every column of S6.
This is a 90 x 6 = 540-fold exact integer check. -/
theorem s6_all_unbiased :
    s6Vectors.all (fun v => isUnbiasedTo v s6) = true := by
  decide

/-- Exact orthogonality predicate. -/
def isOrthoPair (u v : List Q3Q5Int) : Bool :=
  qNormSq (dot u v) == (0, 0)

set_option maxRecDepth 1000000 in
set_option maxHeartbeats 10000000 in
/-- Strong exact obstruction: no two of the 90 vectors are orthogonal.
This is an 8100-fold kernel computation over the exact integer vectors.
The `decide` proof takes about 50 seconds; the heartbeats allowance is
raised locally for this one theorem. -/
theorem s6_no_ortho_exact :
    (List.range 90).all (fun i =>
      (List.range 90).all (fun j =>
        i == j || !(isOrthoPair (List.getD s6Vectors i [])
          (List.getD s6Vectors j [])))) = true := by
  decide

end MUB6
