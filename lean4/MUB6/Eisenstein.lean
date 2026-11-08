/-
  Eisenstein norm obstruction.

  The norm form of the Eisenstein integers is

      N(a + b*omega) = a^2 - a*b + b^2,   omega = (-1 + i*sqrt(3))/2.

  No integer pair (a,b) represents 6.  Equivalently, the norm form is
  never congruent to 2 modulo 4, while 6 = 2 (mod 4).

  This is the arithmetic lemma used to rule out unbiasedness between two
  BH(6,6) bases: the overlap numerator is a sum of sixth roots, hence an
  Eisenstein integer, and unbiasedness would force its norm to be 6.

  This file uses Mathlib for `nlinarith`, `interval_cases` and
  `norm_num`.  The other Lean files in this project are Mathlib-free.
-/

import Mathlib

namespace MUB6

/-- The Eisenstein norm form. -/
def eisNorm (a b : ℤ) : ℤ := a ^ 2 - a * b + b ^ 2

/-- 6 is not represented by the Eisenstein norm form. -/
theorem eisNorm_ne_six (a b : ℤ) : eisNorm a b ≠ 6 := by
  intro h
  -- Complete the square:  4*N = (2a-b)^2 + 3b^2.
  have h4 : (2 * a - b) ^ 2 + 3 * b ^ 2 = 24 := by
    unfold eisNorm at h
    nlinarith
  -- Hence b^2 <= 8, so b is one of -2,-1,0,1,2.
  have hb2 : b ^ 2 ≤ 8 := by
    nlinarith [sq_nonneg (2 * a - b)]
  have hb_bounds : -2 ≤ b ∧ b ≤ 2 := by
    constructor
    · by_contra hlt
      push Not at hlt
      nlinarith [sq_nonneg (b + 3)]
    · by_contra hgt
      push Not at hgt
      nlinarith [sq_nonneg (b - 3)]
  -- From h4 and 3*b^2 >= 0 we also get (2a-b)^2 <= 24, hence
  -- -4 <= 2a-b <= 4 because 2a-b is an integer.
  have hx_sq : (2 * a - b) ^ 2 ≤ 24 := by
    nlinarith [sq_nonneg b]
  have hx_upper : 2 * a - b ≤ 4 := by
    by_contra hgt
    have hx_ge5 : 5 ≤ 2 * a - b := by omega
    nlinarith [sq_nonneg (2 * a - b - 5)]
  have hx_lower : -4 ≤ 2 * a - b := by
    by_contra hlt
    have hx_le_neg5 : 2 * a - b ≤ -5 := by omega
    nlinarith [sq_nonneg (2 * a - b + 5)]
  have hb : b = -2 ∨ b = -1 ∨ b = 0 ∨ b = 1 ∨ b = 2 := by
    omega
  -- For each b, the equation bounds a to [-3,3]; check the finite cases.
  rcases hb with rfl | rfl | rfl | rfl | rfl <;>
    (have ha_lower : -3 ≤ a := by omega
     have ha_upper : a ≤ 3 := by omega
     interval_cases a <;> norm_num at h4)

/-- Modulo 4, the Eisenstein norm form is never 2.  This is the finite
residue check behind the integer statement above. -/
theorem eisNorm_zmod_four_ne_two :
    ∀ a b : ZMod 4, a ^ 2 - a * b + b ^ 2 ≠ (2 : ZMod 4) := by
  decide

end MUB6
