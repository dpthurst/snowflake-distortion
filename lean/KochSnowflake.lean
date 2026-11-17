/-
Copyright: none declared.
Source: `notes/koch-extremals.tex`, Section 3.1.
Scope: this file formalises the exact value of the Koch ratio at the minimising
code pair (5/16, 11/16), namely `2 * 6 ^ (-alpha)`.
It does **not** formalise the global lower bound (that every pair has ratio at
least this value); that requires the interval branch-and-bound certificate.
-/
import Mathlib

noncomputable section

open Complex

namespace KochSnowflake

/-- The sixth root of unity `exp (i*pi/3)`, written in real coordinates. -/
def omega : ℂ := ⟨1 / 2, Real.sqrt 3 / 2⟩

/-- The rotation part `d_k` of the `k`-th Koch similarity, for `k = 0,1,2,3`. -/
def rotation : ℕ → ℂ
  | 0 => 1
  | 1 => omega
  | 2 => 1 - omega
  | _ => 1

/-- The translation part `c_k` of the `k`-th Koch similarity, for `k = 0,1,2,3`. -/
def translation : ℕ → ℂ
  | 0 => 0
  | 1 => 1
  | 2 => 1 + omega
  | _ => 2

/-- The `k`-th Koch map `z ↦ (d_k z + c_k) / 3`. -/
def kochMap (k : ℕ) (z : ℂ) : ℂ := (rotation k * z + translation k) / 3

/-- The coding map at code points `M / 4^k`, for `M < 4^(k+1)`.

The recursion reads the most significant base-4 digit first: at level `k+1` the
leading digit is `M / 4^k` and the remaining code point is `M % 4^k`. -/
def kochPoint : ℕ → ℕ → ℂ
  | _, 0 => 0
  | M, k + 1 => kochMap (M / 4 ^ k) (kochPoint (M % 4 ^ k) k)

/-- The renormalisation exponent `alpha = log 3 / log 4`. -/
def alpha : ℝ := Real.log 3 / Real.log 4

/-- The normalized Koch ratio at the grid pair `(M / 4^k, M' / 4^k)`. -/
def kochRatio (M M' k : ℕ) : ℝ :=
  ‖kochPoint M k - kochPoint M' k‖ /
    |(M : ℝ) / (4 : ℝ) ^ k - (M' : ℝ) / (4 : ℝ) ^ k| ^ alpha

/-- The primitive sixth root has conjugate `1 - omega`. -/
lemma star_omega : star omega = 1 - omega := by
  apply Complex.ext
  · simp [omega]; norm_num
  · simp [omega]

/-- The second Koch rotation is the conjugate of `omega`, as in the source IFS. -/
lemma rotation_two : rotation 2 = star omega := by
  rw [star_omega]
  rfl


/-- The four `kochMap`s are exactly the four similarities of the source IFS. -/
theorem kochMap_matches_source (z : ℂ) :
    kochMap 0 z = z / 3 ∧
    kochMap 1 z = omega * z / 3 + 1 / 3 ∧
    kochMap 2 z = (star omega) * z / 3 + (1 + omega) / 3 ∧
    kochMap 3 z = z / 3 + 2 / 3 := by
  refine ⟨?_, ?_, ?_, ?_⟩
  · simp [kochMap, rotation, translation]
  · simp [kochMap, rotation, translation]; ring
  · rw [kochMap, rotation_two]
    simp [translation]
    ring
  · simp [kochMap, rotation, translation]; ring

/-- The code point `5/16` is `(3 + omega)/9`. -/
lemma kochPoint_five : kochPoint 5 2 = (3 + omega) / 9 := by
  simp only [kochPoint, kochMap, rotation, translation, omega]
  norm_num [Complex.ext_iff]
  ring_nf

/-- The code point `11/16` is `(5 + omega)/9`. -/
lemma kochPoint_eleven : kochPoint 11 2 = (5 + omega) / 9 := by
  simp only [kochPoint, kochMap, rotation, translation, omega]
  norm_num [Complex.ext_iff]
  ring_nf

/-- The chord of the pair `(5/16,11/16)` is exactly `2/9`. -/
lemma kochPoint_five_sub_eleven :
    kochPoint 5 2 - kochPoint 11 2 = -(2 / 9 : ℝ) := by
  rw [kochPoint_five, kochPoint_eleven]
  apply Complex.ext
  · simp [omega]; ring_nf
  · simp [omega]

/-- `4 ^ alpha = 3` for `alpha = log 3 / log 4`. -/
lemma rpow_four : (4 : ℝ) ^ alpha = 3 := by
  rw [alpha, Real.log_div_log]
  exact Real.rpow_logb (by norm_num) (by norm_num) (by norm_num)

/-- `2 ^ alpha = sqrt 3` for `alpha = log 3 / log 4`. -/
lemma rpow_two : (2 : ℝ) ^ alpha = Real.sqrt 3 := by
  have hsq : ((2 : ℝ) ^ alpha) ^ 2 = 3 := by
    calc
      ((2 : ℝ) ^ alpha) ^ 2 = (2 : ℝ) ^ alpha * (2 : ℝ) ^ alpha := by ring
      _ = ((2 : ℝ) * (2 : ℝ)) ^ alpha := by
        rw [Real.mul_rpow (by norm_num) (by norm_num)]
      _ = (4 : ℝ) ^ alpha := by norm_num
      _ = 3 := rpow_four
  have hpos : 0 < (2 : ℝ) ^ alpha := Real.rpow_pos_of_pos (by norm_num) alpha
  have hs : 0 ≤ Real.sqrt 3 := Real.sqrt_nonneg 3
  nlinarith [Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)]

/-- `8 ^ alpha = 3 * sqrt 3`. -/
lemma rpow_eight : (8 : ℝ) ^ alpha = 3 * Real.sqrt 3 := by
  calc
    (8 : ℝ) ^ alpha = ((2 : ℝ) ^ 3) ^ alpha := by norm_num
    _ = (2 : ℝ) ^ (3 * alpha) :=
      (Real.rpow_mul (by norm_num : (0 : ℝ) ≤ 2) 3 alpha).symm
    _ = (2 : ℝ) ^ (alpha * 3) := by rw [mul_comm]
    _ = ((2 : ℝ) ^ alpha) ^ 3 :=
      Real.rpow_mul (by norm_num : (0 : ℝ) ≤ 2) alpha 3
    _ = (Real.sqrt 3) ^ 3 := by simp [rpow_two]
    _ = 3 * Real.sqrt 3 := by
      rw [show (Real.sqrt 3) ^ 3 = (Real.sqrt 3) ^ 2 * Real.sqrt 3 by ring]
      rw [Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)]

/-- The scalar identity behind the endpoint computation:
`(2/9) / (3/8)^alpha = 2 * 6^(-alpha)`. -/
lemma rho_value_algebra :
    (2 / 9 : ℝ) / (3 / 8 : ℝ) ^ alpha = 2 * (6 : ℝ) ^ (-alpha) := by
  have h8 := rpow_eight
  have h2 := rpow_two
  have h3pos : 0 < (3 : ℝ) ^ alpha := Real.rpow_pos_of_pos (by norm_num) alpha
  have hs3pos : 0 < Real.sqrt 3 := Real.sqrt_pos.2 (by norm_num)
  rw [Real.div_rpow (by norm_num : (0 : ℝ) ≤ 3) (by norm_num : (0 : ℝ) ≤ 8), h8]
  rw [Real.rpow_neg (by norm_num : (0 : ℝ) ≤ 6)]
  rw [show (6 : ℝ) = 2 * 3 by norm_num,
    Real.mul_rpow (by norm_num : (0 : ℝ) ≤ 2) (by norm_num : (0 : ℝ) ≤ 3), h2]
  field_simp [h3pos.ne', hs3pos.ne']
  rw [Real.sq_sqrt (show (0 : ℝ) ≤ 3 by norm_num)]
  norm_num

/-- **The value at the minimising configuration.**
At `(x,y) = (5/16, 11/16)`, the Koch ratio is exactly `2 * 6^(-alpha)`. -/
theorem rho_value_at_minimiser :
    kochRatio 5 11 2 = 2 * (6 : ℝ) ^ (-alpha) := by
  unfold kochRatio
  rw [kochPoint_five_sub_eleven]
  have hden :
      |(((5 : ℕ) : ℝ)) / (4 : ℝ) ^ 2 - (((11 : ℕ) : ℝ)) / (4 : ℝ) ^ 2| ^ alpha =
        (3 / 8 : ℝ) ^ alpha := by
    norm_num
  rw [hden]
  norm_num [Complex.normSq, Complex.sq_norm]
  exact rho_value_algebra

end KochSnowflake
