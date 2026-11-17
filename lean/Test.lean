import Mathlib
noncomputable section
open Complex
namespace Koch

def omega : ℂ := ⟨1 / 2, Real.sqrt 3 / 2⟩
def rot : ℕ → ℂ | 0 => 1 | 1 => omega | 2 => 1 - omega | _ => 1
def trans : ℕ → ℂ | 0 => 0 | 1 => 1 | 2 => 1 + omega | _ => 2
def Fq : ℕ → ℕ → ℂ
  | _, 0 => 0
  | M, k + 1 => (rot (M / 4 ^ k) * Fq (M % 4 ^ k) k + trans (M / 4 ^ k)) / 3
noncomputable def alpha : ℝ := Real.log 3 / Real.log 4

lemma Fq_5_2 : Fq 5 2 = (3 + omega) / 9 := by
  simp only [Fq, rot, trans, omega]; norm_num [Complex.ext_iff]; ring_nf
lemma Fq_11_2 : Fq 11 2 = (5 + omega) / 9 := by
  simp only [Fq, rot, trans, omega]; norm_num [Complex.ext_iff]; ring_nf
lemma Fq_diff : Fq 5 2 - Fq 11 2 = -(2/9 : ℝ) := by
  rw [Fq_5_2, Fq_11_2]; apply Complex.ext <;> simp [omega]; ring_nf

lemma rpow_four : (4:ℝ)^alpha = 3 := by
  rw [alpha, Real.log_div_log]
  exact Real.rpow_logb (by norm_num) (by norm_num) (by norm_num)

lemma rpow_two : (2:ℝ)^alpha = Real.sqrt 3 := by
  have hsq : ((2:ℝ)^alpha)^2 = 3 := by
    calc ((2:ℝ)^alpha)^2 = (2:ℝ)^alpha * (2:ℝ)^alpha := by ring
      _ = ((2:ℝ)*(2:ℝ))^alpha := by rw [Real.mul_rpow (by norm_num) (by norm_num)]
      _ = (4:ℝ)^alpha := by norm_num
      _ = 3 := rpow_four
  have hpos : 0 < (2:ℝ)^alpha := Real.rpow_pos_of_pos (by norm_num) alpha
  have hs : 0 ≤ Real.sqrt 3 := Real.sqrt_nonneg 3
  nlinarith [Real.sq_sqrt (show (0:ℝ) ≤ 3 by norm_num)]

lemma rpow_eight : (8:ℝ)^alpha = 3 * Real.sqrt 3 := by
  calc (8:ℝ)^alpha = ((2:ℝ)^3)^alpha := by norm_num
    _ = (2:ℝ)^(3*alpha) := (Real.rpow_mul (by norm_num : (0:ℝ) ≤ 2) 3 alpha).symm
    _ = (2:ℝ)^(alpha*3) := by rw [mul_comm]
    _ = ((2:ℝ)^alpha)^3 := Real.rpow_mul (by norm_num : (0:ℝ) ≤ 2) alpha 3
    _ = (Real.sqrt 3)^3 := by simpa [rpow_two]
    _ = 3 * Real.sqrt 3 := by
        rw [show (Real.sqrt 3)^3 = (Real.sqrt 3)^2 * Real.sqrt 3 by ring]
        rw [Real.sq_sqrt (show (0:ℝ) ≤ 3 by norm_num)]


lemma rho_value_algebra :
    (2/9:ℝ) / (3/8:ℝ)^alpha = 2 * (6:ℝ)^(-alpha) := by
  have h8 := rpow_eight
  have h2 := rpow_two
  have h3pos : 0 < (3:ℝ)^alpha := Real.rpow_pos_of_pos (by norm_num) alpha
  have hs3pos : 0 < Real.sqrt 3 := Real.sqrt_pos.2 (by norm_num)
  rw [Real.div_rpow (by norm_num : (0:ℝ) ≤ 3) (by norm_num : (0:ℝ) ≤ 8), h8]
  rw [Real.rpow_neg (by norm_num : (0:ℝ) ≤ 6)]
  rw [show (6:ℝ) = 2*3 by norm_num, Real.mul_rpow (by norm_num : (0:ℝ) ≤ 2) (by norm_num : (0:ℝ) ≤ 3), h2]
  field_simp [h3pos.ne', hs3pos.ne']
  rw [Real.sq_sqrt (show (0:ℝ) ≤ 3 by norm_num)]
  norm_num


/-- The grid-coded ratio at the code point `M / 4^k` versus `M' / 4^k`. -/
noncomputable def rhoGrid (M M' k : ℕ) : ℝ :=
  ‖Fq M k - Fq M' k‖ /
    |(M:ℝ) / (4:ℝ)^k - (M':ℝ) / (4:ℝ)^k| ^ alpha

/-- The exact value at the first minimising configuration. -/
theorem rho_min_value :
    rhoGrid 5 11 2 = 2 * (6:ℝ)^(-alpha) := by
  unfold rhoGrid
  rw [Fq_diff]
  have hden2 : |(((5:ℕ):ℝ)) / (4:ℝ)^2 - (((11:ℕ):ℝ)) / (4:ℝ)^2| ^ alpha = (3/8:ℝ)^alpha := by norm_num
  rw [hden2]
  norm_num [Complex.normSq, Complex.sq_norm]
  exact rho_value_algebra

end Koch
