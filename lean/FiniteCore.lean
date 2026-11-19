import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum

namespace LargeN

-- Cross-multiplication in the continuum transfer comparison, for d*j = dJ.
theorem continuumCross (b dJ : ℚ) : b + b*dJ = b*(1+dJ) := by ring

-- Centered random-hopping pair integral after removal of the anomalous hopping term.
theorem hoppingCross (h v u : ℚ) (dyson : u*h = 1-v*h^2) :
    h^2*u = h*(1-v*h^2) := by
  calc
    h^2*u = h*(u*h) := by ring
    _ = h*(1-v*h^2) := by rw [dyson]

-- One periodic two-site pair: the general correlation identity sums this square.
theorem overlapSquare (x y z w : ℚ) :
    (x-y)^2 + (z-w)^2 = x^2+z^2+y^2+w^2-2*(x*y+z*w) := by ring

-- The two exponential integrals of the static dot and the r=2 proportionality.
theorem dotStatic : 1 / (2 * (1 + 1 / 3)) = (3 / 8 : ℚ) := by norm_num

theorem dotSlack : (5 : ℚ)*(2+5/144+5/24) < 12 := by norm_num

end LargeN
