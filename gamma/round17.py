"""Supplemental gamma=2 Rayleigh certificate for rounding to 17 places.

Run from this directory with the same dependencies as replay.py, for example:
    uv run --no-project --with numpy==2.5.5 --with scipy --with python-flint==0.9.0 python -u round17.py

The finite-section Rayleigh principle supplies a lower bound for the infinite
operator; two_tight.certify supplies the upper bound.  The eigensolver only
proposes a finite, positive binary64 vector, whose entries are exact dyadics
when converted to Arb.  All arithmetic used to certify the bounds, including
pi and the square root, is performed with Arb.  The operator characterization
and the upper certificate have the analytic assumptions described in RIGOR.md.
"""

import two_tight as two
from flint import arb, arb_poly, ctx
import numpy as np
from scipy.linalg import matmul_toeplitz
from scipy.signal import fftconvolve
from scipy.sparse.linalg import LinearOperator, eigsh


TRIAL_SIZE = 16384
ROUND_LOWER = "0.182726247772329515"
ROUND_UPPER = "0.182726247772329525"


def exact_actions(vectors):
    """Enclose finite gamma=2 operator actions on Arb vectors."""
    n = len(vectors[0])
    if n == 0 or any(len(v) != n for v in vectors):
        raise ValueError("operator vectors must have the same nonzero length")
    # For x_i = v_i/sqrt(2i+1), use the independent reviewer's
    # convolution indexing: with d_0=0 and d_k=1/k^2, the coefficients
    # at n+i and n-i-1 of d * (reverse(x), x) sum the difference
    # and Hankel kernels.  The latter includes its nonzero diagonal.
    roots = [1 / arb(2 * i + 1).sqrt() for i in range(n)]
    kernel = [arb(0)] + [1 / arb(k * k) for k in range(1, 2 * n)]
    mass = []
    harmonic = arb(0)
    for i in range(n):
        harmonic += kernel[i]
        mass.append(2 * harmonic / (2 * i + 1))
    actions = []
    for v in vectors:
        x = [value * root for value, root in zip(v, roots)]
        conv = arb_poly(x[::-1] + x) * arb_poly(kernel)
        actions.append([roots[i] * (conv[n + i] + conv[n - i - 1])
                        - mass[i] * v[i] for i in range(n)])
    return actions


def dot(v, w):
    return sum((a * b for a, b in zip(v, w)), arb(0))


def rayleigh_lower(n=TRIAL_SIZE):
    """Enclose the temperature associated with an n-site Rayleigh quotient."""
    proposed = two.floating_trial(n)
    v = [arb(float(value)) for value in proposed]
    if len(v) != n or not all(value > 0 for value in v):
        raise ArithmeticError("trial vector must have n strictly positive entries")
    action, = exact_actions([v])
    numerator = dot(v, action)
    norm2 = dot(v, v)
    if not norm2 > 0:
        raise ArithmeticError("trial norm is not strictly positive")
    quotient = numerator / norm2
    if not quotient > 0:
        raise ArithmeticError("Rayleigh quotient is not strictly positive")
    return quotient.sqrt() / (2 * arb.pi())


def reviewer_bad_u():
    """Recreate and certify the reviewer's invalid second-eigenvalue bound."""
    ctx.prec = 160
    ctx.threads = 1
    n = two.TRIAL_SIZE
    w = np.arange(1, 2 * n, 2, dtype=float)
    roots = 1 / np.sqrt(w)
    pp = 1 / np.arange(1, 2 * n, dtype=float) ** 2
    col = np.r_[0., pp[:n - 1]]
    mass = 2 * np.r_[0., np.cumsum(pp[:n - 1])] / w

    def apply(v):
        x = roots * v
        return roots * (matmul_toeplitz((col, col), x, check_finite=False)
                        + fftconvolve(pp, x[::-1])[n - 1:2 * n - 1]) - mass * v

    op = LinearOperator((n, n), matvec=apply, dtype=np.float64)
    _, vecs = eigsh(op, k=2, which="LA", v0=np.ones(n), tol=1e-13)
    v, w = [[arb(float(value)) for value in vecs[:, j]] for j in (1, 0)]
    av, aw = exact_actions([v, w])
    norm = dot(w, w)
    if not norm > 0:
        raise ArithmeticError("second trial vector has nonpositive norm")
    q = dot(w, aw) / norm
    bad_u = q.lower() - arb("1e-12")

    # Positive definiteness of G - bad_u I on this two-dimensional span
    # proves lambda_2(G) > bad_u by the min-max principle.  In particular,
    # the determinant check also proves the proposed vectors are independent.
    x = dot(v, av) - bad_u * dot(v, v)
    y = dot(v, aw) - bad_u * dot(v, w)
    z = dot(w, aw) - bad_u * dot(w, w)
    if not (x > 0 and x * z - y ** 2 > 0):
        raise ArithmeticError("reviewer bad-U two-dimensional witness failed")
    return bad_u


def certify():
    ctx.prec = 160
    ctx.threads = 1
    lower = rayleigh_lower()
    _, upper, _ = two.certify()
    # Arb comparisons prove strict separation of entire balls, not merely
    # their displayed midpoints.  Both halfway boundaries are excluded.
    if not (lower > arb(ROUND_LOWER) and lower < upper
            and upper < arb(ROUND_UPPER)):
        raise ArithmeticError("17-place rounding was not certified")
    print("PASS gamma=2 supplemental 16384-site Rayleigh lower:", lower.str(30), flush=True)
    print("PASS gamma=2 existing infinite-operator upper:", upper.str(30), flush=True)
    print("PASS strictly between 17-place rounding boundaries:",
          ROUND_LOWER, ROUND_UPPER, flush=True)
    print("CERTIFIED gamma=2 T_lin/g rounds to 0.18272624777232952", flush=True)
    return lower, upper


if __name__ == "__main__":
    certify()
