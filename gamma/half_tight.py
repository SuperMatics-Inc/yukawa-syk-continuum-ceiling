"""Certified gamma=1/2 linear-instability temperature, k_B T_lin/g.

Run from GAMMA (no cache writes to the workspace):
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 uv run --no-project --with numpy==2.5.5 --with python-flint==0.9.0 python -u half_tight.py

Use the compact self-adjoint KAY operator from RIGOR.md, with w_i=2i+1,
 t_i=i+1/2, p_k=1/sqrt(k) for k>0, and p_0=0. Its matrix is
 G_ij=(p_|i-j|+p_(i+j+1))/sqrt(w_i*w_j) - delta_ij*2*sum(p_1..p_i)/w_i.
Split at N into blocks G=[[A,B],[B*,C]]. The positive trial is h on
0<=i<N, z_m=eta/t_m for N<=m<M, and zero on m>=M. M is implicit.

Lemma (all expressions in the proof are real, before ball evaluation).
Let x_i=h_i/sqrt(w_i)>0, a=sqrt(2)*sum(x_i), and b_m=(B*h)_m.
For m>=N,
  t_m*b_m = sum_i x_i*((1-t_i/t_m)^(-1/2)+(1+t_i/t_m)^(-1/2))/sqrt(2).
The sum of the two powers is 2 at zero and has nonnegative derivative
on [0,1) by direct differentiation. Hence a<=t_m*b_m<=a+e, where
 e=t_N*b_N-a. For any positive lambda use eta=a/lambda. Write
 D=sum_{N<=m<M} 1/t_m^2 and D+=1/N-1/M. Each summand is
 (1/m-1/(m+1))*4*m*(m+1)/(4*m*(m+1)+1), so
 D+*4*N*(N+1)/(4*N*(N+1)+1) <= D <= D+.
The same telescoping bound gives sum_{m>=M} 1/t_m^2 <= 1/M.

For the C block, write C=P-diag(d_n), P entrywise nonnegative and
 d_n=2*sum(p_1..p_n)/w_n. The weighted row of P at n is
 sum_{m>=N} (p_|n-m|+p_(n+m+1))/w_m (weights 1/sqrt(w_n)).
Split at m<=n/2, n/2<m<=2n, m>2n. The first piece is zero if n<2N;
otherwise it is <=sqrt(2/n)*(1/N+log(n/(2*N))) by harmonic integration,
and thus <=(1/N+2/exp(1))/sqrt(N), using max_{u>=1} log(u)/sqrt(u)=2/exp(1).
In the middle, w_m>=n, both distance sums are bounded by sum_{k<=n} p_k
and there are <=2n Hankel terms each <=1/sqrt(n). This gives <=6/sqrt(n).
On the right the sum is <=sqrt(2)*sum_{m>2n} m^(-3/2) <= 2/sqrt(n).
Therefore ||P||<=cP=(8+1/N+2/exp(1))/sqrt(N) by Schur. Also
 0<=d_n<=d=2/sqrt(N), since sum_{k<=n} p_k <= 2*sqrt(n)-1.
In particular ||C||<=cP+d, and for z>=0,
 -d*||z||^2 <= <z,Cz> <= cP*||z||^2, using <z,Pz>>=0.

For a global second-eigenvalue bound restrict to u_0=0. For the same
weights, a symmetric nonnegative P satisfies the weighted row form
bound by 2|u_i u_j| <= (w_i/w_j)^(1/2)|u_i|^2+(w_j/w_i)^(1/2)|u_j|^2.
Thus an upper bound on all signed weighted rows (subtracting d_n)
of the deletion is a form bound and, by min-max, bounds lambda_2.
For n>=L the same row splitting with first column m=1 bounds the row
by [sqrt(2)*(1+log(n/2))+8]/sqrt(n), decreasing for n>=L>=2*exp(1).
For 1<=n<L sum explicitly through J>=2*(L-1), retaining -d_n. The
remaining row is <=sqrt(2)*sum_{m>J} m^(-3/2)<=2*sqrt(2/J).
These form arguments apply on finitely supported vectors and extend by
continuity; the positive Schur weights need not be square-summable.

Finally ||B||^2 <= beta2 := sum_{i<N} (2/w_i)*(1/(N-i)+1/(N-i)^2):
 G_im^2 <= 2/(w_i*(m-i)^2) for m>=N, and summing uses the integral test.
Let E=<h,Ah> and H=||h||^2. The exact Rayleigh quotient q of u=(h,z)
lies between
  (E+(2*a*eta-d*eta^2)*D)/(H+eta^2*D),
  (E+(2*(a+e)*eta+cP*eta^2)*D)/(H+eta^2*D).
The normalized variance ||(G-q)u||^2/||u||^2 is at most
 ||(G-lambda)u||^2/||u||^2, since q minimizes the squared residual.
The head residual norm is <=||(A-lambda)h||+sqrt(beta2)*eta*sqrt(D+).
In the tail, B*h-lambda*z has disjoint active and omitted parts with
squared norms <=e^2*D+ and <=(a+e)^2/M. Add the norm of Cz, at most
 (cP+d)*eta*sqrt(D+), to their joint norm. Square the two block bounds,
add, and divide by H+eta^2*D to get R2, an upper bound on the variance.
If l<=q<=r and b>=lambda_2 with l>b, Temple gives
 l<=rho(G)<=r+R2/(l-b). Indeed the spectral measure of the unit trial
satisfies E_spec[(rho-t)*(t-b)]<=0: the integrand is zero at rho and
nonpositive everywhere else. Rearrangement gives (rho-q)*(q-b)
<=E_spec[t^2]-q^2. The map rho -> rho^2/(2*pi) is increasing for rho>0.

Arb computes the entire certificate, including the finite matvec by
polynomial convolution, with no dense matrix. NumPy FFT proposes only
binary64 trial entries; no floating-point residual is used as evidence.
"""

import time

import numpy as np
from flint import arb, arb_poly, ctx


HEAD = 1 << 17
SUPPORT = 1 << 34
PLACES = 16


def propose_head(n):
    """Shifted power iteration; neither convergence nor accuracy is assumed."""
    w = np.arange(1, 2 * n, 2, dtype=np.float64)
    root_inv = 1 / np.sqrt(w)
    p = np.empty(2 * n)
    p[0] = 0
    p[1:] = 1 / np.sqrt(np.arange(1, 2 * n, dtype=np.float64))
    mass = 2 * np.cumsum(p[:n]) / w
    length = 1 << (4 * n - 2).bit_length()
    kernel = np.fft.rfft(p, length)
    h = root_inv / np.sqrt(np.sum(root_inv**2))
    for _ in range(56):
        x = root_inv * h
        folded = np.concatenate((x[::-1], x))
        conv = np.fft.irfft(kernel * np.fft.rfft(folded, length), length)
        new_h = root_inv * (conv[n:2 * n] + conv[n - 1::-1]) + (1 - mass) * h
        h = new_h / np.sqrt(np.sum(new_h**2))
    if not (np.all(np.isfinite(h)) and np.all(h > 0)):
        raise ArithmeticError('the finite head trial is not strictly positive and finite')
    return [arb(float(entry)) for entry in h]


def head_action(h, p):
    """Exact finite action enclosed by one Arb polynomial multiplication."""
    n = len(h)
    root_inv = [1 / arb(2 * i + 1).sqrt() for i in range(n)]
    x = [entry * root for entry, root in zip(h, root_inv)]
    # Reflect x_j to the indices -j-1. Products with the positive-lag
    # kernel give the negative lags by reflection about n-1/2. At row
    # i the two coefficients add all |i-j| and i+j+1 kernels, with p_0=0.
    conv = arb_poly(x[::-1] + x) * arb_poly(p[:2 * n])
    harmonic = arb(0)
    action = []
    for i in range(n):
        harmonic += p[i]
        action.append(root_inv[i] * (conv[n + i] + conv[n - 1 - i])
                      - 2 * harmonic * h[i] / (2 * i + 1))
    return action, x


def second_upper(p):
    """Signed weighted Schur bound on the compression deleting index zero."""
    first_infinite = 192
    last_column = 384
    b = arb(5) / 4
    root2 = arb(2).sqrt()
    high_rows = (root2 * (1 + (arb(first_infinite) / 2).log()) + 8) / arb(first_infinite).sqrt()
    if not high_rows < b:
        raise ArithmeticError('second spectrum: high-row estimate failed')
    remainder = 2 * (arb(2) / last_column).sqrt()
    weights = [1 / arb(2 * j + 1) for j in range(1, last_column + 1)]
    harmonic = arb(0)
    for i in range(1, first_infinite):
        harmonic += p[i]
        row = sum(((p[abs(i - j)] + p[i + j + 1]) * weight
                   for j, weight in enumerate(weights, 1)), arb(0))
        if not row - 2 * harmonic / (2 * i + 1) + remainder < b:
            raise ArithmeticError(f'second spectrum: signed weighted row {i} not certified')
    return b


def compute(n=HEAD, m=SUPPORT):
    """Return enclosing eigenvalue endpoints and diagnostic upper bounds."""
    if n < 288 or m <= n:
        raise ValueError('require 288 <= head size < support size')
    ctx.prec = 160
    ctx.threads = 1
    p = [arb(0)] + [1 / arb(k).sqrt() for k in range(1, 2 * n + 1)]
    b = second_upper(p)
    h = propose_head(n)
    ah, x = head_action(h, p)
    head_norm2 = sum((v**2 for v in h), arb(0))
    head_energy = sum((v * av for v, av in zip(h, ah)), arb(0))
    shift = (head_energy / head_norm2).mid()
    if not shift > 0:
        raise ArithmeticError('nonpositive shift for the tail trial')
    head_residual = sum(((av - shift * v)**2 for v, av in zip(h, ah)), arb(0)).sqrt()

    root2 = arb(2).sqrt()
    a = root2 * sum(x, arb(0))
    first_tail_root = (arb(n) + arb(1) / 2).sqrt()
    e = sum((v * (first_tail_root * (p[n - i] + p[n + i + 1]) - 2)
             for i, v in enumerate(x)), arb(0)) / root2
    eta = a / shift
    mass_upper = arb(m - n) / (arb(n) * m)
    scale = 4 * n * (n + 1)
    mass = (mass_upper * scale / (scale + 1)).union(mass_upper)
    norm2 = head_norm2 + eta**2 * mass
    tail_norm = eta * mass_upper.sqrt()

    p_bound = (8 + 1 / arb(n) + 2 / arb(1).exp()) / arb(n).sqrt()
    d_bound = 2 / arb(n).sqrt()
    beta2 = sum((arb(2 * (k + 1)) / ((2 * i + 1) * k**2)
                 for i, k in enumerate(range(n, 0, -1))), arb(0))
    q_lo = ((head_energy + (2 * a * eta - d_bound * eta**2) * mass) / norm2).lower()
    q_hi = ((head_energy + (2 * (a + e) * eta + p_bound * eta**2) * mass) / norm2).upper()
    if not (0 < q_lo < q_hi and b < q_lo):
        raise ArithmeticError('Temple separation or Rayleigh endpoint ordering failed')
    residual_head = head_residual + beta2.sqrt() * tail_norm
    residual_tail = (e**2 * mass_upper + (a + e)**2 / m).sqrt() + (p_bound + d_bound) * tail_norm
    variance_upper = (residual_head**2 + residual_tail**2) / norm2
    rho_hi = (q_hi + variance_upper / (q_lo - b)).upper()
    return q_lo, rho_hi, b, q_hi - q_lo, variance_upper


def decimal_units(endpoint, outward_up):
    scale = 10**PLACES
    value = endpoint * scale
    integer = (value.upper().ceil() if outward_up else value.lower().floor()).unique_fmpz()
    if integer is None:
        raise ArithmeticError('outward decimal rounding did not isolate an integer')
    return int(integer)


def format_units(units):
    integer, fraction = divmod(units, 10**PLACES)
    return f'{integer}.{fraction:0{PLACES}d}'


def main():
    start = time.perf_counter()
    rho_lo, rho_hi, b, rayleigh_width, variance = compute()
    factor = 2 * arb.pi()
    lo = decimal_units(rho_lo**2 / factor, False)
    hi = decimal_units(rho_hi**2 / factor, True)
    if not (0 < lo < hi and (hi - lo) * 10**6 <= lo):
        raise ArithmeticError('outward printed endpoints fail the 1e-6 relative-width target')
    print(f'gamma=1/2 head={HEAD} support={SUPPORT}  k_B T_lin/g in [{format_units(lo)}, {format_units(hi)}]')
    print(f'relative_width={((arb(hi - lo) / lo)).str(10)}  lambda_2<={b.str(8)}')
    print(f'rayleigh_enclosure_width={rayleigh_width.str(10)}  variance_upper<={variance.upper().str(10)}')
    print(f'elapsed={time.perf_counter() - start:.2f}s')


if __name__ == '__main__':
    main()
