"""Exact arithmetic checks; not a replacement for the analytic proofs in RIGOR.md."""
import argparse
from fractions import Fraction as Q

MUTANTS = ('spectral-half', 'bubble-sign', 'bubble-half', 'bubble-shift',
           'hopping-sign', 'log-linear', 'damping-sign', 'damping-half',
           'omit-transfer-index', 'plateau-no-2', 'dot-static-half', 'dot-bound-10')


def arctan(x, n):
    s = sum(((-1)**k * x**(2*k+1) / (2*k+1) for k in range(n)), Q(0))
    t = s + (-1)**n * x**(2*n+1) / (2*n+1)
    return min(s, t), max(s, t)


def log1p(x, n=20):
    # log(1+x)=2 atanh(x/(2+x)); bound the positive remainder geometrically.
    y = x / (2+x)
    s = 2 * sum((y**(2*k+1) / (2*k+1) for k in range(n)), Q(0))
    return s, s + 2*y**(2*n+1) / ((2*n+1)*(1-y*y))


def check(m):
    shi = Q(1, 2) if m == 'spectral-half' else Q(1)
    modes = ((Q(1, 3), Q(1)), (Q(2, 3), Q(3)))
    rho = tuple((shi*w/(2*o), o) for w, o in modes)
    assert sum(2*o*w for w, o in rho) == 1
    for s in (Q(1), Q(4), Q(25)):
        d = sum(2*o*w/(s+o**2) for w, o in rho)
        assert d == sum(w/(s+o**2) for w, o in modes)
        assert 0 < s*d < 1
        assert d <= 1/(s+1/sum(w/o**2 for w, o in modes))

    # Two equally weighted energy levels, f(-1)=3/4 and f(1)=1/4.
    spin_coefficient = -2 if m == 'bubble-sign' else 1 if m == 'bubble-half' else 2
    assert spin_coefficient*Q(1, 4)*2*Q(1, 2)/4 == Q(1, 8)
    assert spin_coefficient*Q(1, 4)*2*Q(1, 2)/(4+9) == Q(1, 26)
    v = (Q(0), Q(-1, 3), Q(-1), Q(1), Q(1, 3), Q(0))
    for j in (1, 2, 3):
        shifted = v[j:] + v[:j]
        diff = sum((x-y)**2 for x, y in zip(v, shifted))
        assert sum(x*x for x in v)-sum(x*y for x,y in zip(v, shifted)) == diff/(1 if m == 'bubble-shift' else 2)

    h, v = Q(1, 3), Q(2)
    u = 1/h-v*h
    assert h**2/(1+(1 if m == 'hopping-sign' else -1)*v*h**2) == h/u

    a5, b5 = arctan(Q(1, 5), 26)
    a239, b239 = arctan(Q(1, 239), 10)
    p_lo, p_hi = 16*a5-4*b239, 16*b5-4*a239
    assert Q('3.141592653589793238462643383279') < p_lo < p_hi < Q('3.141592653589793238462643383280')
    # Test exact transfer, enclosing log(1+x), without substituting its strong-coupling limit.
    # x uses pi in the Matsubara frequency AND in the damping coefficient.
    for temp in (Q(1, 37), Q(1, 9), Q(2)):
        for lam in (Q(1, 31), Q(7), Q(1000)):
            e, wd = Q(1), Q(1, 5)
            for j in (1, 2, 11):
                nus = (2*p_lo*temp*j, 2*p_hi*temp*j)
                denom = tuple((nu/wd)**2 + 2*p*lam*nu/e for nu,p in zip(nus,(p_lo,p_hi)))
                lo, hi = log1p(1/denom[1])[0], log1p(1/denom[0])[1]
                d = e*temp/(lam*wd**2)
                k = 1 if m == 'omit-transfer-index' else j
                factor = -1 if m == 'damping-sign' else Q(1, 2) if m == 'damping-half' else Q(1)
                for p, nu, den in zip((p_lo,p_hi), nus, denom):
                    f = e/(2*p*nu)
                    assert lam/den == f/(1+factor*d*k)
                hi = 1/denom[0] if m == 'log-linear' else hi
                assert 0 < lam*lo <= lam*hi < lam/denom[0]
                assert lam/denom[0] <= e/(2*p_lo*nus[0]*(1+d))
                assert lam*hi < lam*wd**2/nus[0]**2

    # Input enclosures from GAMMA/RESULT.md, conditional on its tail proof.
    c1_lo, c1_hi = Q('0.255836478838462'), Q('0.255836590500742')
    m2 = 1 if m == 'plateau-no-2' else 2
    lo, hi = c1_lo/(m2*p_hi), c1_hi/(m2*p_lo)
    assert Q('0.0407176402') < lo <= hi < Q('0.0407176581')
    c2_lo = Q('0.1827262477723295078215179')
    c2_hi = Q('0.1827262477723295176168720')
    # Integral of t(x)^2: after x=2 sinh(v), integrate exp(-v)+exp(-3v).
    integral = (Q(1)+Q(1,3))/(2 if m == 'dot-static-half' else 1)
    dot_squared_factor = 1/(2*integral)  # r=2
    assert dot_squared_factor == Q(3,8)
    assert Q('0.1118965174138946')**2 < c2_lo**2*dot_squared_factor < c2_hi**2*dot_squared_factor < Q('0.1118965174138948')**2
    gamma = 10 if m == 'dot-bound-10' else 5
    assert gamma*(2+Q(gamma,144)+Q(gamma,24)) < 12
    print('PASS: f-sum, bubble sign/factor/shift, hopping, exact log, gamma=1 division, dot factors')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mutant', choices=MUTANTS)
    check(parser.parse_args().mutant)
