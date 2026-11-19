"""Preliminary Arb finite-block/infinite-tail bounds, NOT the requested narrow certificate.
Run: uv run --no-project --with numpy --with scipy --with python-flint python preliminary.py
"""
import math
import numpy as np
from flint import arb, arb_mat, ctx, fmpq

ctx.prec = 160
N = 192

def q(a, b=1):
    return arb(fmpq(a, b))

def pw(a, x):
    return (arb(a).log()*x).exp()

def tail_power(a, s):
    return pw(a, -s) + pw(a, 1-s)/(s-1)

def tail_block(g, N):
    # Weighted Schur bound; proof in RIGOR.md.
    nn = arb(N)
    a = pw(2, g)*pw(nn, -g)*(1+1/(g*arb(1).exp()))
    c = pw(nn, -g)/g
    if g < 1:
        b = 2/nn+(2/(1-g)+2)*pw(nn, -g)
    elif g == 1:
        b = (2*(1+nn.log())+2)/nn
    else:
        b = 2*g.zeta()/nn+2*pw(nn, -g)
    return a+b+c

def compute(gn, gd, tail_factor=1, verbose=True):
    g = q(gn, gd)
    coef = [arb(0)] + [pw(k, -g) for k in range(1, 2*N+1)]
    root = [1/arb(2*n+1).sqrt() for n in range(N)]
    accum = arb(0)
    M = arb_mat(N, N)
    F = np.empty((N, N))
    for i in range(N):
        if i: accum += coef[i]
        for j in range(i, N):
            v = coef[i+j+1]*root[i]*root[j]
            if i != j: v += coef[j-i]*root[i]*root[j]
            else: v -= 2*accum/(2*i+1)
            M[i,j] = M[j,i] = v
            F[i,j] = F[j,i] = float(v.mid())
    eigenvalues, vectors = np.linalg.eigh(F)
    if vectors[0,-1] < 0:
        vectors[:,-1] *= -1
    if not np.all(vectors[:,-1] > 0):
        raise ArithmeticError('test vector is not positive')
    v = arb_mat([[float(z)] for z in vectors[:,-1]])
    lower = (v.transpose()*M*v)[0,0]/(v.transpose()*v)[0,0]
    # Cholesky residual in Frobenius norm bounds block spectrum.
    s = float(eigenvalues[-1] + 1e-8)
    tau = 5e-9
    L = arb_mat(np.linalg.cholesky((s-tau)*np.eye(N)-F).tolist())
    E = -M - L*L.transpose()
    for i in range(N): E[i,i] += arb(s)-arb(tau)
    fro2 = sum((E[i,j]*E[i,j] for i in range(N) for j in range(N)), arb(0))
    if not fro2 < arb(tau)**2: raise ArithmeticError('Cholesky residual not certified')
    beta2 = sum((2*tail_power(N-i, 2*g+1)/(2*i+1) for i in range(N)), arb(0))
    c = tail_block(g, N)
    if not c < arb(s): raise ArithmeticError('tail norm exceeds block bound')
    upper = (arb(s)+c)/2+(((arb(s)-c)/2)**2+beta2).sqrt()
    # A rank-one deflation bounds the second eigenvalue. Temple's inequality
    # then uses the residual of the trial vector, including its infinite tail.
    vv = v / (v.transpose()*v)[0,0].sqrt()
    P = M - 2*arb(s)*(vv*vv.transpose())
    pf = F - 2*s*np.outer(vectors[:,-1], vectors[:,-1])
    s2 = float(np.linalg.eigvalsh(pf)[-1] + 1e-8)
    L2 = arb_mat(np.linalg.cholesky((s2-tau)*np.eye(N)-pf).tolist())
    E2 = -P - L2*L2.transpose()
    for i in range(N): E2[i,i] += arb(s2)-arb(tau)
    fro22 = sum((E2[i,j]*E2[i,j] for i in range(N) for j in range(N)), arb(0))
    if not fro22 < arb(tau)**2: raise ArithmeticError('deflation residual not certified')
    second = (arb(s2)+c)/2+(((arb(s2)-c)/2)**2+beta2).sqrt()
    if not second < lower: raise ArithmeticError('second-eigenvalue bound is too large')
    r = M*vv-lower*vv
    r2 = sum((r[i,0]**2 for i in range(N)), arb(0))
    rt = sum((abs(vv[i,0])*(2*tail_power(N-i, 2*g+1)/(2*i+1)).sqrt()
              for i in range(N)), arb(0))
    upper_temple = lower + (r2+arb(tail_factor)*rt**2)/(lower-second)
    upper = upper_temple if upper_temple < upper else upper
    ea = sum((q(1,4**j)*(g+2*j+1).zeta() for j in range(80)), arb(0))
    ea += arb(0, (q(1,4**80)*(g+1).zeta()*q(4,3)).upper())
    lo = pw(lower, 1/g)/(2*arb.pi())
    hi = pw(upper, 1/g)/(2*arb.pi())
    pub = pw(ea, 1/g)/(2*arb.pi())
    if verbose:
        print(f'gamma={gn}/{gd} lower={lo.str(18)} upper={hi.str(18)} EA={pub.str(18)}')
        print('  relative_width=', ((hi-lo)/lo).str(8), 'EA_gap_improvement=',
              ((pub-lo)/(hi-lo)).str(8), 'tail=',c.str(8),'beta2=',beta2.str(8))
    return lo, hi

if __name__ == '__main__':
    for gn,gd in [(1,2),(1,1),(3,2),(2,1),(3,1)]:
        N = 1200 if (gn,gd)==(1,1) else 192
        compute(gn,gd)
