"""Independent Hessian-basis trial and infinite-tail enclosure implementation."""
import numpy as np
from scipy.linalg import eigh
from flint import arb, arb_mat, fmpq, ctx

ctx.prec = 160

def power(a, g):
    return (arb(a).log()*g).exp()

def check(gn, gd, n):
    g = arb(fmpq(gn, gd))
    # Generalized eigenproblem K x = rho W x, W = diag(2i+1).
    p = np.arange(1, 2*n+1, dtype=float)**(-gn/gd)
    h = np.cumsum(p[:n])
    K = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            entry = p[i+j]
            if i != j:
                entry += p[j-i-1]
            else:
                entry -= 2*(h[i-1] if i else 0)
            K[i,j] = K[j,i] = entry
    W = np.diag(2*np.arange(n)+1)
    _, eigvec = eigh(K, W, subset_by_index=[n-1,n-1], driver='gvx')
    x = eigvec[:,0]
    if x[0] < 0: x = -x
    if np.any(x <= 0): raise ArithmeticError('nonpositive trial vector')
    # No floating-point arithmetic enters the Rayleigh quotient.
    pp = [arb(0)] + [power(k, -g) for k in range(1,2*n+1)]
    xx = [arb(float(z)) for z in x]
    numerator = arb(0)
    denominator = arb(0)
    hball = arb(0)
    roots = [1/arb(2*i+1).sqrt() for i in range(n)]
    M = arb_mat(n,n)
    for i in range(n):
        if i: hball += pp[i]
        numerator += (pp[2*i+1]-2*hball)*xx[i]**2
        M[i,i] = (pp[2*i+1]-2*hball)*roots[i]**2
        denominator += (2*i+1)*xx[i]**2
        for j in range(i+1, n):
            numerator += 2*(pp[j-i]+pp[i+j+1])*xx[i]*xx[j]
            M[i,j] = M[j,i] = (pp[j-i]+pp[i+j+1])*roots[i]*roots[j]
    rho = numerator/denominator
    # Bound the second full-operator eigenvalue by rank-one deflation, but
    # verify the finite factorization using a maximum absolute row-sum norm.
    y = arb_mat([[xx[i]/roots[i]] for i in range(n)])
    y = y/(y.transpose()*y)[0,0].sqrt()
    shift = 4*float(rho.mid())
    P = M-arb(shift)*(y*y.transpose())
    yf = x*np.sqrt(2*np.arange(n)+1)
    yf = yf/np.linalg.norm(yf)
    Gf = K/np.sqrt(np.outer(2*np.arange(n)+1,2*np.arange(n)+1))
    Pf = Gf-shift*np.outer(yf,yf)
    cutoff = float(np.linalg.eigvalsh(Pf)[-1]+1e-6)
    margin = 5e-7
    chol = arb_mat(np.linalg.cholesky((cutoff-margin)*np.eye(n)-Pf).tolist())
    err = -P-chol*chol.transpose()
    for i in range(n): err[i,i] += arb(cutoff)-arb(margin)
    for i in range(n):
        if not sum((abs(err[i,j]) for j in range(n)),arb(0)) < arb(margin):
            raise ArithmeticError('independent factorization failed')
    t = arb(n)
    a = power(2,g)*power(t,-g)*(1+1/(g*arb(1).exp()))
    if g < 1:
        b = 2/t+(2/(1-g)+2)*power(t,-g)
    elif g == 1:
        b = (4+2*t.log())/t
    else:
        b = 2*g.zeta()/t+2*power(t,-g)
    c = a+b+power(t,-g)/g
    def cross(i):
        u = arb(n-i)
        exponent = 2*g+1
        return 2*(power(u,-exponent)+power(u,1-exponent)/(exponent-1))/(2*i+1)
    coupling = sum((cross(i) for i in range(n)),arb(0))
    second = (cutoff+c)/2+(((cutoff-c)/2)**2+coupling).sqrt()
    if not second < rho: raise ArithmeticError('second eigenvalue not separated')
    residual = M*y-rho*y
    r2 = sum((residual[i,0]**2 for i in range(n)),arb(0))
    rt = sum((abs(y[i,0])*cross(i).sqrt() for i in range(n)),arb(0))
    top = rho+(r2+rt**2)/(rho-second)
    lo = power(rho, 1/g)/(2*arb.pi())
    hi = power(top, 1/g)/(2*arb.pi())
    print(f'gamma={gn}/{gd} n={n} independent_lower={lo.str(18)} independent_upper={hi.str(18)}', flush=True)
    return lo, hi

if __name__ == '__main__':
    for gn, gd, n in [(1,2,384),(1,1,1400),(3,2,384),(2,1,384),(3,1,384)]:
        check(gn, gd, n)
