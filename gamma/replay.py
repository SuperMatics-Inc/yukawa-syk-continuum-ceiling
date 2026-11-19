"""Replay both implementations and explicit adversarial checks.

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --no-project --with numpy --with scipy --with python-flint python replay.py
The tight gamma=1/2 and gamma=2 enclosures are replayed after the baseline checks.
"""
from pathlib import Path

from flint import arb
import preliminary as p
import independent as ind
import half_tight
import two_tight
import round17

# These are the outward endpoints printed in RESULT.md, not rounded
# midpoints of the balls returned by the certificate implementations.
published = {
    '0.5': ('0.8094383586859154', '0.8094386018859510'),
    '1': ('0.255836478838462', '0.255836590500742'),
    '1.5': ('0.199990970667518', '0.199990972152204'),
    '2': ('0.1827262477723295078215179', '0.1827262477723295176168720'),
    '3': ('0.17099873689534872', '0.17099873689534873'),
}
result_text = (Path(__file__).resolve().parent / 'RESULT.md').read_text(encoding='utf-8')
for label, endpoints in published.items():
    rows = [line for line in result_text.splitlines() if line.startswith(f'| {label} |')]
    if len(rows) != 1 or tuple(part.strip() for part in rows[0].split('|')[2:4]) != endpoints:
        raise ArithmeticError(f'RESULT.md printed endpoints changed for gamma={label}')


def check_published(label, lo, hi):
    outer_lo, outer_hi = (arb(x) for x in published[label])
    if not (outer_lo < lo < hi < outer_hi):
        raise ArithmeticError(f'gamma={label} is not strictly inside its printed endpoints')
    print(f'PASS gamma={label} printed endpoint containment')


def reject_low_u_mutation():
    # Inject the reviewer's exact U at the actual Temple assignment, in a
    # fresh copy of the module.
    path = Path(two_tight.__file__)
    source = path.read_text(encoding='utf-8')
    assignment = 'second = gap_bound(proposed, rho)'
    if source.count(assignment) != 1:
        raise ArithmeticError('low-U mutation injection point changed')
    bad_u = round17.reviewer_bad_u()
    source = source.replace(assignment, 'second = bad_u', 1)
    mutant = dict(vars(two_tight))
    mutant['bad_u'] = bad_u
    exec(compile(source, str(path), 'exec'), mutant)
    try:
        mutant['certify']()
    except mutant['_SecondBoundError'] as error:
        if not (error.actual is bad_u and error.actual < error.validated):
            raise ArithmeticError('low-U regression did not reach the Temple guard') from error
    else:
        raise ArithmeticError('injected low-U Temple mutant survived')
    print('REJECT injected low-U Temple mutant')


# The guard must also reject equality and overlapping/indeterminate balls.
for actual, validated in [(arb(1), arb(1)), (arb(1, arb('1e-10')), arb(1))]:
    try:
        two_tight._require_strict_bound(actual, validated)
    except two_tight._SecondBoundError:
        pass
    else:
        raise ArithmeticError('non-strict or indeterminate U comparison accepted')

results = {}
for gn, gd, n0, n1 in [(1,2,192,384),(1,1,1200,1400),(3,2,192,384),(2,1,192,384),(3,1,192,384)]:
    p.N = n0
    lo, hi = p.compute(gn,gd)
    lo2, hi2 = ind.check(gn,gd,n1)
    if not (lo < hi and lo2 < hi2):
        raise ArithmeticError('reversed or indeterminate enclosure endpoints')
    if not (lo < hi2 and lo2 < hi):
        raise ArithmeticError('independent enclosures do not intersect')
    if gn/gd >= 1 and not ((hi-lo)/lo < arb('0.000001') and
                            (hi2-lo2)/lo2 < arb('0.000001')):
        raise ArithmeticError('relative-width requirement failed')
    results[(gn,gd)] = (lo,hi,lo2,hi2)

for label, key in [('1', (1, 1)), ('1.5', (3, 2)), ('3', (3, 1))]:
    check_published(label, *results[key][2:])

# Independently cited exact-arithmetic eigenvalue enclosure (Kondratev,
# as reported in arXiv:2611.20496v1, section 1.4). This absolute anchor
# detects a common incorrect temperature normalization in both programs.
kondratev_lo = arb('1.318139').sqrt() / (2 * arb.pi())
kondratev_hi = arb('1.31815').sqrt() / (2 * arb.pi())
for endpoint in results[(2,1)]:
    if not kondratev_lo < endpoint < kondratev_hi:
        raise ArithmeticError('gamma=2 absolute normalization anchor failed')
print('PASS gamma=2 absolute normalization anchor')
# Wrong temperature normalization, 1/pi instead of 1/(2*pi).
lo,hi,lo2,hi2 = results[(2,1)]
if not 2*lo > hi2: raise ArithmeticError('normalization mutant survived')
print('REJECT normalization mutant')

# Unbalanced removal of the self-energy partner of a regulated n=m term;
# the pure same-frequency term itself cancels when both partners are kept.
# Unit regulator leaves +1/(2n+1) on the diagonal. A two-vector Rayleigh
# lower bound on the mutant must exceed the unmutated independent upper.
rho_mutant = (arb(2) + (arb(3)**(-2)-2)/3 + arb(1)/3 + 2*(1+arb(2)**(-2))/arb(3).sqrt())/2
t_mutant = rho_mutant.sqrt()/(2*arb.pi())
if not t_mutant > hi2: raise ArithmeticError('thermal-cancellation mutant survived')
print('REJECT unbalanced thermal-cancellation mutant')

# Deliberately halve the squared tail residual in the primary Temple bound.
# This creates an invalid upper bound below the independently computed lower.
p.N = 192
_, bad_hi = p.compute(2,1,tail_factor=0.5,verbose=False)
if not bad_hi < lo2: raise ArithmeticError('half-tail mutant survived')
print('REJECT half-tail-residual mutant')
half_lo_rho, half_hi_rho, _, _, _ = half_tight.compute()
half_lo = half_lo_rho**2 / (2 * arb.pi())
half_hi = half_hi_rho**2 / (2 * arb.pi())
if not (half_lo < half_hi and (half_hi-half_lo)/half_lo < arb('1e-6')):
    raise ArithmeticError('tight gamma=1/2 width requirement failed')
if not (half_lo > results[(1,2)][2] and half_hi < results[(1,2)][3]):
    raise ArithmeticError('tight gamma=1/2 is not inside independent enclosure')
print('PASS gamma=1/2 tight interval:', half_lo.str(18), half_hi.str(18))
check_published('0.5', half_lo, half_hi)
two_lo, two_hi, _ = two_tight.certify()
if not (two_lo < two_hi and two_hi-two_lo < arb('9.80e-18')):
    raise ArithmeticError('tight gamma=2 width requirement failed')
if not (two_lo > lo2 and two_hi < hi2):
    raise ArithmeticError('tight gamma=2 is not inside independent enclosure')
check_published('2', two_lo, two_hi)
if not arb(published['2'][1]) - arb(published['2'][0]) < arb('9.80e-18'):
    raise ArithmeticError('printed gamma=2 absolute width requirement failed')
reject_low_u_mutation()
round17.certify()
print('PASS replay: all five relative-width targets and gamma=2 absolute width')
