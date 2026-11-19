# Gamma-model critical temperatures: certified enclosures

**CERTIFIED, conditional on the analytic tail estimates in `RIGOR.md`:** the intervals below enclose $k_BT_c/g$ for the linearized isotropic gamma model. All five have relative width at most $10^{-6}$. Printed endpoints are rounded outward. EA is the upper bound of Elezaby–Abanov, arXiv:2512.20009v2, Eq. 44, evaluated here. The improvement factor is $(\mathrm{EA}-L)/(U-L)$, not a ratio of two published interval widths.

| gamma | L | U | EA upper | Previous literature value | Improvement factor |
|---:|---:|---:|---:|---|---:|
| 0.5 | 0.8094383586859154 | 0.8094386018859510 | 1.412472180812331 | numerical curve in KAY and EA | $2.48\times10^6$ |
| 1 | 0.255836478838462 | 0.255836590500742 | 0.318309886183791 | numerical plots; no precise tabulation verified | $5.59\times10^5$ |
| 1.5 | 0.199990970667518 | 0.199990972152204 | 0.225747978483961 | numerical plots; no precise tabulation verified | $1.73\times10^7$ |
| 2 | 0.1827262477723295078215179 | 0.1827262477723295176168720 | 0.197837817171107 | prior supplied certified interval $[0.182726247772231,0.182726247772398]$ | $1.54\times10^{15}$ |
| 3 | 0.17099873689534872 | 0.17099873689534873 | 0.178899776239800 | numerical plots; no precise tabulation verified | $7.90\times10^{14}$ |

**CERTIFIED:** relative widths from the printed endpoints are approximately $3.01\times10^{-7}$, $4.37\times10^{-7}$, $7.43\times10^{-9}$, $5.37\times10^{-17}$, and $5.85\times10^{-17}$, respectively. The absolute width for gamma=2 is below $9.80\times10^{-18}$. A supplemental 16,384-site Arb Rayleigh lower bound above $0.1827262477723295163523252$ and the certified upper bound together establish that $C_\infty$ rounds to **$0.18272624777232952$ at 17 decimal places**.

**PROVED:** with $g^2=\lambda\Omega^2$ the gamma=2 normalization is the strong-coupling Einstein ceiling. The previously circulated $0.18272624777235$ is a rounded upper-end figure, not the best estimate of the constant. The zero-transfer thermal term cancels between gap and self-energy; neither the nonzero-transfer diagonal Hankel contribution nor the negative mass term is dropped.
 
**NUMERICAL:** A blind independent calculation reported $\gamma=0.5$: $0.80943841333914$ (extrapolated), $\gamma=1$: $0.255836584718537$, $\gamma=1.5$: $0.1999909720339806$, $\gamma=2$: $0.182726247772329537\pm0.000000000000000022$, and $\gamma=3$: $0.1709987368953487281277$. These are cross-checks, not rigorous enclosures; the stated gamma=2 uncertainty overlaps the tighter interval.

**PROVED:** $T_c/g\to\infty$ as $\gamma\downarrow0$ and $T_c/g\to 1/(2\pi)$ as $\gamma\to\infty$. No monotonicity in gamma is claimed.

**[INFERENCE] Literature gate:** the three specified papers were read; searches located the earlier Moon–Chubukov and Wu–Abanov–Chubukov series. KAY and Elezaby–Abanov had numerical results at gamma=0.5; the contribution here is a certified enclosure, not its first numerical evaluation. No certified $10^{-6}$-relative table for the entire requested set was found; this is not an exhaustive priority search. Gamma=2 was already tightly certified in the previously supplied interval, so the new bound is a refinement, not a new physical constant. Kondratev's coarser exact-arithmetic eigenvalue enclosure is described in `RIGOR.md`.

**CERTIFIED, replay:** from this directory run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 uv run --no-project --with numpy==2.5.5 --with scipy --with python-flint==0.9.0 python -u replay.py`. The complete guarded replay of code commit `685537b3fa` exited 0 on the Mac in 13.86 s wall time, including the exact invalid-$U$ injection and the 17-place rounding check; `mac-replay.log` preserves its output. The baseline independent generalized-Hessian enclosures contain the tight half and two enclosures, but are not comparably narrow. The thermal mutant is an independent two-site witness rather than an injected matrix-constructor fault.

**PROVED, finite core only:** `C:/Users/davis/.elan/bin/lake.exe env lean D:/dev/scratch/ceiling-20261118/GAMMA/FiniteCore.lean`, from `C:/dev/.cache/lean/mathlib-v4.35.0`, exited 0 in 77.58 s. Lean checks the finite two-block bound and two-site small-gamma algebra, not the infinite tail. `#print axioms` reports only Mathlib's standard `propext`, `Classical.choice`, and `Quot.sound`.

**[INFERENCE] Publication verdict:** the fully enclosed table and reusable certificates may support a short KAY-program addendum; whether they merit a standalone article needs an exhaustive priority check and scientific review beyond the numerical precision improvement. The second implementation independently corroborates but does not match the narrow widths at gamma=0.5 and 2. No monotonicity proof or gamma=0.25 certificate is claimed. No cloud resources or external actions were used.
