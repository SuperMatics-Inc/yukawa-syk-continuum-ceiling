# Gamma-model certificate: proof map and limitations

## Scope and sources
**PROVED (conditional on the cited KAY operator theorem).** All temperatures here are the normal-state **linear-instability** temperature; identifying it with a thermodynamic phase transition assumes the standard continuous-transition hypothesis. No material measurement is asserted. KAY, [arXiv:2409.00533v2](https://arxiv.org/abs/2409.00533), Eqs. 42–46, proves compact self-adjointness and the eigenvalue characterization; [arXiv:2409.00532v2](https://arxiv.org/abs/2409.00532), Eq. 1 and its discussion, gives the dispersive-phonon normalization and the gamma=2 matching. Elezaby–Abanov [arXiv:2512.20009v2](https://arxiv.org/abs/2512.20009), Eq. 44, gives the comparison upper bound. The three papers were fully read via the legal paperpipe pipeline; validated briefs are in `D:/dev/.cache/paperpipe/papers/doi_10.1007_s10955-025-03446-5/brief.md`, `.../doi_10.1007_s10955-025-03468-z/brief.md`, and `.../doi_10.1103_wpn8-8j2t/brief.md`.

Let `w_n=2n+1`, `H_n=Σ_{k=1}^n k^{-γ}`. For `n,m≥0`,
```
G_nm = [(1-δ_nm)|n-m|^{-γ} + (n+m+1)^{-γ}]/sqrt(w_n w_m)
       - δ_nm 2H_n/w_n,
ρ = sup spec G,   k_B T_lin/g = ρ^{1/γ}/(2π).
```
**PROVED (source equation).** The same-frequency thermal term cancels between the self-energy and the pairing vertex before the gamma kernel is assigned zero on the diagonal. It is *not* the diagonal Hankel term `(2n+1)^{-γ}/w_n`, which stays, nor the finite diagonal self-energy `-2H_n/w_n`, which also stays. The original singular static terms must not be given independent unbalanced values. At gamma=2, `g²=λΩ²` for Einstein phonons; the fixed-g, Ω→0 kernel is the gamma kernel off the canceled diagonal. KAY reports `0.1827262477…` from stabilization of two independent finite-section computations, not a ten-digit rigorous enclosure in the cited papers.

## Infinite tail and finite matrix
**PROVED (analytic reduction; independently reviewed).** A read-only independent mathematical reviewer checked the Schur, cross-block, deflation and Temple inequalities and reported no findings, conditional on KAY, the spectral/min–max theorems and Arb's interval guarantees. For `N≥2`, write `G=[[A,B],[B*,C]]`. Let `P` be `C` with its negative diagonal term discarded. It is entrywise nonnegative, and `C≤P` in quadratic-form sense. The Schur test with positive weights `x_n=1/sqrt(2n+1)` gives `sup spec C ≤ c_N`, where
```
a = 2^γ N^{-γ}(1 + 1/(eγ)),    d = N^{-γ}/γ,
b = 2/N + (2/(1-γ)+2)N^{-γ}                for 0<γ<1;
b = [4+2 log N]/N                            for γ=1;
b = 2ζ(γ)/N + 2N^{-γ}                        for γ>1;
c_N = a+b+d.
```
For a row `n≥N`, split tail columns into `m≤n/2`, `n/2<m≤2n`, and `m>2n`. The left term is bounded by `2^γ n^{-γ}(1+log(n/N))`; its supremum is at most `a`. In the middle, `1/w_m≤1/n`, the difference part sums to at most twice the harmonic sum through `n`, and the Hankel part has at most `2n` terms each at most `n^{-γ}`; the standard integral bounds on harmonic sums give `b`. On the right both kernels are at most `(m/2)^{-γ}`; using `1/w_m≤1/(2m)` and the integral test gives `d`. This is an all-infinite tail estimate, including `γ≤1`; it does not sum the divergent ζ(γ) there.

For the cross block, `k=m-i≥N-i` and `m≥N` imply
```
Σ_{m≥N} |G_im|² ≤ R_i² := 2/(2i+1) * [ (N-i)^(-2γ-1) + (N-i)^(-2γ)/(2γ) ].
||B||² ≤ β² := Σ_{i<N} R_i².
```
The bound uses `(m+i+1)^{-γ}≤(m-i)^{-γ}`, `1/(2m+1)≤1/[2(m-i)]`, and a power-series integral test. All signs are one-sided.

**PROVED (finite algebra + analytic spectral theorem).** Both implementations select a positive finite vector `v`, normalize it, evaluate `q=<v,Av>` with Arb, and take `r²=||(A-q)v||² + (Σ_i |v_i| R_i)²`. The tail part bounds `||B*v||²` by the triangle inequality. To bound the second eigenvalue, compute `P_A=A-d vv*` for a positive rational/binary64 `d`. A numerically proposed Cholesky factor is checked with Arb: `(s-margin)I-P_A=LL*+E` and `||E||<margin`, so `P_A≤sI`. One implementation uses an Arb Frobenius error norm; the independent generalized-Hessian implementation uses the maximum absolute row sum. By min–max, `λ₂(G)≤sup spec(G-d vv*)`; the latter is at most
```
b₂ = (s+c_N)/2 + sqrt(((s-c_N)/2)^2 + β²).
```
If `q>b₂`, the Temple/residual inequality yields `q≤ρ≤q+r²/(q-b₂)`. The finite block is the exact `N×N` principal block, not a finite-cutoff self-energy. The lower bound is the exact Rayleigh quotient of the positive Arb-represented test vector. The conversion `ρ↦ρ^(1/γ)/(2π)` is monotone on positive numbers. `FiniteCore.lean` checks the two-block algebra and the two-site gamma-down-to-zero Rayleigh identity; the infinite-dimensional Schur, min–max, spectral, power-tail and interval-library results are explicit trusted analytic inputs, **not Lean-checked**. `#print axioms` lists only Lean/Mathlib's standard `propext`, `Classical.choice`, `Quot.sound`, no custom axioms or `sorry`.

## Tight gamma=1/2: implicit long trial and explicit tail lemma
**PROVED (conditional on KAY, Arb and standard spectral theorems).** The prior generic bounds remain valid but are too wide at gamma=1/2. `half_tight.py` instead splits `G=[[A,B],[B*,C]]` at `N=2^17`. It uses a positive numerical head `h`, and an exact finite trial tail `z_m=η/(m+1/2)` for `N≤m<M=2^34`, zero for `m≥M`. The support size `M` is implicit. Let `t_i=i+1/2`, `x_i=h_i/sqrt(2i+1)`, `a=sqrt(2) Σ_{i<N} x_i`, and `f_m=(B*h)_m`. For `m≥N`,
```
t_m f_m = Σ_i (x_i/sqrt(2)) [(1-t_i/t_m)^(-1/2)+(1+t_i/t_m)^(-1/2)].
```
The bracket is increasing for `0≤t_i/t_m<1` and equals 2 at zero. Hence, with `e=t_N f_N-a`, the **all-infinite forcing bound** is `a/t_m≤f_m≤(a+e)/t_m`. The code computes `e` with Arb and chooses any rigorously positive numerical shift `λ` and `η=a/λ`; no eigenvector convergence is assumed. Set `D=Σ_{m=N}^{M-1}t_m^-2` and `D+=1/N-1/M`. Since `t_m^-2=(1/m-1/(m+1))*[4m(m+1)/(4m(m+1)+1)]`,
```
D+ * 4N(N+1)/(4N(N+1)+1) ≤ D ≤ D+,
Σ_{m≥M} t_m^-2 ≤ 1/M.
```
For the tail block `C=P-diag(d_n)`, `P` is entrywise nonnegative. Schur weights `1/sqrt(2n+1)` and row ranges `m≤n/2`, `n/2<m≤2n`, `m>2n` give respective upper bounds `(1/N+2/e_Euler)/sqrt(N)`, `6/sqrt(N)`, `2/sqrt(N)`; the first uses `max_{u≥1} log(u)/sqrt(u)=2/e_Euler`. Therefore `||P||≤cP=(8+1/N+2/e_Euler)/sqrt(N)`. From `Σ_{k≤n} k^-1/2≤2sqrt(n)-1`, `0≤d_n≤d=2/sqrt(N)`. Since `z≥0`, `-d||z||²≤<z,Cz>≤cP||z||²` and `||C||≤cP+d`.

**PROVED (second-eigenvalue bound).** Deleting index 0 gives, by min–max, an upper bound on `λ₂(G)`. On this compression, weighted AM–GM bounds its form by the supremum of signed weighted row sums, including the negative diagonal. The code checks rows `1≤n<192` through column 384 with the infinite remainder `2sqrt(2/384)`; for `n≥192` the sum is bounded by the decreasing function `[sqrt(2)(1+log(n/2))+8]/sqrt(n)`. Arb verifies every bound is strictly below `b=5/4`. This uses no square summability of the weights.

**PROVED (interval Rayleigh and residual).** Let `H=||h||²`, `E=<h,Ah>`, and `β²=Σ_{i<N} [2/(2i+1)] [(N-i)^-2+(N-i)^-1]`, which bounds `||B||²`. The tail lemma gives the Rayleigh interval
```
q ∈ [(E+(2aη-dη²)D)/(H+η²D), (E+(2(a+e)η+cPη²)D)/(H+η²D)] =: [l,r].
R_head ≤ ||(A-λ)h|| + βη sqrt(D+).
R_tail ≤ sqrt(e² D+ + (a+e)²/M) + (cP+d)η sqrt(D+).
```
The two terms inside the square root bound disjoint active and omitted tail residuals. `V=(R_head²+R_tail²)/(H+η²D)` bounds the variance about the exact Rayleigh quotient, because that quotient minimizes residual norm over shifts. If `l>b`, the Temple spectral-measure inequality gives `l≤ρ(G)≤r+V/(l-b)`. The full exact finite head action is enclosed with Arb polynomial convolution, while binary64 FFTs only propose `h`. This construction achieves the gamma=1/2 table interval without extrapolation.

**[INFERENCE] Independent review.** A read-only adversarial reviewer independently derived the paired forcing bound, all three Schur constants, signed-row second-eigenvalue bound, interval Rayleigh–Temple step, convolution indices and outward rounding; no blocking finding was reported. The reviewer did not run the code, and the proof still relies on KAY, standard spectral theory and Arb.

## Tight gamma=2: large finite trial with a small certified gap block
**PROVED (conditional on the same analytic inputs).** `two_tight.py` uses a positive `N=8192`-component finite test vector, with Toeplitz and Hankel actions computed by 160-bit Arb polynomial multiplication. It certifies the infinite second-eigenvalue bound using the separately deflated `192×192` principal block, an Arb-checked Cholesky residual, and the generic all-tail Schur/cross-block estimates above; it does not require a dense `8192×8192` factorization. The entire infinite trial residual beyond `N` is bounded by the gamma=2 instance of `R_i²` in Eq. (cross block), and Temple is applied with the exact finite Rayleigh quotient. The independently checked generalized-Hessian interval contains this tighter result but is substantially wider, so it is not an independent verification at the new `10^-17` scale.

**[INFERENCE] Independent review.** A read-only reviewer checked the gamma=2 convolution indexing, interval normalization, infinite boundary and tail estimates, prefix deflation, Cholesky residual and Temple argument, reporting no blocking finding. The reviewer caught an omitted margin in the earlier written proof map, corrected above; the implementation had used the correct shifted identity. The review did not run code and does not replace the shared KAY/Arb assumptions.

**CERTIFIED (17-place rounding, independent review).** `GAMMA-REVIEW/check.py` separately evaluated the 16,384-site Arb Rayleigh quotient, obtaining a lower bound above `0.1827262477723295163523252`; the gamma=2 certified upper lies below `0.182726247772329525`. These bounds lie strictly between the half-unit boundaries `0.182726247772329515` and `0.182726247772329525`, so the value rounds to `0.18272624777232952` to 17 decimal places. `round17.py` reproduces the enlarged finite-block lower bound. The eigensolver only proposes exact dyadic trial-vector entries; Arb evaluates the Rayleigh quotient.

## Arithmetic, comparison, replay, and adversaries
**CERTIFIED (subject to the stated proof and code).** `replay.py` uses 160-bit python-flint Arb throughout proof-relevant operations; floating point only proposes trial vectors, factorizations and spectral cutoffs. Baseline implementations use different finite dimensions and bases: compact symmetrized vs unweighted generalized Hessian, and different factorization error norms. They share the generic tail lemma and Temple theorem. Their enclosures contain the new tighter gamma=1/2 and gamma=2 intervals but are not comparably narrow; those tighter intervals do not have a second comparably precise implementation. All five final intervals have relative width below `10^-6`. The Elezaby–Abanov comparison uses `B_γ=Σ_{j≥0}4^{-j}ζ(γ+2j+1)` from Eq. 44, with the geometric remainder bounded using monotonicity of ζ. The table's improvement is `(EA_upper−our_L)/(our_U−our_L)`, not a ratio of two published intervals.

Replay from the workspace:
```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 uv run --no-project --with numpy==2.5.5 --with scipy --with python-flint==0.9.0 python -u replay.py
```
The script rejects three concrete adversaries: doubling the temperature normalization, cancellation of only one regulated same-frequency partner (the unit static regulator leaves a nonzero diagonal), and halving the squared residual tail in the gamma=2 upper bound. For the last, the mutated upper lies below the independently recomputed finite-section lower. These mutants do not prove that all possible tail, thermal, or normalization errors are caught. `preliminary.compute(..., tail_factor≠1)` exists *only* to create an invalid mutant.

**CERTIFIED (replay gap closure, subject to replay).** The reviewer found that an invalid second-eigenvalue bound could pass the previous interval-overlap gate. At the Temple call, `two_tight.py` now recomputes the all-infinite deflation/Schur/Arb-Cholesky bound from the operator and requires a strict outward inequality for the supplied bound; this recomputes the same analytic algorithm, not an independent implementation of it. `round17.reviewer_bad_u()` reconstructs the reviewer's $q_{\rm second}^{\rm lower}-10^{-12}$, and an Arb positive-definite two-dimensional witness proves it is below $\lambda_2$. The replay injects that exact bound into the actual Temple path and requires rejection there. It also requires the computed endpoints to lie strictly inside all five published decimal pairs and checks the $9.80\times10^{-18}$ gamma=2 absolute width.

**CERTIFIED (review and limitations).** A separate read-only code reviewer found no issue in the finite matrices, normalization formula or Arb bound directions but identified three replay-validation gaps: endpoint ordering, enforcement of the target on both implementations, and a common normalization error. The replay now checks each pair's order and width and requires both gamma=2 enclosures to lie within Kondratev's independently cited exact-arithmetic enclosure after conversion. The corrected replay exited 0. The reviewer did not rerun the corrected code. The thermal mutant is a two-site witness and does not inject a fault into either constructor; the half-tail mutant does exercise the changed certificate path.

**PROVED (limits).** The 2-site gamma→0 generalized Rayleigh quotient is `5/3`; finite-dimensional continuity makes `ρ(γ)>1+ε` for all sufficiently small positive gamma, so `T_lin/g→∞`. The 1-site Rayleigh bound gives `ρ≥1`. The Elezaby–Abanov bound obeys `B_γ≤(4/3)ζ(γ+1)`, hence `1≤(2πT_lin/g)≤B_γ^(1/γ)→1` as γ→∞. Thus the large-gamma limit is `1/(2π)`. These are asymptotic statements, not claims of novelty or a monotonicity theorem.

**PROVED (prior certified comparison).** The published [arXiv:2611.20496v1, §1.4](https://arxiv.org/html/2611.20496v1) attributes the exact-arithmetic enclosure `1.318139 < g(2) < 1.31815` to Kondratev's preprint; this is a relative width about `8.35×10⁻⁶` for the eigenvalue and about `4.17×10⁻⁶` after taking its square root for temperature. The published paper itself uses the eigenvalue's definition and does not claim that its displayed ten-digit decimal is a certified enclosure. This is a secondary-source citation to Kondratev's bound, not a claim that his entire preprint was independently checked here.

## Gaps
**[INFERENCE]** No monotonicity of `T_lin/g` versus gamma was established, and gamma=0.25 was not certified. The previous gamma=2 interval supplied by the program was already certified to `10^-6` relative precision; the new interval is a refinement. Its value is near `0.18272624777232951`, not the earlier rounded upper-end figure `0.18272624777235`. The two tighter scripts are not Lean formalized, nor are their full precision results duplicated by an independent implementation. No cloud resources or outward actions were used.
