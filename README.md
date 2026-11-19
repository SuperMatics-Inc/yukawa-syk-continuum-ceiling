# Replay supplement for a sharp bound on pair formation in the continuum Yukawa–SYK saddle

Davis Rens, SuperMatics Inc. ([ORCID](https://orcid.org/0000-0001-5839-5197)).

This repository accompanies the preprint [A sharp bound on pair formation in the continuum Yukawa–SYK saddle](https://doi.org/10.5281/zenodo.24366469). The archived software release is [10.5281/zenodo.24366471](https://doi.org/10.5281/zenodo.24366471). The DOI links resolve after the records are published. The interval for the gamma-model constant relies on a companion certificate, [10.5281/zenodo.24366341](https://doi.org/10.5281/zenodo.24366341), whose pinned sources are bundled here in `gamma/`.

## Scope

For the regularized infinite-band normal-state equation of Gnezdilov and Boyack, the preprint proves that the linearized pairing scale satisfies $k_B T_\mathrm{lin} < (C_1/2\pi)E$ at every finite coupling, with the limit attained as the coupling goes to infinity. Here $E = k_D^2/(2m_0)$ equals the Fermi energy exactly when $k_D = k_F$, and $C_1$ is the critical-temperature constant of the scalar gamma model at exponent one.

This is a statement about a linear instability at $N = \infty$ in a model. It is not an equilibrium transition temperature, a prediction for any material, or a measurement. It does not control a first-order transition, a competing saddle, the lower band edge, or phase coherence at finite $N$.

| Item | Label | Evidence in this repository |
|---|---|---|
| The strict bound at every finite coupling and its attainment as the coupling goes to infinity | PROVED | `paper/main.tex`, Section 1. The proof uses the scalar transfer-comparison theorem of the cited Kondratev and Rens works as an external input. |
| $0.0407176402 < C_1/(2\pi) < 0.0407176581$ | CERTIFIED, conditional on the analytic tail proof in the companion gamma certificate | `gamma/` (interval-arithmetic replay) and `replay.py` (exact rational division of the imported endpoints $[0.255836478838462, 0.255836590500742]$ for $C_1$). |
| Exact rational identities: f-sum, bubble sign, factor and shift, hopping, exact logarithm, gamma = 1 division, dot factors | CERTIFIED (exact `Fraction` arithmetic) | `replay.py`, with 12 mutants that must fail, run by `check_mutants.py`. |
| Five finite rational identities | CERTIFIED (Lean 4.35.0, Mathlib at the pinned commit) | `lean/FiniteCore.lean`. This is not a formalization of the spectral argument. |
| The plateau near $0.041E$, the three-frequency value $0.038E$, and the dot estimate $0.1827\sqrt{3/8}\,\Omega_0$ | Literature values, not recomputed here | Cited by DOI in the paper. |

The local `replay.py` does not recompute $C_1$. It takes the endpoints above as inputs and divides them outward. The upstream computation is the one replayed in `gamma/`, and its analytic tail proof is written in `gamma/RIGOR.md`.

## Contents

- `paper/`: LaTeX source and bibliography of the preprint.
- `replay.py`: exact-arithmetic checks for the continuum note, with `--mutant NAME` options.
- `check_mutants.py`: runs the baseline replay and every mutant.
- `lean/`: Lake project with `FiniteCore.lean` (five finite rational identities over $\mathbb{Q}$ used in the cross-multiplication, hopping, overlap and dot steps).
- `gamma/`: byte-exact snapshot of the gamma-model certificate sources and the log of one complete replay.
- `RIGOR.md`: the research record behind the preprint, byte-exact, because its hash is pinned.
- `MANIFEST.sha256`: SHA-256 of every pinned file. Check with `sha256sum -c MANIFEST.sha256`.

`RIGOR.md` covers more than the preprint. Only its Section 6 (sharp continuum saturation) and the preprint are claimed in this release. Its other sections (the finite-band case, the dot and the lattice) are an unedited research record, depend on a separate finite-band manuscript that is not part of this release, and contain local workspace paths and arXiv identifiers as written. They are kept because the certificate pins their bytes.

## Replay

Everything except the gamma replay uses only the Python 3 standard library.

```sh
python -B replay.py
python -B check_mutants.py
sha256sum -c MANIFEST.sha256
```

The gamma-model certificate, from `gamma/`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 \
  uv run --no-project --with numpy==2.5.5 --with scipy --with python-flint==0.9.0 python -u replay.py
```

It reads the printed endpoints from `gamma/RESULT.md` and prints PASS lines. `gamma/mac-replay.log` preserves one complete run (13.86 s). `gamma/RESULT.md` and `gamma/RIGOR.md` record the upstream commit and the proof status.

Lean, from `lean/` with elan installed:

```sh
lake exe cache get
lake build
```

The toolchain is `leanprover/lean4:v4.35.0`, and `lake-manifest.json` pins Mathlib to `56a005855ea17e0b0958f73ec801ae140b0a1679`. `FiniteCore.lean` contains no `sorry`, no axiom and no `native_decide`.

The paper, with a TeX Live installation that provides `pdflatex` and `bibtex`, from `paper/`:

```sh
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

Continuous integration (`.github/workflows/ci.yml`) runs the hash check, the rational replay, the mutants, the gamma replay, the Lean build and the paper build on standard GitHub-hosted runners, with no uploaded artifacts and no GitHub cache. The Lean job downloads the public Mathlib build cache.

## Licenses and provenance

The replay code and Lean sources are licensed under Apache-2.0 (see `LICENSE`). The preprint is distributed under CC-BY-4.0 with its Zenodo record. The files in `gamma/` are copies of the gamma-model certificate record, which carries its own Zenodo record.

AI agents contributed proofs, certificates, code, Lean sources and manuscript drafts and edits under the author's direction and review.
