# Chapter 1 — Linear Algebra: Design Spec

**Date:** 2026-06-08
**Branch:** `chapter-1-linear-algebra`
**File to author:** `parts/01-foundations/01-linear-algebra.qmd`
**Canonical anchors:** Strang (2016), Axler (2015), Trefethen & Bau (1997).

## Goal

Replace the scaffolded stub with a complete chapter that reactivates undergraduate
linear algebra and ramps toward the rigor later MMM chapters require. Scope is
**MMM-targeted core**: every topic earns its place by a concrete downstream use in
regression, Bayesian inference, sampling, or optimization. Deliverable this cycle is
a **full chapter draft** — prose, proofs, worked examples, a runnable Python tie-in,
and full C/B/P/A exercises.

The chapter follows the book's fixed template:
Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises (C/B/P/A).

## Through-line (anchor example)

Every abstraction lands on one concrete object: the **design matrix `X` of a
media-mix regression**. Rows are weeks; columns are media-channel spends plus
controls (price, seasonality, intercept); `y` is weekly sales. The object recurs
through the whole chapter so motivation is structural, not bolted on:

- **Vectors/matrices** — a week is a row vector, a channel is a column vector, `X`
  is the data; the matrix–vector product `Xβ` is a *linear combination of the
  channel columns*.
- **Subspaces & rank** — the column space is "all sales patterns the media can
  explain"; rank deficiency *is* multicollinearity (two channels that always move
  together).
- **Projection & least squares** — fitting `ŷ = Xβ̂` is projecting `y` onto the
  column space; the normal equations `XᵀXβ̂ = Xᵀy` are the chapter's keystone.
- **Eigen / spectral & PD** — `XᵀX` is symmetric PSD; its eigenvalues diagnose
  collinearity (near-zero → unstable `β̂`) and later *are* the curvature (Hessian)
  of the regression loss.
- **SVD & Cholesky** — SVD as the numerically honest way to solve least squares;
  Cholesky of `XᵀX` (or a covariance matrix) as the factorization that powers
  Bayesian posteriors and HMC mass matrices downstream.

Collinearity, stable estimation, and covariance structure are the spine — exactly
what Parts II–IV build on.

## Section-by-section content

### Motivation
Open with the design matrix and one question: *"You ran TV and search together
every week. How much did each drive sales — and when can you even tell them
apart?"* That question is collinearity, and answering it honestly motivates the
whole chapter.

### Theory & Proofs
A ladder, each rung a downstream payoff:

1. Vectors, matrices, matrix–vector product as a **linear combination of columns**
   (Strang's framing).
2. Column space, null space, rank; rank deficiency ⟺ collinearity.
3. Orthogonality, projection onto a subspace; **the projection theorem → normal
   equations** (proved).
4. Eigenvalues/eigenvectors; **spectral theorem for symmetric matrices** (real
   symmetric ⇒ orthonormal eigenbasis, real eigenvalues) — proved at the real
   symmetric level.
5. Positive (semi)definiteness; **PD ⟺ all eigenvalues > 0** (proved); the
   quadratic form `βᵀXᵀXβ`.
6. SVD — stated with geometric meaning and its link to the eigendecomposition of
   `XᵀX`; existence **cited** (Trefethen & Bau), not proved.
7. Cholesky factorization of an SPD matrix — stated, with a "why it matters"
   pointer to sampling and Bayesian posteriors.

### Worked Examples
- (a) A 3-week, 2-channel `X`: compute `XᵀX`, solve the normal equations by hand.
- (b) Make the two channels nearly collinear and watch `XᵀX` become near-singular
  (tiny eigenvalue, blown-up `β̂`).

### Code Tie-in
NumPy / SciPy, runnable and self-contained:
- Build a small synthetic `X`.
- Solve least squares three ways: `np.linalg.solve(XᵀX, Xᵀy)`, `np.linalg.lstsq`,
  and via SVD.
- Show they agree on well-conditioned data and diverge in stability as collinearity
  rises; print eigenvalues / condition number.
- `scipy.linalg.cholesky` on `XᵀX` (SPD case), with a one-line teaser that this is
  how correlated draws are generated later.

### Exercises
- **C — Conceptual:** what is the column space of a media matrix; why does
  collinearity hurt estimation.
- **B — By hand:** small projections; eigenvalues of a 2×2 symmetric matrix.
- **P — Prove it:** normal equations minimize SSE; symmetric ⇒ orthogonal
  eigenvectors for distinct eigenvalues; PD ⟹ invertible.
- **A — Applied / code:** condition number vs. estimate variance; Cholesky-sampling
  teaser.

Exercise solutions go in the shared `appendix/solutions.qmd`, gated by the
`show-solutions` metadata flag (no inline links, per the preface).

## Rigor level

Chapter 1 is the *start* of the book's ramp, so it leans **Strang-computational** —
intuition and geometry first — while proving a handful of genuinely illuminating
results at full rigor (projection / normal equations, spectral theorem for
symmetric matrices, PD characterization) and *foreshadowing* Axler's more abstract
treatment for later chapters. Mechanical results (SVD existence) are cited. This
matches the preface's "prove where illuminating, cite where mechanical."

## Conventions & constraints

- Quarto `.qmd`, KaTeX math (HTML) / LaTeX (PDF); number-sections on. Keep math to
  constructs KaTeX renders.
- Citations via `references.bib` keys already present: `@strang2016`, `@axler2015`,
  `@trefethen1997`.
- Code is original, minimal, runnable against pinned `requirements.txt`; no
  proprietary data or priors.
- Match the established voice of `index.qmd` (direct, no hand-waving).

## Out of scope (YAGNI)

Determinants beyond a passing mention, full Gaussian-elimination/LU theory, general
(non-symmetric) eigentheory, complex inner-product spaces, abstract vector-space
axioms as a standalone treatment, and tensor/multilinear algebra. Any of these are
introduced later only where a specific chapter needs them.

## Success criteria

- Chapter renders cleanly (`quarto render`) in HTML and PDF with no math/build
  errors.
- Every theory subsection ties back to the design-matrix through-line.
- The three named proofs are present and correct.
- The code tie-in runs top-to-bottom and demonstrates the agree-then-diverge
  stability story.
- All four exercise tiers are populated, with matching solutions in the appendix.
