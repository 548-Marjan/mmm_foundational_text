# Chapter 3 — Probability & Statistics: Design Spec

**Date:** 2026-06-08
**Branch:** `ch3-probability-statistics`
**File to author:** `parts/01-foundations/03-probability-statistics.qmd`
**Canonical anchors:** Casella & Berger (2002); Wasserman (2004).

## Goal

Replace the scaffolded stub with a complete chapter that reactivates calculus-based
probability/statistics and ramps toward the rigor Parts II–III require. Scope is
**MMM-targeted core**: every topic earns its place by a concrete downstream use in
Bayesian regression (Part II) or sampling (Part III). Deliverable this cycle is a
**full chapter draft** — prose, proofs, worked examples, a runnable Python tie-in,
full C/B/P/A exercises with appendix solutions, and the closing **Summary** section.

The chapter follows the book's fixed template, now ending with a Summary:
Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises (C/B/P/A)
→ Summary.

## Through-line (anchor)

Chapter 1 treated `X`, `y` as fixed numbers; Chapter 2 treated the loss `L(β)` as a
deterministic surface. Chapter 3 makes `y` **random**, anchored on the probabilistic
model behind the regression fit all along:

$$y = X\beta + \varepsilon, \qquad \varepsilon \sim \mathcal{N}(0,\sigma^2 I).$$

Every concept lands on this model:

- **Random variables, expectation, variance** — `ε` and `y` as random vectors;
  `E[y] = X\beta`, `Cov(y) = \sigma^2 I`.
- **Distributions** — the Gaussian as default noise model; Bernoulli/Binomial,
  Poisson, Gamma named where MMM meets them.
- **Joint / conditional / independence** — rows of `y` independent given `X,β`;
  regression *is* conditioning (`y \mid X`).
- **Multivariate normal** — the distribution of `y`; density via `\Sigma^{-1}` and
  `\det\Sigma`; sampling needs the **Cholesky factor from Ch 1** (`y = X\beta + Lz`).
- **Likelihood & MLE** — the keystone: maximizing the Gaussian likelihood of `y`
  **recovers the least-squares problem of Ch 1–2** (OLS = Gaussian MLE); the
  Fisher information is `X^\top X/\sigma^2`, tying precision to Ch 1 conditioning.
- **LLN & CLT** — why sample averages converge (licenses Monte-Carlo in Part III)
  and why the Gaussian is ubiquitous.
- **Bayes' theorem** — flips `y \mid \beta` into `\beta \mid y`, the doorway to Part II.

Spine: *the thing you were minimizing is a negative log-likelihood, and treating
`y` as random is what lets Part II put a posterior on `β`.*

## Section-by-section content

### Motivation
Open on the gap: Ch 1–2 fit `β̂` as if `y` were exact, but weekly sales are noisy.
Driving question: *"If you reran the same campaign, you'd get different sales — so
how much of `β̂` is signal, and how much is luck?"* Answering it requires treating
`y` as a draw from a distribution, and reveals the minimized loss was a likelihood.

### Theory & Proofs
A ladder, each rung paying off on the regression model:

1. Random variables, expectation, variance; linearity of expectation; `Var`, `Cov`,
   and the **covariance matrix** of a random vector. Establish the identity
   `Cov(Ay) = A\,Cov(y)\,A^\top` (**proved**) — the same identity behind Ch 1's
   Cholesky sampling.
2. Key distributions, compact: Gaussian (the star); Bernoulli/Binomial, Poisson,
   Gamma named for where MMM uses them (conversions, counts, positive priors).
   Stated, not belabored.
3. Joint, marginal, conditional distributions; independence; conditioning as the
   meaning of regression. Establish `E[y]=X\beta`, `Cov(y)=\sigma^2 I`.
4. **The multivariate normal** — density, mean vector & covariance matrix, the role
   of `\Sigma^{-1}` and `\det\Sigma`; the affine property `y = \mu + Lz` with
   `z\sim\mathcal N(0,I)` (Cholesky from Ch 1). **Proved:** an affine transform of a
   Gaussian is Gaussian (mean/covariance computation, stated rigorously).
5. **Likelihood, log-likelihood, MLE** — the keystone. **Proved:** for
   `y\sim\mathcal N(X\beta,\sigma^2 I)`, the MLE `\hat\beta` solves the normal
   equations — **Gaussian MLE = least squares** — and the Fisher information is
   `X^\top X/\sigma^2`, tying estimate precision back to Ch 1's conditioning.
6. **LLN & CLT** — weak LLN **proved via Chebyshev**; CLT **stated and cited**
   to [@casella2002; @wasserman2004]. Framed as the license for Monte-Carlo
   estimation in Part III.
7. **Bayes' theorem** — conditional-probability form and the
   `posterior \propto likelihood \times prior` form; one line connecting `\beta\mid y`
   to Part II. Stated.

### Worked Examples
- (a) For a tiny 3-week model, compute `E[y]`, `Cov(y)`, and the Gaussian
  log-likelihood, and show maximizing it gives the **same `\hat\beta`** as Ch 1's
  normal equations.
- (b) A CLT / standard-error demo by hand: the sampling distribution of a sample
  mean, showing the `\sqrt n` shrinkage of its standard error.

### Code Tie-in
NumPy/SciPy, runnable and self-contained:
- Simulate `y = X\beta + \varepsilon` many times; show the empirical sampling
  distribution of `\hat\beta` is centered at `\beta` with covariance
  `\sigma^2 (X^\top X)^{-1}` (the inverse Fisher information).
- Draw correlated noise via the Cholesky factor (link to Ch 1).
- Overlay a Gaussian to illustrate the CLT for one coefficient.

### Exercises
- **C — Conceptual:** what `Cov(y)=\sigma^2 I` means; why OLS being the Gaussian
  MLE matters; what Bayes' theorem buys you over a point estimate.
- **B — By hand:** compute `E`/`Var`/`Cov` for small cases; a 2×2 multivariate
  normal density; a small Bayes-rule update.
- **P — Prove it:** `Cov(Ay)=A\,Cov(y)\,A^\top`; Gaussian MLE solves the normal
  equations; weak LLN via Chebyshev.
- **A — Applied / code:** Monte-Carlo the sampling distribution of `\hat\beta` and
  compare to `\sigma^2(X^\top X)^{-1}`; a CLT convergence experiment.

Exercise solutions go in the shared `appendix/solutions.qmd`, gated by the
`show-solutions` metadata flag (no inline links, per the preface).

### Summary
Per the convention added for all chapters: a **Key concepts** recap and a
**Key identities** list (inline math only, to preserve display-math parity),
covering the covariance-transform identity, `E[y]`/`Cov(y)`, the multivariate
normal density, Gaussian MLE = OLS, Fisher information `X^\top X/\sigma^2`, the
weak LLN, and Bayes' rule.

## Rigor level

Further up the book's ramp than Ch 1 — **Casella & Berger-style measured
precision**, still applied. **Prove in full** the load-bearing results
(`Cov(Ay)=A\,Cov(y)\,A^\top`; affine-Gaussian; **Gaussian MLE = least squares**;
**weak LLN via Chebyshev**). **State and cite** the mechanical/heavy results (CLT;
full distribution catalog; measure-theoretic underpinnings). Keep probability
**elementary, not measure-theoretic** (densities and sums, per a calculus-based
prerequisite) but precise about independence, conditioning, and covariance matrices.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned`
  inside `$$`, never bare `\begin{align}`; keep `$$` delimiters on their own lines
  so display-math counts stay even.
- Citations via existing `references.bib` keys: `@casella2002`, `@wasserman2004`
  (and `@strang2016` where the Cholesky/covariance link to Ch 1 is invoked).
- Code original, minimal, runnable against pinned `requirements.txt` (numpy, scipy);
  no proprietary data or priors. Verify all numeric claims with NumPy before commit.
- Match the established voice of Ch 1–2 and the cross-chapter references
  (`X^\top X`, `\hat\beta`, `\kappa(X)`, the loss `L`, `\nabla^2 L = 2X^\top X`).

## Out of scope (YAGNI)

Moment-generating functions as a primary tool, order statistics, sufficiency &
completeness, full transformation theory, measure theory, hypothesis-testing
machinery (p-values, confidence-interval theory beyond the sampling-distribution
demo), and the generalized-linear-model link functions (named only). Conjugate
priors and posterior computation are Part II; Markov chains and MCMC are Part III.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI
  (`.github/workflows/render.yml`) on the PR; Quarto is not installed locally.
- Every theory subsection ties back to the `y = X\beta + \varepsilon` through-line.
- The four named proofs are present and correct.
- The code tie-in runs top-to-bottom and demonstrates `\hat\beta`'s sampling
  distribution matching `\sigma^2(X^\top X)^{-1}`.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–2.
