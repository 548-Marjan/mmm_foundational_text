# Part II, Ch 3 — Hierarchical Models: Design Spec

**Date:** 2026-06-09
**Branch:** `ch6-hierarchical-models`
**File to author:** `parts/02-regression-bayes/03-hierarchical-regression.qmd`
**Canonical anchors:** Gelman & Hill (2006) *Data Analysis Using Regression and Multilevel/Hierarchical Models* (primary); McElreath (2020) *Statistical Rethinking*; Gelman et al. *BDA3* (2013); Hoff (2009).

## Goal

Replace the scaffolded stub with the complete final chapter of Part II. It turns Chapter 5's
single-population Bayesian regression into a **multilevel** model: per-geography channel
coefficients that **borrow strength** from one another through partial pooling. Scope is
**MMM-targeted core** — every topic earns its place by a downstream MMM use (geo-level MMM is
standard industry practice) or by setting up **Part III's MCMC**. Deliverable this cycle is a
**full chapter draft**: prose, proofs, worked examples, a runnable Python tie-in, full C/B/P/A
exercises with appendix solutions, and the closing **Summary** (auto-included convention).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

**Retitle:** the stub's H1 is `# Bayesian & Hierarchical Regression`; since Chapter 5 already
owns the Bayesian foundations, retitle to `# Hierarchical Models` and update the anchors line
to Gelman & Hill; McElreath.

## Through-line (anchor)

A geo-level MMM run across many DMAs (designated market areas). The recurring model is sales
in geography $g$, week $t$: $y_g = X_g\beta_g + \varepsilon_g$, $\varepsilon_g\sim\mathcal N(0,\sigma^2 I)$,
with **per-geo** coefficients $\beta_g$. Three stances frame the whole chapter:

1. **Complete pooling** — one national $\beta$ for all geos. Ignores real regional
   heterogeneity; biased.
2. **No pooling** — a separate $\beta_g$ fit per geo. Unbiased but wildly high-variance for
   small DMAs: a geo with 8 weeks of data returns a garbage TV coefficient (the "Des Moines
   TV = 9.2" problem).
3. **Partial pooling** — $\beta_g\sim\mathcal N(\mu,\tau^2)$, the principled compromise that
   **borrows strength** across geos. Each geo's estimate is shrunk toward the population mean
   $\mu$ by an amount set by its own precision against the population precision.

Partial pooling **is Chapter 5's "precisions add," now operating across groups** — except that
the prior mean $\mu$ and between-group variance $\tau^2$ are themselves **learned** from all the
geos. The keystone is the shrinkage formula; James–Stein is the frequentist proof that the
shrinkage strictly reduces total error; and the chapter ends by deriving the **Gibbs full
conditionals** of the hierarchical model — the explicit, concrete handoff to Part III.

## Section-by-section content

### Motivation
Open on the geo-MMM dilemma: run MMM across 50 DMAs and you face complete pooling (one national
model that hides that TV works in some regions and not others) versus no pooling (50 separate
models whose small-DMA coefficients are noise). Driving question (verbatim emphasis):
*"Des Moines came back with a TV coefficient of 9.2 — absurd on 10 weeks of data. Do you trust
it, throw it out, or shrink it toward the national average — and by exactly how much?"* Partial
pooling answers "shrink, and here is precisely how much," and it falls straight out of the
Chapter 5 machinery. Name the payoffs: principled shrinkage (the keystone), the James–Stein
guarantee that it beats per-geo fits, a learned pooling strength $\tau^2$, and — because the
joint posterior is non-conjugate — the Gibbs sampler that motivates Part III.

### Theory & Proofs
A ladder, each rung an MMM payoff:

1. **Three stances: complete, no, and partial pooling.** Formalize the two-level model. For
   clarity start with one scalar parameter per geo: group estimate $\bar y_g$ with sampling
   variance $\sigma^2/n_g$, and the population model $\beta_g\sim\mathcal N(\mu,\tau^2)$. State
   complete pooling as $\tau^2\to0$ (all $\beta_g=\mu$) and no pooling as $\tau^2\to\infty$
   (each $\beta_g$ free). Set up the bias (complete) vs variance (no pooling) tension.
2. **Partial pooling = precision-weighted shrinkage** (**KEYSTONE proof A**). For the Gaussian
   hierarchical model with known $\sigma^2,\mu,\tau^2$, prove (this is Chapter 5 Rung 3 applied
   per group, completing the square) that the posterior mean is
   $$\hat\theta_g = \lambda_g\,\bar y_g + (1-\lambda_g)\,\mu,\qquad
     \lambda_g = \frac{n_g/\sigma^2}{\,n_g/\sigma^2 + 1/\tau^2\,}.$$
   Read it off: $\lambda_g$ is the **shrinkage weight**; large $n_g$ → $\lambda_g\to1$ (trust the
   geo's own data, no pooling); small $n_g$ → $\lambda_g\to0$ (shrink to the grand mean, complete
   pooling). Partial pooling interpolates the two stances, per geo, by precision.
3. **James–Stein: shrinkage dominates** (**KEYSTONE proof B**). State the James–Stein estimator
   for $G\ge3$ group means and prove it has strictly smaller **total** mean-squared risk than the
   per-group MLE $\{\bar y_g\}$ under the canonical equal-variance Gaussian setup. Use Stein's
   unbiased risk estimate (SURE) / Stein's lemma to compute the risk and exhibit the strict
   improvement. Connect the empirical-Bayes shrinkage factor to $\lambda_g$ of Rung 2: shrinkage
   is not a Bayesian indulgence — it provably lowers total error, which is why pooling wins.
4. **Where $\mu$ and $\tau^2$ come from.** Partial pooling needs the hyperparameters. Two routes:
   **empirical Bayes** (plug in point estimates of $\mu,\tau^2$ from the marginal likelihood) and
   **full Bayes** (hyperpriors on $\mu,\tau^2$). The **variance-components** reading: $\tau^2$ is
   the pooling dial, estimated from between-group spread, so the data decide how much to pool.
   State the partial-pooling caveat that empirical Bayes understates uncertainty by ignoring
   hyperparameter error — full Bayes (and MCMC) fixes this.
5. **The vectorized geo-MMM model.** Lift to the real design: $\beta_g\in\mathbb R^p$ per geo,
   $\beta_g\sim\mathcal N(\mu_\beta,\Sigma_\beta)$, with the Chapter 5 Gaussian likelihood per geo.
   Pooling now happens per coefficient (each channel partially pooled across geos). One paragraph:
   the identical machinery applies to **any** grouping — product segments, time blocks — not just
   geography; geography is the running example, not a limitation.
6. **The joint posterior is non-conjugate — but the full conditionals are** (**DERIVATION → Part
   III bridge**). Write the joint posterior $p(\beta_{1:G},\mu_\beta,\Sigma_\beta,\sigma^2\mid y)$
   and note it has **no closed form** jointly (Chapter 5's conjugacy wall, now structural in the
   hierarchy). Then derive each **full conditional** and show each is conjugate:
   - $\beta_g\mid\text{rest}$: Gaussian — exactly Chapter 5 Rung 4 (Bayesian linear regression with
     prior $\mathcal N(\mu_\beta,\Sigma_\beta)$). Give $\Sigma_g=(X_g^\top X_g/\sigma^2+\Sigma_\beta^{-1})^{-1}$,
     mean $\Sigma_g(X_g^\top y_g/\sigma^2+\Sigma_\beta^{-1}\mu_\beta)$.
   - $\mu_\beta\mid\text{rest}$: Gaussian — Chapter 5 Rung 3 with the $\beta_g$'s as "data."
   - $\Sigma_\beta\mid\text{rest}$: Inverse-Wishart (or Inverse-Gamma in the scalar/diagonal case);
     **state** the conjugate form.
   - $\sigma^2\mid\text{rest}$: Inverse-Gamma; **state** the conjugate form.
   Conclude: cycling through these conditionals — sampling each in turn — **is Gibbs sampling**, the
   simplest MCMC and the subject of Part III. "Everything is in hand; Part III shows how iterating
   these conjugate conditionals draws from the joint posterior we cannot write down." Cite
   `[@gelman2013]`, `[@gelmanhill2006]`. End the Theory section on this baton-pass.

### Worked Examples
- (a) **Shrinkage by hand**: three geos with different sample sizes (e.g. $n=(4, 16, 64)$),
  shared $\sigma^2$, given $\mu$ and $\tau^2$. Compute each $\lambda_g$ and the partial-pooled
  estimate $\hat\theta_g=\lambda_g\bar y_g+(1-\lambda_g)\mu$, showing the small-$n$ geo shrinks
  hardest toward $\mu$ and the large-$n$ geo barely moves. All numbers NumPy-verified.
- (b) **James–Stein beats the MLE**: a small fixed example with known truth (e.g. $G=6$ group
  means), compute total squared error of the raw MLE vs the James–Stein estimate, and show JS wins.
  NumPy-verified.

### Code Tie-in
NumPy/SciPy, runnable and self-contained, headless (`MPLBACKEND=Agg python3`):
- Simulate a geo-MMM: $G$ geos with **unequal** week counts $n_g$, true $\beta_g\sim\mathcal N(\mu,\tau^2)$.
- Fit three ways: **complete pooling** (one OLS on stacked data), **no pooling** (per-geo OLS),
  and **partial pooling** (closed-form empirical-Bayes shrinkage with $\tau^2$ estimated from the
  between-group spread).
- The classic **shrinkage plot**: per-geo no-pooling estimates and partial-pooled estimates, with
  connecting lines showing small/noisy geos pulled hardest toward the grand mean.
- Report **total RMSE vs truth** for the three methods; partial pooling wins.
- A compact **hand-rolled Gibbs sampler** (a dozen lines) for the scalar hierarchical model that
  cycles the Rung 6 conditionals and recovers the closed-form/EB shrinkage answer — a concrete
  preview of Part III. Verify the Gibbs posterior mean matches the closed form.

### Exercises
- **C — Conceptual:** the three pooling stances and when each is right; what $\lambda_g$ means and
  what drives it; why $\tau^2$ is the pooling dial; why each full conditional being conjugate is
  exactly what makes Gibbs (Part III) possible.
- **B — By hand:** compute a shrinkage factor $\lambda_g$ and partial-pooled estimate; compute a
  James–Stein estimate and its total-error improvement; compute the parameters of a $\beta_g$ or
  $\mu$ full conditional numerically from given quantities.
- **P — Prove it:** the partial-pooling posterior-mean shrinkage formula (complete the square);
  James–Stein dominance (SURE / Stein's lemma risk inequality for $G\ge3$); the $\beta_g$ full
  conditional is Gaussian with the stated mean and covariance.
- **A — Applied / code:** the three-way pooling simulation with the shrinkage plot and RMSE
  comparison; a sweep of $\tau^2$ (or of the small geo's $n_g$) showing shrinkage strength change;
  implement the Gibbs sampler and show it matches the closed-form posterior.

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions` (no inline links),
heading `## Chapter 6 — Hierarchical Models` to match siblings (em-dash; `### C — …`, `### B — …`,
`### P — …`, `### A — …`).

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only) covering the three pooling
stances, the shrinkage factor $\lambda_g$ and the partial-pooling posterior mean, James–Stein
dominance, $\tau^2$ as the learned pooling dial, the vectorized hierarchical regression model, and
the conjugate full conditionals that make Gibbs/Part III possible.

## Rigor level

Closing chapter of Part II — **Gelman–Hill / BDA3 style**: clean analytical results, every result
cashed out on the geo-MMM model, no measure theory. **Prove in full**: (1) the partial-pooling
posterior-mean shrinkage formula by completing the square (leaning on Chapter 5 Rung 3); (2)
**James–Stein dominance** for $G\ge3$ via Stein's unbiased risk estimate; (3) the $\beta_g$ full
conditional is Gaussian (Chapter 5 Rung 4 reused). **Derive** the joint posterior's other full
conditionals and **state** their conjugate families (Gaussian for $\mu_\beta$, Inverse-Gamma /
Inverse-Wishart for the variances). Anchor to the design matrix and the media-mix story; keep the
James–Stein treatment at the canonical equal-variance Gaussian setup (cite the general case).

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned` inside `$$`, never
  bare `\begin{align}`; `$$` delimiters on their own lines so display-math counts stay even.
- Citations via existing `references.bib` keys: `@gelmanhill2006`, `@mcelreath2020`, `@gelman2013`,
  `@hoff2009` (and `@casella2002` for distribution facts if needed). **Do not invent keys.**
- Code original, minimal, runnable against pinned `requirements.txt` (numpy, scipy, matplotlib);
  verify every numeric claim with NumPy before commit; the cell must run headless under
  `MPLBACKEND=Agg python3` and end figures with `plt.show()`.
- Match the established voice of Ch 1–5 and reuse their symbols: $X^\top X$, $\hat\beta$,
  $\mathcal N(\mu,\Sigma)$, the Chapter 5 posterior $\Sigma_{\text{post}},\mu_{\text{post}}$, the
  "precisions add" motif, ridge $(X^\top X+\lambda I)^{-1}X^\top y$. Reuse Chapter 5 results
  explicitly (Rung 3 for $\mu$, Rung 4 for $\beta_g$) rather than re-deriving from scratch.
- Git identity in the fresh worktree: `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"`.

## Out of scope (YAGNI)

General GLMM / non-Gaussian likelihoods (logistic, Poisson multilevel); crossed and non-nested
random effects (only nested geo grouping here); REML estimation theory and its asymptotics;
deep Inverse-Wishart theory (state the conjugate form, do not develop it); MCMC convergence /
mixing theory and any actual sampler beyond the illustrative scalar Gibbs preview (that is Part
III); empirical-Bayes asymptotic optimality; non-linear adstock/saturation hierarchy (Part V).
The chapter keeps $\Sigma_\beta$ scalar/diagonal ($\tau^2$ per coefficient) for the closed forms,
naming the full-covariance extension once.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR (Quarto not
  installed locally).
- Every theory subsection ties back to the three-stances / shrinkage / Part III through-line.
- The three named proofs (partial-pooling shrinkage formula, James–Stein dominance, $\beta_g$
  Gaussian full conditional) are present and correct, plus the derivation of the remaining
  conjugate full conditionals.
- The code tie-in runs top-to-bottom headless: it shows the shrinkage plot, partial pooling's
  lowest RMSE, and a Gibbs sampler whose posterior mean matches the closed form.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–5; H1 retitled to `# Hierarchical Models`.
