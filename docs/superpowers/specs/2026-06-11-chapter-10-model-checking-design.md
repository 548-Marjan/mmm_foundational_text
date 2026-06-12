# Chapter 10 — Model Checking & Calibration: Design Spec

## Placement & Role

- **File:** `parts/03-sampling/04-model-checking.qmd`
- **Part:** III — Computation: Sampling (this is the **final** chapter of Part III)
- **Predecessors it builds on:** Ch6 (hierarchical models / partial pooling), Ch8 (MCMC, ESS via autocorrelation), Ch9 (HMC/NUTS, divergences, Neal's funnel, non-centered reparameterization)
- **Hand-off it pays off:** Ch9 closed by promising that "the formal machinery for detecting these pathologies — $\hat{R}$, effective sample size, and divergence counts — is the subject of Chapter 10."
- **Hands off to:** Part IV — Optimization (a *fitted, checked* posterior is the input to budget optimization).

## Driving Questions (the spine)

The chapter is organized around two questions a practitioner must answer before trusting any posterior:

1. **"Did the sampler converge?"** — computational diagnostics: did the Markov chain actually explore the posterior, and is the Monte Carlo error small enough?
2. **"Is the model any good?"** — calibration & checking, in two sub-answers: **(B1) self-consistency** — even with perfect sampling, does the model reproduce the observed data and are its credible intervals honest (PPC, SBC)? and **(B2) out-of-sample prediction** — does it predict held-out data well, and better than an alternative (ELPD, WAIC, PSIS-LOO)? Framing the two halves explicitly keeps the checking/comparison seam a feature rather than a blur.

The keystone result is **Simulation-Based Calibration (SBC)** rank-uniformity (Talts et al. 2018), which formally connects "correct posterior computation" to a checkable, distribution-free signal.

## Scope — Comprehensive (nothing important omitted)

### Part A — "Did the sampler converge?" (computational diagnostics)

1. **Monte Carlo standard error (MCSE).** Why a posterior mean estimate has its own uncertainty; $\text{MCSE} = \text{sd}/\sqrt{N_\text{eff}}$; relation to the CLT for ergodic averages (Ch8).
2. **$\hat{R}$ (Gelman–Rubin).** Between-chain variance $B$ and within-chain variance $W$; the pooled variance estimator $\widehat{\text{var}}^+ = \frac{N-1}{N}W + \frac{1}{N}B$; $\hat{R} = \sqrt{\widehat{\text{var}}^+/W}$. **Full proof** that $\hat{R}\to 1$ as chains converge (B/N → between-chain variance of means → 0 contribution at stationarity; both estimators consistent for the marginal variance). **Split-$\hat{R}$** (halving chains to catch trends) and **rank-normalized $\hat{R}$** (Vehtari et al. 2021) to handle heavy tails / non-finite variance.
3. **Effective sample size (ESS).** $N_\text{eff} = N/(1 + 2\sum_{k\ge1}\rho_k)$ from the autocorrelation sum. **Full derivation** of the variance-inflation factor (extends Ch8's autocorrelation result). **Bulk-ESS** (rank-normalized, for central tendency) vs **tail-ESS** (for interval estimates) — why a single ESS is insufficient.
4. **HMC-specific signals (from Ch9).** Divergence count as a non-ignorable bias signal; BFMI (energy Bayesian fraction of missing information); tree-depth saturation in NUTS. These are *targeted* diagnostics that $\hat{R}$/ESS can miss.

### Part B — "Is the model any good?" (calibration & checking)

5. **Prior predictive checks.** Draw $\theta\sim p(\theta)$, $\tilde{y}\sim p(y\mid\theta)$; do simulated datasets look plausible? Catches absurd priors before fitting.
6. **Posterior predictive checks (PPC).** Replicated data $y^\text{rep}\sim p(y^\text{rep}\mid y)=\int p(y^\text{rep}\mid\theta)p(\theta\mid y)\,d\theta$; test quantities $T(y,\theta)$; **posterior predictive p-value** $p_B = \Pr(T(y^\text{rep},\theta)\ge T(y,\theta)\mid y)$. **Proof** that for a well-calibrated test statistic the p-value is $\text{Uniform}(0,1)$ (probability-integral-transform argument), and the standard caveat that PPC p-values are conservative because the data are used twice.
7. **Simulation-Based Calibration (SBC) — KEYSTONE.** The construction: draw $\theta_0\sim p(\theta)$, $y\sim p(y\mid\theta_0)$, then $\{\theta_1,\dots,\theta_L\}\sim p(\theta\mid y)$; compute the rank $r=\#\{\theta_\ell < \theta_0\}$. **Full proof** of the rank-uniformity theorem: because the *data-averaged posterior equals the prior*, $\int p(\theta\mid y)p(y)\,dy = p(\theta)$, the rank of the prior draw among posterior draws is discrete-uniform on $\{0,\dots,L\}$. Interpreting rank-histogram shapes (∪ = under-dispersed/overconfident posterior; ∩ = over-dispersed; left/right slope = biased). Tie-back to Ch9: a centered funnel fails SBC; the non-centered form passes.
8. **Out-of-sample predictive accuracy (model comparison).** Expected log pointwise predictive density (**ELPD**); **WAIC** $= \text{lppd} - p_\text{WAIC}$ with the effective-parameter penalty; **PSIS-LOO** (Pareto-smoothed importance sampling leave-one-out) with the **Pareto-$\hat{k}$** reliability diagnostic (Vehtari, Gelman & Gabry 2017). State the WAIC↔LOO asymptotic equivalence; emphasize these answer "which model predicts better," distinct from "is this model calibrated."

### Synthesis

9. **A model-checking workflow / checklist for an MMM.** Order of operations (prior predictive → fit → convergence diagnostics → posterior predictive → SBC on the pipeline → LOO for comparison); what each catches; what a green/red light means for downstream budget optimization (Part IV). MMM relevance is carried **in prose** throughout (e.g., $\hat{R}$ on adstock/saturation/ROAS parameters, PPC on weekly sales, SBC on the adstock→saturation→regression pipeline).

## Proof Inventory (the rigor the chapter is graded on)

Three results get the full theorem→proof (`\blacksquare`) treatment:

- **P1 — $\hat{R}\to 1$** via the $B$/$W$ variance decomposition.
- **P2 — ESS variance-inflation factor** $1 + 2\sum_k\rho_k$ from the autocorrelated-mean variance.
- **P3 — SBC rank-uniformity** (keystone), from the data-averaged-posterior identity.

A fourth, lighter result is the **posterior-predictive p-value uniformity** (PIT argument), stated and proved compactly. WAIC/LOO equivalence and PSIS theory are **stated precisely** with references, not proved (they are asymptotic/empirical-process results beyond the book's scope).

## Running Example — Eight Schools, Adapted to a Geo-Level MMM

The chapter teaches the **canonical eight-schools hierarchical model** by name — the model on which essentially the entire MCMC-diagnostics literature (the funnel, $\hat{R}$, SBC) was developed — and then **immediately adapts it to the geo-level MMM** so the reader sees the marketing structure underneath. The same model carries every worked example and the code tie-in.

- **Canonical model (eight schools, @gelman2013):** $y_j\mid\theta_j\sim\mathcal{N}(\theta_j,\sigma_j^2)$, $\theta_j\mid\mu,\tau\sim\mathcal{N}(\mu,\tau^2)$, $j=1,\dots,J$, with weak priors on $\mu,\tau$. Introduced first because it is *the* canonical funnel-prone hierarchy: it exhibits divergences, $\hat{R}$/ESS degradation, and SBC/PPC behavior all at once, and it connects directly to Ch6 (partial pooling) and Ch9 (funnel + non-centering).
- **Adapted to MMM (the quick translation):** the same likelihood with marketing nouns — $\theta_g$ is a **geo-level effect** (e.g. incremental ROAS or baseline sales for region $g$), partially pooled toward a national mean $\mu$ with cross-geo SD $\tau$, observed as a noisy estimate $y_g$ with standard error $\sigma_g$ from that geo's regression. The funnel pathology becomes a concrete MMM trap: **few geos $\Rightarrow$ $\tau$ weakly identified $\Rightarrow$ funnel $\Rightarrow$ divergences** — exactly what bites a geo-MMM fit on a handful of regions. This is presented as a one-paragraph "same math, marketing labels" bridge right after the canonical model, then used interchangeably.
- **Worked Examples (exact arithmetic a reader can verify):**
  - (a) Compute $\hat{R}$ by hand from a small set of toy chains (e.g., 4 chains × a few draws) with $B$, $W$, $\widehat{\text{var}}^+$, $\hat{R}$ shown to 3–4 decimals.
  - (b) Compute an ESS / variance-inflation factor by hand from a short autocorrelation sequence.
  - (c) Compute one SBC rank and read a small rank histogram; show why a ∪-shape means overconfidence.
- **Code Tie-in (single runnable `{python}` cell, NumPy/Matplotlib only, headless):**
  - Fit the centered and non-centered hierarchical-normal with the **from-scratch HMC/leapfrog sampler built in Ch9** (so divergences are genuine, not a proxy) — **no PyMC/Stan dependency**, consistent with the book's from-scratch ethos.
  - Compute split-$\hat{R}$, bulk/tail-ESS, MCSE, and the genuine divergence count (centered vs non-centered) from the HMC draws.
  - Posterior predictive check on the $y_j$ with a test quantity; report a Bayesian p-value.
  - Run a small SBC loop and plot the rank histogram (centered fails, non-centered ≈ uniform).
  - All numeric assertions NumPy-verified before authoring; figures end with `plt.show()`.

## Template Conformance

Six required H2 sections in order: **Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary** (Summary auto-included; inline math only). Exercises in four tiers **C / B / P / A**, self-contained (no inline solution links). Appendix solutions go in `appendix/solutions.qmd` under a new `## Chapter 10 — Model Checking & Calibration` block (em-dash), gated by `::: {.content-visible when-meta="show-solutions"}`, inserted before the file's final closing `:::`.

## KaTeX / Build Conventions

- `aligned` inside `$$…$$`; never bare `\begin{align}`. `$$` delimiters on their own lines. Even `$` count.
- Single `{python}` cell runs top-to-bottom under `MPLBACKEND=Agg python3`.
- Real gate: CI `quarto render` (HTML + PDF) green on the PR. Quarto is not installed locally.

## Citations

- Existing: `@talts2018` (SBC), `@gelman2013` (BDA3 — PPC, p-values, WAIC, LOO), `@gelmanhill2006`, `@mcelreath2020`, `@betancourt2017`, `@hoffman2014`.
- **To add to `references.bib`:** `@vehtari2021` (rank-normalization, folding, localization — improved $\hat{R}$ and bulk/tail ESS) and `@vehtari2017` (PSIS-LOO / WAIC practical Bayesian model evaluation).

## Out of Scope (deliberately)

- Variational inference diagnostics (the book is MCMC-centric through Part III).
- Cross-validation schemes beyond LOO (k-fold, time-series CV) — mentioned in one sentence as MMM-relevant future reading, not developed.
- Stan/PyMC tooling — diagnostics are computed from raw draws to keep the chapter self-contained.
