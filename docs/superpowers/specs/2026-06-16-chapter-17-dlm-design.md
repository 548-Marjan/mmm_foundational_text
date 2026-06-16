# Chapter 17 — Dynamic Linear Models & State-Space MMM: Design Spec

**Date:** 2026-06-16
**Status:** Approved (brainstorm complete). Next: implementation plan via `writing-plans`.
**Part:** V — Marketing Mix Modeling: Modeling (third chapter).
**File:** `parts/05-mmm-modeling/03-dlm.qmd` (currently a template stub).
**Canonical anchors:** `@west1997` (West & Harrison, *Bayesian Forecasting and Dynamic Models*),
`@durbin2012` (Durbin & Koopman, *Time Series Analysis by State Space Methods*) — both to be added
to `references.bib`. Library-agnostic: name no MMM/PPL library; `numpy` is fine.

This chapter implements two locked decisions:
- **Reorg plan (`docs/notes/2026-06-12-part-v-vi-reorg-plan.md`, Ch. 16 section):** present **classical
  DLM first, then Bayesian DLM**; surface the three payoffs (rhymes Ch4→Ch5; the conjugacy-vs-MCMC
  boundary; the clean downstream split — classical Kalman → BSTS in Ch. 20, Bayesian DLM → prior store
  in Ch. 22).
- **DLM-vs-BTVC decision (`docs/notes/2026-06-14-dlm-vs-btvc-structure-decision.md`):** **DLM is the
  required core baseline; BTVC is one flagged, skippable Advanced/Extensions endpoint.** The escalation
  is static (Ch. 15–16) → dynamic via DLM → dynamic via BTVC. Push BTVC-specific mechanics (SVI,
  folded-normal hierarchy) to an appendix; keep the core arc independent of BTVC.

---

## 1. Role in the book

Chapter 16 fit an MMM with **constant** coefficients. This chapter lets them move. Real channel
effectiveness drifts — media efficiency decays, creative wears out, the market shifts — so a constant
$\beta$ is a falsifiable null, not a law. The dynamic linear model (DLM) is the state-space machinery
that makes coefficients time-varying while keeping inference tractable.

It is the time-series payoff the book has pointed at since Chapter 7, which foreshadowed that AR(1)
persistence and adstock carryover — both first-order recursions — would be "unified by state-space
models in Part V." It rhymes Chapter 4 → Chapter 5 (frequentist → Bayesian) at the time-series level,
and it is the book's **cleanest conjugacy-vs-MCMC boundary**: with Gaussian errors and *known*
variances the Kalman recursions are the exact closed-form posterior (no sampling); make the variances
unknown and you fall back to MCMC — retroactively sharpening Part III.

The chapter follows the fixed template: Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary, with gated appendix solutions (and a BTVC-mechanics appendix block).

## 2. Scope decisions (locked)

- **Classical first, then Bayesian** (reorg plan).
- **Proof slate:** Kalman-filter derivation (keystone) + FFBS correctness + the conjugacy-vs-MCMC
  boundary; the RTS smoother is stated and derived briefly, with the full backward recursion set as a
  P-tier exercise. (Approved in brainstorm.)
- **BTVC is a flagged Advanced/Extensions section**, conceptual in the main text (kernel-smoothed knots
  $\beta = Kb$ with $J \ll T$; folded-normal positivity); SVI/Orbit mechanics go to a dedicated
  appendix block, clearly marked skippable. The core arc does not depend on it.
- **Library-agnostic:** the Kalman filter, RTS smoother, FFBS, and the Gibbs sampler are implemented
  from scratch in `numpy`.

## 3. Theory rungs

1. **Static to dynamic.** Why $\beta_t$ drifts; the constant-coefficient model as the null this chapter
   relaxes (parallels the i.i.d.→AR(1) relaxation of Chapter 7, now at the coefficient level).
2. **State-space form.** Observation $y_t = F_t^\top \theta_t + v_t$, $v_t \sim \mathcal N(0, V_t)$;
   evolution (state) $\theta_t = G_t \theta_{t-1} + w_t$, $w_t \sim \mathcal N(0, W_t)$. Two anchors:
   the **local level** ($F=G=1$, a random walk plus noise) and the **regression DLM** ($F_t = x_t$ the
   channel inputs, $G=I$, so $\beta_t$ random-walks — the time-varying-coefficient MMM).
3. **The Kalman filter.** Forward predict→update recursion giving the exact Gaussian filtering posterior
   $p(\theta_t \mid y_{1:t})$.
4. **The RTS smoother.** Backward recursion for $p(\theta_t \mid y_{1:T})$ (full-sample); stated and
   sketched, full derivation as exercise.
5. **Bayesian DLM.** Priors on $(\theta_0, V, W)$; the full posterior over states *and* variances via
   **FFBS** (the Kalman filter run forward, then sampled backward) inside a Gibbs sampler with
   inverse-gamma variance updates.
6. **The conjugacy-vs-MCMC boundary.** Known $(V,W)$ ⇒ exact closed-form posterior + marginal
   likelihood by prediction-error decomposition; unknown $(V,W)$ ⇒ non-conjugate, requires MCMC.
7. **BTVC (Advanced/Extensions, flagged).** The production reframing that resolves the DLM dilemma:
   coefficients as $J \ll T$ kernel-smoothed latent knots $\beta = Kb$ (parameter count decoupled from
   series length), folded-normal two-layer hierarchy for positivity + shrinkage. Conceptual only;
   mechanics in the appendix.

## 4. Proof slate

| # | Result | Statement / method |
|---|---|---|
| **P1 (keystone)** | Kalman filter = exact Gaussian conditioning | By induction on $t$. Assume $\theta_{t-1}\mid y_{1:t-1}\sim\mathcal N(m_{t-1},C_{t-1})$. **Predict:** $\theta_t\mid y_{1:t-1}\sim\mathcal N(a_t,R_t)$ with $a_t=G_t m_{t-1}$, $R_t=G_t C_{t-1}G_t^\top+W_t$. The joint $(\theta_t,y_t)\mid y_{1:t-1}$ is Gaussian with $y_t$-mean $f_t=F_t^\top a_t$, variance $Q_t=F_t^\top R_t F_t+V_t$, cross-covariance $R_t F_t$. **Update** by the Gaussian-conditioning formula (Ch. 3/5): $\theta_t\mid y_{1:t}\sim\mathcal N(m_t,C_t)$, $m_t=a_t+A_t(y_t-f_t)$, $C_t=R_t-A_tQ_tA_t^\top$, Kalman gain $A_t=R_tF_tQ_t^{-1}$. |
| **P2** | FFBS correctness | The joint state posterior factors by the Markov property: $p(\theta_{0:T}\mid y_{1:T})=p(\theta_T\mid y_{1:T})\prod_{t=0}^{T-1}p(\theta_t\mid\theta_{t+1},y_{1:t})$. The forward filter supplies $p(\theta_t\mid y_{1:t})$; each backward conditional $p(\theta_t\mid\theta_{t+1},y_{1:t})$ is Gaussian (condition the filtered $\theta_t$ on the evolution equation linking it to $\theta_{t+1}$). Sampling $\theta_T$ then each $\theta_t$ backward draws an exact joint posterior path — "filter forward, sample backward." |
| **P3 (boundary)** | Known variances ⇒ closed form; unknown ⇒ MCMC | With $(V,W)$ known, P1 already gives the exact posterior with no sampling, and the marginal likelihood factors as $p(y_{1:T})=\prod_t \mathcal N(y_t;f_t,Q_t)$ (prediction-error / innovations decomposition) — usable for MLE of $(V,W)$. Treating $(V,W)$ as unknown random variables makes the posterior $p(\theta_{0:T},V,W\mid y)$ non-conjugate jointly; a Gibbs sampler restores tractability — FFBS draws $\theta_{0:T}\mid V,W,y$ (Gaussian, P2) and conjugate inverse-gamma draws give $V,W\mid\theta_{0:T},y$. This is the exact frontier where conjugacy ends and MCMC begins. |

## 5. Worked examples + anchor numbers (NumPy-verified in the plan)

- **WE1 — local-level steady-state Kalman gain.** Local level with $V=W=1$. The filtering-variance
  Riccati recursion $C = (C+W) - (C+W)^2/(C+W+V)$ has fixed point $R = C+W$ solving $R^2 - R - 1 = 0$,
  i.e. $R = \varphi = \tfrac{1+\sqrt5}{2}$, giving steady-state Kalman gain
  $K^\star = R/(R+1) = \tfrac{\sqrt5-1}{2} \approx 0.6180$ (the golden ratio) — and steady-state
  filtered variance $C^\star = 0.6180$. Show one predict–update step by hand and that iterating
  converges to $K^\star$. The steady-state filter is an EWMA with smoothing constant $K^\star$.
- **WE2 — regression-DLM recovery.** A single time-varying coefficient $\beta_t$ following a random
  walk; the Kalman filter tracks the drift, and the RTS smoother is strictly tighter (smoother RMSE <
  filter RMSE, since it uses future data). Recover the **time-varying channel ROI**.
- **WE3 — the boundary in numbers.** With known $(V,W)$ the prediction-error decomposition gives the
  exact log marginal likelihood; with unknown $(V,W)$ the Gibbs/FFBS sampler recovers $(V,W)$ near
  their true values — the same model, the two inference regimes.

## 6. Code tie-in

A single runnable `{python}` cell implementing from scratch (numpy only): `kalman_filter`,
`rts_smoother`, `ffbs`, and `gibbs_dlm`. It (a) runs the local level at $V=W=1$ and asserts the
steady-state gain converges to $(\sqrt5-1)/2 \approx 0.618$; (b) simulates a regression DLM with a
drifting $\beta_t$, filters and smooths, and asserts smoother RMSE < filter RMSE; (c) runs the Gibbs/FFBS
sampler with unknown variances and asserts it recovers $(V,W)$ within tolerance; (d) plots the true vs
filtered vs smoothed $\beta_t$ (the time-varying coefficient / ROI). Runs headless under
`MPLBACKEND=Agg`; all pinned numbers asserted.

## 7. BTVC appendix block

A dedicated, clearly-marked **Advanced/Extensions** appendix subsection (skippable) giving the BTVC
mechanics kept out of the main text: the kernel-smoothing map $\beta = Kb$ with $J \ll T$ knots, the
folded-normal two-layer hierarchy for positivity and shrinkage, and a one-paragraph account of why
stochastic variational inference (SVI) is used at scale. Conceptual and self-contained; no library
named. This satisfies the reorg's deferred "mark the Advanced/Extensions track" item.

## 8. Exercises, Summary, housekeeping

- **Exercises** C/B/P/A. Include the **RTS smoother backward recursion as a P-tier proof**; a B-tier
  by-hand Kalman step; an A-tier extending the sampler or comparing filter vs smoother. Appendix
  solutions appended to `appendix/solutions.qmd` as `## Chapter 17 — Dynamic Linear Models &
  State-Space MMM`, gated by `::: {.content-visible when-meta="show-solutions"}`, before the final `:::`.
- **Summary** auto-included; **Key concepts** and **Key identities** both **bulleted** (inline math
  only).
- **Bib:** add `@west1997`, `@durbin2012`. **Anchor line:** replace the stub's
  `*Canonical anchors: West & Harrison 1997; Durbin & Koopman.*` with
  `*Canonical anchors: @west1997; @durbin2012.*`.

## 9. Conventions (enforced per task)

- KaTeX: `aligned` inside `$$ … $$` (never bare `\begin{align}`); `$$` on their own lines; even
  `$`-count per line. **Never** `\begin{psmallmatrix}` (undefined in CI LuaLaTeX/PDF). Matrices use
  `bmatrix`/`pmatrix`.
- Build-based verification: the code cell runs headless; the real gate is CI `quarto render`
  (HTML + PDF) — watch **both** jobs to green before reporting done.
- Cross-references use the current (post-renumber) chapter numbers: Gaussian conditioning Ch. 3/5,
  Time Series Ch. 7, Gibbs Ch. 6/9, DGP Ch. 15, Building & Fitting Ch. 16; forward to BSTS Ch. 20 and
  the prior store Ch. 22.

## 10. Critical files

- `parts/05-mmm-modeling/03-dlm.qmd` — the chapter (replace stub; keep H1; fix anchor line).
- `appendix/solutions.qmd` — append the Chapter 17 solutions block before the final `:::`.
- `references.bib` — add `@west1997`, `@durbin2012`.
- Read-only exemplars: `parts/05-mmm-modeling/02-building-fitting.qmd` (predecessor; the static fit this
  chapter makes dynamic), `parts/02-regression-bayes/04-time-series.qmd` (AR(1)/state-space foreshadow),
  `parts/02-regression-bayes/02-bayesian-inference.qmd` (Gaussian conditioning / MAP voice).
