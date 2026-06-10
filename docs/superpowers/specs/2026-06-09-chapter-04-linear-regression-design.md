# Part II, Ch 1 — Linear Regression: Design Spec

**Date:** 2026-06-09
**Branch:** `ch4-linear-regression`
**File to author:** `parts/02-regression-bayes/01-linear-regression.qmd`
**Canonical anchors:** Hastie, Tibshirani & Friedman, *ESL* (2009); Casella & Berger (2002) for distributional mechanics.

## Goal

Replace the scaffolded stub with a complete chapter that turns the OLS estimator of
Part I into applied statistical methodology — inference, diagnostics, and
regularization — and bridges to the Bayesian Inference chapter that follows. Scope
is **MMM-targeted core**: every topic earns its place by a downstream use in MMM
practice or by setting up Bayesian/hierarchical regression. Deliverable this cycle is
a **full chapter draft** — prose, proofs, worked examples, a runnable Python tie-in,
full C/B/P/A exercises with appendix solutions, and the closing **Summary** (now an
auto-included convention for every chapter).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

## Through-line (anchor)

Continues the Ch 3 model `y = X\beta + \varepsilon`, `\varepsilon\sim\mathcal N(0,\sigma^2 I)`.
Ch 3 asked "what is `\hat\beta`?" (OLS = Gaussian MLE, covariance `\sigma^2(X^\top X)^{-1}`).
This chapter asks the analyst's three follow-ups:

1. **Is `\hat\beta` the best we can do?** → Gauss–Markov: OLS is BLUE.
2. **How sure are we?** → inference: standard errors, t-tests, confidence intervals.
3. **What happens when media channels collide?** → multicollinearity inflates
   `(X^\top X)^{-1}` (large `\kappa(X)` from Ch 1), so standard errors explode.

The third problem is the hinge. The **bias–variance tradeoff** shows accepting a
little bias can slash variance, and **ridge regression**
`\hat\beta_{\text{ridge}} = (X^\top X + \lambda I)^{-1}X^\top y` does exactly that,
re-conditioning the `X^\top X` collinearity made near-singular. Closing reveal:
ridge is the **MAP estimate under a Gaussian prior** `\beta\sim\mathcal N(0,\tau^2 I)`
with `\lambda=\sigma^2/\tau^2` — penalization *is* a prior, the doorway to the
Bayesian chapter.

Spine: OLS is optimal and quantifiable → collinearity (the MMM disease) breaks it →
regularization fixes it → regularization is secretly Bayesian.

## Section-by-section content

### Motivation
Ch 3 gave `\hat\beta` and its spread `\sigma^2(X^\top X)^{-1}`, but an analyst must
*act*: report a number with error bars, decide whether a channel "worked," and cope
when TV and search are nearly collinear. Driving question: *"You fit the model and
TV's coefficient is 2.3. Is that real, how sure are you, and why did the number
swing wildly when you added search?"* The chapter turns the estimator into inference,
then fixes what collinearity breaks.

### Theory & Proofs
A ladder, each rung an MMM payoff:

1. OLS recap + **finite-sample sampling distribution** `\hat\beta\sim\mathcal N(\beta,\sigma^2(X^\top X)^{-1})`
   under the Gaussian model (from Ch 3); **unbiasedness** `\mathbb{E}[\hat\beta]=\beta`
   (**proved**, via `\hat\beta=\beta+(X^\top X)^{-1}X^\top\varepsilon`).
2. **Gauss–Markov theorem** — OLS is the Best Linear Unbiased Estimator
   (**proved**); requires only `\mathbb{E}[\varepsilon]=0`, `\operatorname{Cov}(\varepsilon)=\sigma^2 I`
   (no normality).
3. **Estimating `\sigma^2`** — unbiased `\hat\sigma^2 = \text{RSS}/(n-p)`; the role of
   `n-p` degrees of freedom (stated, with intuition; full unbiasedness proof cited).
4. **Inference** — standard error `\operatorname{se}(\hat\beta_j)=\hat\sigma\sqrt{[(X^\top X)^{-1}]_{jj}}`;
   the t-statistic and **confidence interval** `\hat\beta_j \pm t_{n-p,1-\alpha/2}\,\operatorname{se}(\hat\beta_j)`
   (derived from the sampling distribution); one line on the F-test for joint
   significance. **State and cite** that `\text{RSS}/\sigma^2\sim\chi^2_{n-p}`
   independent of `\hat\beta` (hence exact t) — `[@hastie2009; @casella2002]`, no
   Cochran proof.
5. **Diagnostics & goodness of fit** — the **hat matrix** `H=X(X^\top X)^{-1}X^\top`
   as the Ch 1 projection (idempotent, symmetric); residuals; `R^2 = 1 - \text{RSS}/\text{TSS}`;
   the four assumptions (linearity, independence, homoscedasticity, normality) and
   what each one buys.
6. **Multicollinearity & the bias–variance tradeoff** — variance inflation: small
   eigenvalues / large `\kappa(X)` blow up `(X^\top X)^{-1}`; the **VIF**
   `\text{VIF}_j = 1/(1-R_j^2)`. Then the **bias–variance decomposition** of
   mean-squared error (**proved**): `\text{MSE} = \text{bias}^2 + \text{variance}`,
   motivating trading bias for variance.
7. **Ridge regression** — `\hat\beta_{\text{ridge}}=(X^\top X+\lambda I)^{-1}X^\top y`;
   shrinks coefficients and lifts every eigenvalue of `X^\top X` by `\lambda`
   (re-conditioning); **keystone (proved): ridge = MAP under `\beta\sim\mathcal N(0,\tau^2 I)`**
   with `\lambda=\sigma^2/\tau^2`. One line naming lasso / elastic-net and pointing
   to the Bayesian chapter.

### Worked Examples
- (a) Fit a tiny model: compute `\hat\beta`, `\hat\sigma^2`, one standard error, and a
  95% confidence interval by hand.
- (b) Two near-collinear channels: show the OLS standard error is huge and the VIF
  large, then apply ridge and watch the coefficients stabilize.

### Code Tie-in
NumPy/SciPy, runnable and self-contained:
- Fit OLS on synthetic media data; print coefficients with standard errors and 95% CIs.
- Crank channel correlation up; tabulate VIF and `\operatorname{se}(\hat\beta)`
  exploding.
- Fit ridge across a `\lambda` grid; plot coefficient paths shrinking/stabilizing;
  confirm `\hat\beta_{\text{ridge}}\to\hat\beta_{\text{OLS}}` as `\lambda\to 0`.

### Exercises
- **C — Conceptual:** what a confidence interval does/doesn't claim; why collinearity
  inflates standard errors; what ridge buys and what it costs.
- **B — By hand:** compute a CI; compute a VIF for a 2-predictor case; compute
  `\hat\beta_{\text{ridge}}` for a small `\lambda`.
- **P — Prove it:** OLS unbiased; Gauss–Markov (BLUE); ridge = Gaussian-prior MAP.
- **A — Applied / code:** coefficient-path plot vs `\lambda`; a bias–variance
  simulation (MSE of OLS vs ridge as collinearity rises).

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions`
(no inline links).

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only), covering the
sampling distribution, Gauss–Markov/BLUE, `\hat\sigma^2` and degrees of freedom, the
CI formula, the hat matrix and `R^2`, VIF, the bias–variance decomposition, and
ridge = `(X^\top X+\lambda I)^{-1}X^\top y` = Gaussian-prior MAP.

## Rigor level

First chapter of Part II — **ESL-style**: precise theorem statements with clean
linear-algebra proofs, always cashed out in the regression application. **Prove in
full**: OLS unbiasedness, Gauss–Markov (BLUE), the bias–variance decomposition, and
the keystone ridge = Gaussian-prior MAP. **State and cite** the distributional
mechanics (`\hat\sigma^2 \perp \hat\beta`, `\text{RSS}/\sigma^2\sim\chi^2_{n-p}`,
exact t/F) — `[@hastie2009; @casella2002]`, no Cochran's-theorem proof. Keep
everything anchored to the design matrix and the media-mix collinearity story.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned`
  inside `$$`, never bare `\begin{align}`; `$$` delimiters on their own lines so
  display-math counts stay even.
- Citations via existing `references.bib` keys: `@hastie2009`, `@casella2002`
  (and `@strang2016` where the Ch 1 projection/conditioning link is invoked).
- Code original, minimal, runnable against pinned `requirements.txt`
  (numpy, scipy, matplotlib); verify every numeric claim with NumPy before commit;
  the cell must run headless (`MPLBACKEND=Agg python3`).
- Match the established voice of Ch 1–3 and reuse their symbols: `X^\top X`,
  `\hat\beta`, `\kappa(X)`, `\sigma^2(X^\top X)^{-1}`, the loss `L`, the Gaussian
  likelihood `\ell`.

## Out of scope (YAGNI)

Generalized least squares / weighted least squares beyond a mention, the full
GLM/link-function family, lasso coordinate descent and LARS (lasso named only),
cross-validation theory (used at most as a one-line motivation in an exercise),
robust/quantile regression, and time-series error structure (autocorrelation is
named as an assumption to relax, not developed). Bayesian posterior computation is
the next chapter; hierarchical/multilevel models are the one after.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR
  (Quarto not installed locally).
- Every theory subsection ties back to the `y=X\beta+\varepsilon` / collinearity
  through-line.
- The four named proofs (unbiasedness, Gauss–Markov, bias–variance, ridge = MAP) are
  present and correct.
- The code tie-in runs top-to-bottom (headless) and demonstrates SE/VIF inflation
  under collinearity and ridge coefficient-path stabilization.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–3.
