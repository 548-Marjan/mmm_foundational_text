# Part II, Ch 2 — Bayesian Inference: Design Spec

**Date:** 2026-06-10
**Branch:** `ch5-bayesian-inference`
**File to author:** `parts/02-regression-bayes/02-bayesian-inference.qmd`
**Canonical anchors:** Gelman et al., *BDA3* (2013); Hoff (2009).

## Goal

Replace the scaffolded stub with a complete chapter that turns Chapter 4's ridge
point estimate into a full posterior distribution, develops the conjugate machinery
of Bayesian inference, and motivates the MCMC of Part III. Scope is **MMM-targeted
core**: every topic earns its place by a downstream use in Bayesian MMM or by setting
up Part III (sampling) and the hierarchical chapter. Deliverable this cycle is a
**full chapter draft** — prose, proofs, worked examples, a runnable Python tie-in,
full C/B/P/A exercises with appendix solutions, and the closing **Summary**
(auto-included convention for every chapter).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

## Through-line (anchor)

Chapter 4 ended on a cliffhanger: ridge is the **MAP estimate** — the single highest
point of a Gaussian posterior under the prior `\beta\sim\mathcal N(0,\tau^2 I)`. This
chapter delivers the **whole posterior** that point sat on top of. The anchor stays
the regression model `y = X\beta + \varepsilon`, but `\beta` is now genuinely random,
governed by a prior. Spine:

1. **The Bayesian object is a distribution, not a point.** prior `p(\beta)` →
   likelihood `p(y\mid\beta)` (Ch 3) → posterior `p(\beta\mid y)\propto p(y\mid\beta)p(\beta)`
   → posterior predictive. Estimation becomes conditioning.
2. **Conjugacy makes it tractable.** Beta–Binomial warm-up (does a channel convert?),
   then the Gaussian conjugate.
3. **Bayesian linear regression is the payoff.** Gaussian prior + Gaussian likelihood
   give `\beta\mid y\sim\mathcal N(\mu_{\text{post}},\Sigma_{\text{post}})`, whose mean
   **is the Chapter 4 ridge estimator** — the MAP was just the tip. Hence **credible
   intervals**, the honest answer to Ch 4's "how sure are you?", contrasted with
   frequentist confidence intervals.
4. **Prediction with honest uncertainty.** The posterior predictive propagates
   parameter uncertainty into sales forecasts — what budget planning needs.
5. **Priors are modeling choices.** Weakly-informative priors; the MMM-critical
   encoding of sign constraints (non-negative media effects).
6. **The wall.** Conjugacy is a luxury of Gaussian-linear models; the real MMM
   (adstock, saturation) is non-conjugate, the posterior has no closed form — which is
   precisely why Part III builds MCMC.

## Section-by-section content

### Motivation
Open on Ch 4's cliffhanger: ridge gave one number, the posterior mode, discarding the
rest of the distribution. Driving question (verbatim emphasis): *"Ridge handed you
`\hat\beta_{\text{TV}} = 1.8`. But how plausible is 1.5? Is 0 ruled out? And what will
next quarter's sales actually be?"* Only the full posterior answers these — and getting
it is just conditioning, which Chapter 3 already taught.

### Theory & Proofs
A ladder, each rung an MMM payoff:

1. **The Bayesian framework** — prior, likelihood, posterior
   `p(\beta\mid y)\propto p(y\mid\beta)p(\beta)`, the evidence
   `p(y)=\int p(y\mid\beta)p(\beta)\,d\beta` as normalizer, and the **posterior
   predictive** `p(\tilde y\mid y)=\int p(\tilde y\mid\beta)p(\beta\mid y)\,d\beta`.
   Estimation = conditioning.
2. **Conjugacy & the Beta–Binomial** (**proof**) — define a conjugate prior; derive the
   Beta–Binomial posterior `\text{Beta}(\alpha+k,\ \beta+n-k)` by kernel-matching;
   interpret the prior parameters as pseudo-counts. The clean 1-D warm-up.
3. **The Gaussian conjugate** (**proof**) — Normal prior `\mathcal N(\mu_0,\tau^2)` +
   Normal likelihood (known variance `\sigma^2`, `n` observations) → Normal posterior;
   **precisions add**, and the posterior mean is a precision-weighted average of prior
   mean and sample mean — proved by completing the square.
4. **Bayesian linear regression** (**KEYSTONE proof**) — prior `\beta\sim\mathcal N(0,\tau^2 I)`,
   Gaussian likelihood → `\beta\mid y\sim\mathcal N(\mu_{\text{post}},\Sigma_{\text{post}})`
   with `\Sigma_{\text{post}}=(X^\top X/\sigma^2 + I/\tau^2)^{-1}` and
   `\mu_{\text{post}}=(X^\top X+\lambda I)^{-1}X^\top y`, `\lambda=\sigma^2/\tau^2`.
   **The posterior mean is exactly the Chapter 4 ridge estimator.**
5. **Credible intervals & posterior summaries** — the credible interval and its honest
   interpretation (a direct probability statement about `\beta`), contrasted sharply
   with the frequentist confidence interval; posterior mean / median / MAP as point
   summaries.
6. **The posterior predictive** (**derivation**) — for Bayesian linear regression,
   `\tilde y\mid y\sim\mathcal N(\tilde x^\top\mu_{\text{post}},\ \tilde x^\top\Sigma_{\text{post}}\tilde x + \sigma^2)`,
   decomposing forecast uncertainty into **parameter uncertainty + irreducible noise**.
7. **Priors in practice & the wall** — informative vs weakly-informative priors; the
   MMM-critical encoding of **sign constraints** (non-negative media effects via a
   half-normal prior), which *breaks* Gaussian conjugacy; adstock/saturation making the
   posterior intractable, so closed forms die and **Part III's MCMC** becomes necessary.
   State and cite the general no-closed-form fact `[@gelman2013]`.

### Worked Examples
- (a) **Beta–Binomial by hand**: prior `\text{Beta}(2,2)`, observe 8 conversions in
  `n=10`, posterior `\text{Beta}(10,4)`; posterior mean `10/14\approx0.714` vs MLE
  `0.8`, showing shrinkage toward the prior mean.
- (b) **Bayesian linear regression on a tiny design**: compute `\Sigma_{\text{post}}`,
  `\mu_{\text{post}}`, confirm `\mu_{\text{post}}` equals the ridge estimate, and give a
  95% credible interval for a coefficient — contrasted with Ch 4's confidence interval.

### Code Tie-in
NumPy/SciPy, runnable and self-contained:
- Beta–Binomial: prior → posterior updating as data arrives (plot prior/likelihood/posterior,
  and posterior sharpening with more observations).
- Bayesian linear regression in closed form: verify `\mu_{\text{post}}` equals
  `(X^\top X+\lambda I)^{-1}X^\top y` (ridge); draw posterior coefficient samples; plot
  credible intervals; show a posterior-predictive band; sweep prior strength `\tau^2`
  to interpolate from OLS (vague prior, `\tau^2\to\infty`) to heavy shrinkage.

### Exercises
- **C — Conceptual:** prior/posterior/predictive interpretation; credible vs confidence
  interval; why a non-negativity prior breaks conjugacy and forces computation.
- **B — By hand:** a Beta–Binomial update; a Normal–Normal posterior mean (precision-weighted);
  compute `\mu_{\text{post}}` and show it equals ridge for a small design.
- **P — Prove it:** Beta–Binomial conjugacy; the Gaussian (Normal–Normal) posterior; the
  Bayesian-linear-regression posterior with `\mu_{\text{post}}` = ridge.
- **A — Applied / code:** Beta–Binomial sequential updating; Bayesian linear regression
  posterior + predictive with credible bands; a prior-sensitivity sweep over `\tau^2`.

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions`
(no inline links). Heading format `## Chapter 5 — Bayesian Inference` to match siblings.

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only), covering Bayes'
update, conjugacy, Beta–Binomial, the Gaussian conjugate (precisions add), the Bayesian
linear-regression posterior (`\Sigma_{\text{post}}`, `\mu_{\text{post}}` = ridge),
credible intervals, the posterior predictive variance decomposition, and the conjugacy
wall pointing to MCMC.

## Rigor level

Conceptual heart of Part II — **BDA3 / Hoff-style**: clean analytical conjugate
derivations, every result cashed out on the regression/MMM model, no measure theory.
**Prove in full**: Beta–Binomial conjugacy (kernel-matching), the Gaussian conjugate
posterior (completing the square, precisions add), and the keystone Bayesian
linear-regression posterior with mean = ridge. **Derive** the Gaussian posterior
predictive (a Gaussian integral); **state and cite** that without conjugacy the posterior
and predictive have no closed form, handing the baton to Part III `[@gelman2013]`. Keep
priors honest and practical (weakly-informative defaults, pseudo-count reading, the
non-negativity/sign-constraint point that breaks conjugacy). Anchor to the design matrix
and the media-mix story; make the credible-vs-confidence contrast concrete on a channel
coefficient.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned` inside
  `$$`, never bare `\begin{align}`; `$$` delimiters on their own lines so display-math
  counts stay even.
- Citations via existing `references.bib` keys: `@gelman2013`, `@hoff2009` (and
  `@hastie2009` where the Ch 4 ridge link is invoked, `@casella2002` for distributions).
- Code original, minimal, runnable against pinned `requirements.txt`
  (numpy, scipy, matplotlib); verify every numeric claim with NumPy before commit; the
  cell must run headless (`MPLBACKEND=Agg python3`).
- Match the established voice of Ch 1–4 and reuse their symbols: `X^\top X`, `\hat\beta`,
  ridge `(X^\top X+\lambda I)^{-1}X^\top y`, `\lambda=\sigma^2/\tau^2`, the Gaussian
  likelihood `\ell`, `\mathcal N(\mu,\Sigma)`, the Cholesky factor `L`.
- Git identity in a fresh worktree: `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"`.

## Out of scope (YAGNI)

The full exponential-family / sufficiency theory, a broad conjugate catalog beyond
Beta–Binomial and the Gaussian (Gamma–Poisson, Dirichlet–Multinomial named at most once),
Normal–Inverse–Gamma for jointly unknown `\sigma^2` (we take `\sigma^2` known/fixed for
the closed forms, noting the extension), Bayesian decision theory (loss/Bayes risk),
Bayes factors and formal model comparison, empirical Bayes, and any actual MCMC
algorithm (that is Part III). Hierarchical/multilevel priors are the next chapter.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR
  (Quarto not installed locally).
- Every theory subsection ties back to the posterior / ridge-mode through-line.
- The three named conjugate proofs (Beta–Binomial, Gaussian conjugate, Bayesian linear
  regression with mean = ridge) are present and correct, plus the predictive derivation.
- The code tie-in runs top-to-bottom (headless) and verifies `\mu_{\text{post}}` equals
  the ridge solution and shows credible/predictive bands.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–4.
