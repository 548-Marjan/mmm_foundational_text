# Chapter 15 — Building & Fitting an MMM: Design Spec

**Date:** 2026-06-15
**Status:** Approved (brainstorm complete). Next: implementation plan via `writing-plans`.
**Part:** V — Marketing Mix Modeling: Modeling (second chapter).
**File:** `parts/05-mmm-modeling/02-building-fitting.qmd` (currently a template stub).
**Canonical anchors:** `@jin2017` (the MMM DGP), `@gelman2013` (Bayesian workflow). Cross-refs
to Ch5 (MAP/priors), Ch9 (NUTS), Ch10 (diagnostics), Ch14 (the DGP and its wound).
**House rule:** name **no** probabilistic-programming or MMM library in the text. Fitting is done
with `scipy.optimize` (a general scientific tool, fine to name) and a hand-rolled Metropolis
sampler reusing Part III's machinery.

---

## 1. Role in the book

Chapter 15 delivers the preface's promise to "build a full Bayesian MMM end-to-end." It inverts
Chapter 14: where Ch14 defined the data-generating process and simulated a synthetic sales series,
Ch15 takes that series and **recovers the parameters by Bayesian inference** — then confronts the
identification wound Ch14 ended on, in practice. It is largely a *synthesis* chapter: it reuses
Part II (MAP, priors — Ch5), Part III (sampling and the convergence diagnostics — Ch9, Ch10), and
the optimization/curvature tools of Ch2 and Ch11. Little new machinery; much integration.

The controlled payoff: because we built the data, we know the ground truth. Fitting can be judged
by **parameter recovery**, not just fit — and the chapter shows recovery succeeding when spend
varies enough and failing (a wide posterior ridge) when it does not, handing off to Part VI.

The chapter follows the fixed template: Motivation → Theory & Proofs → Worked Examples → Code
Tie-in → Exercises (C/B/P/A) → Summary (auto-included), with gated appendix solutions.

## 2. Scope decisions (locked in brainstorming)

- **Theorem weight:** two real proofs plus one short rung (below). Not a pure-workflow chapter.
- **Fitting approach in code:** **MAP + Laplace via `scipy.optimize`, plus a short hand-rolled
  Metropolis** — no PPL named. The Hessian at the mode exhibits the identification ridge as a
  near-zero eigenvalue (numerical face of the Fisher-information proof); Metropolis shows the
  posterior is wide along that eigenvector.
- **Model scope:** single-series (national) MMM. Hierarchical/geo pooling is Ch6's machinery and is
  mentioned only in passing; time-varying coefficients are Ch16 (DLM). Keep Ch15 to the
  static single-series Bayesian fit — well-bounded.

## 3. Theory rungs

1. **From DGP to likelihood.** Ch14's generative model becomes a Gaussian likelihood
   $y_t \sim \mathcal{N}\big(\mu_t(\theta),\,\sigma^2\big)$, $\mu_t(\theta)=\text{baseline}_t+\sum_i \beta_i\,g(a_\lambda(x_{i,t}); K_i, n_i)$. Because $g$ is a rational power function of $(K,n)$, the
   posterior is non-conjugate — closed-form (Ch5) is unavailable. Short rung.
2. **Priors.** Weakly-informative priors carrying domain knowledge and regularizing the wound:
   $\lambda \sim \text{Beta}(a,b)$ on $[0,1)$; $K \sim$ positive (half-normal or log-normal);
   $n \sim$ positive; $\beta_i \sim$ positive half-normal (effects are non-negative); $\sigma \sim$
   half-normal. State the logic: weak data ⇒ priors must carry identification (bridge from Ch14).
3. **Posterior geometry and sampling.** Ch14's identification ridge reappears as strong posterior
   correlation; the geometry is curved and, in the confined-spend case, nearly flat along the ridge.
   NUTS (Ch9) is the production tool for this curvature; divergences (Ch10) flag where step size
   fails. (We *use* NUTS conceptually and demonstrate with MAP+Laplace+Metropolis in code.)
4. **The fitting workflow.** Prior predictive check → fit → convergence diagnostics (Ch10: $\hat R$,
   ESS, divergences) → posterior predictive check → **parameter recovery vs known ground truth**.
5. **Identification in practice.** Recovery succeeds with varied spend (P1 full-rank Fisher); the
   ridge stays wide with confined spend (P1 singular Fisher; P2 prior dominates). Motivates Part VI.

## 4. Proof slate (two proofs + one rung)

| # | Result | Statement |
|---|---|---|
| **Rung (non-conjugacy)** | No closed-form posterior | The likelihood factor $g(a_\lambda(x); K, n) = (a^n)/(K^n + a^n)$ is a rational power function of the parameters $(K,n,\lambda)$; multiplying by any standard prior yields a kernel that is not a recognized distribution in $(K,n,\lambda)$, so the posterior has no closed form and must be approximated (MAP/Laplace) or sampled (MCMC). |
| **P1 (keystone)** | Local identifiability via Fisher information | For the regression $y_t \sim \mathcal N(\mu_t(\theta),\sigma^2)$, the parameter $\theta$ is **locally identified** at $\theta_0$ iff the Fisher information $I(\theta_0) = \tfrac{1}{\sigma^2}\sum_t \nabla_\theta\mu_t\,\nabla_\theta\mu_t^\top$ is nonsingular. Proof: $I$ is the Gauss–Newton Hessian of the negative log-likelihood; it is singular iff the gradient vectors $\{\nabla_\theta\mu_t\}$ span a subspace of dimension $< \dim\theta$, i.e. iff a nonzero direction $d$ has $\nabla_\theta\mu_t^\top d = 0$ for all $t$ — a direction along which $\mu$ (hence the likelihood) is flat to first order. Spend confined to a single operating point makes every $\nabla_\theta\mu_t$ collinear (recovering Ch14's P4); spend spanning enough distinct operating points that the gradients are linearly independent makes $I$ nonsingular. **Answers Ch14's open question — how much spend variation is enough: enough distinct operating points that the parameter-gradient vectors are independent.** |
| **P2** | Prior dominates the unidentified ridge | Suppose the log-likelihood is flat along a direction $d$ at $\theta_0$: $\ell(\theta_0 + s\,d) = \ell(\theta_0)$ for $s$ in a neighborhood. Then along that ridge the posterior satisfies $p(\theta_0 + s\,d \mid y) \propto p(\theta_0 + s\,d)$ — the likelihood contributes a constant, so the posterior's variation along $d$ is exactly the (renormalized) prior. Priors therefore carry the identification load on unidentified directions. Callback: this is the Bayesian face of Ch5's ridge regularization. |

## 5. Worked examples + anchor numbers

All numerics NumPy-verified in the plan.

- **WE1 — Fisher singular vs nonsingular (the keystone in arithmetic).** Single channel,
  $K=3$, $n=2$, $\beta=4$; parameters of interest $(\beta, K)$. Gradients
  $\partial\mu/\partial\beta = g(x;K,n)$ and $\partial\mu/\partial K = \beta\,g_K(x;K,n)$ with
  $g_K = -\,nK^{n-1}x^n/(K^n+x^n)^2$. At a single operating point $x_0=3$:
  $\nabla\mu = (0.5,\,-\tfrac23)$ for every $t$ ⇒ $\det I = 0$ (rank 1, the $\beta$–$K$ ridge).
  Add a second point $x_1=6$ (where $g=0.8$, $\nabla\mu=(0.8,\,-0.4267)$): the two gradients are
  not parallel, $\det(J^\top J) = 0.1024 > 0$ ⇒ identified.
- **WE2 — parameter recovery.** Fit a wide-spend synthetic series; report MAP / posterior means
  with intervals bracketing the true $(\beta, K, n, \lambda)$. (Exact recovered values pinned in
  the plan after the fit is run.)
- **WE3 — the ridge made numerical.** Confined-spend series: the Hessian at the MAP has a near-zero
  smallest eigenvalue; the posterior (Metropolis) is wide along the corresponding $\beta$–$K$
  eigenvector while the data fit remains good. The wound, quantified.

## 6. Code tie-in

A single runnable `{python}` cell. It (a) regenerates a single-channel synthetic series from known
truth (reusing Ch14's `geometric_adstock` and `hill`); (b) defines the log-posterior (Gaussian
likelihood + the priors of Rung 2); (c) **Scenario A — wide spend:** `scipy.optimize.minimize` finds
the MAP recovering the truth, the Hessian is well-conditioned (all eigenvalues sizable), and a short
hand-rolled Metropolis sampler returns tight posteriors bracketing truth; (d) **Scenario B — confined
spend:** the MAP still fits, but the Hessian has a tiny smallest eigenvalue and Metropolis shows the
posterior spread wide along the ridge eigenvector. Plots: $\beta$–$K$ posterior scatter (blob in A vs
ridge in B) and the Hessian eigenvalue spectra. Asserts pin: $\det I = 0$ at one operating point and
$0.1024$ at $\{3,6\}$; MAP recovers truth in A to a tolerance; the smallest Hessian eigenvalue in B is
far smaller than in A. Runs headless under `MPLBACKEND=Agg`.

## 7. Exercises, Summary, housekeeping

- **Exercises** C/B/P/A (self-contained); appendix solutions appended to `appendix/solutions.qmd`
  as `## Chapter 15 — Building & Fitting an MMM`, gated by
  `::: {.content-visible when-meta="show-solutions"}`, before the file's final closing `:::`.
- **Summary** auto-included; **Key concepts** and **Key identities** both **bulleted** (inline math
  only) — per the formatting-consistency rule.
- **Anchor line:** the stub has no anchor line; add one directly under the H1 reading exactly
  `*Canonical anchors: @jin2017; @gelman2013.*` (names no library).

## 8. Conventions (enforced per task)

- KaTeX: `aligned` inside `$$ … $$` (never bare `\begin{align}`); `$$` delimiters on their own lines;
  even `$`-count per line. **Never** `\begin{psmallmatrix}` (PDF/LuaLaTeX undefined; passes
  KaTeX/HTML but breaks the PDF build) — use `bmatrix`/`pmatrix`/`smallmatrix`.
- Build-based verification: the code cell runs headless; the real gate is CI `quarto render`
  (HTML + PDF) on the PR — watch **both** jobs to green before reporting done.
- Library-agnostic: no PPL/MMM library named in text; `numpy`/`scipy`/`matplotlib` are fine.

## 9. Critical files

- `parts/05-mmm-modeling/02-building-fitting.qmd` — the chapter (replace stub body; keep H1
  `# Building & Fitting an MMM`; set the anchor line).
- `appendix/solutions.qmd` — append the Chapter 15 solutions block before the final `:::`.
- Read-only exemplars: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (immediate predecessor; reuse its
  `geometric_adstock`/`hill`), `parts/02-regression-bayes/02-bayesian-inference.qmd` (Ch5 MAP/priors
  voice), `parts/03-sampling/04-model-checking.qmd` (Ch10 diagnostics vocabulary).
