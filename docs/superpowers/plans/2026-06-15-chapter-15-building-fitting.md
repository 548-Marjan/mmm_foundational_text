# Chapter 15 — Building & Fitting an MMM Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 15, which inverts Ch14's data-generating process by Bayesian inference — recovering the response-curve parameters from the synthetic series, proving when that is possible (Fisher-information identifiability) and how priors carry the load when it is not (prior-dominates-the-ridge), and demonstrating both in code via MAP + Laplace + a hand-rolled Metropolis sampler.

**Architecture:** Replace the template stub `parts/05-mmm-modeling/02-building-fitting.qmd` with a full chapter following the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary). Append a gated solutions block to `appendix/solutions.qmd`. No bundled cross-file fixes this chapter (only the in-file anchor line is added).

**Tech Stack:** Quarto (`.qmd`), KaTeX math, one `{python}` cell (numpy + scipy.optimize + matplotlib, run headless under `MPLBACKEND=Agg`). Real verification gate is CI `quarto render` (HTML + PDF).

---

## Conventions (enforce in EVERY task)

- **KaTeX:** display math uses `$$ … $$` with delimiters **on their own lines**; multi-line uses `aligned` **inside** `$$ … $$` — never bare `\begin{align}`. Every line has an **even** count of `$`. Inline math uses single `$`.
- **PDF safety:** **never** `\begin{psmallmatrix}` (undefined in CI LuaLaTeX — renders in KaTeX so HTML passes but the PDF stage fails). Use `bmatrix`, `pmatrix`, or `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]`.
- **Library-agnostic:** name **no** probabilistic-programming or MMM library in prose. `numpy`/`scipy`/`matplotlib` are fine. Anchors: `@jin2017`, `@gelman2013`.
- **Notation (fixed):** $\theta=(\beta,K,n,\lambda,\ldots)$ parameters; $\mu_t(\theta)$ mean function; $g(x;K,n)=x^n/(K^n+x^n)$ Hill; $a_\lambda$ adstock; $I(\theta)$ Fisher information; channel index $i$; $T$ series length.
- **Verification:** the single code cell runs top-to-bottom under `MPLBACKEND=Agg python3` with all asserts passing; figures end `plt.show()`. Pinned numbers below are NumPy-verified.
- **Commits:** one per task, prefix `feat(ch15): …`. Run from the worktree root.

## Verified anchor numbers (do not change)

- **Fisher (WE1), single channel $K=3, n=2, \beta=4$, parameters $(\beta,K)$:** gradient
  $\nabla\mu(x) = \big(g(x;K,n),\ \beta\,g_K(x;K,n)\big)$ with $g_K = -nK^{n-1}x^n/(K^n+x^n)^2$.
  At $x=3$: $\nabla\mu=(0.5,\,-\tfrac23)$ (since $g_K(3;3,2)=-\tfrac16$). At $x=6$: $g(6;3,2)=0.8$,
  $\nabla\mu=(0.8,\,-0.42667)$. One operating point $\{3\}$: $\det I = 0$. Two points $\{3,6\}$:
  $\det(J^\top J) = 0.1024$.
- **Recovery / ridge (illustrative, fixed seeds — assert inequalities, not these exact values):**
  Scenario A (wide spend) MAP $\approx(\beta\,3.60,\,K\,2.65,\,n\,2.50,\,\lambda\,0.47)$ vs truth
  $(4,3,2,0.5)$; Scenario B (confined spend) MAP $K\approx4.88$. Smallest $\beta$–$K$ Hessian
  eigenvalue: A $\approx123$, B $\approx0.69$ (ratio $\approx178$). Posterior $K$ std: A $\approx0.18$,
  B $\approx1.58$. $|K_A-3|\approx0.35 < |K_B-3|\approx1.88$.
- **Robust asserts (seed-independent, large margins):** (1) $\det I = 0$ at one point and
  $\approx0.1024$ at $\{3,6\}$; (2) smallest Hessian eigenvalue A $> 20\times$ B; (3)
  $|K_A-3| < |K_B-3|$; (4) posterior $K$ std B $>$ A.

## File structure

- **`parts/05-mmm-modeling/02-building-fitting.qmd`** — the chapter. Replace stub body; keep H1 `# Building & Fitting an MMM`; add the anchor line.
- **`appendix/solutions.qmd`** — append `## Chapter 15 — Building & Fitting an MMM` solutions block (gated) before the file's final closing `:::`.
- Read-only exemplars: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (predecessor; reuse `geometric_adstock`/`hill`), `parts/02-regression-bayes/02-bayesian-inference.qmd` (Ch5 MAP/prior voice), `parts/03-sampling/04-model-checking.qmd` (Ch10 diagnostics vocabulary).

---

### Task 1: Front matter (anchor line) and Motivation

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd`

- [ ] **Step 1: Fix front matter.** Keep H1 `# Building & Fitting an MMM`. Immediately under it add the anchor line (the stub has none):

```
*Canonical anchors: @jin2017; @gelman2013.*
```

Delete the `::: {.callout-note … Stub …}` block.

- [ ] **Step 2: Write `## Motivation`** (3–5 paragraphs). Beats:
  - Ch14 built the data-generating process and simulated a synthetic sales series from known parameters. This chapter runs the machine backward: given only spend and sales, recover the adstock and saturation parameters — the inverse problem.
  - This is where the whole stack converges: the likelihood comes from Ch14's structural equation, the priors and MAP from Ch5, the sampler from Ch8–9, the convergence and predictive checks from Ch10. Little new theory; the chapter is the integration the preface promised — a Bayesian MMM end-to-end.
  - The controlled advantage of synthetic data: we know the truth, so success is *parameter recovery*, not merely good fit. A model can fit the in-sample series beautifully and still have the wrong curve.
  - Plant the tension carried from Ch14's wound: recovery will succeed when spend varies enough and fail — a wide posterior ridge — when it does not. This chapter makes that precise (when is the curve identified?) and shows priors carrying the load on the unidentified directions. The failure mode is the motivation for Part VI.

- [ ] **Step 3: Verify & commit.** `grep -c '\$' …` even; `grep -c 'begin{align}' …` → 0. Confirm no library named: `grep -ciE 'pymc|stan|orbit|lightweight_mmm' …` → 0.
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): front matter and Motivation"
```

---

### Task 2: Theory rungs 1–2 (likelihood, non-conjugacy, priors)

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd`

- [ ] **Step 1: Open `## Theory & Proofs`; write Rung 1 — from DGP to likelihood.** Ch14's generative model becomes a Gaussian likelihood. With mean function
$$
\mu_t(\theta) = \text{baseline}_t + \sum_{i=1}^{C} \beta_i\, g\big(a_\lambda(x_{i,t}); K_i, n_i\big),
$$
the observations are $y_t \mid \theta \sim \mathcal N(\mu_t(\theta), \sigma^2)$, and the log-likelihood is $\ell(\theta) = -\tfrac{1}{2\sigma^2}\sum_t (y_t - \mu_t(\theta))^2 - T\log\sigma + \text{const}$. State that, unlike Ch5's linear-Gaussian model, $\mu_t$ is **nonlinear** in $(K,n,\lambda)$.

- [ ] **Step 2: Write the non-conjugacy rung-proof.** Short proposition: the posterior has no closed form. Argument: the likelihood's dependence on $(K,n,\lambda)$ enters through $g(a; K,n)=a^n/(K^n+a^n)$, a rational power function; the product of this with any standard prior on $(K,n,\lambda)$ is not the kernel of a recognized distribution (no conjugate family closes under it), so the normalizing constant is not available analytically. Conclusion: approximate the posterior (MAP + Laplace) or sample it (MCMC). One paragraph; not a heavy proof.

- [ ] **Step 3: Write Rung 2 — priors.** Motivate weakly-informative priors that encode domain knowledge and regularize: $\lambda \sim \text{Beta}(2,2)$ on $[0,1)$ (carryover between none and total); $K \sim \text{HalfNormal}$ (positive half-saturation, scaled to plausible spend); $n \sim \text{HalfNormal}$ (positive shape, prior mass near $1$–$2$); $\beta_i \sim \text{HalfNormal}$ (effects non-negative — advertising does not reduce sales); $\sigma \sim \text{HalfNormal}$. State the throughline to be proved later: when data are weak (the wound), the prior is what makes the posterior proper and concentrated.

- [ ] **Step 4: Verify & commit.** Even `$`-count; no `begin{align}`/`psmallmatrix`.
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): likelihood, non-conjugacy, and priors"
```

---

### Task 3: Rung 3 + Proof P1 (Fisher-information identifiability, keystone)

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd`

- [ ] **Step 1: Write Rung 3 — posterior geometry and sampling.** The identification ridge of Ch14 reappears as strong posterior correlation and, in the confined-spend case, a nearly flat direction. This curved, ill-conditioned geometry is exactly what gradient-based samplers handle: NUTS (Ch9) adapts to it, and divergences (Ch10) flag where a fixed step size fails. State that the chapter will *demonstrate* the geometry concretely via the mode and its curvature (MAP + Laplace) plus a Metropolis sampler, rather than re-deriving NUTS.

- [ ] **Step 2: Write Proof P1 — local identifiability via Fisher information (keystone).** State the proposition and prove it in full prose:
  - **Setup.** For $y_t \sim \mathcal N(\mu_t(\theta),\sigma^2)$ iid, the Fisher information is
  $$
  I(\theta) = \frac{1}{\sigma^2}\sum_{t} \nabla_\theta\mu_t(\theta)\,\nabla_\theta\mu_t(\theta)^\top.
  $$
  Derive it: $\nabla_\theta\ell = \tfrac{1}{\sigma^2}\sum_t (y_t-\mu_t)\nabla_\theta\mu_t$; the Hessian is $\nabla^2_\theta\ell = \tfrac{1}{\sigma^2}\sum_t\big[(y_t-\mu_t)\nabla^2\mu_t - \nabla\mu_t\nabla\mu_t^\top\big]$; taking $-\mathbb E[\cdot]$ and using $\mathbb E[y_t-\mu_t]=0$ kills the first term, leaving $I(\theta)$ above.
  - **Identifiability criterion.** $\theta_0$ is locally identified iff $I(\theta_0)$ is nonsingular. Proof of the key direction: if $I(\theta_0)$ is singular, pick $d\neq0$ with $I(\theta_0)d=0$; then $0 = d^\top I d = \tfrac{1}{\sigma^2}\sum_t (\nabla\mu_t^\top d)^2$, forcing $\nabla\mu_t^\top d = 0$ for every $t$. So moving $\theta$ along $d$ leaves every $\mu_t$ unchanged to first order — the mean function, hence the likelihood, is flat along $d$, and $\theta_0$ is not locally identified. Conversely a nonsingular $I$ admits no such flat direction.
  - **Apply to the MMM.** $\nabla_\theta\mu_t$ depends on the operating point $x_t$ through the transforms. If spend is confined to one operating point $x_0$, every $\nabla\mu_t$ equals one vector $v$, so $I = (T/\sigma^2)\,vv^\top$ has rank $1 < \dim\theta$: singular — **recovering Ch14's P4**. If spend spans enough distinct operating points that the gradient vectors are linearly independent, $I$ is full rank and $\theta_0$ is identified.
  - **Answer Ch14's open question.** "How much variation is enough?" — enough distinct operating points that the parameter-gradient vectors $\{\nabla_\theta\mu_t\}$ span the parameter space. For the two-parameter $(\beta,K)$ slice this means at least two operating points with non-parallel gradients (the WE1 arithmetic).

- [ ] **Step 3: Verify & commit.**
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): posterior geometry and Fisher-information identifiability (P1)"
```

---

### Task 4: Rungs 4–5 + Proof P2 (prior dominates the ridge)

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd`

- [ ] **Step 1: Write Rung 4 — the fitting workflow.** The disciplined sequence: prior predictive check (do prior draws produce plausible sales?) → fit (MAP/Laplace or MCMC) → convergence diagnostics ($\hat R$, ESS, divergences — Ch10) → posterior predictive check (does the fitted model reproduce held-out features?) → **parameter recovery** against the known truth (the synthetic-data privilege). One paragraph per stage, concise.

- [ ] **Step 2: Write Proof P2 — prior dominates the unidentified ridge.** Proposition and proof:
  - Suppose the log-likelihood is flat along a direction $d$ at $\theta_0$: $\ell(\theta_0 + s\,d) = \ell(\theta_0)$ for $s$ in a neighborhood (an unidentified ridge, as P1 produces under confined spend). Then $L(\theta_0+s\,d) = e^{\ell(\theta_0)} =: L_0$ is constant.
  - The posterior is $p(\theta\mid y) \propto p(\theta)\,L(\theta)$, so along the ridge $p(\theta_0+s\,d\mid y) \propto p(\theta_0+s\,d)\,L_0 \propto p(\theta_0+s\,d)$.
  - Normalizing along the ridge, the conditional posterior density of $s$ is $p(\theta_0+s\,d)\big/\!\int p(\theta_0+s'd)\,ds'$ — exactly the prior restricted to the ridge. The likelihood contributes nothing along $d$; the prior alone determines the posterior there.
  - **Interpretation.** On identified directions the data update the prior; on unidentified directions the posterior *is* the prior. Priors therefore carry the identification load — the Bayesian face of the ridge regularization in Ch5, and the reason a well-chosen prior keeps the MMM posterior proper even at the wound.

- [ ] **Step 3: Write Rung 5 — identification in practice.** Tie the two proofs to what the code will show: with wide-ranging spend the Fisher information is well-conditioned and recovery succeeds; with confined spend a near-zero eigenvalue appears (the ridge), recovery of the shape parameters fails, and the posterior there is the prior. This residual non-identification — present even with perfect sampling — is what Part VI's experiments and calibration exist to break.

- [ ] **Step 4: Verify & commit.**
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): fitting workflow, prior-dominates-the-ridge (P2), identification in practice"
```

---

### Task 5: Worked Examples

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd` (replace `## Worked Examples` stub)

- [ ] **Step 1: WE1 — Fisher singular vs nonsingular (keystone in arithmetic).** Single channel, $K=3$, $n=2$, $\beta=4$, parameters $(\beta,K)$. Gradients $\partial\mu/\partial\beta = g(x;K,n)$, $\partial\mu/\partial K = \beta\,g_K(x;K,n)$ with $g_K = -nK^{n-1}x^n/(K^n+x^n)^2$. At $x_0=3$: $g(3)=0.5$, $g_K(3)=-\tfrac16$, so $\nabla\mu=(0.5,-\tfrac23)$. With all data at $x_0$, $I \propto \nabla\mu\,\nabla\mu^\top$ has $\det = 0$ — the $\beta$–$K$ ridge. Add $x_1=6$: $g(6)=0.8$, $\nabla\mu=(0.8,-0.42667)$; the two gradients are not parallel, so $\det(J^\top J)=0.1024>0$ — identified. State the lesson: identifiability is about the *spread of operating points*, not the number of weeks.

- [ ] **Step 2: WE2 — parameter recovery (wide spend).** Describe fitting the wide-spend synthetic series and recovering $(\beta,K,n,\lambda)\approx(3.60,2.65,2.50,0.47)$ against truth $(4,3,2,0.5)$, with credible intervals bracketing the truth. Emphasize: recovery, not just fit — the synthetic-data privilege.

- [ ] **Step 3: WE3 — the ridge made numerical (confined spend).** Same model, spend confined near $x_0=3$. The MAP still fits the data, but $K$ is badly recovered ($\approx4.9$ vs $3$), the $\beta$–$K$ Hessian has a tiny smallest eigenvalue ($\approx0.69$ vs $\approx123$ for wide spend), and the posterior $K$ is wide ($\text{std}\approx1.6$ vs $\approx0.18$). Quantifies the wound: good fit, unrecoverable shape.

- [ ] **Step 4: Verify & commit.** Even `$`-count; confirm `grep -c 'Worked numerical example' …` → 0.
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): worked examples (Fisher, recovery, the ridge)"
```

---

### Task 6: Code Tie-in — fit the MMM, exhibit the ridge

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd` (replace `## Code Tie-in` stub)

- [ ] **Step 1: Write the section intro** (2–3 sentences): we fit the single-channel model two ways — wide-ranging spend and confined spend — using `scipy.optimize` for the MAP and a short hand-rolled Metropolis sampler for the posterior, and read the identification ridge off the Hessian's smallest eigenvalue. Fixed seeds make it reproducible.

- [ ] **Step 2: Write the single `{python}` cell.** Use exactly this (verified headless and reproducible):

```{python}
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def geometric_adstock(x, lam):
    out = np.empty_like(x, dtype=float); p = 0.0
    for t, xt in enumerate(x):
        p = xt + lam * p; out[t] = p
    return out

def hill(x, K, n):
    x = np.asarray(x, dtype=float); return x**n / (K**n + x**n)

def hill_K(x, K, n):
    x = np.asarray(x, dtype=float); return -n * K**(n-1) * x**n / (K**n + x**n)**2

# --- WE1: Fisher information is singular at one operating point, full-rank at two ---
beta, K, n = 4.0, 3.0, 2.0
def grad_mu(x): return np.array([hill(x, K, n), beta * hill_K(x, K, n)])
J1 = np.array([grad_mu(3.0)])
J2 = np.array([grad_mu(3.0), grad_mu(6.0)])
assert np.isclose(np.linalg.det(J1.T @ J1), 0.0)
assert np.isclose(np.linalg.det(J2.T @ J2), 0.1024)

# --- synthetic data from known truth ---
gen = np.random.default_rng(15)
K_t, n_t, lam_t, beta_t, b0_t, sig_t = 3.0, 2.0, 0.5, 4.0, 10.0, 0.3
T = 104; t = np.arange(T)
def make(x):
    return b0_t + beta_t * hill(geometric_adstock(x, lam_t), K_t, n_t) + gen.normal(0, sig_t, T)
xA = np.clip(3.0 + 2.5 * np.sin(2*np.pi*t/52) + gen.normal(0, 1.2, T), 0.05, None)  # wide spend
xB = np.clip(3.0 + gen.normal(0, 0.15, T), 0.05, None)                              # confined spend
yA, yB = make(xA), make(xB)

# --- negative log-posterior (Gaussian likelihood + weakly-informative priors) ---
def neg_log_post(th, x, y):
    beta, K, n, lam, b0, ls = th
    if not (0 < lam < 1 and K > 0 and n > 0 and beta > 0):
        return 1e12
    sig = np.exp(ls)
    mu = b0 + beta * hill(geometric_adstock(x, lam), K, n)
    ll = -0.5*np.sum((y - mu)**2)/sig**2 - T*ls
    lp = (-0.5*(beta/5)**2 - 0.5*(K/5)**2 - 0.5*(n/3)**2
          + np.log(lam) + np.log(1-lam)            # Beta(2,2)
          - 0.5*((b0-10)/5)**2 - 0.5*(sig/2)**2)
    return -(ll + lp)

start = np.array([3.0, 3.0, 2.0, 0.5, 10.0, np.log(0.5)])
opt = dict(maxiter=20000, xatol=1e-7, fatol=1e-7)
mapA = minimize(neg_log_post, start, args=(xA, yA), method="Nelder-Mead", options=opt).x
mapB = minimize(neg_log_post, start, args=(xB, yB), method="Nelder-Mead", options=opt).x

# --- beta-K Hessian (Laplace curvature) at each MAP ---
def hess_betaK(th, x, y, h=1e-3):
    H = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            pp, pm, mp, mm = (th.copy() for _ in range(4))
            pp[a]+=h; pp[b]+=h; pm[a]+=h; pm[b]-=h; mp[a]-=h; mp[b]+=h; mm[a]-=h; mm[b]-=h
            H[a, b] = (neg_log_post(pp,x,y)-neg_log_post(pm,x,y)
                       -neg_log_post(mp,x,y)+neg_log_post(mm,x,y))/(4*h*h)
    return H
eigA = np.linalg.eigvalsh(hess_betaK(mapA, xA, yA))
eigB = np.linalg.eigvalsh(hess_betaK(mapB, xB, yB))

# --- short hand-rolled Metropolis for the posterior (fixed seed) ---
def metropolis(x, y, start, nsteps=40000):
    mc = np.random.default_rng(0); th = start.copy(); cur = neg_log_post(th, x, y); keep = []
    sc = np.array([0.15, 0.2, 0.1, 0.02, 0.15, 0.05])
    for i in range(nsteps):
        prop = th + mc.normal(0, sc); p = neg_log_post(prop, x, y)
        if np.log(mc.uniform()) < -(p - cur):
            th, cur = prop, p
        if i > 10000 and i % 5 == 0:
            keep.append(th.copy())
    return np.array(keep)
postA = metropolis(xA, yA, mapA)
postB = metropolis(xB, yB, mapB)

# --- the story, asserted (robust, seed-independent margins) ---
assert eigA.min() > 20 * eigB.min()                 # confined spend -> near-flat ridge
assert abs(mapA[1] - 3.0) < abs(mapB[1] - 3.0)      # wide spend recovers K; confined does not
assert postB[:, 1].std() > postA[:, 1].std()        # posterior K wider under the wound
print(f"MAP wide   : beta={mapA[0]:.2f} K={mapA[1]:.2f} n={mapA[2]:.2f} lam={mapA[3]:.2f}")
print(f"MAP confined: beta={mapB[0]:.2f} K={mapB[1]:.2f} n={mapB[2]:.2f} lam={mapB[3]:.2f}")
print(f"smallest beta-K eig  wide={eigA.min():.1f}  confined={eigB.min():.2f}")
print(f"posterior K std  wide={postA[:,1].std():.2f}  confined={postB[:,1].std():.2f}")

# --- figures ---
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].scatter(postA[:, 0], postA[:, 1], s=4, alpha=0.3, label="wide spend")
ax[0].scatter(postB[:, 0], postB[:, 1], s=4, alpha=0.3, label="confined spend")
ax[0].axhline(3.0, ls="--", color="grey"); ax[0].axvline(4.0, ls="--", color="grey")
ax[0].set_xlabel(r"$\beta$"); ax[0].set_ylabel(r"$K$"); ax[0].set_title("posterior: blob vs ridge"); ax[0].legend()
ax[1].bar(["wide\nmin", "wide\nmax", "confined\nmin", "confined\nmax"],
          [eigA.min(), eigA.max(), eigB.min(), eigB.max()])
ax[1].set_yscale("log"); ax[1].set_title(r"$\beta$–$K$ Hessian eigenvalues")
plt.tight_layout()
plt.show()
```

- [ ] **Step 3: Run it headless.** Extract the cell to `/tmp/ch15_cell.py` and run `MPLBACKEND=Agg python3 /tmp/ch15_cell.py`. Expected: the four print lines, no assertion error, exit 0 (a benign Agg `plt.show()` warning is fine).

- [ ] **Step 4: Commit.**
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): code tie-in fitting the MMM and exhibiting the ridge"
```

---

### Task 7: Exercises (C / B / P / A)

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd` (replace `## Exercises` stub)

> **Controller note:** controller writes Exercises directly. Self-contained — no inline solution links. 2–3 per tier, using only this chapter's concepts.

- [ ] **Step 1: `### C — Conceptual / Reading Comprehension`.** e.g.: why does a good in-sample fit not imply the right response curve? In words, what does a singular Fisher information mean geometrically? Why is a non-negative prior on $\beta$ both domain knowledge and regularization?

- [ ] **Step 2: `### B — By Hand`.** e.g.: for $K=2, n=2, \beta=3$ compute the $(\beta,K)$ gradient at $x=2$ and at $x=4$ and show the two are not parallel (so two operating points identify $(\beta,K)$); show that all-equal spend makes the $2\times2$ Fisher determinant zero; given a flat ridge direction $d$, state the posterior along $d$.

- [ ] **Step 3: `### P — Prove It`.** e.g.: derive $I(\theta)=\tfrac1{\sigma^2}\sum_t\nabla\mu_t\nabla\mu_t^\top$ from the Gaussian log-likelihood; prove $I$ singular $\Rightarrow$ a likelihood-flat direction exists; prove that if $\ell$ is flat along $d$ then the posterior along $d$ equals the renormalized prior.

- [ ] **Step 4: `### A — Applied / Code`.** e.g.: rerun the fit at an intermediate spend spread and watch the smallest Hessian eigenvalue grow; add a prior-predictive check (draw parameters from the priors, simulate sales, plot); vary the prior on $K$ and show the confined-spend posterior tracks the prior (P2 empirically).

- [ ] **Step 5: Commit.**
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): exercises (C/B/P/A)"
```

---

### Task 8: Summary (auto-included)

**Files:** Modify `parts/05-mmm-modeling/02-building-fitting.qmd` (add `## Summary` after Exercises)

> **Controller note:** controller writes Summary directly. **Key concepts** and **Key identities** BOTH bulleted (inline math only) — never a run-on paragraph.

- [ ] **Step 1: Write `## Summary`** with a 1–2 sentence lead, then:
  - `**Key concepts**` bullets: the inverse problem (recover parameters, judged by recovery not fit); the MMM likelihood is non-conjugate ⇒ MAP/MCMC; weakly-informative priors as domain knowledge + regularizer; Fisher information and local identifiability (the spread of operating points, not the week count); the prior carries unidentified directions; the fitting workflow (prior predictive → fit → diagnostics → posterior predictive → recovery); residual non-identification motivates Part VI.
  - `**Key identities**` bullets: likelihood $y_t\sim\mathcal N(\mu_t(\theta),\sigma^2)$; mean $\mu_t=\text{baseline}_t+\sum_i\beta_i g(a_\lambda(x_{i,t});K_i,n_i)$; Fisher $I(\theta)=\tfrac1{\sigma^2}\sum_t\nabla\mu_t\nabla\mu_t^\top$; identifiability $\Leftrightarrow I$ nonsingular; ridge gradient $\nabla\mu=(g,\ \beta g_K)$ with $g_K=-nK^{n-1}x^n/(K^n+x^n)^2$; prior-dominates-ridge $p(\theta_0+sd\mid y)\propto p(\theta_0+sd)$.
  - Closing sentence handing to Ch16 (letting the coefficients move over time — dynamic models).

- [ ] **Step 2: Verify (bulleted, even `$`) & commit.**
```
git add parts/05-mmm-modeling/02-building-fitting.qmd
git commit -m "feat(ch15): summary"
```

---

### Task 9: Appendix solutions block

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Find the insertion point.** The file is one `::: {.content-visible when-meta="show-solutions"}` wrapper; its final closing `:::` is the last line. Insert the new block immediately before that final `:::`, after the Chapter 14 block, matching the existing per-chapter `## Chapter N — Title` structure.

- [ ] **Step 2: Write `## Chapter 15 — Building & Fitting an MMM`** with full solutions for every C/B/P/A exercise from Task 7. For B, show arithmetic: $K=2,n=2,\beta=3$ gradient at $x=2$: $g(2;2,2)=0.5$, $g_K(2;2,2)=-2\cdot2\cdot4/8^2=-0.25$, $\nabla\mu=(0.5,-0.75)$; at $x=4$: $g(4;2,2)=16/20=0.8$, $g_K(4;2,2)=-2\cdot2\cdot16/20^2=-0.16$, $\nabla\mu=(0.8,-0.48)$; ratios $0.5/0.8\neq -0.75/-0.48$, not parallel. For P, give the derivations from Tasks 3–4. For A, give a short headless-verified `{python}` cell (run it before committing). Provide the prose for the ridge-tracks-prior demonstration.

- [ ] **Step 3: Verify.** Balanced `:::` for the new content (no new fence needed — it lives inside the existing wrapper); even `$`-count on added lines (`awk` scan from the `## Chapter 15` line to end → no odd lines); any A code cell runs headless. Commit.
```
git add appendix/solutions.qmd
git commit -m "feat(ch15): appendix solutions"
```

---

### Task 10: Final review, headless re-run, PR

**Files:** all of the above.

- [ ] **Step 1: Structure/KaTeX lint** on `parts/05-mmm-modeling/02-building-fitting.qmd`: H1 present; the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises`, `## Summary`); even `$`-count; zero `begin{align}`; zero `psmallmatrix`; zero PPL/MMM library names.
```bash
f=parts/05-mmm-modeling/02-building-fitting.qmd
grep -c '\$' "$f"; grep -c 'begin{align}' "$f"; grep -c 'psmallmatrix' "$f"; grep -ciE 'pymc|stan|orbit|lightweight_mmm' "$f"
grep -nE '^#( |$)|^## ' "$f"
awk '{n=gsub(/\$/,"$"); if(n%2) print "ODD line "NR}' "$f"
```
Expected: even; 0; 0; 0; H1 then the six `##` headings in order (code-comment `#` lines inside the fence are fine); no odd lines.

- [ ] **Step 2: Re-run the chapter code cell headless** under `MPLBACKEND=Agg python3` — exit 0, all asserts pass.

- [ ] **Step 3: Push and open PR.**
```bash
git push -u origin worktree-ch15-building-fitting
gh pr create --base main --head worktree-ch15-building-fitting \
  --title "feat: Chapter 15 — Building & Fitting an MMM" \
  --body "Part V fitting chapter. Inverts Ch14's DGP by Bayesian inference: non-conjugacy rung, Fisher-information local-identifiability keystone (answering Ch14's 'how much variation is enough'), and prior-dominates-the-ridge. Fits via scipy MAP + Laplace + a hand-rolled Metropolis (no PPL named); the Hessian's near-zero eigenvalue exhibits the identification ridge; parameter recovery vs known ground truth. CI quarto render (HTML+PDF) is the gate."
```

- [ ] **Step 4: Watch CI to green.** Watch BOTH the HTML and PDF render jobs of `render.yml` on the PR to a green conclusion before reporting done. The `deploy` job is skipped on PR branches. If the PDF stage fails on a LaTeX/KaTeX gap, fix and re-push (HTML-side lint cannot catch PDF-only gaps). Then hand off to the user to merge.

---

## Self-Review (controller, before dispatching Task 1)

**Spec coverage:** Motivation (T1) ✓; likelihood + non-conjugacy rung (T2) ✓; priors (T2) ✓; posterior geometry/sampling (T3) ✓; **P1** Fisher identifiability keystone (T3) ✓; fitting workflow (T4) ✓; **P2** prior-dominates-ridge (T4) ✓; identification in practice (T4) ✓; WE1 Fisher / WE2 recovery / WE3 ridge (T5) ✓; two-scenario code with MAP+Laplace+Metropolis, ridge via eigenvalue (T6) ✓; C/B/P/A exercises (T7) ✓; bulleted Summary (T8) ✓; gated appendix solutions (T9) ✓; anchor line added, no library named (T1, T10) ✓; CI HTML+PDF gate (T10) ✓.

**Placeholder scan:** no TBD/TODO; all asserted numbers NumPy-verified; the code cell is complete and was run headless; recovery point-values are explicitly illustrative with robust inequality asserts (not fragile equality asserts).

**Consistency:** notation fixed once and reused ($\theta,\mu_t,g,a_\lambda,I,g_K$); `geometric_adstock`/`hill`/`hill_K`/`neg_log_post`/`hess_betaK`/`metropolis` defined once in T6 and referenced consistently; the gradient $\nabla\mu=(g,\beta g_K)$ and $g_K$ formula match between P1, WE1, T6, and the Summary.
