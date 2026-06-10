# Part II, Ch 1 — Linear Regression Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stub `parts/02-regression-bayes/01-linear-regression.qmd` with a complete chapter on OLS inference, diagnostics, and ridge regression — with four proofs (OLS unbiasedness, Gauss–Markov, bias–variance, ridge = Gaussian-prior MAP), worked examples, a runnable NumPy/SciPy tie-in, C/B/P/A exercises with appendix solutions, and a closing Summary.

**Architecture:** A single Quarto `.qmd` authored section by section against the fixed template, continuing the Ch 3 model `y = Xβ + ε` and the media-mix collinearity through-line. Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions`. Verification is build-based, not unit-test-based.

**Tech Stack:** Quarto (book project), KaTeX math, Python 3.11 with NumPy / SciPy / matplotlib (pinned in `requirements.txt`), `references.bib` (keys `@hastie2009`, `@casella2002`, `@strang2016`).

---

## Verification model (read first)

Prose + math + one code cell — no pytest suite. Each task is verified by:

- **PY** — extract the chapter's Python into `/tmp/ch04_codecell.py` and run
  `MPLBACKEND=Agg python3 /tmp/ch04_codecell.py` → clean exit + described output.
- **MATH** — `grep -c '\$\$' parts/02-regression-bayes/01-linear-regression.qmd`
  returns an **even** number; no bare `\begin{align}` (use `aligned` inside `$$`);
  `$$` delimiters on their own lines.
- **RENDER** — full `quarto render` runs in CI (`.github/workflows/render.yml`) on
  the PR; Quarto is **not installed locally**.

Commit after every task. Work ONLY in the worktree (below).

## Worktree (all work happens here)

`/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch4-linear-regression`
Branch `ch4-linear-regression` (own HEAD). Run all git/file commands here. Do NOT
touch `/Users/jameshenson/Documents/mmm_foundational_text`. Do NOT create/switch
branches. If git complains about author identity, run once:
`git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"`.

## File structure

- **Modify:** `parts/02-regression-bayes/01-linear-regression.qmd` — the whole
  chapter (stub now; replace section bodies, preserving H2 headings; ADD `## Summary`
  after `## Exercises`).
- **Modify:** `appendix/solutions.qmd` — append a "Part II · Linear Regression"
  solutions section inside the existing `.content-visible when-meta="show-solutions"`
  block (after the last existing chapter section, before the closing `:::`).
- **Reference (do not edit):** `references.bib`, `index.qmd`, the Ch 1/3 files.

Heading skeleton: `# Linear Regression` → `## Motivation` → `## Theory & Proofs` →
`## Worked Examples` → `## Code Tie-in` → `## Exercises` (`### C…/B…/P…/A…`) →
**`## Summary`**.

## Notation conventions (use consistently)

- Model `$y = X\beta + \varepsilon$`, `$\varepsilon\sim\mathcal N(0,\sigma^2 I)$`,
  `$X\in\mathbb{R}^{n\times p}$`, OLS `$\hat\beta=(X^\top X)^{-1}X^\top y$`.
- Residual sum of squares `$\text{RSS}=\lVert y-X\hat\beta\rVert^2$`, total
  `$\text{TSS}=\lVert y-\bar y\rVert^2$`, estimator `$\hat\sigma^2=\text{RSS}/(n-p)$`.
- Standard error `$\operatorname{se}(\hat\beta_j)=\hat\sigma\sqrt{[(X^\top X)^{-1}]_{jj}}$`;
  hat matrix `$H=X(X^\top X)^{-1}X^\top$`; VIF `$\text{VIF}_j=1/(1-R_j^2)$`.
- Ridge `$\hat\beta_{\text{ridge}}=(X^\top X+\lambda I)^{-1}X^\top y$`; Gaussian prior
  `$\beta\sim\mathcal N(0,\tau^2 I)$`, `$\lambda=\sigma^2/\tau^2$`.
- Reuse Ch 1/3 symbols: `$\kappa(X)$`, `$\sigma^2(X^\top X)^{-1}$`, the likelihood `$\ell$`.

---

### Task 1: Motivation

**Files:** Modify `parts/02-regression-bayes/01-linear-regression.qmd` (remove the
`.callout-note` "Stub" block; replace the `## Motivation` body).

- [ ] **Step 1: Remove stub callout, write Motivation (~250–350 words).**
Delete the `::: {.callout-note ...} **Stub.** ... :::` block (keep `# Linear
Regression` and the `*Canonical anchors: ...*` line). Under `## Motivation`:
- Recall Ch 3: we have `$\hat\beta$` and its spread `$\sigma^2(X^\top X)^{-1}$`, but
  an analyst must *act* — report a number with error bars, judge whether a channel
  "worked," and cope when TV and search nearly coincide.
- Driving question, verbatim emphasis: *"You fit the model and TV's coefficient is
  2.3. Is that real, how sure are you, and why did the number swing wildly when you
  added search?"*
- Preview the arc: OLS is optimal (Gauss–Markov) and quantifiable (standard errors,
  confidence intervals); collinearity breaks the quantification; the bias–variance
  tradeoff and ridge regression fix it; and ridge turns out to be a Gaussian prior in
  disguise — the bridge to the Bayesian chapter. Cite `[@hastie2009]`.

- [ ] **Step 2: MATH check** — `grep -c '\$\$' …01-linear-regression.qmd` even;
`grep -n 'Stub' …` nothing.

- [ ] **Step 3: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): write Motivation, drop stub callout"
```

---

### Task 2: Theory rungs 1–2 (sampling distribution + unbiasedness; Gauss–Markov)

**Files:** Modify the chapter (`## Theory & Proofs`).

- [ ] **Step 1: Rung 1 — OLS sampling distribution & unbiasedness [PROOF] (~300 words).**
`###` heading e.g. "The OLS estimator and its sampling distribution". Recall
`$\hat\beta=(X^\top X)^{-1}X^\top y$`. Substitute the model to get the key decomposition
$$
\hat\beta = (X^\top X)^{-1}X^\top(X\beta+\varepsilon) = \beta + (X^\top X)^{-1}X^\top\varepsilon .
$$
**Prove unbiasedness:** $\mathbb{E}[\hat\beta]=\beta+(X^\top X)^{-1}X^\top\mathbb{E}[\varepsilon]=\beta$.
Then the covariance (rung-1 identity of Ch 3): $\operatorname{Cov}(\hat\beta)=(X^\top X)^{-1}X^\top(\sigma^2 I)X(X^\top X)^{-1}=\sigma^2(X^\top X)^{-1}$.
Conclude that, since `$\hat\beta$` is an affine function of the Gaussian `$\varepsilon$`,
$$
\hat\beta \sim \mathcal N\big(\beta,\ \sigma^2(X^\top X)^{-1}\big),
$$
the finite-sample sampling distribution that all inference rests on. End with `$\blacksquare$`.

- [ ] **Step 2: Rung 2 — the Gauss–Markov theorem [PROOF] (~340 words).**
`###` heading e.g. "Gauss–Markov: OLS is BLUE". State:
> **Theorem (Gauss–Markov).** Assume $\mathbb{E}[\varepsilon]=0$ and $\operatorname{Cov}(\varepsilon)=\sigma^2 I$ (no normality needed). Among all linear unbiased estimators $\tilde\beta=Cy$, the OLS estimator has the smallest covariance: $\operatorname{Cov}(\tilde\beta)-\operatorname{Cov}(\hat\beta)\succeq 0$.

Proof: any linear estimator is $\tilde\beta=Cy$; unbiasedness $\mathbb{E}[Cy]=CX\beta=\beta$
for all $\beta$ forces $CX=I$. Write $C=(X^\top X)^{-1}X^\top + D$; then $CX=I$ implies
$DX=0$. Compute
$$
\operatorname{Cov}(\tilde\beta)=\sigma^2 CC^\top=\sigma^2\big[(X^\top X)^{-1}+DD^\top\big],
$$
because the cross terms $(X^\top X)^{-1}X^\top D^\top$ and its transpose vanish (since
$DX=0\Rightarrow X^\top D^\top=0$). Hence
$\operatorname{Cov}(\tilde\beta)-\operatorname{Cov}(\hat\beta)=\sigma^2 DD^\top\succeq0$,
with equality iff $D=0$, i.e. $\tilde\beta=\hat\beta$. End with `$\blacksquare$`. One
sentence: BLUE is why OLS is the default — but "best *unbiased*" leaves room to beat
it by allowing bias, which rungs 6–7 exploit.

- [ ] **Step 3: MATH check** — `$$` count even.

- [ ] **Step 4: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): theory rungs 1-2, sampling distribution + Gauss-Markov"
```

---

### Task 3: Theory rungs 3–4 (estimating σ²; inference)

**Files:** Modify the chapter (`## Theory & Proofs`).

- [ ] **Step 1: Rung 3 — estimating the noise variance (~230 words).**
`###` heading e.g. "Estimating the noise level". Define the residual vector
$\hat\varepsilon=y-X\hat\beta$ and $\text{RSS}=\lVert\hat\varepsilon\rVert^2$. State the
unbiased estimator
$$
\hat\sigma^2 = \frac{\text{RSS}}{n-p},
$$
and give the intuition for the $n-p$ **degrees of freedom**: the residual is confined
to the $(n-p)$-dimensional space orthogonal to the column space (the fitted values
use up $p$ directions), so dividing by $n-p$ rather than $n$ corrects the downward
bias. State (cite, no proof) that $\text{RSS}/\sigma^2\sim\chi^2_{n-p}$ and is
independent of $\hat\beta$ `[@hastie2009; @casella2002]` — the fact that makes the
next rung's $t$-statistic exact.

- [ ] **Step 2: Rung 4 — inference: standard errors, t-tests, confidence intervals (~320 words).**
`###` heading e.g. "Standard errors, t-tests, and confidence intervals". From rung 1,
$\hat\beta_j\sim\mathcal N(\beta_j,\ \sigma^2[(X^\top X)^{-1}]_{jj})$. Define the
**standard error** $\operatorname{se}(\hat\beta_j)=\hat\sigma\sqrt{[(X^\top X)^{-1}]_{jj}}$.
Because $\hat\sigma^2$ replaces $\sigma^2$ (rung 3) and $\text{RSS}/\sigma^2\sim\chi^2_{n-p}$
independent of $\hat\beta$, the standardized quantity follows a Student $t$:
$$
\frac{\hat\beta_j-\beta_j}{\operatorname{se}(\hat\beta_j)}\sim t_{n-p}.
$$
Derive the **confidence interval** by inverting this pivot:
$$
\hat\beta_j \pm t_{n-p,\,1-\alpha/2}\,\operatorname{se}(\hat\beta_j),
$$
and explain the $t$-test of $H_0:\beta_j=0$ (is the channel's effect distinguishable
from noise?). One sentence on the **F-test** for joint significance of several
coefficients (cite `[@hastie2009]`). Tie to the Motivation: the CI is the precise
answer to "how sure are you."

- [ ] **Step 3: MATH check** — `$$` count even.

- [ ] **Step 4: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): theory rungs 3-4, sigma-hat + inference"
```

---

### Task 4: Theory rungs 5–6 (diagnostics; multicollinearity & bias–variance)

**Files:** Modify the chapter (`## Theory & Proofs`).

- [ ] **Step 1: Rung 5 — diagnostics and goodness of fit (~300 words).**
`###` heading e.g. "Fitted values, residuals, and goodness of fit". Define the **hat
matrix** $H=X(X^\top X)^{-1}X^\top$ so $\hat y=Hy$; note it is the orthogonal
projection onto the column space (Ch 1) and is symmetric and idempotent
($H^2=H$, $H^\top=H$) — state this, one line of justification. Define $R^2 = 1 -
\text{RSS}/\text{TSS}$ with $\text{TSS}=\lVert y-\bar y\rVert^2$, the fraction of
variance explained. List the **four assumptions** — linearity, independence of
errors, homoscedasticity ($\operatorname{Cov}(\varepsilon)=\sigma^2 I$), normality —
and what each buys: linearity+zero-mean for unbiasedness, homoscedasticity for
Gauss–Markov and the simple SE formula, normality for exact $t$/$F$ inference. One
line: residual plots are how you check them.

- [ ] **Step 2: Rung 6 — multicollinearity and the bias–variance tradeoff [PROOF] (~340 words).**
`###` heading e.g. "Multicollinearity and the bias–variance tradeoff". Explain that
$\operatorname{Cov}(\hat\beta)=\sigma^2(X^\top X)^{-1}$ blows up when $X^\top X$ has
small eigenvalues — collinear channels, large $\kappa(X)$ (Ch 1). Define the
**variance inflation factor** $\text{VIF}_j = 1/(1-R_j^2)$, where $R_j^2$ is from
regressing predictor $j$ on the others; a VIF of 100 means that coefficient's
variance is 100× what it would be under orthogonal predictors. Then state and **prove**
the **bias–variance decomposition** for any estimator $\hat\theta$ of a scalar
$\theta$:
> **Proposition.** $\operatorname{MSE}(\hat\theta)=\mathbb{E}[(\hat\theta-\theta)^2]=\operatorname{Var}(\hat\theta)+\big(\mathbb{E}[\hat\theta]-\theta\big)^2$ (variance + bias squared).

Proof: add and subtract $\mathbb{E}[\hat\theta]$ inside the square,
$$
\mathbb{E}[(\hat\theta-\theta)^2]=\mathbb{E}\big[(\hat\theta-\mathbb{E}\hat\theta)^2\big]+2(\mathbb{E}\hat\theta-\theta)\,\mathbb{E}[\hat\theta-\mathbb{E}\hat\theta]+(\mathbb{E}\hat\theta-\theta)^2,
$$
and the cross term vanishes because $\mathbb{E}[\hat\theta-\mathbb{E}\hat\theta]=0$,
leaving variance + bias$^2$. End with `$\blacksquare$`. Punchline: Gauss–Markov makes
OLS best among *unbiased* estimators, but a biased estimator with much smaller
variance can have lower MSE — which is exactly what ridge will do.

- [ ] **Step 3: MATH check** — `$$` count even.

- [ ] **Step 4: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): theory rungs 5-6, diagnostics + bias-variance proof"
```

---

### Task 5: Theory rung 7 (ridge regression; ridge = Gaussian-prior MAP) [KEYSTONE]

**Files:** Modify the chapter (`## Theory & Proofs`).

- [ ] **Step 1: Rung 7 — ridge regression and its Bayesian reading [PROOF] (~360 words).**
`###` heading e.g. "Ridge regression: a penalty that is a prior". Define ridge as the
penalized least-squares problem
$$
\hat\beta_{\text{ridge}} = \arg\min_\beta\ \lVert y-X\beta\rVert^2 + \lambda\lVert\beta\rVert^2,
$$
and derive its closed form by setting the gradient to zero:
$$
-2X^\top(y-X\beta)+2\lambda\beta=0 \ \Longrightarrow\ \hat\beta_{\text{ridge}}=(X^\top X+\lambda I)^{-1}X^\top y .
$$
Explain the conditioning fix: adding $\lambda I$ lifts every eigenvalue of $X^\top X$
by $\lambda$, so the inverse is well-behaved even when $X^\top X$ is near-singular —
ridge directly cures the collinearity disease, at the cost of shrinking coefficients
toward zero (bias). State and **prove** the keystone bridge:
> **Theorem (ridge = Gaussian-prior MAP).** Under the model $y\sim\mathcal N(X\beta,\sigma^2 I)$ with prior $\beta\sim\mathcal N(0,\tau^2 I)$, the posterior mode (MAP estimate) of $\beta$ is the ridge estimator with $\lambda=\sigma^2/\tau^2$.

Proof: by Bayes (Ch 3), $p(\beta\mid y)\propto p(y\mid\beta)\,p(\beta)$, so the log
posterior is
$$
\log p(\beta\mid y) = -\frac{1}{2\sigma^2}\lVert y-X\beta\rVert^2 - \frac{1}{2\tau^2}\lVert\beta\rVert^2 + \text{const}.
$$
Maximizing it equals minimizing $\lVert y-X\beta\rVert^2 + \frac{\sigma^2}{\tau^2}\lVert\beta\rVert^2$,
which is the ridge objective with $\lambda=\sigma^2/\tau^2$; its minimizer is
$(X^\top X+\lambda I)^{-1}X^\top y$. End with `$\blacksquare$`. Close: a penalty *is* a
prior — a strong prior ($\tau^2$ small) means heavy shrinkage ($\lambda$ large). One
line naming **lasso** ($\ell_1$ penalty / Laplace prior) and **elastic-net**, and
pointing to the Bayesian chapter, which replaces this single point estimate with the
full posterior.

- [ ] **Step 2: MATH check** — `$$` count even.

- [ ] **Step 3: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): theory rung 7, ridge regression = Gaussian-prior MAP"
```

---

### Task 6: Worked Examples

**Files:** Modify the chapter (`## Worked Examples`). VERIFY arithmetic in Python before committing.

- [ ] **Step 1: Example (a) — a standard error and a 95% CI by hand.**
`###` heading. Use a small design (intercept + one channel), e.g.
$$
X=\begin{bmatrix}1&1\\1&2\\1&3\\1&4\end{bmatrix},\qquad y=\begin{bmatrix}2\\2\\4\\5\end{bmatrix}.
$$
Compute $X^\top X$, $(X^\top X)^{-1}$, $\hat\beta$, $\text{RSS}$, $\hat\sigma^2=\text{RSS}/(n-p)$
with $n=4,p=2$, then $\operatorname{se}(\hat\beta_1)$ and the 95% CI
$\hat\beta_1\pm t_{2,0.975}\operatorname{se}(\hat\beta_1)$ (note $t_{2,0.975}\approx4.303$).
VERIFY every number in Python (below) and match the prose.

- [ ] **Step 2: Example (b) — collinearity inflates the SE; ridge stabilizes.**
`###` heading. Construct an intercept + two near-duplicate channels (e.g. second and
third columns differing by a small perturbation). Compute and report: the OLS
coefficients with their large standard errors, the VIF (or $\kappa(X)$), then the
ridge estimate for a modest $\lambda$ showing the coefficients pulled back to sane,
stable values. Use Python-computed numbers; round and state rounding.

- [ ] **Step 3: PY check** — write `/tmp/ch04_examples_check.py` reproducing (a)
($\hat\beta$, $\hat\sigma^2$, se, CI) and (b) (OLS se/VIF vs ridge coefficients);
run `python3 /tmp/ch04_examples_check.py`; make prose match; paste key values into the report.

- [ ] **Step 4: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): worked examples (SE/CI by hand; collinearity vs ridge)"
```

---

### Task 7: Code Tie-in (runnable NumPy/SciPy/matplotlib cell)

**Files:** Modify the chapter (`## Code Tie-in`).

- [ ] **Step 1: Write a lead-in + this Quarto Python cell.**

````markdown
```{python}
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

rng = np.random.default_rng(1)
n, sigma = 80, 1.0

def make_data(corr):
    """Synthetic media data: intercept, TV, and a search channel corr-correlated with TV."""
    tv = rng.normal(size=n)
    search = corr * tv + np.sqrt(1 - corr**2) * rng.normal(size=n)
    X = np.column_stack([np.ones(n), tv, search])
    beta = np.array([5.0, 2.0, 1.0])
    y = X @ beta + sigma * rng.normal(size=n)
    return X, y

def ols_fit(X, y):
    n, p = X.shape
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    sigma2 = resid @ resid / (n - p)
    se = np.sqrt(np.diag(sigma2 * XtX_inv))
    return beta, se, n - p

def vif(X, j):
    """Variance inflation factor for column j: regress it on the other columns."""
    xj = X[:, j]
    others = np.delete(X, j, axis=1)
    coef, *_ = np.linalg.lstsq(others, xj, rcond=None)
    resid = xj - others @ coef
    r2 = 1 - (resid @ resid) / ((xj - xj.mean()) @ (xj - xj.mean()))
    return 1.0 / (1.0 - r2)

# OLS coefficients with standard errors and 95% CIs (mild collinearity)
X, y = make_data(corr=0.3)
beta, se, df = ols_fit(X, y)
tcrit = stats.t.ppf(0.975, df)
for j, name in enumerate(["intercept", "tv", "search"]):
    lo, hi = beta[j] - tcrit * se[j], beta[j] + tcrit * se[j]
    print(f"{name:9s} beta={beta[j]:6.3f}  se={se[j]:.3f}  95% CI=[{lo:.3f}, {hi:.3f}]")

# Standard error and VIF inflation as the two channels become collinear
print("\n corr   se(tv)   VIF(tv)   cond(X)")
for corr in [0.0, 0.9, 0.99, 0.999]:
    X, y = make_data(corr)
    _, se, _ = ols_fit(X, y)
    print(f"{corr:5.3f}  {se[1]:7.3f}  {vif(X, 1):8.1f}  {np.linalg.cond(X):8.1f}")

# Ridge coefficient paths: shrink and stabilize the collinear fit
X, y = make_data(corr=0.99)
p = X.shape[1]
lambdas = np.logspace(-3, 3, 60)
paths = np.array([np.linalg.solve(X.T @ X + lam * np.eye(p), X.T @ y) for lam in lambdas])
beta_ols, _, _ = ols_fit(X, y)
print("\nridge -> OLS as lambda -> 0:", np.allclose(paths[0], beta_ols, atol=1e-2))

plt.semilogx(lambdas, paths[:, 1], label="tv")
plt.semilogx(lambdas, paths[:, 2], label="search")
plt.axhline(0, color="0.7", lw=0.8)
plt.xlabel(r"$\lambda$"); plt.ylabel("ridge coefficient"); plt.legend()
plt.title("Ridge coefficient paths (collinear channels)")
plt.show()
```
````

After the cell, write 2–3 sentences reading the output: at mild correlation the CIs
are tight; as `corr → 1` the TV standard error and VIF explode while `cond(X)` grows
(Ch 1 conditioning); ridge pulls the wild collinear coefficients back toward each
other and toward zero as `λ` grows, and reduces to OLS as `λ → 0`.

- [ ] **Step 2: PY check** — copy the cell body to `/tmp/ch04_codecell.py`; run
`MPLBACKEND=Agg python3 /tmp/ch04_codecell.py`; expect exit 0; the CI block prints
finite intervals; the inflation table shows `se(tv)`, `VIF(tv)`, and `cond(X)` all
rising monotonically with `corr`; and `ridge -> OLS as lambda -> 0: True`. Paste exact
stdout into the report.

- [ ] **Step 3: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): runnable code tie-in (SE/CI, VIF inflation, ridge paths)"
```

---

### Task 8: Exercises (C / B / P / A)

**Files:** Modify the chapter (`## Exercises`).

- [ ] **Step 1: Populate the four tiers (2–3 each).**
- **C — Conceptual:** (C1) What does a 95% confidence interval for $\beta_{\text{TV}}$
  claim, and what does it *not* claim? (C2) Explain why two collinear channels produce
  huge standard errors even when the overall fit ($R^2$) is excellent. (C3) Ridge
  introduces bias on purpose — what does it buy, and in what sense is choosing $\lambda$
  the same as choosing a prior?
- **B — By hand:** (B1) For $X=\begin{bmatrix}1&1\\1&2\\1&3\\1&4\end{bmatrix}$,
  $y=(2,2,4,5)^\top$, compute $\hat\beta$, $\hat\sigma^2$, and $\operatorname{se}(\hat\beta_1)$.
  (B2) Two predictors with sample correlation $r=0.95$: compute the VIF.
  (B3) For $X^\top X=\begin{bmatrix}4&2\\2&4\end{bmatrix}$, $X^\top y=(6,4)^\top$,
  compute $\hat\beta_{\text{ridge}}$ at $\lambda=2$.
- **P — Prove it:** (P1) Prove $\hat\beta$ is unbiased, $\mathbb{E}[\hat\beta]=\beta$.
  (P2) Prove the Gauss–Markov theorem (OLS is BLUE). (P3) Prove that the ridge
  estimator is the MAP estimate under $y\sim\mathcal N(X\beta,\sigma^2 I)$,
  $\beta\sim\mathcal N(0,\tau^2 I)$, identifying $\lambda=\sigma^2/\tau^2$.
- **A — Applied / code:** (A1) Plot OLS-coefficient sampling variability: simulate many
  datasets at correlation $\rho\in\{0,0.9,0.99\}$ and show the spread of $\hat\beta_{\text{search}}$
  widening; overlay the theoretical $\operatorname{se}$. (A2) Bias–variance simulation:
  for a fixed collinear design, compare the MSE of OLS vs ridge (over a $\lambda$ grid)
  in recovering the true $\beta$, and identify a $\lambda$ where ridge beats OLS.

No inline solution links.

- [ ] **Step 2: MATH check** — `$$` count even.

- [ ] **Step 3: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): exercises across C/B/P/A tiers"
```

---

### Task 9: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`. VERIFY numeric answers (B1, B2, B3) in Python.

- [ ] **Step 1: Inspect existing structure.** Read `appendix/solutions.qmd`: it uses
`::: {.content-visible when-meta="show-solutions"}` with per-chapter `##` sections
(Ch 1–3 already present). Add a `## Part II · Linear Regression` section INSIDE that
visible block, AFTER the last existing chapter section and BEFORE the closing `:::`,
with `###` tier subheadings (C/B/P/A). Mirror the existing idiom exactly; do not alter
the gating divs or other chapters' solutions.

- [ ] **Step 2: Write solutions.**
- C1–C3: 2–4 sentence model answers.
- B1: $\hat\beta=(X^\top X)^{-1}X^\top y$, $\hat\sigma^2=\text{RSS}/2$, and
  $\operatorname{se}(\hat\beta_1)=\hat\sigma\sqrt{[(X^\top X)^{-1}]_{22}}$ — compute
  explicitly (VERIFY in Python).
- B2: $\text{VIF}=1/(1-0.95^2)=1/(1-0.9025)\approx10.26$ (VERIFY).
- B3: $\hat\beta_{\text{ridge}}=(X^\top X+2I)^{-1}X^\top y=\begin{bmatrix}6&2\\2&6\end{bmatrix}^{-1}\begin{bmatrix}6\\4\end{bmatrix}$
  — compute explicitly (VERIFY; $=(0.875,\ 0.375)$ after checking).
- P1–P3: full self-contained proofs (mirror rungs 1, 2, 7).
- A1: expected result (empirical spread of $\hat\beta_{\text{search}}$ widens, tracking
  $\sigma\sqrt{[(X^\top X)^{-1}]_{jj}}$); short snippet sketch. A2: a $\lambda>0$ exists
  where ridge MSE < OLS MSE under collinearity; snippet sketch + the qualitative U-shape.

- [ ] **Step 3: PY check** — `/tmp/ch04_sol_check.py` computes B1 (se), B2 (VIF), B3
(ridge); run `python3 /tmp/ch04_sol_check.py`; make prose match; paste values into the report.

- [ ] **Step 4: Commit**
```
git add appendix/solutions.qmd
git commit -m "feat(ch04): appendix solutions for B/P/A exercises"
```

---

### Task 10: Summary

**Files:** Modify the chapter (ADD `## Summary` after `## Exercises`).

- [ ] **Step 1: Write the Summary** (Ch 1/3 format: **Key concepts** bullets +
**Key identities** bullets, inline math only to keep display-math parity). Cover:
- Concepts: OLS sampling distribution and unbiasedness; Gauss–Markov/BLUE; $\hat\sigma^2$
  and $n-p$ degrees of freedom; standard errors → $t$ → confidence intervals; the hat
  matrix and $R^2$; collinearity inflates variance (VIF, $\kappa$); the bias–variance
  tradeoff; ridge re-conditions $X^\top X$ and is a Gaussian prior in disguise.
- Identities (inline): $\hat\beta=\beta+(X^\top X)^{-1}X^\top\varepsilon$;
  $\hat\beta\sim\mathcal N(\beta,\sigma^2(X^\top X)^{-1})$;
  $\hat\sigma^2=\text{RSS}/(n-p)$;
  $\operatorname{se}(\hat\beta_j)=\hat\sigma\sqrt{[(X^\top X)^{-1}]_{jj}}$;
  CI $\hat\beta_j\pm t_{n-p,1-\alpha/2}\operatorname{se}(\hat\beta_j)$;
  $H=X(X^\top X)^{-1}X^\top$, $R^2=1-\text{RSS}/\text{TSS}$;
  $\text{VIF}_j=1/(1-R_j^2)$; $\operatorname{MSE}=\text{bias}^2+\text{variance}$;
  $\hat\beta_{\text{ridge}}=(X^\top X+\lambda I)^{-1}X^\top y$, $\lambda=\sigma^2/\tau^2$.

- [ ] **Step 2: MATH check** — `$$` count even; `## Summary` is the LAST `##` heading.

- [ ] **Step 3: Commit**
```
git add parts/02-regression-bayes/01-linear-regression.qmd
git commit -m "feat(ch04): Summary section"
```

---

### Task 11: Integration sweep

**Files:** Possibly modify the chapter and/or `appendix/solutions.qmd` (fixes only).

- [ ] **Step 1: Stub/notation/cross-ref sweep.** Read the whole chapter. Confirm:
notation matches the conventions table; cross-references (Ch 1 `$X^\top X$`/projection/
`$\kappa$`, Ch 3 sampling covariance/Bayes) resolve; the four proofs are present; all
four exercise tiers + Summary populated.
Run `grep -ni 'stub\|_why this matters\|_worked numerical\|_minimal, runnable\|_definitions ->\|TODO\|TBD' parts/02-regression-bayes/01-linear-regression.qmd` → expect no matches.

- [ ] **Step 2: Final MATH check** — `grep -c '\$\$' …` even; no bare `\begin{align}`.

- [ ] **Step 3: Final PY check** — re-extract and run the code cell
(`MPLBACKEND=Agg python3 /tmp/ch04_codecell.py`, exit 0).

- [ ] **Step 4: RENDER gate** — note in the PR that local Quarto is unavailable and CI
is the render authority.

- [ ] **Step 5: Commit any fixes**
```
git add -A
git commit -m "fix(ch04): integration sweep — notation, refs, cleanup"
```

---

## Self-Review (completed by plan author)

**Spec coverage:** Through-line (OLS → collinearity → ridge → Bayes) → Tasks 1–7.
Theory rungs 1–7 → Tasks 2–5. Four named proofs — unbiasedness (Task 2), Gauss–Markov
(Task 2), bias–variance (Task 4), ridge = Gaussian-prior MAP (Task 5) — all present.
Inference (SE/t/CI) → Task 3. Diagnostics/hat matrix/$R^2$/assumptions → Task 4.
VIF/$\kappa$ → Tasks 4 & 7. Worked examples (a)/(b) → Task 6. NumPy/SciPy tie-in with
SE/CI, VIF inflation, ridge paths → Task 7. C/B/P/A exercises → Task 8. Appendix
solutions gated by `show-solutions` → Task 9. Summary → Task 10. Render + notation
success criteria → Task 11. Distributional mechanics ($\chi^2_{n-p}$, exact $t$/$F$)
stated and cited, not proved → Tasks 3–4. Out-of-scope items never introduced.

**Placeholder scan:** No TBD/TODO/"handle edge cases"; the code cell is given in full;
proof contents stated explicitly, not deferred.

**Type/notation consistency:** `make_data`, `ols_fit`, `vif`, `beta`, `se`, `df`,
`tcrit`, `paths`, `lambdas`, `beta_ols` consistent across Tasks 6/7/9; math symbols
($\hat\beta$, $X^\top X$, $\hat\sigma^2$, $\operatorname{se}$, $H$, $R^2$, $\text{VIF}$,
$\hat\beta_{\text{ridge}}$, $\lambda=\sigma^2/\tau^2$) used uniformly per the
conventions table.
