# Chapter 3 — Probability & Statistics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stub `parts/01-foundations/03-probability-statistics.qmd` with a complete, MMM-targeted probability/statistics chapter (7-rung theory ladder with four proofs incl. Gaussian MLE = OLS, worked examples, a runnable NumPy/SciPy tie-in, C/B/P/A exercises with appendix solutions, and a closing Summary).

**Architecture:** A single Quarto `.qmd` authored section by section against the book's fixed template, anchored on the regression as a probabilistic model `y = Xβ + ε, ε ~ N(0, σ²I)`. Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions`. Verification is build-based, not unit-test-based.

**Tech Stack:** Quarto (book project), KaTeX math, Python 3.11 with NumPy / SciPy / matplotlib (pinned in `requirements.txt`), `references.bib` (keys `@casella2002`, `@wasserman2004`, `@strang2016`).

---

## Verification model (read first)

Prose + math + one code cell — no pytest suite. Each task is verified by:

- **PY** — extract the chapter's Python into `/tmp/ch03_codecell.py` and run
  `MPLBACKEND=Agg python3 /tmp/ch03_codecell.py` → expect clean exit + described output.
- **MATH** — every `$$` balanced: `grep -c '\$\$' parts/01-foundations/03-probability-statistics.qmd` returns an **even** number; no bare `\begin{align}` (use `aligned` inside `$$`); `$$` delimiters on their own lines.
- **RENDER** — full `quarto render` runs in CI (`.github/workflows/render.yml`) on the PR; Quarto is **not installed locally** (`brew install quarto` if a local check is wanted).

Commit after every task. Work ONLY in the worktree (below); never touch the main checkout.

## Worktree (all work happens here)

`/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch3-probability-statistics`
Branch `ch3-probability-statistics` (own HEAD). Run all git/file commands here. Do NOT touch `/Users/jameshenson/Documents/mmm_foundational_text`. Do NOT create/switch branches.

## File structure

- **Modify:** `parts/01-foundations/03-probability-statistics.qmd` — the whole chapter
  (currently a stub; replace section bodies, preserving the H2 headings; ADD a
  `## Summary` after `## Exercises`).
- **Modify:** `appendix/solutions.qmd` — append a "Chapter 3 — Probability &
  Statistics" section inside the existing `.content-visible when-meta="show-solutions"`
  block (after the Chapter 2 solutions, before the block's closing `:::`).
- **Reference (do not edit):** `references.bib`, `index.qmd`, the Ch 1/Ch 2 files
  (for voice + cross-references).

Heading skeleton to preserve/extend in the chapter file:
`# Probability & Statistics` → `## Motivation` → `## Theory & Proofs` →
`## Worked Examples` → `## Code Tie-in` → `## Exercises` (`### C…/B…/P…/A…`) →
**`## Summary`** (new).

## Notation conventions (use consistently)

- Model: `$y = X\beta + \varepsilon$`, `$\varepsilon \sim \mathcal{N}(0,\sigma^2 I)$`,
  `$X \in \mathbb{R}^{n\times p}$`, `$\hat\beta$` the estimate.
- Expectation `$\mathbb{E}[\cdot]$`, variance `$\operatorname{Var}(\cdot)$`,
  covariance matrix `$\operatorname{Cov}(\cdot)$`. Gaussian `$\mathcal{N}(\mu,\Sigma)$`.
- Likelihood `$\mathcal{L}(\beta,\sigma^2)$`, log-likelihood `$\ell$`, Fisher
  information `$\mathcal{I}$`. Cholesky `$\Sigma = LL^\top$` (from Ch 1).
- Reuse Ch 1/2 symbols where relevant: `$X^\top X$`, `$\kappa(X)$`, `$\nabla^2 L = 2X^\top X$`.

---

### Task 1: Motivation

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (remove the
`.callout-note` "Stub" block at top; replace the `## Motivation` body).

- [ ] **Step 1: Remove stub callout, write Motivation**

Delete the `::: {.callout-note ...} **Stub.** ... :::` block (keep `# Probability &
Statistics` and the `*Canonical anchors: ...*` line). Under `## Motivation`, write
~250–350 words that:
- Name the gap: Ch 1–2 fit `$\hat\beta$` as if `$y$` were exact, but weekly sales
  are noisy.
- Pose the driving question verbatim in emphasis: *"If you reran the same campaign,
  you'd get different sales — so how much of $\hat\beta$ is signal, and how much is
  luck?"*
- State the resolution: treat `$y$` as a draw from a distribution
  `$y = X\beta + \varepsilon$`, `$\varepsilon\sim\mathcal N(0,\sigma^2 I)$`; this
  reveals the minimized loss `$L$` of Ch 2 was a (negative log-) likelihood.
- Map the ladder (random vectors & covariance → distributions → conditioning →
  multivariate normal → likelihood/MLE → LLN/CLT → Bayes) and forward-point to
  Part II (posterior on `$\beta$`) and Part III (Monte Carlo). Cite
  `[@casella2002]`.

- [ ] **Step 2: MATH check** — `grep -c '\$\$' parts/01-foundations/03-probability-statistics.qmd` → even; `grep -n 'Stub' ...` → nothing.

- [ ] **Step 3: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): write Motivation, drop stub callout"
```

---

### Task 2: Theory rungs 1–3 (random vectors & covariance → distributions → conditioning)

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 1 — random variables, expectation, variance, covariance matrix [PROOF] (~300 words).**
`###` heading e.g. "Random variables, expectation, and covariance". Define a random
variable, `$\mathbb{E}$`, `$\operatorname{Var}$`; linearity of expectation. Define
the **covariance matrix** of a random vector `$w\in\mathbb{R}^m$`:
`$\operatorname{Cov}(w) = \mathbb{E}[(w-\mathbb{E}w)(w-\mathbb{E}w)^\top]$`
(symmetric PSD). State and **prove** the transformation identity:
> **Proposition.** For a fixed matrix `$A$` and vector `$b$`,
> `$\operatorname{Cov}(Aw + b) = A\,\operatorname{Cov}(w)\,A^\top$`.

Proof: let `$\mu=\mathbb{E}w$`; then `$\mathbb{E}[Aw+b]=A\mu+b$`, and
`$\operatorname{Cov}(Aw+b)=\mathbb{E}[A(w-\mu)(w-\mu)^\top A^\top]
=A\,\mathbb{E}[(w-\mu)(w-\mu)^\top]\,A^\top=A\operatorname{Cov}(w)A^\top$`. Note this
is the identity that justified Ch 1's Cholesky sampling (`$\operatorname{Cov}(Lz)=LL^\top$`).

- [ ] **Step 2: Rung 2 — key distributions, compact (~220 words).**
`###` heading e.g. "The distributions we will need". Present the **Gaussian**
`$\mathcal N(\mu,\sigma^2)$` density as the star (state the 1-D density). Then a
compact paragraph naming **Bernoulli/Binomial** (a converted / not-converted user,
counts of conversions), **Poisson** (event counts like store visits), and **Gamma**
(positive quantities; later a prior for variance/spend coefficients) — one line each
on where MMM meets them. State, don't belabor; cite `[@casella2002]`.

- [ ] **Step 3: Rung 3 — joint, marginal, conditional, independence (~280 words).**
`###` heading e.g. "Joint, conditional, and the regression model". Define joint /
marginal / conditional densities and **independence** (`$p(u,v)=p(u)p(v)$`).
Explain regression as conditioning: we model `$y \mid X,\beta$`. Establish, for the
model `$y = X\beta+\varepsilon$` with `$\varepsilon\sim\mathcal N(0,\sigma^2 I)$` and
independent rows:
$$
\mathbb{E}[y] = X\beta, \qquad \operatorname{Cov}(y) = \sigma^2 I,
$$
using linearity and the rung-1 identity (`$\operatorname{Cov}(X\beta+\varepsilon)=\operatorname{Cov}(\varepsilon)=\sigma^2 I$`).

- [ ] **Step 4: MATH check** — `$$` count even.

- [ ] **Step 5: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): theory rungs 1-3 incl. covariance-transform proof"
```

---

### Task 3: Theory rungs 4–5 (multivariate normal → likelihood & MLE)

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 4 — the multivariate normal [PROOF] (~320 words).**
`###` heading e.g. "The multivariate normal". State the density for
`$w\sim\mathcal N(\mu,\Sigma)$`, `$\Sigma\succ0$`:
$$
p(w) = (2\pi)^{-m/2}\,(\det\Sigma)^{-1/2}\,\exp\!\Big(-\tfrac12 (w-\mu)^\top \Sigma^{-1}(w-\mu)\Big).
$$
Explain the role of `$\Sigma^{-1}$` (the quadratic form = squared Mahalanobis
distance) and `$\det\Sigma$` (normalizing volume). State the **affine property** and
**prove** it at the mean/covariance level:
> **Proposition.** If `$z\sim\mathcal N(0,I)$` and `$w = \mu + Lz$`, then
> `$\mathbb{E}[w]=\mu$` and `$\operatorname{Cov}(w)=LL^\top$`; in particular with
> `$\Sigma=LL^\top$` (Cholesky, Ch 1), `$w\sim\mathcal N(\mu,\Sigma)$`.

Proof: `$\mathbb{E}[w]=\mu+L\,\mathbb{E}[z]=\mu$`; by the rung-1 identity
`$\operatorname{Cov}(w)=\operatorname{Cov}(Lz)=L\,\operatorname{Cov}(z)\,L^\top=LIL^\top=LL^\top$`.
That an affine image of a Gaussian is Gaussian (not merely matched in moments) is
standard; cite `[@casella2002]`. Tie back: this is exactly how Ch 1's Cholesky factor
manufactures correlated draws, and it gives the distribution of `$y$`:
`$y\sim\mathcal N(X\beta,\sigma^2 I)$`.

- [ ] **Step 2: Rung 5 — likelihood, log-likelihood, MLE = least squares [PROOF] (~340 words).**
`###` heading e.g. "Likelihood and the Gaussian MLE". Define the **likelihood**
`$\mathcal L(\beta,\sigma^2)=p(y\mid\beta,\sigma^2)$` and log-likelihood `$\ell$`.
For `$y\sim\mathcal N(X\beta,\sigma^2 I)$`:
$$
\ell(\beta,\sigma^2) = -\frac{n}{2}\log(2\pi\sigma^2) - \frac{1}{2\sigma^2}\lVert y - X\beta\rVert^2 .
$$
State and **prove** the keystone:
> **Theorem (Gaussian MLE = least squares).** For fixed `$\sigma^2>0$`, the maximizer
> of `$\ell$` over `$\beta$` is the least-squares solution: any `$\hat\beta$`
> maximizing `$\ell$` solves the normal equations `$X^\top X\hat\beta = X^\top y$`.

Proof: only the `$-\frac{1}{2\sigma^2}\lVert y-X\beta\rVert^2$` term depends on
`$\beta$`; maximizing `$\ell$` in `$\beta$` is therefore minimizing
`$\lVert y-X\beta\rVert^2$` — the exact Ch 1–2 least-squares objective — whose
minimizer solves the normal equations (cite the Ch 1 result / Ch 2 gradient
`$\nabla_\beta\ell = \frac{1}{\sigma^2}X^\top(y-X\beta)$`, set to 0). Then state the
**Fisher information** `$\mathcal I(\beta) = X^\top X/\sigma^2$` (the negative
expected Hessian of `$\ell$` in `$\beta$`), so the MLE's covariance scales like
`$\sigma^2 (X^\top X)^{-1}$` — collinearity (small eigenvalues of `$X^\top X$`,
Ch 1) means a high-variance, imprecise estimate. This is the probabilistic restatement
of Ch 1's conditioning story.

- [ ] **Step 3: MATH check** — `$$` count even.

- [ ] **Step 4: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): theory rungs 4-5, multivariate normal + Gaussian-MLE proof"
```

---

### Task 4: Theory rungs 6–7 (LLN & CLT → Bayes)

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 6 — LLN (proved) & CLT (cited) (~300 words).**
`###` heading e.g. "Law of large numbers and the central limit theorem". State
Markov/Chebyshev briefly, then state and **prove the weak LLN** via Chebyshev:
> **Theorem (Weak LLN).** For i.i.d. `$X_1,\dots,X_n$` with mean `$\mu$` and finite
> variance `$\sigma^2$`, the sample mean `$\bar X_n$` satisfies, for any `$\epsilon>0$`,
> `$\Pr(|\bar X_n - \mu| \ge \epsilon) \to 0$` as `$n\to\infty$`.

Proof: `$\mathbb{E}[\bar X_n]=\mu$` and `$\operatorname{Var}(\bar X_n)=\sigma^2/n$`;
Chebyshev gives `$\Pr(|\bar X_n-\mu|\ge\epsilon)\le \sigma^2/(n\epsilon^2)\to 0$`.
Then **state the CLT** (`$\sqrt n(\bar X_n-\mu)\xrightarrow{d}\mathcal N(0,\sigma^2)$`)
and **cite** `[@casella2002; @wasserman2004]` (no proof). Frame: the LLN licenses
Monte-Carlo estimates (averages of draws converge to expectations) and the CLT both
explains the Gaussian noise model and quantifies estimator uncertainty — the engine
of Part III.

- [ ] **Step 2: Rung 7 — Bayes' theorem (~220 words).**
`###` heading e.g. "Bayes' theorem". State the conditional-probability form
`$p(\beta\mid y) = \dfrac{p(y\mid\beta)\,p(\beta)}{p(y)}$` and the proportional form
$$
p(\beta \mid y) \;\propto\; p(y\mid\beta)\,p(\beta)
\;=\; \underbrace{\mathcal L(\beta)}_{\text{likelihood (rung 5)}}\;\times\;\underbrace{p(\beta)}_{\text{prior}} .
$$
Contrast with MLE: where MLE returns a single `$\hat\beta$`, Bayes returns a whole
**posterior** distribution over `$\beta$`, directly answering the Motivation's
"signal vs luck". One line: computing this posterior — conjugacy, then MCMC — is
Parts II and III. Stated, not proved.

- [ ] **Step 3: MATH check** — `$$` count even.

- [ ] **Step 4: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): theory rungs 6-7, weak LLN proof + CLT + Bayes"
```

---

### Task 5: Worked Examples

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Worked Examples`). VERIFY arithmetic in Python before committing.

- [ ] **Step 1: Example (a) — log-likelihood maximized = normal-equations `$\hat\beta$`.**
`###` heading. Reuse a tiny 3-week, 2-column design (e.g. the Ch 2 worked example
`$X=\begin{bmatrix}1&0\\1&1\\0&1\end{bmatrix}$`, `$y=(1,3,2)^\top$`). Compute
`$\mathbb{E}[y]=X\beta$` symbolically and `$\operatorname{Cov}(y)=\sigma^2 I_3$`,
write the Gaussian log-likelihood `$\ell(\beta)$`, take `$\nabla_\beta\ell = \frac{1}{\sigma^2}X^\top(y-X\beta)=0$`,
and show it yields the SAME normal equations / `$\hat\beta=(1,2)^\top$` as Ch 2.
Emphasize: the noise scale `$\sigma^2$` drops out of the `$\beta$`-MLE.

- [ ] **Step 2: Example (b) — sampling distribution of a mean (CLT / standard error).**
`###` heading. Take `$n$` i.i.d. weekly observations with variance `$\sigma^2$`; show
`$\operatorname{Var}(\bar X_n)=\sigma^2/n$`, so the standard error is
`$\sigma/\sqrt n$` — quadrupling weeks halves the error. Give a concrete numeric
instance (e.g. `$\sigma=2$`, `$n=4\to16\to64$` giving SE `$1\to0.5\to0.25$`) and one
sentence connecting `$\sqrt n$` shrinkage to the CLT and to "more weeks of data → a
sharper `$\hat\beta$`".

- [ ] **Step 3: PY check** — write `/tmp/ch03_examples_check.py` reproducing (a)'s
`$\hat\beta$` (via `np.linalg.solve(X.T@X, X.T@y)`) and (b)'s SE values; run
`python3 /tmp/ch03_examples_check.py`; expect `$\hat\beta=[1,2]$` and SE `[1,0.5,0.25]`.
Fix prose to match.

- [ ] **Step 4: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): worked examples (Gaussian MLE; sampling-distribution SE)"
```

---

### Task 6: Code Tie-in (runnable NumPy/SciPy/matplotlib cell)

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Code Tie-in`).

- [ ] **Step 1: Write a lead-in + this Quarto Python cell.**

````markdown
```{python}
import numpy as np
from scipy.linalg import cholesky
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)
n, sigma = 60, 0.5
X = np.column_stack([np.ones(n), rng.normal(size=n), rng.normal(size=n)])
beta_true = np.array([1.0, 2.0, -0.5])
XtX_inv = np.linalg.inv(X.T @ X)

# Sampling distribution of beta_hat over many simulated datasets y = X beta + noise.
B = 5000
beta_hats = np.empty((B, 3))
for b in range(B):
    y = X @ beta_true + rng.normal(scale=sigma, size=n)
    beta_hats[b] = np.linalg.lstsq(X, y, rcond=None)[0]

emp_mean = beta_hats.mean(axis=0)
emp_cov = np.cov(beta_hats.T)
theory_cov = sigma**2 * XtX_inv
print("empirical mean of beta_hat:", np.round(emp_mean, 3), " (true:", beta_true, ")")
print("max |empirical cov - sigma^2 (X^T X)^-1|:",
      f"{np.max(np.abs(emp_cov - theory_cov)):.2e}")

# Correlated draws via the Cholesky factor (Chapter 1): Cov(Lz) = L L^T.
Sigma = np.array([[1.0, 0.6], [0.6, 1.0]])
L = cholesky(Sigma, lower=True)
draws = L @ rng.standard_normal((2, 10000))
print("Cholesky draws empirical cov:\n", np.round(np.cov(draws), 2))

# CLT in action: the sampling distribution of one coefficient is Gaussian.
plt.hist(beta_hats[:, 2], bins=40, density=True, alpha=0.65)
plt.xlabel(r"$\hat\beta_{\mathrm{search}}$"); plt.ylabel("density")
plt.title("Sampling distribution of one coefficient")
plt.show()
```
````

After the cell, write 2–3 sentences reading the output: the empirical mean of
`$\hat\beta$` sits on `$\beta_{\text{true}}$` (unbiased), its empirical covariance
matches `$\sigma^2(X^\top X)^{-1}$` to a couple of decimals (the inverse Fisher
information of rung 5), the Cholesky draws reproduce the target correlation `$0.6$`
(Ch 1), and the histogram is visibly bell-shaped (the CLT).

- [ ] **Step 2: PY check** — copy the cell body to `/tmp/ch03_codecell.py`; run
`MPLBACKEND=Agg python3 /tmp/ch03_codecell.py`; expect exit 0, `emp_mean ≈ [1, 2, -0.5]`,
`max |emp cov - theory|` small (`< 1e-2`), Cholesky cov ≈ `[[1,0.6],[0.6,1]]`. Paste
exact stdout into the report.

- [ ] **Step 3: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): runnable code tie-in (sampling distribution, Cholesky, CLT)"
```

---

### Task 7: Exercises (C / B / P / A)

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (`## Exercises`).

- [ ] **Step 1: Populate the four tiers (2–3 each).**
- **C — Conceptual:** (C1) In words, what does `$\operatorname{Cov}(y)=\sigma^2 I$`
  assert about the weeks of data? (C2) Why does "OLS is the Gaussian MLE" matter —
  what does it let you say about `$\hat\beta$` that Ch 1's geometry alone could not?
  (C3) What does a posterior `$p(\beta\mid y)$` give you that a single MLE `$\hat\beta$`
  does not?
- **B — By hand:** (B1) For `$w=(w_1,w_2)$` with `$\operatorname{Cov}(w)=\begin{bmatrix}2&1\\1&2\end{bmatrix}$`
  and `$A=\begin{bmatrix}1&1\\1&-1\end{bmatrix}$`, compute `$\operatorname{Cov}(Aw)$`.
  (B2) Write the bivariate-normal density for `$\mu=0$`,
  `$\Sigma=\begin{bmatrix}1&0.5\\0.5&1\end{bmatrix}$` (give `$\Sigma^{-1}$` and
  `$\det\Sigma$` explicitly). (B3) A single Bayes update: prior
  `$\Pr(\text{high season})=0.3$`, likelihoods `$\Pr(\text{spike}\mid\text{high})=0.8$`,
  `$\Pr(\text{spike}\mid\text{low})=0.2$`; compute `$\Pr(\text{high}\mid\text{spike})$`.
- **P — Prove it:** (P1) Prove `$\operatorname{Cov}(Aw+b)=A\operatorname{Cov}(w)A^\top$`.
  (P2) Prove that for `$y\sim\mathcal N(X\beta,\sigma^2 I)$` the MLE of `$\beta$`
  solves the normal equations. (P3) Prove the weak LLN via Chebyshev.
- **A — Applied / code:** (A1) Monte-Carlo: simulate `$y=X\beta+\varepsilon$` `$B$`
  times, collect `$\hat\beta$`, and verify the empirical covariance matches
  `$\sigma^2(X^\top X)^{-1}$`. (A2) CLT experiment: for sample sizes
  `$n\in\{1,5,30,200\}$`, draw many sample means of a non-normal variable (e.g.
  Exponential) and show the histograms approach a Gaussian as `$n$` grows.

No inline solution links (per the preface).

- [ ] **Step 2: MATH check** — `$$` count even.

- [ ] **Step 3: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): exercises across C/B/P/A tiers"
```

---

### Task 8: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`. VERIFY numeric answers (B1, B3, A1) in Python.

- [ ] **Step 1: Inspect existing structure.** Read `appendix/solutions.qmd`: it uses
`::: {.content-visible when-meta="show-solutions"}` with per-chapter `## Chapter N — …`
sections (Ch 1, Ch 2 already present). Add `## Chapter 3 — Probability & Statistics`
INSIDE that visible block, AFTER the Chapter 2 section and BEFORE the closing `:::`,
with `###` tier subheadings (C/B/P/A). Match the existing idiom exactly; do not alter
the gating divs or other chapters' solutions.

- [ ] **Step 2: Write solutions.**
- C1–C3: 2–4 sentence model answers.
- B1: `$\operatorname{Cov}(Aw)=A\Sigma A^\top$` with `$\Sigma=\begin{bmatrix}2&1\\1&2\end{bmatrix}$`,
  `$A=\begin{bmatrix}1&1\\1&-1\end{bmatrix}$` → compute explicitly (VERIFY in Python).
- B2: give `$\Sigma^{-1}=\frac{1}{0.75}\begin{bmatrix}1&-0.5\\-0.5&1\end{bmatrix}$`,
  `$\det\Sigma=0.75$`, and the density.
- B3: Bayes update → `$\Pr(\text{high}\mid\text{spike}) = \frac{0.8\cdot0.3}{0.8\cdot0.3+0.2\cdot0.7}$`
  (VERIFY in Python; = `$0.24/0.38 \approx 0.632$`).
- P1–P3: full self-contained proofs (mirror rungs 1, 5, 6).
- A1: expected result (empirical cov ≈ `$\sigma^2(X^\top X)^{-1}$`); short snippet
  sketch + qualitative trend. A2: histograms converge to Gaussian as `$n$` grows.

- [ ] **Step 3: PY check** — `/tmp/ch03_sol_check.py` computes B1 `$\operatorname{Cov}(Aw)$`
and B3 posterior; run `python3 /tmp/ch03_sol_check.py`; make prose match. Expect B3 ≈ `0.6316`.

- [ ] **Step 4: Commit**
```
git add appendix/solutions.qmd
git commit -m "feat(ch03): appendix solutions for B/P/A exercises"
```

---

### Task 9: Summary

**Files:** Modify `parts/01-foundations/03-probability-statistics.qmd` (ADD `## Summary` after `## Exercises`).

- [ ] **Step 1: Write the Summary section** (matching the Ch 1/Ch 2 format: a
**Key concepts** bulleted recap + a **Key identities** bulleted list, inline math
only to keep display-math parity). Cover:
- Concepts: `$y$` is random (`$\mathbb{E}[y]=X\beta$`, `$\operatorname{Cov}(y)=\sigma^2 I$`);
  covariance transforms as `$A(\cdot)A^\top$`; the multivariate normal and Cholesky
  sampling; **OLS = Gaussian MLE**; Fisher information ties precision to `$X^\top X$`
  conditioning; LLN/CLT license Monte Carlo and explain the Gaussian; Bayes turns the
  likelihood into a posterior on `$\beta$`.
- Identities (inline): `$\operatorname{Cov}(Aw+b)=A\operatorname{Cov}(w)A^\top$`;
  `$\mathbb{E}[y]=X\beta,\ \operatorname{Cov}(y)=\sigma^2 I$`; MVN density;
  `$w=\mu+Lz \Rightarrow w\sim\mathcal N(\mu,LL^\top)$`;
  `$\ell(\beta)=-\frac{1}{2\sigma^2}\lVert y-X\beta\rVert^2+\text{const}$` and MLE
  solves `$X^\top X\hat\beta=X^\top y$`; `$\mathcal I(\beta)=X^\top X/\sigma^2$`;
  weak LLN (`$\operatorname{Var}(\bar X_n)=\sigma^2/n$`); Bayes
  `$p(\beta\mid y)\propto p(y\mid\beta)p(\beta)$`.

- [ ] **Step 2: MATH check** — `$$` count even; `## Summary` is the LAST `##` heading.

- [ ] **Step 3: Commit**
```
git add parts/01-foundations/03-probability-statistics.qmd
git commit -m "feat(ch03): Summary section"
```

---

### Task 10: Integration sweep

**Files:** Possibly modify the chapter and/or `appendix/solutions.qmd` (fixes only).

- [ ] **Step 1: Stub/notation/cross-ref sweep.** Read the whole chapter. Confirm:
notation matches the conventions table; rung cross-references (to Ch 1 `$X^\top X$`/
Cholesky/`$\kappa$`, Ch 2 loss/gradient) resolve; the four proofs are present; all
four exercise tiers + Summary are populated.
Run `grep -ni 'stub\|_why this matters\|_worked numerical\|_minimal, runnable\|_definitions ->\|TODO\|TBD' parts/01-foundations/03-probability-statistics.qmd` → expect no matches.

- [ ] **Step 2: Final MATH check** — `grep -c '\$\$' …` even; no bare `\begin{align}`;
spot-check display blocks for KaTeX validity.

- [ ] **Step 3: Final PY check** — re-extract and run the code cell
(`MPLBACKEND=Agg python3 /tmp/ch03_codecell.py`, exit 0).

- [ ] **Step 4: RENDER gate** — note in the PR that local Quarto is unavailable and
CI (`.github/workflows/render.yml`) is the render authority.

- [ ] **Step 5: Commit any fixes**
```
git add -A
git commit -m "fix(ch03): integration sweep — notation, refs, cleanup"
```

---

## Self-Review (completed by plan author)

**Spec coverage:** Through-line (`y=Xβ+ε`) → Tasks 1–6. Theory rungs 1–7 → Tasks 2–4.
Four named proofs — `Cov(Ay)=A Cov(y)Aᵀ` (Task 2), affine-Gaussian (Task 3 rung 4),
Gaussian MLE = OLS (Task 3 rung 5), weak LLN via Chebyshev (Task 4) — all present.
Multivariate normal + Cholesky link → Task 3. Fisher information `XᵀX/σ²` → Task 3.
Worked examples (a)/(b) → Task 5. NumPy/SciPy tie-in with sampling distribution +
Cholesky + CLT → Task 6. C/B/P/A exercises → Task 7. Appendix solutions gated by
`show-solutions` → Task 8. Summary section → Task 9. Render + notation success
criteria → Task 10. CLT cited not proved → Task 4. Out-of-scope items never
introduced. All spec sections map to a task.

**Placeholder scan:** No TBD/TODO/"handle edge cases"; the code cell is given in full;
proof contents stated explicitly, not deferred.

**Type/notation consistency:** `make`-free direct NumPy; `X`, `beta_true`, `beta_hats`,
`XtX_inv`, `emp_cov`, `theory_cov`, `Sigma`, `L` consistent across Tasks 5/6/8;
math symbols (`$\mathbb{E}$`, `$\operatorname{Cov}$`, `$\mathcal N$`, `$\ell$`,
`$\mathcal I$`, `$X^\top X$`, `$\hat\beta$`) used uniformly per the conventions table.
