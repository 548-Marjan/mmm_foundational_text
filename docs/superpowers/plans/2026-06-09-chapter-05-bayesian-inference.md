# Chapter 5 — Bayesian Inference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/02-regression-bayes/02-bayesian-inference.qmd` with a complete, MMM-anchored Bayesian Inference chapter that turns Chapter 4's ridge point estimate into a full posterior distribution, and add matching appendix solutions.

**Architecture:** One Quarto `.qmd` chapter file authored section-by-section against the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus a `## Chapter 5 — Bayesian Inference` block appended to the shared `appendix/solutions.qmd` (gated by `show-solutions`). Each section is a task; the through-line is Chapter 4's "ridge is the posterior mode — here is the whole posterior." The single runnable Python cell is verified headless before commit; the real render gate is CI (`quarto render`), Quarto not being installed locally.

**Tech Stack:** Quarto book, KaTeX (HTML) / LaTeX (PDF) math, Python code cell (numpy, scipy, matplotlib pinned in `requirements.txt`). Canonical references via `references.bib` keys: `@gelman2013`, `@hoff2009`, `@hastie2009`, `@casella2002`.

---

## Conventions (apply to every task)

- **Math/KaTeX:** Use `aligned` inside `$$ ... $$`, never bare `\begin{align}`. Keep `$$` delimiters on their own lines so display-math counts stay even. Inline math with single `$`. Reuse established symbols: `X^\top X`, `\hat\beta`, ridge `(X^\top X+\lambda I)^{-1}X^\top y`, `\lambda=\sigma^2/\tau^2`, `\mathcal N(\mu,\Sigma)`, likelihood `\ell`.
- **Voice:** Match Chapters 1–4 (see `parts/02-regression-bayes/01-linear-regression.qmd`) — rigorous, prose-driven, each theory subsection ("rung") ending by tying back to the through-line / an MMM payoff. End every full proof with `$\blacksquare$`.
- **Citations:** `[@gelman2013]`, `[@hoff2009]`, `[@hastie2009]`, `[@casella2002]` only (these keys exist in `references.bib`).
- **Verification:** Verify every numeric claim against NumPy before committing. The code cell must run top-to-bottom under `MPLBACKEND=Agg python3`. No pytest exists; the build is the test.
- **Commits:** One commit per task, message prefix `feat(ch05): ...`. Run `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"` once if commits fail with "Author identity unknown".
- **No inline solution links** in exercises; solutions live only in the appendix, gated by `::: {.content-visible when-meta="show-solutions"}`.

## Verified anchor numbers (use these exact values)

- **Worked (a) Beta–Binomial:** prior `Beta(2,2)`, observe `k=8` successes in `n=10` → posterior `Beta(10,4)`; posterior mean `10/14 ≈ 0.714286`; MLE `0.8`. Posterior mean sits between prior mean `0.5` and MLE `0.8` (shrinkage).
- **Worked (b) Bayesian linear regression** on `X = [[1,1],[1,2],[1,3],[1,4]]`, `y = (2,2,4,5)^\top`, with `\sigma^2 = 1`, `\tau^2 = 1` (so `\lambda = \sigma^2/\tau^2 = 1`):
  - `X^\top X = [[4,10],[10,30]]`, `X^\top y = (13,38)^\top`.
  - `\Sigma_{\text{post}} = (X^\top X/\sigma^2 + I/\tau^2)^{-1} = [[0.563636, -0.181818], [-0.181818, 0.090909]]`.
  - `\mu_{\text{post}} = \Sigma_{\text{post}}\,X^\top y/\sigma^2 = (0.418182, 1.090909)^\top`, **exactly equal** to ridge `(X^\top X+\lambda I)^{-1}X^\top y`.
  - Slope posterior: mean `1.090909`, sd `\sqrt{0.090909} = 0.301511`; 95% **credible** interval `1.090909 ± 1.96·0.301511 = [0.500, 1.682]`.
  - Contrast: OLS on the same data gives slope `1.1` exactly, and Chapter 4's Worked Example reports the wider 95% **confidence** interval `[-0.038, 2.238]` (t₂ quantile 4.303). The Bayesian interval is tighter because the prior contributes information.

---

## File Structure

- **Modify (replace body):** `parts/02-regression-bayes/02-bayesian-inference.qmd` — the chapter. Keep the existing `# Bayesian Inference` H1 and the `*Canonical anchors: ...*` line; replace the stub callout and empty section bodies with full content.
- **Modify (append):** `appendix/solutions.qmd` — insert a `## Chapter 5 — Bayesian Inference` block immediately before the closing `:::` of the `show-solutions` visible div (currently line 481). Match the format of the existing `## Chapter 4 — Linear Regression` block.

---

### Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd:1-34` (replace stub callout through the empty `## Motivation` body)

- [ ] **Step 1: Replace the stub header and write the Motivation section**

Keep lines 1–3 (`# Bayesian Inference`, blank, `*Canonical anchors: Gelman et al. (BDA3); Hoff.*`). Delete the `::: {.callout-note ...}` stub block (lines 5–7). Replace the placeholder section bodies as the following tasks dictate; in this task, write `## Motivation`.

The Motivation must:
- Open on Chapter 4's cliffhanger: ridge returned a single number — the posterior **mode** — and threw away the rest of the distribution it sat on top of.
- Pose the driving question verbatim (emphasized): *"Ridge handed you $\hat\beta_{\text{TV}} = 1.8$. But how plausible is 1.5? Is 0 ruled out? And what will next quarter's sales actually be?"*
- Make the thesis explicit: only the **full posterior** answers these, and obtaining it is just **conditioning** — which Chapter 3's Bayes' theorem already taught. Name the payoffs the chapter will deliver: credible intervals (the honest answer to "how sure are you?"), the posterior predictive (honest forecasts for budget planning), priors as explicit modeling choices (including non-negative media effects), and the wall where conjugacy fails and Part III's MCMC becomes necessary.
- ~3 paragraphs, matching the density of Chapter 4's Motivation.

- [ ] **Step 2: Commit**

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): front matter and Motivation"
```

---

### Task 2: Theory & Proofs — rungs 1–4 (framework + three conjugate proofs)

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (the `## Theory & Proofs` section, first four rungs)

Write `## Theory & Proofs` with these subsections (use `###` headings, prose-led, each ending on the through-line). Math must be correct and KaTeX-safe.

- [ ] **Step 1: Rung 1 — The Bayesian framework**

Define prior `p(\beta)`, likelihood `p(y\mid\beta)` (the Chapter 3 Gaussian likelihood), and the posterior via Bayes' theorem:
$$
p(\beta\mid y) = \frac{p(y\mid\beta)\,p(\beta)}{p(y)} \propto p(y\mid\beta)\,p(\beta).
$$
Introduce the evidence (normalizer) `p(y)=\int p(y\mid\beta)p(\beta)\,d\beta` and the **posterior predictive** `p(\tilde y\mid y)=\int p(\tilde y\mid\beta)\,p(\beta\mid y)\,d\beta`. Land the slogan: estimation becomes **conditioning** — the posterior is the whole object, a point estimate is just a summary of it.

- [ ] **Step 2: Rung 2 — Conjugacy & the Beta–Binomial (full proof)**

Define a **conjugate** prior (prior and posterior in the same family). Set up the channel-conversion question: `k` conversions in `n` trials, likelihood `\propto \theta^k (1-\theta)^{n-k}`; prior `\theta\sim\text{Beta}(\alpha,\beta)\propto \theta^{\alpha-1}(1-\theta)^{\beta-1}`. Prove by **kernel-matching** that the posterior is
$$
\theta\mid k \sim \text{Beta}(\alpha+k,\ \beta+n-k),
$$
i.e. multiply, collect exponents, recognize the Beta kernel; the normalizer is forced. Interpret `\alpha,\beta` as **pseudo-counts** (prior successes/failures). End with `$\blacksquare$`. This is the clean 1-D warm-up.

- [ ] **Step 3: Rung 3 — The Gaussian conjugate (full proof)**

Normal prior `\mu\sim\mathcal N(\mu_0,\tau^2)` + Normal likelihood with **known** variance `\sigma^2` and `n` observations with sample mean `\bar y`. Prove by **completing the square** in `\mu` that the posterior is Gaussian with
$$
\frac{1}{\sigma_{\text{post}}^2} = \frac{1}{\tau^2} + \frac{n}{\sigma^2},
\qquad
\mu_{\text{post}} = \sigma_{\text{post}}^2\left(\frac{\mu_0}{\tau^2} + \frac{n\bar y}{\sigma^2}\right).
$$
State the two readings: **precisions add**, and the posterior mean is a **precision-weighted average** of prior mean and sample mean. End with `$\blacksquare$`. Note this is the scalar template for the keystone.

- [ ] **Step 4: Rung 4 — Bayesian linear regression (KEYSTONE, full proof)**

Prior `\beta\sim\mathcal N(0,\tau^2 I)`, Gaussian likelihood `y\mid\beta\sim\mathcal N(X\beta,\sigma^2 I)`. Prove (completing the square in the vector `\beta`, matching a Gaussian kernel) that
$$
\beta\mid y \sim \mathcal N(\mu_{\text{post}},\ \Sigma_{\text{post}}),
\qquad
\Sigma_{\text{post}} = \left(\frac{X^\top X}{\sigma^2} + \frac{I}{\tau^2}\right)^{-1},
\qquad
\mu_{\text{post}} = \Sigma_{\text{post}}\,\frac{X^\top y}{\sigma^2}.
$$
Then show algebraically that `\mu_{\text{post}} = (X^\top X + \lambda I)^{-1}X^\top y` with `\lambda=\sigma^2/\tau^2` — **exactly the Chapter 4 ridge estimator**. Emphasize the payoff: ridge's MAP was the *mode* of this Gaussian (for a Gaussian the mode equals the mean), and now we have the entire distribution it crowned. End with `$\blacksquare$`. Cite `[@hoff2009]` for the conjugate-regression derivation.

- [ ] **Step 5: Verify the keystone numerics against NumPy, then commit**

Confirm on `X=[[1,1],[1,2],[1,3],[1,4]]`, `y=(2,2,4,5)`, `\sigma^2=\tau^2=1` that `\mu_{\text{post}}=(0.418182,1.090909)` equals `np.linalg.solve(XtX+np.eye(2), Xty)`:

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
X=np.array([[1,1],[1,2],[1,3],[1,4]],float); y=np.array([2,2,4,5],float)
XtX=X.T@X; Xty=X.T@y
S=np.linalg.inv(XtX+np.eye(2)); mu=S@Xty
print(np.round(mu,6), np.allclose(mu, np.linalg.solve(XtX+np.eye(2),Xty)))"
```
Expected: `[0.418182 1.090909] True`.

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Theory rungs 1-4 — framework and three conjugate proofs"
```

---

### Task 3: Theory & Proofs — rungs 5–7 (credible intervals, posterior predictive, priors & the wall)

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (append three subsections to `## Theory & Proofs`)

- [ ] **Step 1: Rung 5 — Credible intervals & posterior summaries**

Define the **credible interval**: a central `1-\alpha` interval of the posterior, e.g. `\mu_{\text{post},j} \pm z_{1-\alpha/2}\sqrt{[\Sigma_{\text{post}}]_{jj}}` for the Gaussian posterior. Give its **honest interpretation** — a direct probability statement `P(\beta_j\in[\,\cdot\,]\mid y)=1-\alpha` about the parameter — and contrast it **sharply** with the frequentist confidence interval of Chapter 4 (whose 95% refers to the procedure's long-run coverage, not to a probability about the fixed `\beta_j`). Name the point summaries: posterior mean, median, MAP (mode); for a Gaussian posterior all three coincide, which is why ridge's mode equals the posterior mean.

- [ ] **Step 2: Rung 6 — The posterior predictive (derivation)**

Derive, for Bayesian linear regression, that a new observation at `\tilde x` has
$$
\tilde y\mid y \sim \mathcal N\big(\tilde x^\top\mu_{\text{post}},\ \tilde x^\top\Sigma_{\text{post}}\,\tilde x + \sigma^2\big).
$$
Show the derivation as a Gaussian integral / linear-Gaussian convolution: `\tilde y = \tilde x^\top\beta + \tilde\varepsilon` with `\beta\mid y\sim\mathcal N(\mu_{\text{post}},\Sigma_{\text{post}})` and independent `\tilde\varepsilon\sim\mathcal N(0,\sigma^2)`; the mean and variance follow from linearity and independence. Decompose the predictive variance into **parameter uncertainty** `\tilde x^\top\Sigma_{\text{post}}\tilde x` **+ irreducible noise** `\sigma^2`. Tie to MMM: budget planning needs the full forecast band, not a point.

- [ ] **Step 3: Rung 7 — Priors in practice & the wall**

Discuss informative vs **weakly-informative** priors (the pseudo-count reading from rung 2; vague `\tau^2\to\infty` recovers OLS). Make the MMM-critical point: media effects are **non-negative**, encoded by a **half-normal / truncated** prior on the coefficients — which is **not** conjugate to the Gaussian likelihood, so the tidy closed form **breaks**. Add that real MMM transforms (adstock, saturation) make the likelihood nonlinear in the parameters, so the posterior has **no closed form** at all. State and cite the general fact `[@gelman2013]` that without conjugacy the posterior and predictive are analytically intractable — which is exactly why **Part III builds MCMC**. End the section on this baton-pass.

- [ ] **Step 4: Commit**

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Theory rungs 5-7 — credible intervals, posterior predictive, priors and the wall"
```

---

### Task 4: Worked Examples

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (the `## Worked Examples` section)

- [ ] **Step 1: Example (a) — Beta–Binomial by hand**

Prior `\text{Beta}(2,2)` (mean `0.5`), observe `k=8` conversions in `n=10`. Posterior `\text{Beta}(2+8,\,2+2)=\text{Beta}(10,4)`. Posterior mean `\frac{10}{10+4}=\frac{10}{14}\approx 0.714`, versus the MLE `\hat\theta=8/10=0.8`. Walk the arithmetic and interpret: the posterior mean is pulled from the MLE toward the prior mean `0.5` — **shrinkage** driven by the prior pseudo-counts. Note the posterior is sharper than the prior (data added information).

- [ ] **Step 2: Example (b) — Bayesian linear regression on a tiny design**

Use `X=\begin{bmatrix}1&1\\1&2\\1&3\\1&4\end{bmatrix}`, `y=(2,2,4,5)^\top`, `\sigma^2=1`, `\tau^2=1` (so `\lambda=1`). Show:
$$
X^\top X = \begin{bmatrix}4&10\\10&30\end{bmatrix},\qquad X^\top y=\begin{bmatrix}13\\38\end{bmatrix},
$$
$$
\Sigma_{\text{post}} = (X^\top X + I)^{-1} = \begin{bmatrix}0.5636&-0.1818\\-0.1818&0.0909\end{bmatrix},\qquad
\mu_{\text{post}} = \begin{bmatrix}0.4182\\1.0909\end{bmatrix}.
$$
Confirm `\mu_{\text{post}}` equals the ridge estimate `(X^\top X+I)^{-1}X^\top y`. Give the 95% **credible** interval for the slope: mean `1.0909`, sd `\sqrt{0.0909}=0.3015`, so `1.0909\pm1.96\times0.3015=[0.500,\,1.682]`. Contrast explicitly with Chapter 4's 95% **confidence** interval on the *same* `X,y` (OLS slope `1.1`, interval `[-0.038,\,2.238]` with `t_{2,0.975}=4.303`): the Bayesian interval is tighter and is a direct probability statement, because the prior `\mathcal N(0,\tau^2 I)` supplied information that the frequentist analysis lacked.

- [ ] **Step 3: Verify (b) numerics and commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
X=np.array([[1,1],[1,2],[1,3],[1,4]],float); y=np.array([2,2,4,5],float)
XtX=X.T@X; Xty=X.T@y
S=np.linalg.inv(XtX+np.eye(2)); mu=S@Xty
sd=np.sqrt(S[1,1])
print('mu',np.round(mu,4),'Sigma',np.round(S,4).tolist())
print('CI',np.round(mu[1]+np.array([-1,1])*1.959964*sd,3))
print('ridge ok',np.allclose(mu,np.linalg.solve(XtX+np.eye(2),Xty)))"
```
Expected: `mu [0.4182 1.0909]`, `CI [0.5 1.682]`, `ridge ok True`.

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Worked Examples — Beta-Binomial and Bayesian linear regression"
```

---

### Task 5: Code Tie-in

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (the `## Code Tie-in` section)

- [ ] **Step 1: Write a single self-contained ```` ```{python} ```` cell**

Model the structure and prose-wrapping on Chapter 4's Code Tie-in (`parts/02-regression-bayes/01-linear-regression.qmd:339-409`). The cell (numpy/scipy/matplotlib, seeded `rng = np.random.default_rng(...)`) must, with helper functions:
1. **Beta–Binomial updating:** start from `Beta(2,2)`, feed conversions sequentially (e.g. blocks of data), and plot prior → posterior densities sharpening as `n` grows (use `scipy.stats.beta`). Print the closed-form posterior `Beta(\alpha+k,\beta+n-k)` parameters and confirm the posterior mean moves from prior `0.5` toward the empirical rate.
2. **Bayesian linear regression in closed form** on synthetic media data `y=X\beta+\varepsilon` (intercept + 2 channels): compute `\Sigma_{\text{post}}` and `\mu_{\text{post}}`, and **assert** `np.allclose(mu_post, np.linalg.solve(XtX + lam*np.eye(p), Xty))` (i.e. `\mu_{\text{post}}` = ridge), printing the confirmation.
3. **Posterior samples & credible band:** draw `\beta` samples from `\mathcal N(\mu_{\text{post}},\Sigma_{\text{post}})` (via `np.random.multivariate_normal` or a Cholesky factor `L`), and plot per-coefficient 95% credible intervals.
4. **Posterior-predictive band:** over a grid of `\tilde x`, plot the predictive mean `\tilde x^\top\mu_{\text{post}}` with the `\pm1.96\sqrt{\tilde x^\top\Sigma_{\text{post}}\tilde x + \sigma^2}` band.
5. **Prior-strength sweep:** sweep `\tau^2` (small → large) and show `\mu_{\text{post}}` interpolating from heavy shrinkage toward the OLS solution as `\tau^2\to\infty` (`\lambda\to 0`); print the endpoints.

End every figure path with `plt.show()`. Follow the cell with a prose paragraph reading off the printed numbers (matching Chapter 4's pattern of quoting concrete values).

- [ ] **Step 2: Extract the cell to a temp file and run it headless**

Save the cell body to `/tmp/ch5_code.py` and run:

```bash
MPLBACKEND=Agg python3 /tmp/ch5_code.py
```
Expected: runs top-to-bottom with no error; prints the `\mu_{\text{post}}` = ridge confirmation (`True`) and the `\tau^2` sweep endpoints. If any printed value is quoted in the surrounding prose, make the prose match the actual output.

- [ ] **Step 3: Commit**

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Code Tie-in — Beta-Binomial updating and Bayesian linear regression"
```

---

### Task 6: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (the `## Exercises` section with the four `###` tier headings)

Match Chapter 4's exercise style and heading text exactly: `### C -- Conceptual / Reading Comprehension`, `### B -- By Hand`, `### P -- Prove It`, `### A -- Applied / Code`. Each problem fully self-contained (state all data/matrices inline). **No links to solutions.**

- [ ] **Step 1: Write tier C (Conceptual), 3 problems**

- **C1.** Interpret prior / posterior / posterior predictive; explain why the posterior is "the whole object" and a point estimate is a summary.
- **C2.** Credible vs confidence interval: state precisely what each "95%" claims, and why only the credible interval is a probability statement about `\beta`. (Mirror the contrast in Worked (b).)
- **C3.** Explain why a non-negativity (half-normal) prior on a media coefficient **breaks** Gaussian conjugacy and forces computation (MCMC), connecting to the Part III hand-off.

- [ ] **Step 2: Write tier B (By Hand), 3 problems**

- **B1.** A Beta–Binomial update: given a stated prior `\text{Beta}(\alpha,\beta)` and `k`/`n` (choose clean numbers, e.g. `\text{Beta}(3,1)`, `k=5`, `n=8` → `\text{Beta}(8,4)`), report the posterior and its mean, comparing to the MLE.
- **B2.** A Normal–Normal posterior: given `\mu_0,\tau^2,\sigma^2,n,\bar y`, compute the precision-weighted posterior mean and posterior variance (precisions add). Use clean numbers.
- **B3.** Compute `\mu_{\text{post}}=(X^\top X+\lambda I)^{-1}X^\top y` for a small design and show it equals the ridge estimate; reuse a `2\times2` `X^\top X` so the inverse is by-hand (e.g. the Chapter 4 B3 matrices `X^\top X=[[4,2],[2,4]]`, `X^\top y=(6,4)^\top`, `\lambda=2`).

- [ ] **Step 3: Write tier P (Prove It), 3 problems**

- **P1.** Prove Beta–Binomial conjugacy (kernel-matching), deriving `\text{Beta}(\alpha+k,\beta+n-k)`.
- **P2.** Prove the Normal–Normal posterior (complete the square): precisions add and the mean is precision-weighted.
- **P3.** Prove the Bayesian-linear-regression posterior `\mathcal N(\mu_{\text{post}},\Sigma_{\text{post}})` and that `\mu_{\text{post}}=(X^\top X+\lambda I)^{-1}X^\top y` with `\lambda=\sigma^2/\tau^2` (= ridge).

- [ ] **Step 4: Write tier A (Applied / Code), 2 problems**

- **A1.** Beta–Binomial **sequential updating**: stream data in, plot the posterior sharpening, and confirm the closed-form posterior matches a fine-grid numerical posterior (likelihood × prior, normalized).
- **A2.** Bayesian linear regression: compute the posterior, draw samples, plot coefficient **credible** bands and a **posterior-predictive** band, and run a **prior-sensitivity sweep** over `\tau^2` showing the interpolation from shrinkage to OLS.

- [ ] **Step 5: Commit**

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Exercises — C/B/P/A tiers"
```

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/02-regression-bayes/02-bayesian-inference.qmd` (append `## Summary` as the final section)

- [ ] **Step 1: Write `## Summary`**

Match Chapter 4's Summary format exactly: a one-line lead, then **Key concepts** (bulleted, bold lead-ins) and **Key identities** (bulleted, **inline math only** — single `$`, no display blocks). Cover:
- **Key concepts:** Bayes' update (estimation = conditioning); conjugacy and pseudo-counts; Beta–Binomial; the Gaussian conjugate (precisions add, precision-weighted mean); Bayesian linear regression posterior (mean = ridge, so ridge was the mode); credible vs confidence interval; posterior predictive = parameter uncertainty + noise; the conjugacy wall → MCMC.
- **Key identities (inline):**
  - Bayes: $p(\beta\mid y)\propto p(y\mid\beta)\,p(\beta)$.
  - Beta–Binomial: prior $\text{Beta}(\alpha,\beta)$ + $k$/$n$ → $\text{Beta}(\alpha+k,\beta+n-k)$.
  - Gaussian conjugate: $1/\sigma_{\text{post}}^2 = 1/\tau^2 + n/\sigma^2$ (precisions add).
  - Posterior covariance: $\Sigma_{\text{post}} = (X^\top X/\sigma^2 + I/\tau^2)^{-1}$.
  - Posterior mean = ridge: $\mu_{\text{post}} = (X^\top X+\lambda I)^{-1}X^\top y$, $\lambda=\sigma^2/\tau^2$.
  - Credible interval: $\mu_{\text{post},j}\pm z_{1-\alpha/2}\sqrt{[\Sigma_{\text{post}}]_{jj}}$.
  - Posterior predictive: $\tilde y\mid y\sim\mathcal N(\tilde x^\top\mu_{\text{post}},\ \tilde x^\top\Sigma_{\text{post}}\tilde x + \sigma^2)$.

- [ ] **Step 2: Commit**

```bash
git add parts/02-regression-bayes/02-bayesian-inference.qmd
git commit -m "feat(ch05): Summary — key concepts and identities"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (insert before the closing `:::` at line 481)

- [ ] **Step 1: Append the Chapter 5 solutions block**

Immediately before the final `:::` that closes the `show-solutions` visible div, add a `## Chapter 5 — Bayesian Inference` heading followed by full solutions to **every** exercise from Task 6 (C1–C3, B1–B3, P1–P3, A1–A2). Match the format of the existing `## Chapter 4 — Linear Regression` block (`appendix/solutions.qmd:317-479`): `### C — Conceptual / Reading Comprehension` etc. (note the em-dash `—` used in appendix headings), bold problem labels (`**C1.**`), full worked math for B and P, and a reference Python sketch for the A problems. Each proof ends with `$\blacksquare$`. Verify every numeric answer (e.g. B1 posterior mean `8/12≈0.667`; B3 `\mu_{\text{post}}=(0.875,0.375)` reusing the Chapter 4 B3 arithmetic) against NumPy.

- [ ] **Step 2: Verify by-hand answers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
# B3 reused: XtX=[[4,2],[2,4]], Xty=(6,4), lambda=2
XtX=np.array([[4,2],[2,4]],float); Xty=np.array([6,4],float)
print('B3 mu_post', np.linalg.solve(XtX+2*np.eye(2), Xty))"
```
Expected: `B3 mu_post [0.875 0.375]`.

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch05): appendix solutions for Bayesian Inference exercises"
```

---

### Task 9: Final whole-chapter review

**Files:**
- Review: `parts/02-regression-bayes/02-bayesian-inference.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: KaTeX / structure lint**

Confirm: no bare `\begin{align}` (search for it — must be zero); every `$$` is on its own line; display-math `$$` count is even; the template headings appear in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises`, `## Summary`); all citations use existing keys (`@gelman2013`, `@hoff2009`, `@hastie2009`, `@casella2002`).

```bash
cd /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch5-bayesian-inference
grep -n 'begin{align}' parts/02-regression-bayes/02-bayesian-inference.qmd || echo "no bare align — good"
grep -c '^\$\$' parts/02-regression-bayes/02-bayesian-inference.qmd   # expect even
grep -nE '^\## ' parts/02-regression-bayes/02-bayesian-inference.qmd  # expect the 6 template headings in order
```

- [ ] **Step 2: Re-run the code cell headless one final time**

```bash
MPLBACKEND=Agg python3 /tmp/ch5_code.py && echo "CODE CELL OK"
```
Expected: `CODE CELL OK` with no traceback.

- [ ] **Step 3: Through-line check (manual)**

Confirm every Theory subsection ends by tying back to the posterior / ridge-mode through-line, and that the three named conjugate proofs + the predictive derivation are all present and correct. Fix any gaps inline and commit if changes were made.

---

## Self-Review (completed by plan author)

- **Spec coverage:** Motivation (T1), 7 theory rungs incl. 3 conjugate proofs + predictive derivation (T2–T3), both worked examples (T4), code tie-in with the 5 required behaviors (T5), C/B/P/A exercises (T6), auto-included Summary (T7), appendix solutions (T8), render/lint gate (T9). All spec success criteria mapped.
- **Placeholder scan:** none — every task carries the actual content, formulas, and verified anchor numbers.
- **Consistency:** symbols (`\mu_{\text{post}}`, `\Sigma_{\text{post}}`, `\lambda=\sigma^2/\tau^2`) identical across tasks; the keystone `\mu_{\text{post}} = (X^\top X+\lambda I)^{-1}X^\top y` = ridge is stated identically in T2, T4, T5, T7. The worked-(b) numbers are NumPy-verified and reused consistently.
- **Note:** the real render gate is CI `quarto render` (HTML + PDF) on the PR — Quarto is not installed locally, so local verification is limited to running the Python cell and the lint checks above.
