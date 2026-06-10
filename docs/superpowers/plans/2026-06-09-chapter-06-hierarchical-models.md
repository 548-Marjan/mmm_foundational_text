# Chapter 6 — Hierarchical Models Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/02-regression-bayes/03-hierarchical-regression.qmd` with the complete final chapter of Part II — geo-level partial pooling that turns Chapter 5's single-population Bayesian regression into a multilevel model — and add matching appendix solutions.

**Architecture:** One Quarto `.qmd` chapter authored section-by-section against the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus a `## Chapter 6 — Hierarchical Models` block appended to the shared `appendix/solutions.qmd` (gated by `show-solutions`). The through-line is the geo-MMM three-stances story (complete / no / partial pooling); the chapter reuses Chapter 5's Rung 3 (Gaussian conjugate) and Rung 4 (Bayesian linear regression) results directly and ends by deriving the Gibbs full conditionals that hand off to Part III. The single runnable Python cell is verified headless before commit; the real render gate is CI (`quarto render`), Quarto not being installed locally.

**Tech Stack:** Quarto book, KaTeX (HTML) / LaTeX (PDF) math, Python code cell (numpy, scipy, matplotlib pinned in `requirements.txt`). Canonical references via `references.bib` keys: `@gelmanhill2006`, `@mcelreath2020`, `@gelman2013`, `@hoff2009`.

---

## Conventions (apply to every task)

- **Math/KaTeX:** Use `aligned` inside `$$ ... $$`, never bare `\begin{align}`. Keep `$$` delimiters on their own lines so display-math counts stay even. Inline math with single `$`. Reuse established symbols: `X^\top X`, `\hat\beta`, `\mathcal N(\mu,\Sigma)`, Chapter 5's `\Sigma_{\text{post}}`/`\mu_{\text{post}}`, ridge `(X^\top X+\lambda I)^{-1}X^\top y`, the "precisions add" motif.
- **Voice:** Match Chapters 1–5 (read `parts/02-regression-bayes/02-bayesian-inference.qmd` for the immediate predecessor) — rigorous, prose-led, each theory subsection ("rung") ending by tying back to the through-line / an MMM payoff. End every full proof with `$\blacksquare$`; use `> **Theorem (name).** ...` blockquotes for named results.
- **Reuse Chapter 5 explicitly:** Rung 3 (Gaussian conjugate, "precisions add") and Rung 4 (Bayesian linear regression posterior). Cite/refer rather than re-deriving from scratch.
- **Citations:** `[@gelmanhill2006]`, `[@mcelreath2020]`, `[@gelman2013]`, `[@hoff2009]` only (these keys exist in `references.bib`). Do not invent keys.
- **Verification:** Verify every numeric claim against NumPy before committing. The code cell must run top-to-bottom under `MPLBACKEND=Agg python3` and end figures with `plt.show()`. No pytest exists; the build is the test.
- **Commits:** One commit per task, message prefix `feat(ch06): ...`. Run `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"` once if commits fail with "Author identity unknown".
- **No inline solution links** in exercises; solutions live only in the appendix, gated by `::: {.content-visible when-meta="show-solutions"}`.

## Verified anchor numbers (use these exact values)

- **Worked (a) shrinkage by hand:** three geos, shared `\sigma^2 = 4`, `\tau^2 = 1`, grand mean `\mu = 2`; sample sizes `n = (4, 16, 64)`, group means `\bar y = (6.0, 3.0, 2.2)`.
  - Shrinkage weights `\lambda_g = (n_g/\sigma^2)/(n_g/\sigma^2 + 1/\tau^2)`: `\lambda = (0.5,\ 0.8,\ 0.9412)`.
  - Partial-pooled estimates `\hat\theta_g = \lambda_g\bar y_g + (1-\lambda_g)\mu`: `\hat\theta = (4.0,\ 2.8,\ 2.188)`.
  - The small geo (`n=4`, `\bar y=6`) moves `2.0` toward `\mu=2`; the large geo (`n=64`) barely moves (`0.012`). Small samples shrink hardest.
- **Worked (b) James–Stein beats the MLE:** `G = 6` group means, `\sigma^2 = 1`, one observation each. True means `\theta = (0.5,-0.3,0.8,-0.6,0.2,0.4)`, observed `z = (1.4,-1.1,1.9,-0.2,1.3,-0.7)`.
  - `\lVert z\rVert^2 = 9.0`; James–Stein shrink factor `1 - (G-2)\sigma^2/\lVert z\rVert^2 = 1 - 4/9 = 5/9 \approx 0.5556`.
  - `\hat\theta^{\text{JS}} = (5/9)z = (0.778,-0.611,1.056,-0.111,0.722,-0.389)`.
  - Total squared error vs truth: MLE `\lVert z-\theta\rVert^2 = 5.24` versus JS `\lVert\hat\theta^{\text{JS}}-\theta\rVert^2 = 1.373`. **JS wins decisively.**
- **James–Stein SURE risk identity (for the proof, Rung 3 / P2):** for `z\sim\mathcal N(\theta,\sigma^2 I_G)` and `\delta(z)=z+g(z)`, Stein's lemma gives `\mathbb E\lVert\delta-\theta\rVert^2 = G\sigma^2 + \mathbb E[\,2\sigma^2\nabla\!\cdot g(z) + \lVert g(z)\rVert^2\,]`. With `g(z) = -\tfrac{(G-2)\sigma^2}{\lVert z\rVert^2}z`: `\lVert g\rVert^2 = (G-2)^2\sigma^4/\lVert z\rVert^2`, `\nabla\!\cdot g = -(G-2)^2\sigma^2/\lVert z\rVert^2`, so the bracket equals `-(G-2)^2\sigma^4/\lVert z\rVert^2` and `\mathbb E\lVert\delta-\theta\rVert^2 = G\sigma^2 - (G-2)^2\sigma^4\,\mathbb E[1/\lVert z\rVert^2] < G\sigma^2` for `G\ge3`.

---

## File Structure

- **Modify (replace body):** `parts/02-regression-bayes/03-hierarchical-regression.qmd` — the chapter. **Retitle the H1** from `# Bayesian & Hierarchical Regression` to `# Hierarchical Models`, and set the anchors line to `*Canonical anchors: Gelman & Hill; McElreath.*`. Replace the stub callout and all empty section bodies with full content.
- **Modify (append):** `appendix/solutions.qmd` — insert a `## Chapter 6 — Hierarchical Models` block immediately before the closing `:::` of the `show-solutions` visible div (currently the last line, line 709). Match the format of the existing `## Chapter 5 — Bayesian Inference` block (`appendix/solutions.qmd:481-708`).

---

### Task 1: Retitle + front matter + Motivation

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd:1-11`

- [ ] **Step 1: Retitle and write the Motivation**

Replace line 1 `# Bayesian & Hierarchical Regression` with `# Hierarchical Models`. Replace line 3 with `*Canonical anchors: Gelman & Hill; McElreath.*`. Delete the `::: {.callout-note ...}` stub block (lines 5–7). Replace the `## Motivation` placeholder body with a full section; leave the other section headings/placeholders for later tasks.

The Motivation (~3 paragraphs, Chapter-5 density) must:
- Open on the geo-MMM dilemma: running MMM across ~50 DMAs forces a choice between **complete pooling** (one national model that hides regional heterogeneity) and **no pooling** (50 separate models whose small-DMA coefficients are pure noise).
- Pose the driving question verbatim (emphasized): *"Des Moines came back with a TV coefficient of 9.2 — absurd on 10 weeks of data. Do you trust it, throw it out, or shrink it toward the national average — and by exactly how much?"*
- State the thesis: **partial pooling** ($\beta_g\sim\mathcal N(\mu,\tau^2)$) is the principled middle that borrows strength across geos, and it falls straight out of Chapter 5's "precisions add" — except $\mu$ and $\tau^2$ are now learned from all the geos. Name the payoffs: the shrinkage formula (keystone), the James–Stein guarantee that pooling beats per-geo fits, a learned pooling strength $\tau^2$, and the Gibbs sampler (non-conjugate joint posterior) that motivates Part III.
- Anchor to $y_g = X_g\beta_g + \varepsilon_g$ and the media-mix story.

- [ ] **Step 2: Commit**

```bash
cd /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch6-hierarchical-models
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): retitle, front matter, and Motivation"
```

---

### Task 2: Theory & Proofs — rungs 1–3 (three stances + two keystone proofs)

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (the `## Theory & Proofs` section, first three rungs)

Write `## Theory & Proofs` with a 1–2 sentence ladder lead, then three `###` subsections. **Stop before rung 4** (next task). Prose-led, each rung closes on the through-line.

- [ ] **Step 1: Rung 1 — Three stances: complete, no, and partial pooling**

Formalize the two-level model. Start scalar for clarity: each geo has a group estimate $\bar y_g$ with sampling variance $\sigma^2/n_g$, and the population model $\beta_g\sim\mathcal N(\mu,\tau^2)$ (write $\theta_g$ for the scalar per-geo parameter). Frame the three stances: **complete pooling** = $\tau^2\to0$ (all $\theta_g=\mu$, low variance / high bias), **no pooling** = $\tau^2\to\infty$ (each $\theta_g=\bar y_g$ free, unbiased / high variance, ruinous for small $n_g$), **partial pooling** = finite $\tau^2$ (the compromise). Tie the bias–variance tension back to Chapter 4's MSE decomposition.

- [ ] **Step 2: Rung 2 — Partial pooling = precision-weighted shrinkage (KEYSTONE proof A)**

State and prove (this is Chapter 5 Rung 3 applied per group; complete the square, or cite Rung 3 and specialize):
> **Theorem (partial-pooling shrinkage).** For the Gaussian hierarchical model with known $\sigma^2,\mu,\tau^2$ and group estimate $\bar y_g$ of sampling variance $\sigma^2/n_g$, the posterior mean of $\theta_g$ is $\hat\theta_g = \lambda_g\bar y_g + (1-\lambda_g)\mu$ with $\lambda_g = \dfrac{n_g/\sigma^2}{n_g/\sigma^2 + 1/\tau^2}$.

Proof: the likelihood $\bar y_g\mid\theta_g\sim\mathcal N(\theta_g,\sigma^2/n_g)$ and prior $\theta_g\sim\mathcal N(\mu,\tau^2)$ are conjugate (Chapter 5 Rung 3 with data precision $n_g/\sigma^2$ and prior precision $1/\tau^2$); precisions add and the posterior mean is the precision-weighted average, which rearranges to the stated $\lambda_g$ form. End `$\blacksquare$`. Read off the limits: $n_g$ large → $\lambda_g\to1$ (no pooling), $n_g$ small → $\lambda_g\to0$ (complete pooling). Partial pooling interpolates the two stances **per geo, by precision** — the answer to the Motivation's "by exactly how much?"

- [ ] **Step 3: Rung 3 — James–Stein: shrinkage dominates (KEYSTONE proof B)**

Canonical equal-variance setup: $z\sim\mathcal N(\theta,\sigma^2 I_G)$ (one observation per group, $\sigma^2$ known). State the James–Stein estimator $\hat\theta^{\text{JS}} = \big(1 - \tfrac{(G-2)\sigma^2}{\lVert z\rVert^2}\big)z$ (shrinking toward $0$; note the practical version shrinks toward the grand mean $\bar z$).
> **Theorem (James–Stein dominance).** For $G\ge3$, $\hat\theta^{\text{JS}}$ has strictly smaller total mean-squared risk $\mathbb E\lVert\hat\theta-\theta\rVert^2$ than the MLE $z$, for every $\theta$.

Prove it via **Stein's unbiased risk estimate**. State Stein's lemma identity (cite `[@casella2002]` or `[@hoff2009]`): for $\delta(z)=z+g(z)$ with $g$ weakly differentiable, $\mathbb E\lVert\delta-\theta\rVert^2 = G\sigma^2 + \mathbb E[\,2\sigma^2\nabla\!\cdot g + \lVert g\rVert^2\,]$. Take $g(z) = -\tfrac{(G-2)\sigma^2}{\lVert z\rVert^2}z$ and compute (show the divergence step $\nabla\!\cdot(z/\lVert z\rVert^2) = (G-2)/\lVert z\rVert^2$):
$$
\mathbb E\lVert\hat\theta^{\text{JS}}-\theta\rVert^2 = G\sigma^2 - (G-2)^2\sigma^4\,\mathbb E\!\left[\frac{1}{\lVert z\rVert^2}\right] < G\sigma^2,
$$
the MLE risk, for $G\ge3$. End `$\blacksquare$`. Connect: the empirical-Bayes James–Stein shrinkage factor is exactly the data-estimated form of Rung 2's $\lambda_g$ — shrinkage is not a Bayesian indulgence; it provably lowers total error, which is why pooling wins. Cite `[@gelmanhill2006]`.

- [ ] **Step 4: Verify the Rung-3 worked numerics referenced later, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
z=np.array([1.4,-1.1,1.9,-0.2,1.3,-0.7]); th=np.array([0.5,-0.3,0.8,-0.6,0.2,0.4]); G=6
js=(1-(G-2)/(z@z))*z
print('||z||^2',z@z,'factor',round(1-(G-2)/(z@z),4))
print('MLE SE',round(((z-th)**2).sum(),4),'JS SE',round(((js-th)**2).sum(),4))"
```
Expected: `||z||^2 9.0 factor 0.5556`, `MLE SE 5.24 JS SE 1.3733`.

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Theory rungs 1-3 — pooling stances, shrinkage and James-Stein proofs"
```

---

### Task 3: Theory & Proofs — rungs 4–6 (hyperparameters, vectorized model, Gibbs bridge)

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (append three subsections to `## Theory & Proofs`, before `## Worked Examples`)

- [ ] **Step 1: Rung 4 — Where $\mu$ and $\tau^2$ come from**

Partial pooling needs the hyperparameters. Two routes: **empirical Bayes** (plug in point estimates of $\mu,\tau^2$ from the marginal/between-group spread) and **full Bayes** (hyperpriors on $\mu,\tau^2$, integrated out). The **variance-components** reading: $\tau^2$ is the pooling dial — $\tau^2\to0$ forces complete pooling, $\tau^2\to\infty$ forces no pooling — and it is **estimated from the data** (the between-geo spread of the $\bar y_g$ net of sampling noise), so the amount of shrinkage is *learned*, not chosen. State the empirical-Bayes caveat: plugging in point hyperparameters understates uncertainty by ignoring hyperparameter error; full Bayes (hence MCMC) repairs this. Cite `[@gelmanhill2006]`.

- [ ] **Step 2: Rung 5 — The vectorized geo-MMM model**

Lift from scalar to the real design: per geo $\beta_g\in\mathbb R^p$, likelihood $y_g\mid\beta_g\sim\mathcal N(X_g\beta_g,\sigma^2 I)$ (Chapter 5), population prior $\beta_g\sim\mathcal N(\mu_\beta,\Sigma_\beta)$. Pooling now happens **per coefficient** (each channel's effect partially pooled across geos). Keep $\Sigma_\beta$ diagonal ($\tau^2$ per coefficient) for the closed forms, naming the full-covariance extension once. One paragraph: the identical machinery applies to **any** grouping — product segments, time blocks — geography is the running example, not a limitation.

- [ ] **Step 3: Rung 6 — Non-conjugate joint posterior, conjugate full conditionals (DERIVATION → Part III bridge)**

Write the joint posterior $p(\beta_{1:G},\mu_\beta,\Sigma_\beta,\sigma^2\mid y) \propto \big[\prod_g p(y_g\mid\beta_g,\sigma^2)\,p(\beta_g\mid\mu_\beta,\Sigma_\beta)\big]\,p(\mu_\beta)p(\Sigma_\beta)p(\sigma^2)$ and note it has **no closed form** jointly (Chapter 5's conjugacy wall, now structural in the hierarchy). Then derive each **full conditional** and identify its conjugate family:
- $\beta_g\mid\text{rest}\sim\mathcal N(m_g,V_g)$ with $V_g = (X_g^\top X_g/\sigma^2 + \Sigma_\beta^{-1})^{-1}$, $m_g = V_g(X_g^\top y_g/\sigma^2 + \Sigma_\beta^{-1}\mu_\beta)$ — **exactly Chapter 5 Rung 4** with prior mean $\mu_\beta$ (show this; this is the proof in P3 / chapter).
- $\mu_\beta\mid\text{rest}$: Gaussian — Chapter 5 Rung 3 with the $G$ vectors $\beta_g$ as the "data" (precisions add over geos).
- $\Sigma_\beta\mid\text{rest}$: Inverse-Wishart (Inverse-Gamma per coordinate in the diagonal case) — **state** the conjugate update from $\sum_g(\beta_g-\mu_\beta)(\beta_g-\mu_\beta)^\top$.
- $\sigma^2\mid\text{rest}$: Inverse-Gamma — **state** the update from the pooled residual sum of squares $\sum_g\lVert y_g-X_g\beta_g\rVert^2$.

Conclude: sampling each conditional in turn, cycling, **is Gibbs sampling** — the simplest MCMC and the subject of Part III. "Every conditional is a Chapter 5 result; Part III shows how iterating them draws from the joint posterior we cannot write down." End the Theory section on this baton-pass. Cite `[@gelman2013]`, `[@gelmanhill2006]`.

- [ ] **Step 4: Commit**

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Theory rungs 4-6 — hyperparameters, vectorized model, Gibbs conditionals"
```

---

### Task 4: Worked Examples

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (the `## Worked Examples` section)

Match Chapter 5's worked-example style: tiny problems carried through by hand, **bold step mini-headings**, every number shown, closing interpretation.

- [ ] **Step 1: Example (a) — Shrinkage by hand**

Three geos, shared $\sigma^2 = 4$, $\tau^2 = 1$, grand mean $\mu = 2$; sample sizes $n = (4,16,64)$, group means $\bar y = (6.0, 3.0, 2.2)$. Compute each $\lambda_g = (n_g/\sigma^2)/(n_g/\sigma^2 + 1/\tau^2) = (0.5, 0.8, 0.9412)$ and each $\hat\theta_g = \lambda_g\bar y_g + (1-\lambda_g)\mu = (4.0, 2.8, 2.188)$. Interpret: the small geo ($n=4$, the "Des Moines" case) is pulled $2.0$ units toward the national mean; the large geo ($n=64$) barely moves ($0.012$). Shrinkage is automatic and proportional to how little each geo's data can be trusted.

- [ ] **Step 2: Example (b) — James–Stein beats the MLE**

$G = 6$ group means, $\sigma^2 = 1$, true $\theta = (0.5,-0.3,0.8,-0.6,0.2,0.4)$, observed $z = (1.4,-1.1,1.9,-0.2,1.3,-0.7)$. Compute $\lVert z\rVert^2 = 9.0$, the shrink factor $1 - (G-2)\sigma^2/\lVert z\rVert^2 = 1 - 4/9 = 5/9 \approx 0.556$, and $\hat\theta^{\text{JS}} = (5/9)z = (0.778,-0.611,1.056,-0.111,0.722,-0.389)$. Compare total squared error vs the truth: MLE $\lVert z-\theta\rVert^2 = 5.24$ versus JS $\lVert\hat\theta^{\text{JS}}-\theta\rVert^2 = 1.373$. JS reduces total error by ~74% on this realization — a concrete instance of the Rung 3 dominance theorem.

- [ ] **Step 3: Verify both, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
ns=np.array([4,16,64.]); yb=np.array([6,3,2.2]); s2,t2,mu=4.,1.,2.
lam=(ns/s2)/(ns/s2+1/t2); th=lam*yb+(1-lam)*mu
print('(a)',np.round(lam,4),np.round(th,4))
z=np.array([1.4,-1.1,1.9,-0.2,1.3,-0.7]); tr=np.array([.5,-.3,.8,-.6,.2,.4]); G=6
js=(1-(G-2)/(z@z))*z
print('(b)',z@z,round(5/9,4),round(((z-tr)**2).sum(),4),round(((js-tr)**2).sum(),4))"
```
Expected: `(a) [0.5 0.8 0.9412] [4. 2.8 2.1882]`, `(b) 9.0 0.5556 5.24 1.3733`.

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Worked Examples — shrinkage and James-Stein"
```

---

### Task 5: Code Tie-in

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (the `## Code Tie-in` section)

- [ ] **Step 1: Write a single self-contained ```` ```{python} ```` cell**

Model the structure and prose-wrapping on Chapter 5's Code Tie-in. Intro paragraph, one ```` ```{python} ```` block (numpy/scipy/matplotlib, seeded `rng = np.random.default_rng(...)`), closing paragraph quoting actual printed numbers. The cell must, with small helper functions:
1. **Simulate a geo-MMM:** $G$ geos (e.g. 20) with **unequal** week counts $n_g$ (some small, some large), true per-geo scalar effects $\theta_g\sim\mathcal N(\mu,\tau^2)$ with known $\mu,\tau^2$; generate group means $\bar y_g\sim\mathcal N(\theta_g,\sigma^2/n_g)$.
2. **Fit three ways:** **complete pooling** ($\hat\theta_g=$ grand weighted mean for all $g$), **no pooling** ($\hat\theta_g=\bar y_g$), **partial pooling** (empirical-Bayes: estimate $\mu$ by the precision-weighted mean and $\tau^2$ by a moment estimator on the between-geo spread net of sampling variance — clamp at 0 — then apply $\lambda_g$ shrinkage). Print the three total RMSEs vs the true $\theta_g$; partial pooling should be lowest.
3. **Shrinkage plot:** plot the no-pooling estimates and the partial-pooled estimates on a number line / two-column plot with connecting segments, showing small-$n$ geos pulled hardest toward the grand mean. Mark $\mu$.
4. **Gibbs preview (exact check):** with $\mu,\tau^2,\sigma^2$ held at their **known** values, run a short Gibbs/Monte-Carlo over the scalar conditionals $\theta_g\mid\text{rest}\sim\mathcal N(\lambda_g\bar y_g+(1-\lambda_g)\mu,\ \lambda_g\sigma^2/n_g)$, average the draws, and **assert** the Gibbs posterior means match the closed-form shrinkage estimates: `np.allclose(gibbs_mean, lambda*ybar+(1-lambda)*mu, atol=...)`. Print the confirmation. (Note in prose: the full Part III sampler additionally samples $\mu,\tau^2,\sigma^2$.)

End figures with `plt.show()`. The closing paragraph quotes the three RMSEs and the Gibbs-match confirmation.

- [ ] **Step 2: Extract the cell and run headless**

Save the cell body to `/tmp/ch6_code.py` and run:
```bash
MPLBACKEND=Agg python3 /tmp/ch6_code.py
```
Expected: runs top-to-bottom, prints partial-pooling RMSE lowest and the Gibbs-match `True`. Make the prose match the actual printed numbers.

- [ ] **Step 3: Commit**

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Code Tie-in — three-way pooling, shrinkage plot, Gibbs preview"
```

---

### Task 6: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (the `## Exercises` section with the four `###` tier headings)

Match Chapter 5's exercise style and heading text exactly: `### C -- Conceptual / Reading Comprehension`, `### B -- By Hand`, `### P -- Prove It`, `### A -- Applied / Code`. Each problem fully self-contained (state all data inline). **No links to solutions.**

- [ ] **Step 1: Tier C (3 problems)**
- **C1.** The three pooling stances — define complete / no / partial pooling, state when each is appropriate, and explain how $\tau^2\to0$ and $\tau^2\to\infty$ recover the two extremes.
- **C2.** Interpret the shrinkage weight $\lambda_g$: what drives it toward 0 vs 1, why small-$n$ geos shrink hardest, and why this is the principled answer to "trust it, throw it out, or shrink it?"
- **C3.** Why each full conditional being conjugate (Gaussian for $\beta_g$ and $\mu_\beta$, Inverse-Gamma/Wishart for the variances) is exactly what makes Gibbs sampling — and hence Part III — possible, even though the joint posterior has no closed form.

- [ ] **Step 2: Tier B (3 problems)**
- **B1.** Given $\sigma^2$, $\tau^2$, $\mu$, and a geo's $n_g$, $\bar y_g$ (clean numbers, e.g. $\sigma^2=9,\tau^2=1,\mu=0,n_g=9,\bar y_g=6$ → $\lambda=0.5$, $\hat\theta=3$), compute $\lambda_g$ and the partial-pooled estimate.
- **B2.** Given a small $z$ vector ($G\ge3$, clean $\lVert z\rVert^2$) and $\sigma^2$, compute the James–Stein estimate and its shrink factor, and (given a truth $\theta$) the MLE vs JS total squared error.
- **B3.** Given $X_g^\top X_g$, $X_g^\top y_g$, $\sigma^2$, $\mu_\beta$, $\Sigma_\beta$ (small $2\times2$), compute the $\beta_g$ full-conditional mean and covariance $V_g=(X_g^\top X_g/\sigma^2+\Sigma_\beta^{-1})^{-1}$, $m_g=V_g(X_g^\top y_g/\sigma^2+\Sigma_\beta^{-1}\mu_\beta)$.

- [ ] **Step 3: Tier P (3 problems)**
- **P1.** Prove the partial-pooling posterior-mean shrinkage formula $\hat\theta_g=\lambda_g\bar y_g+(1-\lambda_g)\mu$ with the stated $\lambda_g$ (complete the square / cite Chapter 5 Rung 3).
- **P2.** Prove James–Stein dominance for $G\ge3$ via Stein's unbiased risk estimate (the divergence computation giving $G\sigma^2-(G-2)^2\sigma^4\,\mathbb E[1/\lVert z\rVert^2]$).
- **P3.** Prove the $\beta_g$ full conditional is Gaussian with the stated $V_g,m_g$ (Chapter 5 Rung 4 reused with prior mean $\mu_\beta$).

- [ ] **Step 4: Tier A (2 problems)**
- **A1.** The three-way pooling simulation: simulate a geo-MMM with unequal $n_g$, fit complete/no/partial pooling, produce the shrinkage plot, and compare total RMSE; then sweep $\tau^2$ (or the small geo's $n_g$) and show the shrinkage strength change.
- **A2.** Implement the scalar hierarchical **Gibbs sampler** (cycle the Rung 6 conditionals; with known hyperparameters, or sampling $\mu,\tau^2$ too), and show its posterior means match the closed-form partial-pooling estimates.

- [ ] **Step 5: Commit**

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Exercises — C/B/P/A tiers"
```

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (append `## Summary` as the final section)

- [ ] **Step 1: Write `## Summary`**

Match Chapter 5's Summary format exactly: one-line lead, then **Key concepts** (bold lead-ins) and **Key identities** (**inline math only**). Cover:
- **Key concepts:** the three pooling stances; partial pooling = precision-weighted shrinkage; the shrinkage weight $\lambda_g$ driven by $n_g$; James–Stein dominance (shrinkage provably beats the MLE for $G\ge3$); $\tau^2$ as the learned pooling dial; the vectorized hierarchical regression; conjugate full conditionals → Gibbs → Part III.
- **Key identities (inline):**
  - Shrinkage weight: $\lambda_g = (n_g/\sigma^2)/(n_g/\sigma^2 + 1/\tau^2)$.
  - Partial-pooled mean: $\hat\theta_g = \lambda_g\bar y_g + (1-\lambda_g)\mu$.
  - James–Stein: $\hat\theta^{\text{JS}} = (1 - (G-2)\sigma^2/\lVert z\rVert^2)\,z$, risk $G\sigma^2 - (G-2)^2\sigma^4\,\mathbb E[1/\lVert z\rVert^2] < G\sigma^2$.
  - $\beta_g$ full conditional: $V_g=(X_g^\top X_g/\sigma^2+\Sigma_\beta^{-1})^{-1}$, $m_g=V_g(X_g^\top y_g/\sigma^2+\Sigma_\beta^{-1}\mu_\beta)$.

- [ ] **Step 2: Commit**

```bash
git add parts/02-regression-bayes/03-hierarchical-regression.qmd
git commit -m "feat(ch06): Summary — key concepts and identities"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (insert before the closing `:::` at line 709)

- [ ] **Step 1: Append the Chapter 6 solutions block**

Immediately before the final `:::` that closes the `show-solutions` visible div, add a `## Chapter 6 — Hierarchical Models` heading followed by full solutions to every exercise from Task 6 (C1–C3, B1–B3, P1–P3, A1–A2). Match the format of the existing `## Chapter 5 — Bayesian Inference` block (`appendix/solutions.qmd:481-708`): em-dash sub-headings `### C — …` etc., bold labels (`**C1.**`), full worked math for B and P, plain ```` ```python ```` reference sketches (NOT `{python}`) plus a qualitative takeaway for A. Each proof ends `$\blacksquare$`. Verify every numeric answer against NumPy (B1 $\lambda=0.5,\hat\theta=3$; B2 the JS numbers; B3 the $2\times2$ full-conditional).

- [ ] **Step 2: Verify by-hand answers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
# B1
s2,t2,mu,n,yb=9.,1.,0.,9,6.; lam=(n/s2)/(n/s2+1/t2); print('B1',lam, lam*yb+(1-lam)*mu)
# B3 example: XtXg=[[4,1],[1,3]], Xtyg=[5,2], sigma2=1, mu_beta=[0,0], Sigma_beta=diag(1,1)
XtX=np.array([[4,1],[1,3.]]); Xty=np.array([5,2.]); Sinv=np.eye(2)
V=np.linalg.inv(XtX/1.0+Sinv); m=V@(Xty/1.0+Sinv@np.zeros(2))
print('B3 V',np.round(V,4).tolist(),'m',np.round(m,4))"
```
Expected: `B1 0.5 3.0`; B3 prints a valid SPD $V$ and mean $m$ (use whatever clean $X_g^\top X_g,X_g^\top y_g$ the solution states — keep the chapter exercise and solution numbers identical).

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch06): appendix solutions for Hierarchical Models exercises"
```

---

### Task 9: Final whole-chapter review

**Files:**
- Review: `parts/02-regression-bayes/03-hierarchical-regression.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: KaTeX / structure lint**

```bash
cd /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch6-hierarchical-models
grep -n 'begin{align}' parts/02-regression-bayes/03-hierarchical-regression.qmd || echo "no bare align — good"
grep -c '^\$\$' parts/02-regression-bayes/03-hierarchical-regression.qmd   # expect even
grep -nE '^\## ' parts/02-regression-bayes/03-hierarchical-regression.qmd  # expect 6 template headings in order
grep -nE '^# ' parts/02-regression-bayes/03-hierarchical-regression.qmd    # expect '# Hierarchical Models'
grep -c '^:::$' appendix/solutions.qmd   # expect 2 (file still closes cleanly)
```
Confirm no bare `\begin{align}`, even `$$` count, H1 retitled, the six template headings in order, citations only the four valid keys.

- [ ] **Step 2: Re-run the code cell headless**

```bash
MPLBACKEND=Agg python3 /tmp/ch6_code.py && echo "CODE CELL OK"
```
Expected: `CODE CELL OK`, partial-pooling RMSE lowest, Gibbs match `True`.

- [ ] **Step 3: Through-line check (manual)**

Confirm every theory subsection ends on the three-stances / shrinkage / Part III through-line, the three named proofs + the full-conditional derivations are present and correct, and the chapter reuses Chapter 5 Rung 3/Rung 4 by reference. Fix gaps inline and commit if changed.

---

## Self-Review (completed by plan author)

- **Spec coverage:** retitle + Motivation (T1); rungs 1–3 with both keystone proofs (T2); rungs 4–6 incl. Gibbs full-conditional derivation (T3); both worked examples (T4); code tie-in with three-way pooling + shrinkage plot + Gibbs exact check (T5); C/B/P/A exercises (T6); auto-included Summary (T7); appendix solutions (T8); render/lint gate (T9). All spec success criteria mapped.
- **Placeholder scan:** none — every task carries the actual formulas and NumPy-verified anchor numbers.
- **Consistency:** symbols ($\lambda_g$, $\hat\theta_g$, $\mu_\beta$, $\Sigma_\beta$, $V_g$, $m_g$) identical across tasks; the shrinkage formula and the $\beta_g$ full conditional are stated identically in T2/T3/T6/T7; worked-example numbers are NumPy-verified and reused. The $\beta_g$ full conditional ($V_g,m_g$) is defined in T3 and reused unchanged in T6/T7/T8.
- **Note:** the real render gate is CI `quarto render` (HTML + PDF) on the PR — Quarto is not installed locally, so local verification is the Python cell + lint checks above.
