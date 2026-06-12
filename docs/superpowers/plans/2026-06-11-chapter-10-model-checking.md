# Chapter 10 — Model Checking & Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 10 — Model Checking & Calibration, the final chapter of Part III: comprehensive MCMC convergence diagnostics ($\hat{R}$, ESS, MCSE, divergences) and model calibration/checking (prior/posterior predictive checks, SBC keystone, ELPD/WAIC/PSIS-LOO), anchored on the canonical eight-schools model adapted to a geo-level MMM, shipped as a PR whose CI `quarto render` passes.

**Architecture:** Single `{python}` cell, NumPy/Matplotlib only, reusing the from-scratch HMC/leapfrog sampler from Ch9 (no PyMC/Stan). Three full proofs ($\hat{R}\to1$, ESS variance-inflation, SBC rank-uniformity) plus a compact PIT proof. Fixed chapter template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary).

**Tech Stack:** Quarto `.qmd`, KaTeX math, Python 3 (NumPy, Matplotlib, `MPLBACKEND=Agg`), `references.bib`.

---

## Conventions (enforce in every task)

- **Worktree:** all work in `~/.config/superpowers/worktrees/mmm_foundational_text/ch10-model-checking`. **Use explicit `git -C <worktree>` for every git command** (a bare `cd` in a compound command may not apply and would send commits to main, where a hook blocks them). Identity already set (`jlh530i` / `jlh530i@gmail.com`).
- **KaTeX:** use `aligned` inside `$$…$$`; **never** bare `\begin{align}`. `$$` delimiters on their own lines. Even `$` count per file. Inline math only in Summary.
- **Chapter file:** `parts/03-sampling/04-model-checking.qmd`. Keep H1 `# Model Checking & Calibration`. Replace the stub body; keep the canonical-anchors italic line.
- **Citations:** `@gelman2013`, `@talts2018`, `@gelmanhill2006`, `@mcelreath2020`, `@betancourt2017`, `@hoffman2014`, plus new `@vehtari2021`, `@vehtari2017` (added in Task 9).
- **Verification:** the single `{python}` cell must run top-to-bottom under `MPLBACKEND=Agg python3`; figures end `plt.show()`. Real gate: CI `quarto render` (HTML+PDF) green on the PR (Quarto not installed locally).
- **Voice/structure exemplars (read-only):** `parts/03-sampling/03-hmc-nuts.qmd` (Ch9, immediate predecessor), `parts/02-regression-bayes/03-hierarchical-regression.qmd` (Ch6, the hierarchy this adapts).

## Verified numeric anchors (NumPy-checked)

- **$\hat{R}$ worked:** 4 chains, length $N=4$, values
  chain1 `[3,4,5,4]`, chain2 `[4,5,6,5]`, chain3 `[2,3,4,3]`, chain4 `[5,6,7,6]`.
  Chain means `[4,5,3,6]`, grand mean `4.5`. $B = \frac{N}{M-1}\sum(\bar\theta_m-\bar\theta)^2 = 6.6667$, $W = \frac1M\sum s_m^2 = 0.6667$, $\widehat{\text{var}}^+ = \frac{N-1}{N}W+\frac1N B = 2.1667$, $\hat{R}=\sqrt{2.1667/0.6667}=\mathbf{1.8028}$ (≫ 1.01 → not converged).
- **ESS worked:** AR(1) with $\phi=0.5$ → $\rho_k=0.5^k$; exact integrated factor $1+2\sum_{k\ge1}\phi^k = 1+2\cdot\frac{\phi}{1-\phi}= \mathbf{3.0}$; truncated at 3 lags $=2.75$. With $N=1000$: $N_\text{eff}\approx 333$.
- **SBC single rank worked:** $\theta_0=0.7$, sorted posterior draws `[-0.2,0.1,0.3,0.55,0.8,0.9,1.1,1.3,1.6]` ($L=9$) → rank $=\#\{\text{draws}<\theta_0\}=\mathbf{4}$ (of 9; 10 possible bins $0..9$).
- **SBC uniformity check:** conjugate $\theta\sim\mathcal N(0,1)$, $y\mid\theta\sim\mathcal N(\theta,1)$, posterior $\mathcal N(y/2,1/2)$; over 60k sims with $L=9$, all 10 rank bins ≈ 6000 (range ~5892–6110) → uniform. (Implementer re-derives; exact counts seed-dependent.)
- **Eight-schools / geo-MMM HMC (Code Tie-in):** centered parameterization yields many divergences and degraded tail-ESS / elevated $\hat{R}$ on $\tau$; non-centered yields ≈0 divergences and $\hat{R}<1.01$. **Assert robust relationships** (as Ch9 did with `ess_hmc > 10*ess_rwm`): `div_centered > 10 * div_noncentered` and `rhat_tau_noncentered < 1.01`. Implementer pins exact counts via NumPy before authoring prose.

---

## Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (replace stub body below the H1 + anchors line; write the Motivation H2)

- [ ] **Step 1: Keep H1 and canonical-anchors line; write Motivation.**

Keep:
```markdown
# Model Checking & Calibration

*Canonical anchors: Talts et al.; BDA3.*
```

Write a `## Motivation` section (~350–500 words) covering:
- The closing promise of Ch9 paid off verbatim in spirit: HMC gives you draws, but draws are worthless until you answer two questions. Quote/paraphrase the Ch9 hand-off ("$\hat{R}$, effective sample size, and divergence counts — is the subject of Chapter 10").
- **The two driving questions:** (1) *Did the sampler converge?* — is the chain a faithful sample from the posterior and is Monte Carlo error small? (2) *Is the model any good?* — split into **B1 self-consistency** (does it reproduce the data; are intervals honest?) and **B2 out-of-sample prediction** (does it predict held-out data, and beat alternatives?).
- Why both matter in an MMM: a converged-but-wrong model gives confident, precise, and useless budget recommendations to Part IV; an uncalibrated model's credible intervals on ROAS are dishonest.
- Introduce the running example: name the **canonical eight-schools model** (the model the diagnostics literature was built on) and promise the quick MMM adaptation in the Theory section.
- One-paragraph roadmap of the chapter's two parts.

- [ ] **Step 2: Lint + commit.**

Check even `$` count, no bare `\begin{align}`. Then:
```bash
WT=~/.config/superpowers/worktrees/mmm_foundational_text/ch10-model-checking
git -C "$WT" add parts/03-sampling/04-model-checking.qmd
git -C "$WT" commit -m "feat(ch10): Motivation — two questions for trusting a posterior"
```

---

## Task 2: Theory Part A — MCSE and $\hat{R}$ (with full proof)

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (begin `## Theory & Proofs`; rungs 1–2)

- [ ] **Step 1: Write the canonical model + MMM adaptation bridge (opens Theory).**

Before the diagnostics, state the running model once:
- **Eight schools (@gelman2013):** $y_j\mid\theta_j\sim\mathcal N(\theta_j,\sigma_j^2)$, $\theta_j\mid\mu,\tau\sim\mathcal N(\mu,\tau^2)$, $j=1..J$, weak priors on $\mu,\tau$.
- **MMM adaptation (one paragraph, "same math, marketing labels"):** $\theta_g$ = geo-level effect (incremental ROAS / baseline sales for region $g$), pooled toward national mean $\mu$ with cross-geo SD $\tau$; $y_g$ is that geo's noisy regression estimate with standard error $\sigma_g$. Note the funnel trap: **few geos ⇒ $\tau$ weakly identified ⇒ funnel ⇒ divergences** (forward-reference Ch9).

- [ ] **Step 2: Rung 1 — Monte Carlo standard error (MCSE).**

- Posterior summaries are themselves estimated from finite draws and carry sampling error.
- For an ergodic average (Ch8 CLT), $\widehat{\mathbb E}[\theta]=\frac1N\sum\theta^{(i)}$ has $\text{MCSE}=\sigma_\theta/\sqrt{N_\text{eff}}$, where $N_\text{eff}$ is the effective sample size (defined next rung). Contrast $N$ vs $N_\text{eff}$: autocorrelation inflates the error.
- Practical rule: report MCSE alongside any posterior summary; an MMM ROAS reported to 3 sig-figs with MCSE in the 2nd is false precision.

- [ ] **Step 3: Rung 2 — $\hat{R}$ (Gelman–Rubin) with full proof.**

Define for $M$ chains of length $N$:
$$
\bar\theta_m=\tfrac1N\textstyle\sum_i\theta_m^{(i)},\quad
\bar\theta=\tfrac1M\textstyle\sum_m\bar\theta_m,\quad
B=\tfrac{N}{M-1}\textstyle\sum_m(\bar\theta_m-\bar\theta)^2,\quad
W=\tfrac1M\textstyle\sum_m s_m^2 .
$$
$$
\widehat{\operatorname{var}}^+(\theta)=\tfrac{N-1}{N}W+\tfrac1N B,\qquad
\hat{R}=\sqrt{\widehat{\operatorname{var}}^+/W}.
$$

**Theorem.** If the chains are stationary and have converged to a common target with finite variance, then $\hat{R}\to 1$ as $N\to\infty$; if the chains have not mixed (different stationary locations), $\hat R>1$.

**Proof (full, ends `\blacksquare`).** Use `aligned` inside `$$`:
- At convergence each chain samples the same marginal with variance $\operatorname{var}(\theta)$. $W$ is the mean within-chain sample variance, a consistent estimator of $\operatorname{var}(\theta)$, so $W\to\operatorname{var}(\theta)$.
- $B/N=\frac{1}{M-1}\sum(\bar\theta_m-\bar\theta)^2$ estimates the variance of the chain means. For converged chains the means are i.i.d. with variance $\operatorname{var}(\theta)/N$, so $B\to\operatorname{var}(\theta)$ as well and $B/N\to 0$.
- Hence $\widehat{\operatorname{var}}^+=\frac{N-1}{N}W+\frac1N B\to\operatorname{var}(\theta)$ and $\hat R=\sqrt{\widehat{\operatorname{var}}^+/W}\to\sqrt{\operatorname{var}(\theta)/\operatorname{var}(\theta)}=1$.
- If chains have not mixed, the means differ systematically: $B$ does not shrink ($B/N\not\to0$ relative scale), $\widehat{\operatorname{var}}^+>W$ (the overdispersed-initialization estimator overestimates), so $\hat R>1$. $\blacksquare$

Note threshold: modern practice $\hat R<1.01$.

- [ ] **Step 4: Rung 2 continued — split-$\hat{R}$ and rank-normalized $\hat{R}$.**

- **Split-$\hat{R}$:** halve each chain and treat halves as separate chains; catches within-chain trends/non-stationarity a whole-chain $\hat R$ misses. (@vehtari2021)
- **Rank-normalized $\hat{R}$:** replace draws by normal scores of their ranks before computing $\hat R$; makes the diagnostic valid for heavy-tailed / infinite-variance targets where the variance-based $\hat R$ is undefined. (@vehtari2021)

- [ ] **Step 5: Lint + commit.**

Even `$`/`$$` counts, no bare align, `\blacksquare` present. Then `git -C "$WT" add … && git -C "$WT" commit -m "feat(ch10): Theory rungs 1-2 — MCSE and R-hat (proof)"`.

---

## Task 3: Theory Part A — ESS (full derivation) and HMC signals

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (rungs 3–4)

- [ ] **Step 1: Rung 3 — Effective sample size with full derivation.**

**Theorem (variance inflation).** For a stationary sequence with lag-$k$ autocorrelation $\rho_k$, the variance of the mean of $N$ draws is
$$
\operatorname{var}(\bar\theta)=\frac{\sigma^2}{N}\Big(1+2\textstyle\sum_{k=1}^{N-1}(1-\tfrac{k}{N})\rho_k\Big)\xrightarrow{N\to\infty}\frac{\sigma^2}{N}\Big(1+2\textstyle\sum_{k\ge1}\rho_k\Big),
$$
so the **effective sample size** is $N_\text{eff}=\dfrac{N}{1+2\sum_{k\ge1}\rho_k}$.

**Proof (full, ends `\blacksquare`).** Expand $\operatorname{var}(\bar\theta)=\frac1{N^2}\sum_{i,j}\operatorname{cov}(\theta^{(i)},\theta^{(j)})$, use stationarity $\operatorname{cov}=\sigma^2\rho_{|i-j|}$, count the $N$ diagonal and $2(N-k)$ off-diagonal terms at lag $k$, factor $\sigma^2/N$. (This extends Ch8's autocorrelation result — reference it.)

Then:
- **Integrated autocorrelation time** $\tau_\text{int}=1+2\sum\rho_k$; $N_\text{eff}=N/\tau_\text{int}$. Worked AR(1) sanity in prose: $\phi=0.5\Rightarrow\tau_\text{int}=3$.
- **Bulk-ESS vs tail-ESS (@vehtari2021):** bulk-ESS (rank-normalized) governs central estimates like the posterior mean; tail-ESS (from 5%/95% indicator series) governs interval endpoints. A chain can have fine bulk-ESS but terrible tail-ESS → unreliable credible intervals. Why one number is not enough.
- Connect back: $\text{MCSE}=\sigma/\sqrt{N_\text{eff}}$ from Rung 1 now fully defined.

- [ ] **Step 2: Rung 4 — HMC-specific signals (from Ch9).**

$\hat R$ and ESS are necessary but not sufficient — they can all look fine while a region of the posterior is never visited. Targeted HMC diagnostics:
- **Divergences:** a divergent transition (energy error blows up, Ch9) flags a region of high curvature the integrator cannot resolve; even a few divergences signal biased exploration (the funnel). Not a convergence statistic — a *geometry* alarm.
- **BFMI** (energy Bayesian fraction of missing information): low BFMI ⇒ the momentum resampling cannot move between energy levels efficiently ⇒ slow exploration of the tails.
- **Tree-depth saturation** (NUTS): hitting `max_treedepth` routinely ⇒ trajectories truncated before the U-turn ⇒ inefficiency (not bias).
- Synthesis: divergences = bias alarm; BFMI / tree-depth = efficiency alarms.

- [ ] **Step 3: Lint + commit.** As before; commit `-m "feat(ch10): Theory rungs 3-4 — ESS (derivation) and HMC diagnostics"`.

---

## Task 4: Theory Part B1 — Prior/Posterior predictive checks (+ PIT proof)

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (rungs 5–6)

- [ ] **Step 1: Rung 5 — Prior predictive checks.**

- Draw $\theta\sim p(\theta)$, $\tilde y\sim p(y\mid\theta)$: do simulated datasets look remotely plausible (right order of magnitude, sign, range)? Catches absurd priors before any fitting.
- MMM example in prose: a prior that implies a single channel can drive 10× total sales is caught here, not after a week of sampling.

- [ ] **Step 2: Rung 6 — Posterior predictive checks (PPC) with Bayesian p-value + PIT proof.**

- Replicated data: $p(y^\text{rep}\mid y)=\int p(y^\text{rep}\mid\theta)\,p(\theta\mid y)\,d\theta$. Draw $\theta^{(s)}\sim$ posterior, then $y^{\text{rep},(s)}\sim p(\cdot\mid\theta^{(s)})$.
- Test quantity $T(y,\theta)$ (e.g. max, min, a discrepancy); **posterior predictive p-value** $p_B=\Pr\!\big(T(y^\text{rep},\theta)\ge T(y,\theta)\mid y\big)$, estimated as the fraction of replicates exceeding the observed statistic.
- **Proposition (PIT / uniformity).** If $T(y)$ is a continuous statistic and the model is correct (data truly drawn from the model), then the predictive p-value (equivalently the PIT value $F(T(y))$) is $\text{Uniform}(0,1)$.
- **Proof (compact, ends `\blacksquare`).** For a continuous CDF $F$, the random variable $U=F(T)$ satisfies $\Pr(U\le u)=\Pr(T\le F^{-1}(u))=F(F^{-1}(u))=u$, i.e. $U\sim\text{Uniform}(0,1)$ (the probability integral transform). Hence under a correct model the p-value is uniform; systematic deviation from uniform flags misfit. $\blacksquare$
- **Caveat:** PPC p-values using the data twice (to fit and to test) are *conservative* — they tend toward 0.5; SBC (next rung) repairs this by working in the prior-predictive space. Forward-reference.

- [ ] **Step 3: Lint + commit.** Commit `-m "feat(ch10): Theory rungs 5-6 — prior/posterior predictive checks (PIT proof)"`.

---

## Task 5: Theory Part B1 — SBC keystone (full proof)

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (rung 7, keystone)

- [ ] **Step 1: Rung 7 — Simulation-Based Calibration, the keystone.**

Construction (one parameter $\theta$ for clarity, vector case analogous):
1. Draw $\theta_0\sim p(\theta)$ (a "ground truth" from the prior).
2. Draw a dataset $y\sim p(y\mid\theta_0)$.
3. Draw $L$ posterior samples $\theta_1,\dots,\theta_L\sim p(\theta\mid y)$ (using the *same* inference you want to validate).
4. Compute the rank $r=\#\{\ell:\theta_\ell<\theta_0\}\in\{0,1,\dots,L\}$.

**Theorem (rank uniformity, @talts2018).** If the posterior samples are drawn from the true posterior $p(\theta\mid y)$, then over the joint draw of $(\theta_0,y,\theta_{1:L})$ the rank $r$ is discrete-uniform on $\{0,\dots,L\}$.

**Proof (full, ends `\blacksquare`).** Use `aligned`:
- Key identity — the **data-averaged posterior equals the prior**:
  $$
  \int p(\theta\mid y)\,p(y)\,dy=\int p(\theta\mid y)\!\left(\int p(y\mid\theta')p(\theta')\,d\theta'\right)dy=p(\theta),
  $$
  because $p(\theta\mid y)p(y)=p(\theta,y)$ and $\int p(\theta,y)\,dy=p(\theta)$.
- Therefore $\theta_0$ (a prior draw) and each $\theta_\ell$ (a draw from $p(\theta\mid y)$ with $y$ from $\theta_0$) are, marginally over the construction, draws from the *same* distribution — the prior. So $\theta_0,\theta_1,\dots,\theta_L$ are $L+1$ exchangeable (i.i.d.-in-distribution) draws.
- The rank of one exchangeable value among $L+1$ is uniform on $\{0,\dots,L\}$ (no value is special). $\blacksquare$

- [ ] **Step 2: Reading rank histograms.**

Aggregate ranks over many simulations and histogram into $L+1$ bins:
- **Uniform** ⇒ calibrated inference.
- **∪-shaped** (mass at extremes) ⇒ posterior too narrow / overconfident (true value lands in the tails too often).
- **∩-shaped** (mass in middle) ⇒ posterior too wide / under-confident.
- **Left/right slope** ⇒ posterior biased low/high.
- **Tie-back to Ch9:** a *centered* funnel parameterization fails SBC (biased/∪ ranks on $\tau$) because the sampler misses the neck; the *non-centered* form passes. SBC is how you'd *catch* the Ch9 pathology automatically.

- [ ] **Step 3: Lint + commit.** Commit `-m "feat(ch10): Theory rung 7 — SBC rank-uniformity (keystone proof)"`.

---

## Task 6: Theory Part B2 — Out-of-sample prediction + synthesis checklist

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (rungs 8–9)

- [ ] **Step 1: Rung 8 — ELPD, WAIC, PSIS-LOO (stated precisely, not proved).**

- The goal of B2: estimate **expected log pointwise predictive density** for new data, $\text{elpd}=\sum_i\int p_t(\tilde y_i)\log p(\tilde y_i\mid y)\,d\tilde y_i$ — a proper scoring rule rewarding calibration *and* sharpness.
- **lppd** (in-sample) $=\sum_i\log\frac1S\sum_s p(y_i\mid\theta^{(s)})$ overfits; correct it.
- **WAIC** $=\text{lppd}-p_\text{WAIC}$, with effective-parameter penalty $p_\text{WAIC}=\sum_i\operatorname{var}_s\big(\log p(y_i\mid\theta^{(s)})\big)$. (@gelman2013, @vehtari2017)
- **PSIS-LOO:** leave-one-out elpd via importance sampling with weights $\propto 1/p(y_i\mid\theta^{(s)})$, Pareto-smoothed for stability; the **Pareto-$\hat k$** per point is a *reliability diagnostic* ($\hat k>0.7$ ⇒ estimate untrustworthy for that point). (@vehtari2017)
- State the **WAIC↔LOO asymptotic equivalence**; emphasize these answer "which model predicts better" (B2), distinct from "is this model self-consistent" (B1). MMM use in prose: comparing a model with vs without a competitor-spend control, or two adstock forms.

- [ ] **Step 2: Rung 9 — Synthesis: a model-checking workflow for an MMM.**

Ordered checklist with what each step catches:
1. **Prior predictive** → absurd priors.
2. **Fit**, then **convergence diagnostics** ($\hat R<1.01$ split + rank-normalized, bulk/tail-ESS adequate, MCSE small, **zero divergences**) → sampler trustworthy.
3. **Posterior predictive checks** → structural misfit (e.g. can't reproduce seasonal sales peaks).
4. **SBC on the pipeline** → mis-calibration / coding bugs in the inference.
5. **PSIS-LOO** → compare candidate specifications.
- Green/red light meaning for Part IV: only a model that clears 1–4 should feed budget optimization; a red light on divergences or SBC means the ROAS posterior is untrustworthy and optimization would amplify the error.

- [ ] **Step 3: Lint + commit.** Commit `-m "feat(ch10): Theory rungs 8-9 — out-of-sample prediction and MMM checklist"`.

---

## Task 7: Worked Examples

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (`## Worked Examples`)

- [ ] **Step 1: Example (a) — $\hat{R}$ by hand.**

Use the verified anchor: 4 chains length 4, values `[3,4,5,4] / [4,5,6,5] / [2,3,4,3] / [5,6,7,6]`. Show chain means `[4,5,3,6]`, grand mean `4.5`, then $B=6.6667$, $W=0.6667$, $\widehat{\operatorname{var}}^+=2.1667$, $\hat R=\sqrt{2.1667/0.6667}=1.8028$. Interpret: $\hat R\approx1.80\gg1.01$ — the four chains sit in different places; **not converged**. Note this is the classic signature of an unmixed geo-hierarchy with a poorly identified $\tau$.

- [ ] **Step 2: Example (b) — ESS / variance-inflation by hand.**

AR(1) with $\phi=0.5$, so $\rho_k=0.5^k$. Truncating at 3 lags: $1+2(0.5+0.25+0.125)=2.75$; exact $1+2\cdot\frac{0.5}{0.5}=3.0$. With $N=1000$ draws, $N_\text{eff}\approx1000/3\approx333$ and $\text{MCSE}=\sigma/\sqrt{333}\approx0.055\,\sigma$. Lesson: 1000 correlated draws are worth ~333 independent ones here.

- [ ] **Step 3: Example (c) — one SBC rank + a small histogram read.**

Ground truth $\theta_0=0.7$; $L=9$ posterior draws (sorted) `[-0.2,0.1,0.3,0.55,0.8,0.9,1.1,1.3,1.6]`. Rank $=\#\{<0.7\}=4$ (of 9; bins $0..9$). Then show a sketch of two aggregated histograms: a flat one (calibrated) vs a ∪-shaped one (overconfident), and state the diagnosis for each.

- [ ] **Step 4: Lint + commit.** Commit `-m "feat(ch10): Worked Examples — R-hat, ESS, SBC rank by hand"`.

---

## Task 8: Code Tie-in

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (`## Code Tie-in`, single `{python}` cell)

- [ ] **Step 1: Build the cell (NumPy/Matplotlib only; reuse Ch9 leapfrog HMC).**

The cell, top-to-bottom, on the eight-schools / geo-MMM model:
1. **Model + data:** $J=8$ geos, fixed $\sigma_g$ and observed $y_g$ (use the classic eight-schools values, framed as geo estimates). Define `log_post` and `grad` for both **centered** ($\theta,\mu,\tau$ with $\log\tau$) and **non-centered** ($\tilde\theta,\mu,\log\tau$) parameterizations.
2. **Sampler:** lift the from-scratch leapfrog HMC from Ch9 (resample momentum, $L$ leapfrog steps, Metropolis accept on $\Delta H$, flag a **divergence** when $|\Delta H|>1000$ or energy non-finite). Run multiple chains for each parameterization.
3. **Diagnostics from raw draws (implement, don't import):**
   - `split_rhat(chains)` → report for $\mu$ and $\tau$.
   - `ess(chain)` via autocorrelation sum (bulk) → report.
   - `mcse = sd/sqrt(ess)`.
   - divergence counts per parameterization.
4. **Assertions (robust, seed-independent):**
   - `assert div_centered > 10 * max(div_noncentered, 1)`
   - `assert rhat_tau_noncentered < 1.01`
   - print all diagnostics in a small table.
5. **Posterior predictive check:** draw $y^\text{rep}$ from non-centered posterior, test quantity = sample max (or between-geo SD); compute Bayesian p-value; print it.
6. **SBC loop:** small loop (e.g. a conjugate one-parameter sub-model, or the hierarchical model with modest $L$ and sims for runtime) collecting ranks; plot the rank histogram (expect ≈uniform for correct/non-centered inference).
7. **Figures:** trace-plot or rank histogram + a centered-vs-non-centered comparison; each `plt.show()`.
8. **Adstock coda (short):** compute $\hat R$/ESS on one realistic adstock-decay parameter from a tiny synthetic AR-style fit, captioned "the same diagnostics on a real MMM parameter."

- [ ] **Step 2: Verify headless.**

Extract the cell and run:
```bash
WT=~/.config/superpowers/worktrees/mmm_foundational_text/ch10-model-checking
cd "$WT" && MPLBACKEND=Agg python3 /tmp/ch10_cell.py
```
Expected: assertions pass, diagnostics print, Bayesian p-value in $(0,1)$, SBC histogram roughly flat, no errors. Iterate until clean. **Pin the actual printed divergence counts / $\hat R$ / ESS into the surrounding prose.**

- [ ] **Step 3: Lint + commit.** Even `$`, fenced ```{python}``` cell present. Commit `-m "feat(ch10): Code Tie-in — diagnostics, PPC, SBC on eight-schools/geo-MMM"`.

---

## Task 9: references.bib additions

**Files:**
- Modify: `references.bib` (append two entries)

- [ ] **Step 1: Add the two keys.**

Append:
```bibtex
@article{vehtari2021,
  title={Rank-normalization, folding, and localization: An improved $\widehat{R}$ for assessing convergence of MCMC},
  author={Vehtari, Aki and Gelman, Andrew and Simpson, Daniel and Carpenter, Bob and B{\"u}rkner, Paul-Christian},
  journal={Bayesian Analysis},
  volume={16},
  number={2},
  pages={667--718},
  year={2021}
}

@article{vehtari2017,
  title={Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC},
  author={Vehtari, Aki and Gelman, Andrew and Gabry, Jonah},
  journal={Statistics and Computing},
  volume={27},
  number={5},
  pages={1413--1432},
  year={2017}
}
```

- [ ] **Step 2: Commit.** Commit `-m "feat(ch10): add Vehtari 2021/2017 references"`.

---

## Task 10: Exercises (written directly by the controller, not a subagent)

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (`## Exercises`)

- [ ] **Step 1: Write four tiers, self-contained (no inline solution links).**

- **C — Conceptual (3):** e.g. why $\hat R\approx1$ is necessary but not sufficient; what a ∪-shaped SBC histogram means; why a PPC p-value near 0.5 is not strong evidence of fit (data used twice).
- **B — By hand (3):** compute $\hat R$ from a given 3-chain set; compute $N_\text{eff}$ from given $\rho_k$; compute an SBC rank and diagnose a small histogram. (Use fresh numbers distinct from the Worked Examples; NumPy-verify each before writing.)
- **P — Prove it (2–3):** prove $\hat R\to1$ for converged chains (or the $B/N\to0$ step); prove the variance-inflation $1+2\sum\rho_k$; prove the data-averaged posterior identity underlying SBC.
- **A — Applied/code (2–3):** implement split-$\hat R$ and confirm it flags a trending chain whole-$\hat R$ misses; run a mini-SBC and plot ranks; compute WAIC/LOO-style lppd on a toy model and compare two specs.

Verify every "By hand" numeric answer with NumPy before writing.

- [ ] **Step 2: Lint + commit.** Commit `-m "feat(ch10): Exercises — C/B/P/A tiers"`.

---

## Task 11: Summary (written directly by the controller)

**Files:**
- Modify: `parts/03-sampling/04-model-checking.qmd` (`## Summary`, auto-included; inline math only)

- [ ] **Step 1: Write Key concepts + Key identities.**

- **Key concepts:** the two questions; MCSE; $\hat R$ (split, rank-normalized); bulk/tail ESS; divergences/BFMI/tree-depth; prior & posterior predictive checks; Bayesian p-value; SBC rank-uniformity; ELPD/WAIC/PSIS-LOO + Pareto-$\hat k$; the MMM workflow.
- **Key identities (inline math only):** $\hat R=\sqrt{\widehat{\operatorname{var}}^+/W}$; $\widehat{\operatorname{var}}^+=\frac{N-1}{N}W+\frac1N B$; $N_\text{eff}=N/(1+2\sum_k\rho_k)$; $\text{MCSE}=\sigma/\sqrt{N_\text{eff}}$; $p_B=\Pr(T(y^\text{rep},\theta)\ge T(y,\theta)\mid y)$; SBC rank uniform on $\{0,\dots,L\}$; $\text{WAIC}=\text{lppd}-p_\text{WAIC}$.

- [ ] **Step 2: Lint + commit.** Inline math only (no `$$`). Commit `-m "feat(ch10): Summary — key concepts and identities"`.

---

## Task 12: Appendix solutions block

**Files:**
- Modify: `appendix/solutions.qmd` (insert a Chapter 10 block before the final closing `:::`)

- [ ] **Step 1: Append the solutions block.**

Insert, **before the file's final `:::`**, a block:
```markdown
## Chapter 10 — Model Checking & Calibration

[full solutions to every C/B/P/A exercise, matching Task 10 numbering]
```
The whole appendix content is already wrapped in `::: {.content-visible when-meta="show-solutions"}`; insert inside it (before the last `:::`). Provide complete worked solutions (NumPy-verified for numeric ones). No inline links back.

- [ ] **Step 2: Lint + commit.** Even `$` count in appendix; file still closes with a single `:::`; Ch10 heading present (em-dash). Commit `-m "feat(ch10): appendix solutions for Model Checking & Calibration"`.

---

## Task 13: Final review + PR + CI

- [ ] **Step 1: Chapter lint.** Verify: H1 `# Model Checking & Calibration`; six template H2 headings in order (Motivation, Theory & Proofs, Worked Examples, Code Tie-in, Exercises, Summary); no bare `\begin{align}`; even `$` count; `\blacksquare` proofs present (≥3: $\hat R$, ESS, SBC; plus PIT); citation keys all resolve (`@vehtari2021`/`@vehtari2017` now in bib).
- [ ] **Step 2: Appendix lint.** Even `$` count; closes with single `:::`; `## Chapter 10 — Model Checking & Calibration` present.
- [ ] **Step 3: Re-run the code cell headless** (`MPLBACKEND=Agg python3`) — assertions pass, no errors.
- [ ] **Step 4: Through-line check.** Ch9 hand-off paid off; eight-schools→geo-MMM bridge present; B1/B2 framing consistent; Part IV forward-reference present.
- [ ] **Step 5: Push + PR.**
```bash
WT=~/.config/superpowers/worktrees/mmm_foundational_text/ch10-model-checking
git -C "$WT" push -u origin ch10-model-checking
gh pr create --title "Chapter 10 — Model Checking & Calibration" --body "..."
```
- [ ] **Step 6: Watch CI `quarto render` to a green conclusion** (HTML+PDF; deploy skipped on PR branch). Report PR link + green status. Then the user merges.

---

## Self-Review (controller, before dispatching Task 1)

- **Spec coverage:** MCSE ✓(T2) · $\hat R$+split+rank-norm ✓(T2) · ESS bulk/tail ✓(T3) · divergences/BFMI/tree-depth ✓(T3) · prior predictive ✓(T4) · PPC+p-value+PIT proof ✓(T4) · SBC keystone ✓(T5) · ELPD/WAIC/PSIS-LOO/Pareto-k ✓(T6) · workflow ✓(T6) · eight-schools→geo-MMM ✓(T2/T8) · three `∎` proofs ✓(T2/T3/T5).
- **Placeholder scan:** none (PR body filled at T13).
- **Anchor consistency:** $\hat R=1.8028$, ESS factor $3.0$, SBC rank $4/9$ all NumPy-verified and used identically in Worked Examples (T7); Exercises (T10) use *distinct* fresh numbers.
