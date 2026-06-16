# Time Series, Trend & Seasonality — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write the new foundations chapter "Time Series, Trend & Seasonality" (end of Part II) that pays off the prob-stats i.i.d.-errors IOU — autocorrelation, stationarity, trend/seasonal decomposition, and the Fourier seasonality basis — with two terminology-collision call-outs and a keystone Fourier-orthogonality proof.

**Architecture:** New file `parts/02-regression-bayes/04-time-series.qmd` following the fixed template. Append a gated solutions block to `appendix/solutions.qmd`, add `@hyndman2021` to `references.bib`, and add the file to the Part II `chapters:` list in `_quarto.yml`. **No PR and no book-wide prose renumber in this plan** — both are deferred to the later bundling pass (after PR #37, now merged, is rebased in). Cross-references from this chapter use chapter *names*, never numbers.

**Tech Stack:** Quarto (`.qmd`), KaTeX math, one `{python}` cell (numpy + matplotlib, headless under `MPLBACKEND=Agg`). Real gate is CI `quarto render` (HTML + PDF) at bundle-PR time.

---

## Conventions (enforce in EVERY task)

- **KaTeX:** display `$$ … $$` with delimiters on their own lines; multi-line via `aligned` inside `$$ … $$` — never bare `\begin{align}`. Even `$`-count per line. Inline `$…$`.
- **PDF safety:** never `\begin{psmallmatrix}` (undefined in CI LuaLaTeX — passes KaTeX/HTML, breaks PDF). The PDF gate is unavailable locally this chapter, so hold this strictly.
- **Cross-references by NAME, not number** (renumber-proof): "the MCMC chapter," "the Markov-chains chapter," "the data-generating-process chapter," "Part III," "Part V." Part labels are stable. Do **not** write "Chapter 8" etc.
- **Library-agnostic:** no MMM/PPL library named; `numpy`/`scipy`/`matplotlib` fine.
- **Notation (fixed):** $y_t$ series; $\gamma_k$ autocovariance; $\rho_k = \gamma_k/\gamma_0$ autocorrelation; $\phi$ AR(1) coefficient; $s$ seasonal period; $K$ number of Fourier harmonics; $T$ series length.
- **Verification:** the code cell runs headless with all asserts passing; figures end `plt.show()`. Pinned numbers below are NumPy-verified.
- **Commits:** one per task, prefix `feat(ts): …`, on branch `docs/foundations-additions` (git identity already set).

## Verified anchor numbers (do not change)

- **AR(1)** $\phi=0.7$: $\rho_k=\phi^k$ → $\rho_1=0.7,\ \rho_2=0.49,\ \rho_3=0.343$; stationary variance factor $1/(1-\phi^2)=1.9608$.
- **Fourier orthogonality**, grid $t=0,\dots,s-1$, $s=52$: $\sum_t\cos^2(2\pi j t/s)=\sum_t\sin^2(2\pi j t/s)=s/2=26$; all cross-sums $=0$ (Gram matrix diagonal).
- **Simulated series** (seed 7, $T=312$, trend $0.01t$, seasonal $2\cos\theta+\sin\theta+0.6\cos2\theta$ with $\theta=2\pi t/52$, AR(1) $\phi=0.7,\sigma=0.5$): raw ACF lag1/lag2 $\approx0.958/0.916$ (trend non-stationarity), seasonal spike lag52 $\approx0.73$; after detrend + $K{=}3$ Fourier deseasonalize, residual ACF lag1/2/3 $\approx0.676/0.460/0.334$ (recovers $\phi^k$); $R^2$ by $K{=}1..4$ $\approx0.859/0.905/0.905/0.907$.
- **Robust asserts (seed-independent margins):** Fourier Gram diagonal $=s/2$ and off-diagonal $\approx0$; residual ACF lag-1 within $0.12$ of $\phi=0.7$; $R^2(K{=}2) > R^2(K{=}1)$.

## File structure

- **`parts/02-regression-bayes/04-time-series.qmd`** — the chapter (new). H1 `# Time Series, Trend & Seasonality`.
- **`appendix/solutions.qmd`** — append `## Chapter 7 — Time Series, Trend & Seasonality` block (gated) before the final `:::`.
- **`references.bib`** — add `@hyndman2021`.
- **`_quarto.yml`** — add the new file as the last entry of the Part II `chapters:` list (Task 10).
- Read-only exemplars: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (Part II voice), `parts/03-sampling/02-mcmc.qmd` (chain-ACF collision), `parts/03-sampling/01-markov-chains.qmd` (stationary-distribution collision), `parts/05-mmm-modeling/01-mmm-dgp.qmd` (adstock $\lambda^k$).

---

### Task 1: Front matter, Motivation, and bib entry

**Files:** Create `parts/02-regression-bayes/04-time-series.qmd`; modify `references.bib`

- [ ] **Step 1: Add the bib entry** to `references.bib`:
```bibtex
@book{hyndman2021,
  author    = {Hyndman, Rob J. and Athanasopoulos, George},
  title     = {Forecasting: Principles and Practice},
  edition   = {3rd},
  year      = {2021},
  publisher = {OTexts},
  address   = {Melbourne, Australia},
  url       = {https://otexts.com/fpp3/}
}
```

- [ ] **Step 2: Create the file** with H1 and anchor line:
```
# Time Series, Trend & Seasonality

*Canonical anchors: @jin2017; @hyndman2021.*
```

- [ ] **Step 3: Write `## Motivation`** (3–5 paragraphs). Beats:
  - Every model so far assumed observations are independent — the i.i.d. error of the regression setup. The probability chapter flagged this as temporary: it would be relaxed "when we model autocorrelated or seasonal errors." This chapter relaxes it.
  - Weekly sales are not independent: a high week tends to follow a high week (carryover, momentum), and the calendar imposes a yearly rhythm (holidays, seasons). Dependence across time and seasonality are structure, not noise — and ignoring them biases both estimates and uncertainty.
  - Three tools: the **autocorrelation function** (how a series correlates with its own past), the notion of a **stationary** series (whose statistical character does not drift), and the **decomposition** of a series into trend, season, and remainder — with the seasonal part captured compactly by a **Fourier basis**.
  - Two words collide with earlier uses, and the chapter calls both out explicitly: "autocorrelation" also names a property of an MCMC chain's output (a nuisance to minimize, not signal to model), and "stationary" also names a Markov chain's invariant distribution — a different object entirely.
  - Forward tie: this is the groundwork for the seasonality the MMM data-generating-process chapter builds into its baseline, and the geometric autocorrelation of the simplest time-series model will turn out to be the same algebra as advertising carryover.

- [ ] **Step 4: Verify & commit.** `grep -c '\$' …` even; `grep -c 'begin{align}' …` → 0; confirm no numeric chapter refs: `grep -niE 'chapter [0-9]' parts/02-regression-bayes/04-time-series.qmd` → none.
```
git add parts/02-regression-bayes/04-time-series.qmd references.bib
git commit -m "feat(ts): front matter, Motivation, and bib entry"
```

---

### Task 2: Rungs 1–2 (ACF, stationarity) + collisions + Proof P1

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd`

- [ ] **Step 1: Open `## Theory & Proofs`; write Rung 1 — autocovariance and autocorrelation.** Define a time series $\{y_t\}$; the autocovariance $\gamma_k = \mathrm{Cov}(y_t, y_{t+k})$ and autocorrelation $\rho_k = \gamma_k/\gamma_0$. **Collision call-out #1 (required):** this $\rho_k$ has the same formula as the autocorrelation of an MCMC chain's output (the MCMC chapter, where $n_{\text{eff}} = n/(1+2\sum_k\rho_k)$), but the role is opposite — here $\rho_k$ is **signal in the data we want to model and exploit**; there it is a **sampler nuisance we want to minimize**.

- [ ] **Step 2: Write Rung 2 — covariance-stationarity.** A series is covariance-stationary if $\mathbb E[y_t]=\mu$ is constant, $\mathrm{Var}(y_t)=\gamma_0$ is constant, and $\gamma_k$ depends only on the lag $k$, not on $t$. **Collision call-out #2 (required):** a stationary *series* (moments invariant under a time shift) is a different object from a Markov chain's *stationary distribution* $\pi$ ($\pi P=\pi$, the Markov-chains chapter) — a shared word, unrelated concepts.

- [ ] **Step 3: Write Proof P1 — ACF properties.** State and prove (full prose):
  - **Even:** $\gamma_k = \mathrm{Cov}(y_t,y_{t+k}) = \mathrm{Cov}(y_{t+k},y_t) = \gamma_{-k}$ by stationarity (shift the index), so the autocovariance is symmetric in the lag.
  - **Normalized:** $\rho_0 = \gamma_0/\gamma_0 = 1$.
  - **Bounded:** by Cauchy–Schwarz, $|\gamma_k| = |\mathrm{Cov}(y_t,y_{t+k})| \le \sqrt{\mathrm{Var}(y_t)\mathrm{Var}(y_{t+k})} = \gamma_0$, hence $|\rho_k|\le1$.
  - **Positive-semidefinite:** for any lags $t_1,\dots,t_m$ and weights $a\in\mathbb R^m$, $\mathrm{Var}\big(\sum_i a_i y_{t_i}\big) = \sum_{i,j} a_i a_j \gamma_{|t_i-t_j|} \ge 0$, so the autocovariance matrix $[\gamma_{|t_i-t_j|}]$ is PSD — not every symmetric sequence is a valid autocovariance.

- [ ] **Step 4: Verify & commit.** Even `$`; no `begin{align}`/`psmallmatrix`; no numeric chapter refs.
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): ACF, stationarity, the two collision call-outs, and ACF properties (P1)"
```

---

### Task 3: Rung 3 (AR(1)) + Proof P2 (stationarity & geometric ACF)

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd`

- [ ] **Step 1: Write Rung 3 — the AR(1) model.** $y_t = \phi\,y_{t-1} + \varepsilon_t$ with $\varepsilon_t$ i.i.d. $\mathcal N(0,\sigma^2)$ and $|\phi|<1$. This is the simplest model of autocorrelated errors — the concrete object that "relaxes" the i.i.d. assumption flagged earlier. State that $\phi$ controls persistence: each shock decays by a factor $\phi$ per step.

- [ ] **Step 2: Write Proof P2 — AR(1) stationarity and geometric ACF.** State the proposition ($|\phi|<1 \Rightarrow$ covariance-stationary, mean $0$, $\gamma_0=\sigma^2/(1-\phi^2)$, $\rho_k=\phi^k$) and prove:
  - **MA($\infty$) form.** Unrolling the recursion, $y_t = \sum_{j=0}^{\infty}\phi^j\varepsilon_{t-j}$; the series converges in mean square because $\sum_j\phi^{2j}<\infty$ when $|\phi|<1$.
  - **Mean.** $\mathbb E[y_t] = \sum_j\phi^j\mathbb E[\varepsilon_{t-j}] = 0$.
  - **Variance.** Independence of the $\varepsilon$'s gives $\gamma_0 = \mathrm{Var}(y_t) = \sigma^2\sum_{j=0}^\infty\phi^{2j} = \dfrac{\sigma^2}{1-\phi^2}$.
  - **Autocovariance.** Write $y_{t+k} = \phi^k y_t + \sum_{j=0}^{k-1}\phi^j\varepsilon_{t+k-j}$; the second sum is independent of $y_t$, so $\gamma_k = \mathbb E[y_t y_{t+k}] = \phi^k\gamma_0$, giving $\rho_k = \phi^k$. None of these depend on $t$, so the process is covariance-stationary.
  - **Anchor + cross-link.** At $\phi=0.7$: $\rho_k = 0.7^k$ ($0.7, 0.49, 0.343$), variance factor $1/(1-0.49)\approx1.9608$. Note explicitly that the geometric autocorrelation $\rho_k=\phi^k$ is the same algebra as the geometric **adstock** carryover weights $\lambda^k$ from the data-generating-process chapter — persistence in errors and carryover in advertising share one mathematical form.

- [ ] **Step 3: Verify & commit.**
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): AR(1) stationarity and geometric ACF (P2), the adstock cross-link"
```

---

### Task 4: Rungs 4–5 (decomposition, Fourier) + Proof P3 (keystone)

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd`

- [ ] **Step 1: Write Rung 4 — trend and seasonal decomposition.** The additive model $y_t = T_t + S_t + R_t$: a slow trend $T_t$, a periodic seasonal $S_t$ of period $s$ (e.g. $s=52$ weeks for a yearly cycle), and a stationary remainder $R_t$ (often modeled as AR(1)). The seasonal is identified by a zero-sum constraint over one period, $\sum_{j=0}^{s-1} S_{t+j}=0$, so it carries no trend. Estimating $T_t$ and $S_t$ and subtracting leaves a remainder whose autocorrelation can be examined.

- [ ] **Step 2: Write Rung 5 — Fourier seasonality.** Rather than $s-1$ free seasonal levels, represent the periodic seasonal as a short harmonic sum
$$
S_t = \sum_{j=1}^{K}\Big[a_j\cos\!\big(\tfrac{2\pi j t}{s}\big) + b_j\sin\!\big(\tfrac{2\pi j t}{s}\big)\Big],
$$
$K$ harmonics, $2K$ parameters. A yearly cycle is captured by a handful of harmonics ($K$ small), which is exactly how seasonality enters an MMM as a small set of regressor columns and how the data-generating-process chapter's baseline seasonality is parameterized.

- [ ] **Step 3: Write Proof P3 — Fourier orthogonality (keystone).** State the proposition and prove:
  - **Claim.** On the regular grid $t=0,1,\dots,s-1$, for integers $1\le j,k<s/2$,
  $$
  \sum_{t=0}^{s-1}\cos\!\big(\tfrac{2\pi j t}{s}\big)\cos\!\big(\tfrac{2\pi k t}{s}\big) = \tfrac{s}{2}\,\delta_{jk}, \quad
  \sum_{t=0}^{s-1}\sin\!\big(\tfrac{2\pi j t}{s}\big)\sin\!\big(\tfrac{2\pi k t}{s}\big) = \tfrac{s}{2}\,\delta_{jk}, \quad
  \sum_{t=0}^{s-1}\cos\!\big(\tfrac{2\pi j t}{s}\big)\sin\!\big(\tfrac{2\pi k t}{s}\big) = 0.
  $$
  - **Proof tool.** The root-of-unity sum $\sum_{t=0}^{s-1} e^{2\pi i m t/s} = s$ if $m\equiv0\pmod s$ and $0$ otherwise (geometric series). Write the products of sines/cosines via Euler's formula as combinations of $e^{2\pi i(j\pm k)t/s}$ and apply the identity; the only surviving terms are those with $j=k$, yielding $s/2$.
  - **Consequences.** The harmonic vectors are mutually orthogonal, so $\{1,\cos_j,\sin_j\}$ is an orthogonal basis of the period-$s$ signals: any period-$s$ sequence is represented exactly by $\lfloor s/2\rfloor$ harmonics, a truncation to $K$ harmonics is the least-squares projection onto the leading frequencies, and — because the design columns are orthogonal — the seasonal coefficients $\{a_j,b_j\}$ are estimable independently (no collinearity among Fourier features).

- [ ] **Step 4: Verify & commit.**
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): decomposition, Fourier seasonality, and the orthogonality keystone (P3)"
```

---

### Task 5: Worked Examples

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd` (replace `## Worked Examples` stub)

- [ ] **Step 1: WE1 — ACF of a seasonal series.** Take a clean period-4 pattern $y_t = \cos(2\pi t/4)$ sampled over several periods. Its autocorrelation returns to $+1$ at lags $4,8,\dots$ and dips to $-1$ at lags $2,6,\dots$: the ACF itself is periodic with the series' period, so a peak in the ACF at lag $s$ is the signature of period-$s$ seasonality. Contrast with short-lag dependence (AR(1)), whose ACF decays monotonically without a seasonal peak.

- [ ] **Step 2: WE2 — AR(1) geometric ACF.** $\phi=0.7$: from $\rho_k=\phi^k$, $\rho_1=0.7$, $\rho_2=0.49$, $\rho_3=0.343$ — a geometric decay; the stationary variance is $\gamma_0=\sigma^2/(1-0.49)\approx1.9608\,\sigma^2$. State the cross-link in arithmetic: this is the same geometric sequence as adstock carryover $\lambda^k$ — at $\lambda=0.7$ the carryover weights are $0.7,0.49,0.343,\dots$, identical algebra in a different setting.

- [ ] **Step 3: WE3 — Fourier orthogonality.** On the weekly grid $s=52$: $\sum_{t=0}^{51}\cos^2(2\pi t/52)=26=s/2$; $\sum_{t=0}^{51}\cos(2\pi t/52)\cos(4\pi t/52)=0$ (different harmonics); $\sum_{t=0}^{51}\cos(2\pi t/52)\sin(2\pi t/52)=0$ (cosine vs sine). Conclude the harmonic regressors are orthogonal, so adding harmonics never disturbs the coefficients already estimated.

- [ ] **Step 4: Verify & commit.** Even `$`; `grep -c 'Worked numerical example' …` → 0.
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): worked examples (seasonal ACF, AR(1), Fourier orthogonality)"
```

---

### Task 6: Code Tie-in

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd` (replace `## Code Tie-in` stub)

- [ ] **Step 1: Section intro** (2–3 sentences): simulate a weekly series with trend, yearly Fourier seasonality, and AR(1) noise; read its sample ACF (trend non-stationarity + a seasonal spike), then decompose — fit trend and Fourier seasonal by least squares — and confirm the residual ACF recovers the AR(1) geometric decay. Fixed seed for reproducibility.

- [ ] **Step 2: Write the single `{python}` cell** (verified headless and reproducible):

```{python}
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)
s, T = 52, 312
t = np.arange(T)
phi = 0.7

# --- WE2: AR(1) theory ---
assert np.allclose([phi**k for k in range(4)], [1, 0.7, 0.49, 0.343])
assert np.isclose(1 / (1 - phi**2), 1.9608, atol=1e-4)

# --- Fourier design (intercept + linear trend + K harmonics) ---
def fourier_design(t, s, K):
    cols = [np.ones_like(t, dtype=float), t.astype(float)]
    for j in range(1, K + 1):
        cols.append(np.cos(2 * np.pi * j * t / s))
        cols.append(np.sin(2 * np.pi * j * t / s))
    return np.column_stack(cols)

# --- P3 / WE3: harmonic columns are orthogonal on one period ---
tp = np.arange(s)
H = fourier_design(tp, s, 3)[:, 2:]            # drop intercept + trend, keep harmonics
G = H.T @ H
assert np.allclose(np.diag(G), s / 2)          # each harmonic norm^2 = s/2 = 26
assert np.allclose(G - np.diag(np.diag(G)), 0, atol=1e-8)   # off-diagonals vanish

# --- simulate trend + seasonal + AR(1) noise ---
trend = 0.01 * t
season = 2.0 * np.cos(2*np.pi*t/s) + 1.0 * np.sin(2*np.pi*t/s) + 0.6 * np.cos(4*np.pi*t/s)
ar = np.zeros(T); eps = rng.normal(0, 0.5, T)
for k in range(1, T):
    ar[k] = phi * ar[k-1] + eps[k]
y = trend + season + ar

def acf(x, maxlag):
    x = x - x.mean(); g0 = np.mean(x * x)
    return np.array([np.mean(x[:len(x)-k] * x[k:]) / g0 for k in range(maxlag + 1)])

raw = acf(y, 60)                               # slow decay (trend) + spike at lag 52

# --- decompose: least-squares fit of trend + K=3 Fourier, inspect residual ---
X = fourier_design(t, s, 3)
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
resid = y - X @ beta
res_acf = acf(resid, 10)

# --- the story, asserted ---
assert abs(res_acf[1] - phi) < 0.12            # residual ACF recovers AR(1) phi
def r2(K):
    Xk = fourier_design(t, s, K); b, *_ = np.linalg.lstsq(Xk, y, rcond=None)
    return 1 - np.var(y - Xk @ b) / np.var(y)
assert r2(2) > r2(1)                           # added harmonic raises fit
print(f"raw ACF lag1={raw[1]:.3f} lag2={raw[2]:.3f} seasonal lag52={raw[52]:.3f}")
print(f"residual ACF lag1..3 = {np.round(res_acf[1:4], 3)}  (AR1 truth 0.7, 0.49, 0.343)")
print(f"R^2 by K=1..4: {[round(r2(K), 3) for K in (1, 2, 3, 4)]}")

# --- figures ---
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].stem(range(61), raw)
ax[0].axvline(52, ls="--", color="grey"); ax[0].set_title("raw ACF: trend decay + seasonal spike")
ax[0].set_xlabel("lag")
ax[1].stem(range(11), res_acf)
ax[1].plot(range(11), [phi**k for k in range(11)], "r--", label=r"$\phi^k$")
ax[1].set_title("residual ACF recovers AR(1)"); ax[1].set_xlabel("lag"); ax[1].legend()
ax[2].plot([1, 2, 3, 4], [r2(K) for K in (1, 2, 3, 4)], "o-")
ax[2].set_title(r"fit $R^2$ vs Fourier order $K$"); ax[2].set_xlabel("K"); ax[2].set_ylabel(r"$R^2$")
plt.tight_layout()
plt.show()
```

- [ ] **Step 3: Run headless.** Extract the cell to `/tmp/ts_cell.py`, run `MPLBACKEND=Agg python3 /tmp/ts_cell.py`. Expected: three print lines, no assertion error, exit 0 (benign Agg `plt.show()` warning is fine).

- [ ] **Step 4: Commit.**
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): code tie-in — ACF, decomposition, residual AR(1) recovery, Fourier order"
```

---

### Task 7: Exercises (C / B / P / A)

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd` (replace `## Exercises` stub)

> **Controller note:** controller writes Exercises directly. Self-contained, name-based refs, 2–3 per tier.

- [ ] **Step 1: `### C`.** e.g.: in one sentence each, distinguish the series ACF from an MCMC chain's autocorrelation, and a stationary series from a Markov stationary distribution; why does a peak in the ACF at lag $s$ indicate seasonality?; why are Fourier seasonal features preferable to $s-1$ dummy variables?
- [ ] **Step 2: `### B`.** e.g.: for AR(1) with $\phi=0.5$, compute $\rho_1,\rho_2,\rho_3$ and the variance factor $1/(1-\phi^2)$; verify $\sum_{t=0}^{3}\cos^2(2\pi t/4)=2=s/2$ on the period-4 grid; show $\sum_{t=0}^{3}\cos(2\pi t/4)\sin(2\pi t/4)=0$.
- [ ] **Step 3: `### P`.** e.g.: prove $|\rho_k|\le1$ via Cauchy–Schwarz; prove the AR(1) autocovariance $\gamma_k=\phi^k\gamma_0$ from the MA($\infty$) form; prove $\sum_{t=0}^{s-1}\cos(2\pi j t/s)\sin(2\pi k t/s)=0$ using the root-of-unity sum.
- [ ] **Step 4: `### A`.** e.g.: simulate AR(1) at $\phi\in\{0.3,0.7,0.95\}$ and overlay the sample ACFs against $\phi^k$; fit Fourier seasonality at $K=1,2,3$ to a simulated series and plot residual variance vs $K$; show that orthogonality means the $K=1$ coefficients are unchanged when you add the second harmonic.
- [ ] **Step 5: Commit.**
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): exercises (C/B/P/A)"
```

---

### Task 8: Summary

**Files:** Modify `parts/02-regression-bayes/04-time-series.qmd` (add `## Summary` after Exercises)

> **Controller note:** controller writes Summary directly. **Key concepts** and **Key identities** BOTH bulleted (inline math), per the formatting rule.

- [ ] **Step 1: Write `## Summary`** with a 1–2 sentence lead, then:
  - `**Key concepts**` bullets: autocovariance/autocorrelation and the chain-ACF collision; covariance-stationarity and the Markov-stationary-distribution collision; AR(1) as the minimal autocorrelated-error model; trend/seasonal/remainder decomposition; Fourier seasonality and orthogonality; the adstock cross-link (geometric ACF = geometric carryover).
  - `**Key identities**` bullets: $\rho_k=\gamma_k/\gamma_0$, $|\rho_k|\le1$; covariance-stationarity ($\gamma_k$ depends only on $k$); AR(1) $\gamma_0=\sigma^2/(1-\phi^2)$, $\rho_k=\phi^k$; decomposition $y_t=T_t+S_t+R_t$; Fourier seasonal $S_t=\sum_{j=1}^K[a_j\cos(2\pi jt/s)+b_j\sin(2\pi jt/s)]$; orthogonality $\sum_t\cos^2=\sum_t\sin^2=s/2$, cross-sums $0$.
  - Closing sentence: this groundwork feeds the seasonality baseline of the MMM data-generating-process chapter and, in Part V, the dynamic state-space models that let coefficients move over time.
- [ ] **Step 2: Verify (bulleted, even `$`) & commit.**
```
git add parts/02-regression-bayes/04-time-series.qmd
git commit -m "feat(ts): summary"
```

---

### Task 9: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Insertion point.** The file is one `::: {.content-visible when-meta="show-solutions"}` wrapper; its final closing `:::` is the last line. Insert the new block immediately before that final `:::`, after the last existing chapter block, as `## Chapter 7 — Time Series, Trend & Seasonality` (heading text only — no numeric chapter cross-refs inside).

- [ ] **Step 2: Write full solutions** for every C/B/P/A item from Task 7. Pin arithmetic: AR(1) $\phi=0.5$ → $\rho_1=0.5,\rho_2=0.25,\rho_3=0.125$, variance factor $1/(1-0.25)=1.3333$; period-4 grid $\sum\cos^2=2$, $\sum\cos\sin=0$. For A, give a short headless-verified `{python}` cell (overlay sample ACFs vs $\phi^k$, or residual-variance-vs-$K$); run it before committing.

- [ ] **Step 3: Verify.** Balanced `:::` (no new fence — lives inside the existing wrapper); even `$`-count on added lines (`awk` scan from the `## Chapter 7 — Time` line to end → no odd lines); any A cell runs headless. Commit.
```
git add appendix/solutions.qmd
git commit -m "feat(ts): appendix solutions"
```

---

### Task 10: Add to `_quarto.yml`, final local review

**Files:** Modify `_quarto.yml`

- [ ] **Step 1: Add the chapter to the build.** In `_quarto.yml`, append to the **Part II** `chapters:` list (after `parts/02-regression-bayes/03-hierarchical-regression.qmd`):
```yaml
        - parts/02-regression-bayes/04-time-series.qmd
```
Do **not** touch any other part, and do **not** renumber prose cross-references anywhere — that is the deferred bundling pass.

- [ ] **Step 2: Structure/KaTeX lint** on the chapter:
```bash
f=parts/02-regression-bayes/04-time-series.qmd
grep -c '\$' "$f"; grep -c 'begin{align}' "$f"; grep -c 'psmallmatrix' "$f"; grep -niE 'chapter [0-9]' "$f"
grep -nE '^## ' "$f"
awk '{n=gsub(/\$/,"$"); if(n%2) print "ODD line "NR}' "$f"
```
Expected: even `$`; `0`; `0`; no numeric chapter refs; the six template headings in order; no odd lines.

- [ ] **Step 3: Re-run the code cell headless** under `MPLBACKEND=Agg python3` — exit 0, all asserts pass.

- [ ] **Step 4: Commit. Do NOT open a PR** (the branch is held open for the bundle; the renumber + bundled PR come in the later pass).
```
git add _quarto.yml
git commit -m "feat(ts): add Time Series chapter to Part II in _quarto.yml"
```

- [ ] **Step 5: Report** to the controller that the chapter is complete on `docs/foundations-additions` and that the remaining work (rebase onto merged `main`, the `+1` prose renumber, and the bundled PR with CI HTML+PDF verification) is the next, separate phase.

---

## Self-Review (controller, before dispatching Task 1)

**Spec coverage:** Motivation + IOU payoff (T1) ✓; bib `@hyndman2021` (T1) ✓; ACF + collision #1 (T2) ✓; stationarity + collision #2 (T2) ✓; P1 ACF properties (T2) ✓; AR(1) + P2 geometric ACF + adstock cross-link (T3) ✓; decomposition (T4) ✓; Fourier + P3 orthogonality keystone (T4) ✓; WE1/WE2/WE3 (T5) ✓; code cell with ACF/decomposition/residual-AR(1)/Fourier-order, headless-verified (T6) ✓; C/B/P/A (T7) ✓; bulleted Summary (T8) ✓; gated appendix solutions (T9) ✓; `_quarto.yml` insertion, no renumber, no PR (T10) ✓; cross-refs by name throughout ✓; deferred renumber/PR called out (T10 Step 5) ✓.

**Placeholder scan:** no TBD/TODO; all asserted numbers NumPy-verified; the code cell is complete and run headless; recovery/ACF claims use robust tolerance asserts.

**Consistency:** notation fixed once ($y_t,\gamma_k,\rho_k,\phi,s,K,T$) and reused; `fourier_design`/`acf`/`r2` defined once in T6; the $\rho_k=\phi^k$ identity and the $s/2$ orthogonality value match between P2/P3, WE2/WE3, the code asserts, and the Summary.
