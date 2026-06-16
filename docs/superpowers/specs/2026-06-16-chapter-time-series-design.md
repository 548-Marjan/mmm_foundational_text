# Chapter (new, end of Part II) — Time Series, Trend & Seasonality: Design Spec

**Date:** 2026-06-16
**Status:** Approved (brainstorm complete). Next: implementation plan via `writing-plans`.
**Branch:** `docs/foundations-additions` (the foundations-additions + renumber bundle).
**Placement:** **end of Part II — Regression & Bayesian Inference**, after Hierarchical
Regression. Becomes the new **Chapter 7** once the deferred `+1` renumber is applied.
**File:** `parts/02-regression-bayes/04-time-series.qmd` (new file).
**Canonical anchors:** `@jin2017` (seasonality in MMM), `@hyndman2021` (Forecasting:
Principles and Practice — to be added to `references.bib`). Existing stats refs
`@casella2002` / `@wasserman2004` are available for moment/Cauchy–Schwarz citations.

---

## 1. Role in the book

This chapter pays off the promissory note left in the probability-statistics chapter, whose
linear-regression setup assumes i.i.d. Gaussian errors "one we will relax when we model
autocorrelated or seasonal errors" (`parts/01-foundations/03-probability-statistics.qmd:75`).
It supplies the foundations for data structured in time — dependence across lags, seasonality,
and trend — and sets up two things the rest of the book already leans on: the **seasonality
baseline** the MMM data-generating process uses ($\sin(2\pi t/52)$ in Chapter 14's DGP) and the
**Fourier seasonality features** used to model a yearly cycle in practice (the backend's
`yearly_seasonality` Fourier order).

It is a **foundations** chapter, deliberately elementary. It is **not** the dynamic-linear-model /
state-space chapter (that is a later Part V chapter, Kalman/FFBS); this one is the groundwork:
autocorrelation, stationarity, decomposition, and the Fourier seasonality basis. AR/MA is kept
minimal — AR(1) only — with no ARIMA.

The chapter follows the fixed template: Motivation → Theory & Proofs → Worked Examples → Code
Tie-in → Exercises (C/B/P/A) → Summary (auto-included), with gated appendix solutions.

## 2. Two terminology-collision call-outs (required)

The book has used the words "autocorrelation" and "stationary" in two *other* senses already; this
chapter must explicitly disambiguate both, because the collisions are real and confusing:

1. **Series ACF vs MCMC-chain autocorrelation.** The autocorrelation $\rho_k$ of a time series
   (this chapter) shares its formula with the autocorrelation of an MCMC chain's output (the MCMC
   chapter, $n_{\text{eff}} = n/(1 + 2\sum_k \rho_k)$), but the *role is opposite*: here $\rho_k$ is
   **signal in the data we want to model**; there it is a **sampler nuisance we want to minimize**.
2. **Series stationarity vs Markov stationary distribution.** A covariance-stationary *series*
   (moments invariant under a time shift) is a different object from a Markov chain's *stationary
   distribution* $\pi$ ($\pi P = \pi$, the Markov-chains chapter). A shared word, different concepts.

## 3. Theory rungs

1. **Time series and the i.i.d. break.** Define a series $\{y_t\}$; the **autocovariance**
   $\gamma_k = \mathrm{Cov}(y_t, y_{t+k})$ and **autocorrelation** $\rho_k = \gamma_k/\gamma_0$.
   State call-out #1 here.
2. **Covariance-stationarity.** Constant mean $\mathbb E[y_t] = \mu$, constant variance, and
   autocovariance depending only on the lag $k$ (not on $t$). State call-out #2 here.
3. **AR(1) — the minimal autocorrelated-error model.** $y_t = \phi\,y_{t-1} + \varepsilon_t$,
   $\varepsilon_t$ i.i.d. $\mathcal N(0,\sigma^2)$, $|\phi| < 1$. This is the model that "relaxes"
   the i.i.d. assumption (the IOU), in its simplest form.
4. **Trend and seasonal decomposition.** Additive model $y_t = T_t + S_t + R_t$ (trend +
   seasonal + remainder); the seasonal component is periodic with period $s$ and is identified by a
   zero-sum-over-a-period constraint $\sum_{j=0}^{s-1} S_{t+j} = 0$.
5. **Fourier seasonality.** Represent the periodic seasonal as a finite harmonic sum
   $S_t = \sum_{j=1}^{K} \big[a_j\cos(2\pi j t/s) + b_j\sin(2\pi j t/s)\big]$ — $K$ harmonics,
   $2K$ parameters. This is exactly how a yearly cycle ($s=52$ weeks) enters the model as a small
   set of regressor columns (the Fourier order $K$), and how Chapter 14's baseline seasonality is
   parameterized.

## 4. Proof slate (one keystone + two lighter)

| # | Result | Statement |
|---|---|---|
| **P1 (light)** | ACF properties | The autocovariance is even, $\gamma_k = \gamma_{-k}$; $\rho_0 = 1$ and $|\rho_k| \le 1$ by Cauchy–Schwarz; and the autocovariance sequence is positive-semidefinite (any finite covariance matrix $[\gamma_{|i-j|}]$ is PSD because it is a genuine covariance matrix). |
| **P2** | AR(1) stationarity and geometric ACF | For $\lvert\phi\rvert < 1$, the AR(1) process is covariance-stationary with mean $0$, variance $\gamma_0 = \sigma^2/(1-\phi^2)$, and autocorrelation $\rho_k = \phi^{\,k}$. Proof: solve the recursion as $y_t = \sum_{j\ge0}\phi^j\varepsilon_{t-j}$ (valid since $\lvert\phi\rvert<1$), compute $\gamma_0$ as a geometric series, and $\gamma_k = \phi^k\gamma_0$ from the recursion. **Cross-link:** the geometric ACF $\rho_k=\phi^k$ rhymes with adstock's geometric carryover weights $\lambda^k$ — same algebra, foreshadowing Part V. |
| **P3 (keystone)** | Fourier seasonality representation | On the regular grid $t = 0,1,\dots,s-1$, the vectors $\{\cos(2\pi j t/s)\}$ and $\{\sin(2\pi j t/s)\}$ for $1 \le j < s/2$ are mutually **orthogonal**, with $\sum_{t=0}^{s-1}\cos^2(2\pi j t/s) = \sum_{t=0}^{s-1}\sin^2(2\pi j t/s) = s/2$ and all cross-sums zero. Consequently any period-$s$ sequence is exactly represented by $\lfloor s/2\rfloor$ harmonics, a truncation to $K$ harmonics is the least-squares projection onto the leading frequencies, and orthogonality makes the seasonal-regression coefficients $\{a_j,b_j\}$ estimable independently (the design columns are uncorrelated). |

## 5. Worked examples + anchor numbers (NumPy-verified in the plan)

- **WE1 — ACF of a seasonal series.** Build a small periodic series and show its sample ACF spikes
  at the seasonal lag (e.g. a period-4 pattern has an ACF peak at lag 4), distinguishing a seasonal
  signal from short-lag dependence.
- **WE2 — AR(1) geometric ACF.** $\phi = 0.7$: $\rho_k = 0.7^k$ giving $\rho_1 = 0.7$,
  $\rho_2 = 0.49$, $\rho_3 = 0.343$; stationary variance $\gamma_0 = \sigma^2/(1-0.49) \approx 1.9608\,\sigma^2$.
  Explicitly contrast $\rho_k = \phi^k$ with adstock's $\lambda^k$.
- **WE3 — Fourier orthogonality.** On the $s=52$ weekly grid, verify
  $\sum_t \cos(2\pi t/52)^2 = 26 = s/2$, $\sum_t \cos(2\pi t/52)\cos(4\pi t/52) = 0$, and
  $\sum_t \cos(2\pi t/52)\sin(2\pi t/52) = 0$; show a $K{=}1$ vs $K{=}2$ harmonic fit to a seasonal
  shape and how added harmonics sharpen it.

## 6. Code tie-in

A single runnable `{python}` cell: simulate a weekly series = linear trend + yearly Fourier seasonal
($s=52$) + AR(1) noise ($\phi=0.7$); compute and plot the **sample ACF** (showing the seasonal spike
near lag 52 and the geometric short-lag decay from the AR(1) part); build the **Fourier design matrix**
($\cos,\sin$ at $K$ harmonics), assert its columns are orthogonal ($X^\top X$ diagonal, entries $s/2$),
and least-squares-fit the seasonal, showing residual variance dropping as $K$ grows. Asserts pin:
AR(1) $\rho_k=\phi^k$ to the sample ACF within tolerance; Fourier column orthogonality exact; $R^2$
increasing in $K$. Runs headless under `MPLBACKEND=Agg`.

## 7. Exercises, Summary, housekeeping

- **Exercises** C/B/P/A (self-contained); appendix solutions appended to `appendix/solutions.qmd`
  as `## Chapter 7 — Time Series, Trend & Seasonality` (heading text only; **no chapter *number* in
  cross-references** — see §8), gated by `::: {.content-visible when-meta="show-solutions"}`, before
  the file's final closing `:::`.
- **Summary** auto-included; **Key concepts** and **Key identities** both **bulleted** (inline math
  only) — per the formatting-consistency rule.
- **Bib:** add `@hyndman2021` (Hyndman & Athanasopoulos, *Forecasting: Principles and Practice*) to
  `references.bib`.
- **Anchor line** under the H1: `*Canonical anchors: @jin2017; @hyndman2021.*`

## 8. Sequencing & the deferred renumber (critical)

This chapter is written **now**, but the book-wide `+1` renumber is **deferred** until PR #37
(Chapter 15, Building & Fitting) merges to `main` and `docs/foundations-additions` is rebased onto
it (the handoff's locked decision #3). Therefore:

- **Cross-references from this chapter use chapter *names*, not numbers** ("the MCMC chapter," "the
  Markov-chains chapter," "Part III," "Part V," "the MMM data-generating-process chapter"). This
  makes the new chapter renumber-proof and keeps the global renumber pass from having to edit its
  prose. Part labels (Part II/III/V) are stable and may be used directly.
- **`_quarto.yml`:** the new file is added to the Part II `chapters:` list as the final entry. (This
  alone shifts Quarto's *rendered* numbers for Part III onward; the textual "Chapter N" references
  elsewhere in the book are corrected in the later renumber pass, not in this chapter's work.)
- **CI/PDF verification** of this chapter happens when the **bundled PR** opens (no PR before then,
  per the handoff). Until then, verify locally: KaTeX/structure lint and a headless run of the code
  cell. Hold the PDF-safety conventions strictly (no `psmallmatrix`, etc.) since the PDF gate is not
  available locally.

## 9. Conventions (enforced per task)

- KaTeX: `aligned` inside `$$ … $$` (never bare `\begin{align}`); `$$` delimiters on their own lines;
  even `$`-count per line. **Never** `\begin{psmallmatrix}` (undefined in the CI LuaLaTeX/PDF build).
- Library-agnostic: name no MMM/PPL library in prose; `numpy`/`scipy`/`matplotlib` are fine.
- Build-based verification: the code cell runs headless; the real gate is CI `quarto render`
  (HTML + PDF) at bundle-PR time.

## 10. Critical files

- `parts/02-regression-bayes/04-time-series.qmd` — the chapter (new file).
- `appendix/solutions.qmd` — append the new solutions block before the final `:::`.
- `references.bib` — add `@hyndman2021`.
- `_quarto.yml` — add the new file as the last entry of the Part II `chapters:` list.
- Read-only exemplars: `parts/02-regression-bayes/03-hierarchical-regression.qmd` (Part II voice),
  `parts/03-sampling/02-mcmc.qmd` (the chain-autocorrelation collision source),
  `parts/03-sampling/01-markov-chains.qmd` (the stationary-distribution collision source),
  `parts/05-mmm-modeling/01-mmm-dgp.qmd` (adstock $\lambda^k$ for the AR(1) cross-link, and the
  seasonality baseline this chapter underpins).
