# Chapter 20 — Quasi-Experimental Design Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (streamlined variant — one implementer subagent per task, controller reviews math inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Write Chapter 20, which supplies the interventional identification Chapter 19 proved necessary — organized around the secant/tangent/validation taxonomy, with two keystone proofs (DiD parallel-trends; IV-is-a-secant) and full worked/code treatment of DiD, IV, synthetic control, and ITS/BSTS, all recovering the same lift on one geo panel.

**Architecture:** A single Quarto `.qmd` chapter on the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus an appendix solutions block and two new bib entries. Verification is build-based: NumPy anchor checks + one headless-runnable `{python}` cell + CI `quarto render`. No pytest — the chapter's empirical gates replace TDD, as in Chapters 1–19.

**Tech Stack:** Quarto (`.qmd`, KaTeX), Python (`numpy`, `scipy.optimize`, `matplotlib`), BibTeX.

---

## Spec

`docs/superpowers/specs/2026-06-17-chapter-20-quasi-experimental-design.md` (committed). Read for full rung rationale.

## House rules (enforce every task)

- **Never** name PyMC-Marketing or any MMM/PPL/sampler library. Describe BSTS/CausalImpact as *methods*, never products. `numpy`/`scipy`/`matplotlib` only.
- **Never** `\begin{psmallmatrix}`. Use `bmatrix`/`pmatrix`/`smallmatrix`.
- KaTeX: multi-line math uses `aligned` inside `$ … $`/`$$ … $$`, never bare `\begin{align}`. `$$` on own lines; even `$$` count per file.
- "Key identities" in the Summary must be a **bulleted list**, not a run-on paragraph.
- Cross-reference siblings by name + post-renumber number (below).
- Commit with identity: `git -c user.name='jlh530i' -c user.email='jlh530i@gmail.com' commit -m "..."`.

## Cross-references (post-renumber numbering)

- **Chapter 17 — Dynamic Linear Models:** the local-level Kalman filter reused for the ITS/BSTS counterfactual.
- **Chapter 18 — Budget Optimization:** the response curve $S$ these designs identify is the one the optimizer consumes.
- **Chapter 19 — Causal Inference Foundations:** the confounding these designs defeat; the IV teaser cashed here.
- **Chapter 21 — Advanced Calibration:** ingests each design's estimand as a taxonomy-tagged likelihood factor.

## Files

- **Modify:** `parts/06-causal-calibration/02-quasi-experimental.qmd` — replace stub body; keep H1 `# Quasi-Experimental Design`. The stub has an italic anchors line (`*Canonical anchors: Angrist & Pischke; Abadie, Diamond & Hainmueller (synthetic control); Brodersen et al. (BSTS).*`) — convert to resolving Pandoc citations or keep; remove the stub callout.
- **Modify:** `appendix/solutions.qmd` — append `## Chapter 20 — Quasi-Experimental Design` after the Chapter 19 block, inside the `::: {.content-visible when-meta="show-solutions"}` div (before its closing `:::`).
- **Modify:** `references.bib` — add `@abadie2010`, `@brodersen2015` (`@angrist2009` already present).
- Read-only voice exemplars: `parts/06-causal-calibration/01-causal-foundations.qmd` (Ch. 19, immediate predecessor), `parts/05-mmm-modeling/03-dlm.qmd` (Ch. 17, the Kalman machinery reused).

---

## Verified Anchors (NumPy-checked — use these exact numbers)

**Anchor curve:** $S(x) = 2\sqrt{x}$ (concave; Chapter 18's surrogate). $S'(x) = 1/\sqrt{x}$.

| Quantity | Value |
|---|---|
| $S(4)$, $S(9)$ | $4$, $6$ |
| Lift $S(9)-S(4)$ | $2$ |
| Secant slope $[S(9)-S(4)]/(9-4)$ | $0.4$ |
| MVT point $\xi$ with $S'(\xi)=0.4$ | $\xi = (1/0.4)^2 = 6.25$ |
| Tangents $S'(4)$, $S'(9)$ | $0.5$, $0.3333$ |

**Geo panel (DiD):** baseline demand (region fixed effects) $\alpha_T=10$, $\alpha_C=8$; treated spend $4\to9$, control held at $4$.

| Cell | Value |
|---|---|
| $Y_{T,\text{pre}}=\alpha_T+S(4)$ | $14$ |
| $Y_{T,\text{post}}=\alpha_T+S(9)$ | $16$ |
| $Y_{C,\text{pre}}=Y_{C,\text{post}}=\alpha_C+S(4)$ | $12$ |
| Naive post cross-section $Y_{T,\text{post}}-Y_{C,\text{post}}$ | $4$ (= lift $2$ + baseline gap $2$, **biased**) |
| **DiD** $(16-14)-(12-12)$ | $2$ (= lift, confound differenced out) |

**IV / Wald:** instrument shifts complier spend $4\to9$. $\hat\beta_{\text{IV}}=\Delta Y/\Delta X = (S(9)-S(4))/(9-4)=2/5=0.4$ = secant slope = $S'(6.25)$; **not** $S'(4)=0.5$ nor $S'(9)=0.333$.

**Synthetic control:** treated pre-level $14$; convex weights on controls match it; post gap $=2$.

**ITS/BSTS (Kalman):** treated pre-level $\approx14$; forecast counterfactual $14$; observed post $16$; gap $=2$.

Treat noisy simulated values as "$\approx$"; assert in code with tolerances (e.g. `abs(did-2)<0.1`, `abs(wald-0.4)<0.05`, `abs(gap-2)<0.15`). Seed `np.random.default_rng(20)`.

---

## Task 1 — Front matter + Motivation

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (replace stub `## Motivation`; keep H1; remove stub callout).

- [ ] **Step 1: Write Motivation** (3–5 paragraphs): Chapter 19 proved observational MMM cannot identify the curve when the confounder is unobserved; the escape is to *generate* spend variation independent of demand. Introduce the chapter's organizing question — *"an experiment gives you one number, but a number about what part of the curve?"* — and the secant/tangent/validation taxonomy as the answer. Preview the designs (DiD, IV, synthetic control, ITS/BSTS) and the two surprises: DiD differences out the very confound Chapter 19 flagged, and IV's number is a secant (interval-average slope), not a point marginal return. Name the payoff: each design hands Chapter 21 a calibrated functional of the same curve $S$ the Chapter 18 optimizer needs. Cross-reference Chapters 17 (Kalman), 18 (the curve), 19 (confounding/IV teaser), 21 (calibration).
- [ ] **Step 2: Lint** — even `$$`, no `\begin{align}`/`psmallmatrix`, H1 intact, no banned library names.
- [ ] **Step 3: Commit** — `feat(ch20): motivation`.

## Task 2 — Theory rungs 1–3 (taxonomy; DiD keystone; IV-secant keystone)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 1 — the taxonomy.** Define $S(x)$ and its functionals: **secant** $[S(x_1)-S(x_0)]/(x_1-x_0)$ (or level change $S(x_1)-S(x_0)$); **tangent** $S'(x)$ (or curvature $S''(x)$); **validation** = not a functional of $S$. State the thesis: every design identifies one of these, and which one is what the experiment legitimately says about the curve. Present the taxonomy table (secant: DiD/synthetic control/ITS/geo lift/IV; tangent: RDD/RKD; validation: non-$S$ estimands).
- [ ] **Step 2: Rung 2 — Proof P1 (KEYSTONE): DiD under parallel trends.** Potential-outcomes panel ($T,C$ × pre,post). Parallel-trends assumption $\mathbb E[Y_{T,\text{post}}(0)-Y_{T,\text{pre}}(0)]=\mathbb E[Y_{C,\text{post}}(0)-Y_{C,\text{pre}}(0)]$. **Theorem:** DiD estimand $=$ ATT $\mathbb E[Y_{T,\text{post}}(1)-Y_{T,\text{post}}(0)]$. **Proof:** write observed $=$ potential outcomes, add/subtract $Y_{T,\text{post}}(0)$, cancel the common trend via parallel trends; the region fixed effect (baseline demand = Chapter 19 confounder) differences out. MMM reading: ATT $=$ lift $S(x_1)-S(x_0)$, a **secant**. Anchor: $\alpha_T=10,\alpha_C=8$, $4\to9$: naive $=4$, DiD $=2=S(9)-S(4)$. Multi-line math as `aligned` inside `$$`.
- [ ] **Step 3: Rung 3 — Proof P2 (KEYSTONE): IV is a secant in disguise.** Wald estimator $\hat\beta_{\text{IV}}=\frac{\mathbb E[Y\mid W=1]-\mathbb E[Y\mid W=0]}{\mathbb E[X\mid W=1]-\mathbb E[X\mid W=0]}$; instrument validity (relevance, exclusion, independence, monotonicity) ⇒ LATE over compliers. **Theorem:** for $Y_i(x)=S(x)+\alpha_i$, $\hat\beta_{\text{IV}}=\frac{S(x_1)-S(x_0)}{x_1-x_0}=\frac{1}{x_1-x_0}\int_{x_0}^{x_1}S'(u)\,du$ — the average derivative = secant slope (FTC), $=S'(\xi)$ for interior $\xi$ (MVT), **not** the tangent at the operating point. Anchor: $4\to9$ ⇒ $0.4=S'(6.25)$, $\ne S'(4)=0.5,\ S'(9)=0.333$. Cashes Chapter 19's IV teaser.
- [ ] **Step 4: Lint + commit** — `feat(ch20): theory rungs 1-3 (taxonomy, DiD, IV-secant)`.

## Task 3 — Theory rungs 4–6 (synthetic control; ITS/BSTS; tangent designs + handoff)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (continue `## Theory & Proofs`).

- [ ] **Step 1: Rung 4 — synthetic control (constructed counterfactual).** One treated unit, many controls; convex weights $w\ge0,\sum w=1$ matching the treated **pre-period**; post gap $=$ lift. Identifying intuition: matching the pre-period trajectory matches latent demand factors (factor-model view); note Abadie convex-hull/pre-fit conditions — no formal theorem. Estimand is a **secant**. Anchor: pre-match $14$, post gap $2$.
- [ ] **Step 2: Rung 5 — ITS and BSTS (state-space counterfactual).** ITS fits the treated unit's pre-intervention trajectory and extrapolates it as the counterfactual; the post gap is the lift. The modern construction builds that extrapolation with the **Chapter 17 local-level Kalman filter** — make the reuse explicit (the forecast level *is* the counterfactual; the forecast variance gives the band). Estimand is a **secant**. Describe BSTS/CausalImpact as a method, never a product. Anchor: pre-level $14$, forecast $14$, observed post $16$, gap $2$.
- [ ] **Step 3: Rung 6 — tangent designs (brief) + handoff.** RDD identifies the response **level** at a cutoff, RKD the **kink** (curvature $S''$) — the only true point-local **tangent** designs, obtained as $\delta\to0$. Mention switchback and staggered-adoption DiD briefly. Close: each design hands Chapter 21 a taxonomy-tagged functional of $S$ (secant/tangent/validation) — the schema the calibration likelihood and prior store consume. Loop: intervene → identify → calibrate → re-optimize.
- [ ] **Step 4: Lint + commit** — `feat(ch20): theory rungs 4-6 (synthetic control, ITS/BSTS, tangent designs)`.

## Task 4 — Worked Examples (WE1, WE2, WE3)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (`## Worked Examples`). Use the exact Verified Anchors.

- [ ] **Step 1: WE1 — geo-lift DiD (secant).** Tabulate the 2×2 panel ($Y_T=14\to16$, $Y_C=12\to12$). Show naive post cross-section $16-12=4$ is biased (lift $2$ + baseline gap $2$); DiD $(16-14)-(12-12)=2$ recovers the lift $S(9)-S(4)=2$, differencing out the confound; parallel trends holds by construction; lift is a secant (slope $0.4$).
- [ ] **Step 2: WE2 — IV / Wald (secant in disguise).** Instrument shifts complier spend $4\to9$; $\hat\beta_{\text{IV}}=2/5=0.4$ = secant = $S'(6.25)$, not $S'(4)=0.5$ or $S'(9)=0.333$. Contrast with Chapter 19's confounded naive slope. Lesson: IV gives an interval-average slope; you must know where on the curve it sits to act at an operating point.
- [ ] **Step 3: WE3 — synthetic control.** One treated region, several controls with different baselines; convex weights match the treated pre-level $14$; synthetic stays $\approx14$ post while treated rises to $16$, gap $=2$. Pre-fit is the identification lever; estimand is a secant. (Exact weights pinned in the Code Tie-in.)
- [ ] **Step 4: Lint + commit** — `feat(ch20): worked examples (geo-DiD, IV/Wald, synthetic control)`.

## Task 5 — Code Tie-in (controller writes directly)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (`## Code Tie-in`).

> Controller writes directly, then verifies headless. Single `{python}` cell, `numpy`+`scipy`+`matplotlib`, every claim asserted, figure ends `plt.show()`.

- [ ] **Step 1: Write the cell.** `rng=np.random.default_rng(20)`. Build one geo panel from $S(x)=2\sqrt x$ + region fixed effects + noise. (1) **DiD** 2×2 estimate, assert $\approx2$ and naive $\approx4$; (2) **IV/Wald** via $\mathrm{cov}(W,Y)/\mathrm{cov}(W,X)$, assert $\approx0.4$, print $S'(4),S'(9)$ to show it is neither; (3) **synthetic control** convex weights via `scipy.optimize.minimize` on the simplex (or `nnls` + renormalize) matching treated pre-period, assert post gap $\approx2$; (4) **ITS/BSTS** local-level Kalman (reuse the Chapter 17 recursion: predict/update with state/obs variances) on the treated pre-period, forecast counterfactual, assert post gap $\approx2$ with a band; (5) figure: recovered lifts side by side (all $\approx2$) + a panel drawing the secant of $S$ on $[4,9]$ vs the tangents at $4,9$.
- [ ] **Step 2: Verify headless** — extract the cell, run `MPLBACKEND=Agg python3`; confirm all asserts pass.
- [ ] **Step 3: Lint + commit** — `feat(ch20): code tie-in — DiD, IV, synthetic control, ITS-Kalman`.

## Task 6 — Exercises (C / B / P / A — controller writes directly)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (`## Exercises`). Self-contained; no inline solution links.

- [ ] **Step 1: Write exercises.**
  - **C:** why every QE estimand is a functional of $S$ and why the taxonomy matters; why IV's number is a secant not a point marginal return; what parallel trends assumes and when a geo experiment violates it.
  - **B:** on $S=2\sqrt x$: lift and secant slope for $4\to9$; the MVT point $\xi=6.25$; the DiD 2×2 arithmetic and the naive bias.
  - **P:** (i) prove DiD identifies the ATT under parallel trends; (ii) prove the Wald estimand $=\frac{1}{x_1-x_0}\int_{x_0}^{x_1}S'(u)\,du$ (secant via FTC), stating the instrument-validity conditions.
  - **A:** break parallel trends (control given a trend) and show DiD bias; widen the instrument's complier interval and show the Wald secant drifts; add a second control and show synthetic-control pre-fit improves the counterfactual.
- [ ] **Step 2: Lint + commit** — `feat(ch20): exercises (C/B/P/A)`.

## Task 7 — Summary (controller writes directly)

**Files:** Modify `parts/06-causal-calibration/02-quasi-experimental.qmd` (`## Summary`).

- [ ] **Step 1: Write Summary.** Opening prose + bulleted **Key concepts** (the taxonomy; DiD differences out confounds; IV is a secant; synthetic control / ITS-BSTS as constructed counterfactuals; RDD/RKD as the true tangent designs; every estimand a functional of $S$) + bulleted **Key identities** (inline math, bulleted):
  - secant $[S(x_1)-S(x_0)]/(x_1-x_0)$;
  - DiD $=(\bar Y_{T,\text{post}}-\bar Y_{T,\text{pre}})-(\bar Y_{C,\text{post}}-\bar Y_{C,\text{pre}})$;
  - parallel trends (untreated trends equal);
  - Wald $=\Delta Y/\Delta X=\frac{1}{x_1-x_0}\int_{x_0}^{x_1}S'$;
  - tangent $S'(x)$ vs curvature $S''(x)$.
  Close tying forward to Chapter 21 (ingests the tagged functional) and back to Chapter 18 (the curve).
- [ ] **Step 2: Lint + commit** — `feat(ch20): summary`.

## Task 8 — Appendix solutions

**Files:** Modify `appendix/solutions.qmd` (append after the Chapter 19 block, before the final `:::`).

- [ ] **Step 1: Append** `## Chapter 20 — Quasi-Experimental Design` inside the gated div, in chapter order. Full C/B/P/A: C prose; B the lift/secant/$\xi$/DiD arithmetic; P both proofs in full (DiD parallel-trends; Wald = average-derivative secant via FTC, with the instrument conditions); A expected code behavior (parallel-trends break ⇒ DiD bias; wider complier interval ⇒ secant drift; extra control ⇒ better pre-fit).
- [ ] **Step 2: Lint** (even `$$` in appendix; gated div still closed) **+ commit** — `feat(ch20): appendix solutions`.

## Task 9 — Bibliography + final review + PR

**Files:** Modify `references.bib`; final lint of the chapter.

- [ ] **Step 1: Add bib entries** to `references.bib`:
  ```bibtex
  @article{abadie2010,
    author    = {Abadie, Alberto and Diamond, Alexis and Hainmueller, Jens},
    title     = {Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California's Tobacco Control Program},
    journal   = {Journal of the American Statistical Association},
    volume    = {105}, number = {490}, pages = {493--505}, year = {2010},
    publisher = {Taylor \& Francis}
  }
  @article{brodersen2015,
    author    = {Brodersen, Kay H. and Gallusser, Fabian and Koehler, Jim and Remy, Nicolas and Scott, Steven L.},
    title     = {Inferring Causal Impact Using Bayesian Structural Time-Series Models},
    journal   = {The Annals of Applied Statistics},
    volume    = {9}, number = {1}, pages = {247--274}, year = {2015},
    publisher = {Institute of Mathematical Statistics}
  }
  ```
  Confirm every `@key` cited in the chapter resolves (`@angrist2009`, `@abadie2010`, `@brodersen2015`).
- [ ] **Step 2: Full chapter lint** — six headings present and in order; H1 `# Quasi-Experimental Design` intact; even `$$`; no `\begin{align}`; no `psmallmatrix`; no banned library names (word-boundary grep: `pymc`, `stan`, `orbit`, `numpyro`, `pyro`); all citation keys valid.
- [ ] **Step 3: Re-run the code cell headless** one final time — all asserts pass.
- [ ] **Step 4: Push + open PR** against `main` (`gh pr create`); summarize (taxonomy spine, DiD + IV-secant keystones, the four designs recovering lift $2$ on one panel, ITS reusing Ch. 17 Kalman). Start background CI-render watcher (`gh run watch` to a green conclusion) and merge watcher (`until gh pr view <N> --json state -q .state | grep -qE 'MERGED|CLOSED'; do sleep 120; done`, run_in_background).
- [ ] **Step 5:** Report PR URL + CI status. **User merges** (never self-merge).

---

## Self-Review (controller, before dispatching)

1. **Spec coverage:** Rungs 1–6 → Tasks 2–3; WE1–3 → Task 4; Code (4 designs + taxonomy figure) → Task 5; C/B/P/A → Task 6; Summary → Task 7; appendix → Task 8; bib + CI → Task 9. The two keystones (DiD, IV-secant) are Rungs 2–3 with full proofs; synthetic control / ITS-BSTS are worked (Rungs 4–5 + code), RDD/RKD brief (Rung 6). All spec sections mapped. ✅
2. **Placeholders:** none — every anchor is a NumPy-verified number; every proof has its argument. ✅
3. **Consistency:** the curve $S(x)=2\sqrt x$, the geo panel ($\alpha_T=10,\alpha_C=8$, $4\to9$), and the recovered lift $2$ / secant $0.4$ / $\xi=6.25$ are used identically across Rungs 2–5, WE1–3, and the Code Tie-in; seed 20 throughout. ✅
