# Chapter 19 — Causal Inference Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (streamlined variant — one implementer subagent per task, controller reviews math inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Write Chapter 19, which opens Part VI by explaining the Chapter 18 wound causally — spend is endogenous, so the fitted response coefficient is structurally biased — leading with potential outcomes, bridging to DAGs/do-operator, proving the OVB and backdoor-adjustment theorems, and handing off to the interventional designs of Chapter 20.

**Architecture:** A single Quarto `.qmd` chapter on the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus an appendix solutions block and three new bib entries. Verification is build-based: NumPy anchor checks + one headless-runnable `{python}` cell + CI `quarto render`. No pytest — the chapter's empirical gates replace TDD, exactly as Chapters 1–18.

**Tech Stack:** Quarto (`.qmd`, KaTeX), Python (`numpy`, `numpy.linalg.lstsq`, `matplotlib`), BibTeX.

---

## Spec

`docs/superpowers/specs/2026-06-17-chapter-19-causal-foundations-design.md` (committed). Read it for full rung rationale.

## House rules (enforce every task)

- **Never** name PyMC-Marketing or any MMM/PPL/sampler library. `numpy`/`scipy`/`matplotlib` only.
- **Never** `\begin{psmallmatrix}` (breaks CI PDF). Use `bmatrix`/`pmatrix`/`smallmatrix`.
- KaTeX: multi-line math uses `aligned` inside `$ … $`/`$$ … $$`, never bare `\begin{align}`. `$$` on own lines; even `$$` count per file.
- "Key identities" in the Summary must be a **bulleted list**, not a run-on paragraph.
- Cross-reference siblings by name + post-renumber number (below).
- Commit with identity: `git -c user.name='jlh530i' -c user.email='jlh530i@gmail.com' commit -m "..."`.

## Cross-references (post-renumber numbering)

- **Chapter 4 — Linear Regression:** the endogeneity teaser planted there is cashed out here.
- **Chapter 16 — Building & Fitting an MMM:** the Fisher-information ridge; its causal reading is confounding.
- **Chapter 18 — Budget Optimization:** the wound this chapter explains structurally.
- **Chapter 20 — Quasi-Experimental Designs:** interventional identification (IV teased here).
- **Chapter 21 — Advanced Calibration:** folds interventional estimates into the posterior.

## Files

- **Modify:** `parts/06-causal-calibration/01-causal-foundations.qmd` — replace stub body; keep H1 `# Causal Inference Foundations`. (The stub also has an italic anchors line under H1: `*Canonical anchors: Imbens & Rubin 2015; Angrist & Pischke; Pearl.*` — keep or fold into Motivation; do not leave a dangling stub callout.)
- **Modify:** `appendix/solutions.qmd` — append `## Chapter 19 — Causal Inference Foundations` after the Chapter 18 block, inside the existing `::: {.content-visible when-meta="show-solutions"}` div (before its closing `:::`).
- **Modify:** `references.bib` — add `@imbens2015`, `@angrist2009`, `@pearl2009`.
- Read-only voice exemplars: `parts/05-mmm-modeling/04-budget-optimization.qmd` (Ch. 18, immediate predecessor), `parts/05-mmm-modeling/02-building-fitting.qmd` (Ch. 16, the ridge).

---

## Verified Anchors (NumPy-checked, seed 19, N=200,000 — use these exact numbers)

**Confounded MMM DGP:** demand $Z\sim N(0,1)$ (the confounder); spend $X=\delta Z+\nu$, $\nu\sim N(0,1)$; sales $Y=\beta X+\gamma Z+\varepsilon$, $\varepsilon\sim N(0,1)$. Parameters $\beta=2$, $\gamma=3$, $\delta=1$.

| Quantity | Value |
|---|---|
| True response slope $\beta$ | $2$ |
| Naive OLS slope ($Y\sim X$), analytic plim | $\beta+\gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X) = 2 + 3\cdot\tfrac12 = 3.5$ |
| Naive OLS slope, simulated (seed 19) | $\approx 3.50$ |
| OVB bias | $\gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X) = 3\cdot 0.5 = 1.5$ |
| Adjusted OLS slope ($Y\sim X+Z$) | $\approx 2.00$ (recovers $\beta$) |
| Interventional slope $dE[Y\mid\mathrm{do}(x)]/dx$ | $2$ (cut $Z\to X$) |
| Observational slope $dE[Y\mid x]/dx$ | $\approx 3.5$ |
| Intervention check (regenerate $X=\nu\perp Z$, re-fit naive) | $\approx 2.00$ |

Key derived facts: $\mathrm{cov}(X,Z)=\delta\,\mathrm{var}(Z)=1$; $\mathrm{var}(X)=\delta^2\mathrm{var}(Z)+\mathrm{var}(\nu)=2$; bias $=\gamma\cdot 1/2 = 1.5$. Treat simulated values as "$\approx$" in prose; assert in code with tolerances (e.g. `abs(naive-3.5)<0.02`, `abs(adj-2.0)<0.02`).

---

## Task 1 — Front matter + Motivation

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (replace stub `## Motivation`; keep H1; remove the stub callout).

- [ ] **Step 1: Write Motivation** (3–5 paragraphs): Part V ended on a wound — the optimizer ran over a curve the observational data only weakly pins. This chapter gives the reason: **spend is endogenous**. Managers spend more in periods they expect to sell well, so spend correlates with demand drivers that themselves move sales; the fitted slope confuses "spend works" with "spend happens when demand is high." State the driving question — *"a fitted coefficient is causal only under conditions we must name"* — and preview the two tools: potential outcomes (what *would* sales have been at a different spend?) and the do-operator (intervening vs observing). Name the payoff arc: a concrete confounded MMM where the naive slope overstates the truth by $75\%$, recovered only by controlling the confounder — and the wound is that in real MMM the confounder is usually unobserved, which is why Part VI must turn to experiments. Cross-reference Chapter 4 (endogeneity teaser), Chapter 18 (the wound), and Part VI's path (Chapter 20 designs, Chapter 21 calibration). Keep the canonical-anchors mention (Imbens & Rubin, Angrist & Pischke, Pearl) either as a kept line or folded into prose.
- [ ] **Step 2: Lint** — even `$$`, no `\begin{align}`/`psmallmatrix`, H1 intact, no banned library names.
- [ ] **Step 3: Commit** — `feat(ch19): motivation`.

## Task 2 — Theory rungs 1–3 (potential outcomes; ignorability; OVB keystone)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 1 — the causal question (potential outcomes).** Define $Y_i(x)$ = period $i$'s sales if spend were set to $x$; the dose–response $\mu(x)=\mathbb E[Y_i(x)]$ is the response curve the MMM wants and the Chapter 18 optimizer needs. State the **fundamental problem of causal inference** (only $Y_i(X_i)$ is observed). The fitted MMM regresses observed sales on observed spend; whether its slope equals $d\mu/dx$ depends on how spend was chosen — and spend is an endogenous management decision, not randomized.
- [ ] **Step 2: Rung 2 — ignorability and the observational/interventional gap.** Define unconfoundedness/ignorability $Y_i(x)\perp X_i\mid Z_i$ (among periods with the same demand state $Z$, spend is as-good-as-random) plus positivity/overlap. Under it, $\mu(x)=\mathbb E_Z[\mathbb E[Y\mid X=x,Z]]$ (adjustment). Without it, $\mathbb E[Y\mid X=x]$ measures something else. Name the estimand: the response slope $d\mu/dx$ (the marginal return the optimizer consumes).
- [ ] **Step 3: Rung 3 — Proof P1 (KEYSTONE): omitted-variable / confounding bias.** Linear structural MMM $Y=\beta X+\gamma Z+\varepsilon$, endogenous spend $X=\delta Z+\nu$. **Theorem:** $\operatorname{plim}\hat\beta_{\text{naive}} = \beta + \gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X)$. **Proof:** $\mathrm{cov}(X,Y)=\beta\,\mathrm{var}(X)+\gamma\,\mathrm{cov}(X,Z)+\mathrm{cov}(X,\varepsilon)$; with $\varepsilon\perp X$, divide by $\mathrm{var}(X)$. Bias $=\gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X)$, zero **iff** $\mathrm{cov}(X,Z)=0$ or $\gamma=0$. MMM reading: endogenous spend biases the slope and the bias does **not** shrink with $N$ — the rigorous wound. **Anchor:** $\beta=2,\gamma=3,\delta=1,\mathrm{var}(Z)=\mathrm{var}(\nu)=1\Rightarrow$ bias $=1.5$, $\operatorname{plim}\hat\beta_{\text{naive}}=3.5$ (a $75\%$ overstatement). Present multi-line math as `aligned` inside `$$`.
- [ ] **Step 4: Lint + commit** — `feat(ch19): theory rungs 1-3 (potential outcomes, ignorability, OVB)`.

## Task 3 — Theory rungs 4–6 (SCM/DAG/do-operator; backdoor P2; the wound)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (continue `## Theory & Proofs`).

- [ ] **Step 1: Rung 4 — SCMs, DAGs, the do-operator (bridge).** Introduce the structural causal model (each variable a function of parents + noise) and its DAG. The confounded MMM is the fork $X\leftarrow Z\rightarrow Y$ with $X\rightarrow Y$; $Z$ opens a **backdoor path** $X\leftarrow Z\rightarrow Y$ carrying non-causal association. Define the **do-operator**: $P(Y\mid\mathrm{do}(X=x))$ is the distribution in the mutilated graph with $X$'s incoming arrows cut (truncated factorization). Contrast $P(Y\mid\mathrm{do}(x))$ (cut $Z\to X$) with $P(Y\mid x)$; in the linear-Gaussian MMM the interventional slope is $\beta$ and the observational slope is $\beta+\text{bias}$ — Rung 3, graphical. (If a graph is drawn, use text/ASCII or a simple TikZ-free description — do NOT introduce `psmallmatrix` or fragile LaTeX.)
- [ ] **Step 2: Rung 5 — Proof P2: backdoor adjustment + identification.** **Theorem:** if $Z$ satisfies the backdoor criterion relative to $(X,Y)$ (blocks every backdoor path; no descendant of $X$), then $P(Y\mid\mathrm{do}(X=x))=\sum_z P(Y\mid X=x,Z=z)P(Z=z)$. **Proof** via truncated factorization: the interventional joint deletes $X$'s mechanism, $P(\cdot\mid\mathrm{do}(x))=[\prod_{V\neq X}P(V\mid\mathrm{pa}_V)]_{X=x}$; marginalize $Z$, whose distribution is unchanged by the intervention. **Bridge:** the graphical backdoor criterion equals the potential-outcomes ignorability $Y(x)\perp X\mid Z$ — two languages, one condition. Conclusion: observational MMM recovers the causal curve **iff** a sufficient adjustment set is observed and conditioned on.
- [ ] **Step 3: Rung 6 — the wound restated, and the handoff.** In MMM the sufficient adjustment set is usually **unobserved** (you cannot measure every demand shock, competitor move, price change). Then the backdoor cannot close, ignorability fails, and **no observational estimator recovers $\beta$** — irreducible bias, the causal content of Chapter 16's ridge and Chapter 18's wound. Escape = change the data-generating process: **Chapter 20** supplies interventional identification (geo experiments, holdouts, IV-style natural experiments; IV is the unobserved-confounder escape, teased here), **Chapter 21** calibrates them back in. Name the loop: intervene → identify → calibrate → re-optimize.
- [ ] **Step 4: Lint + commit** — `feat(ch19): theory rungs 4-6 (SCM/do-operator, backdoor adjustment, the wound)`.

## Task 4 — Worked Examples (WE1, WE2, WE3)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (`## Worked Examples`). Use the exact Verified Anchors.

- [ ] **Step 1: WE1 — confounded MMM, adjustment in action.** DGP $Z\sim N(0,1)$, $X=Z+\nu$, $Y=2X+3Z+\varepsilon$. Show (i) naive OLS slope $\approx 3.50$, biased high by $1.5 = \gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X)=3\cdot\tfrac12$; (ii) regressing $Y$ on $X$ **and** observed $Z$ recovers $\hat\beta\approx 2.00$ — backdoor adjustment made arithmetic. The fitted curve is wrong by $75\%$ until $Z$ is controlled.
- [ ] **Step 2: WE2 — the unobserved confounder (the wound).** Same DGP, $Z$ now unobserved; only the naive regression is feasible, stuck at $\approx 3.50$. Quantify the irreducible bias; state every observational estimator from $(X,Y)$ inherits it — Chapter 18's wound in causal language, Chapter 16's ridge its posterior shadow. Tease an instrument as the escape (Chapter 20).
- [ ] **Step 3: WE3 — observation vs intervention (do-operator).** Observational slope $dE[Y\mid x]/dx\approx 3.5$ vs interventional $dE[Y\mid\mathrm{do}(x)]/dx=2$. Demonstrate numerically by severing $Z\to X$: regenerate $X=\nu$ independent of demand (randomized assignment), keep $Y=2X+3Z+\varepsilon$, recover $\hat\beta\approx 2.0$ from the naive regression. Randomizing spend closes the backdoor by construction — the math of "run an experiment."
- [ ] **Step 4: Lint + commit** — `feat(ch19): worked examples (confounded MMM, the wound, do-operator)`.

## Task 5 — Code Tie-in (controller writes directly)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (`## Code Tie-in`).

> Controller writes directly, then verifies headless. Single `{python}` cell, `numpy`+`matplotlib`, every claim asserted, figure ends `plt.show()`.

- [ ] **Step 1: Write the cell.** `rng=np.random.default_rng(19)`, `N` large. (1) simulate confounded DGP ($\beta=2,\gamma=3,\delta=1$); (2) naive OLS via `lstsq` ($Y\sim X$) → $\approx 3.50$, assert analytic OVB $=1.5$; (3) adjusted OLS ($Y\sim X+Z$) → $\approx 2.00$; (4) hidden-confounder note (adjustment unavailable); (5) intervention: regenerate $X=\nu\perp Z$, re-fit naive → $\approx 2.00$; (6) figure: scatter $(X,Y)$ (subsample for clarity) colored by $Z$, with observational fit (slope $3.5$) and true/interventional slope (slope $2.0$) overlaid. Assert all against analytic values.
- [ ] **Step 2: Verify headless** — extract the cell, run `MPLBACKEND=Agg python3`; confirm all asserts pass.
- [ ] **Step 3: Lint + commit** — `feat(ch19): code tie-in — confounding, adjustment, intervention`.

## Task 6 — Exercises (C / B / P / A — controller writes directly)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (`## Exercises`). Self-contained; no inline solution links.

- [ ] **Step 1: Write exercises.**
  - **C:** why a fitted MMM coefficient is not automatically causal; the fundamental problem of causal inference in MMM terms; why confounding bias does not vanish with more data; what ignorability means for spend and why it is usually indefensible observationally.
  - **B:** derive $\operatorname{plim}\hat\beta_{\text{naive}}=\beta+\gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X)$; evaluate the bias for the WE1 parameters; show it is zero when $\mathrm{cov}(X,Z)=0$.
  - **P:** (i) prove the OVB theorem; (ii) prove the backdoor adjustment via truncated factorization, and show ignorability $Y(x)\perp X\mid Z$ ⇔ $Z$ blocks the backdoor in the fork $X\leftarrow Z\rightarrow Y$.
  - **A:** sweep confounding strength $\delta$, plot naive bias vs $\delta$ (linear in $\mathrm{cov}(X,Z)$); add a second observed confounder and show the adjustment set must be complete; simulate an instrument and recover $\beta$ (Chapter 20 teaser).
- [ ] **Step 2: Lint + commit** — `feat(ch19): exercises (C/B/P/A)`.

## Task 7 — Summary (controller writes directly)

**Files:** Modify `parts/06-causal-calibration/01-causal-foundations.qmd` (`## Summary`).

- [ ] **Step 1: Write Summary.** Opening prose + bulleted **Key concepts** (endogeneity of spend; potential outcomes; ignorability; confounding/OVB; SCM/do-operator; backdoor adjustment; the unobserved-confounder wound) + bulleted **Key identities** (inline math, bulleted):
  - dose–response $\mu(x)=\mathbb E[Y(x)]$;
  - ignorability $Y(x)\perp X\mid Z$;
  - OVB $\operatorname{plim}\hat\beta=\beta+\gamma\,\mathrm{cov}(X,Z)/\mathrm{var}(X)$;
  - backdoor adjustment $P(Y\mid\mathrm{do}(x))=\sum_z P(Y\mid x,z)P(z)$;
  - observation-vs-intervention slope gap.
  Close tying forward to Chapter 20 (interventional identification) and Chapter 21 (calibration).
- [ ] **Step 2: Lint + commit** — `feat(ch19): summary`.

## Task 8 — Appendix solutions

**Files:** Modify `appendix/solutions.qmd` (append after the Chapter 18 block, before the final `:::`).

- [ ] **Step 1: Append** `## Chapter 19 — Causal Inference Foundations` inside the gated div, in chapter order. Full C/B/P/A solutions: C prose; B the plim derivation + bias $=1.5$ + the $\mathrm{cov}(X,Z)=0$ case; P both proofs in full (OVB plim; truncated-factorization backdoor + ignorability⇔backdoor equivalence); A expected code behavior (bias linear in $\delta$; complete adjustment set; instrument recovers $\beta$).
- [ ] **Step 2: Lint** (even `$$` in appendix; gated div still closed) **+ commit** — `feat(ch19): appendix solutions`.

## Task 9 — Bibliography + final review + PR

**Files:** Modify `references.bib`; final lint of the chapter.

- [ ] **Step 1: Add bib entries** to `references.bib`:
  ```bibtex
  @book{imbens2015,
    author    = {Imbens, Guido W. and Rubin, Donald B.},
    title     = {Causal Inference for Statistics, Social, and Biomedical Sciences: An Introduction},
    year      = {2015},
    publisher = {Cambridge University Press},
    address   = {Cambridge}
  }
  @book{angrist2009,
    author    = {Angrist, Joshua D. and Pischke, J{\"o}rn-Steffen},
    title     = {Mostly Harmless Econometrics: An Empiricist's Companion},
    year      = {2009},
    publisher = {Princeton University Press},
    address   = {Princeton, NJ}
  }
  @book{pearl2009,
    author    = {Pearl, Judea},
    title     = {Causality: Models, Reasoning, and Inference},
    edition   = {2nd},
    year      = {2009},
    publisher = {Cambridge University Press},
    address   = {Cambridge}
  }
  ```
  Confirm every `@key` cited in the chapter resolves in `references.bib`.
- [ ] **Step 2: Full chapter lint** — six headings present and in order (Motivation, Theory & Proofs, Worked Examples, Code Tie-in, Exercises, Summary); H1 `# Causal Inference Foundations` intact; even `$$`; no `\begin{align}`; no `psmallmatrix`; no banned library names (word-boundary grep: `pymc`, `stan`, `orbit`, `numpyro`, `pyro`); all citation keys valid.
- [ ] **Step 3: Re-run the code cell headless** one final time — all asserts pass.
- [ ] **Step 4: Push + open PR** against `main` (`gh pr create`); summarize (potential-outcomes lead, OVB + backdoor keystones, the confounded-MMM anchor recovering $\beta=2$, the wound → Chapter 20/21). Start background CI-render watcher (poll `gh run watch` to a green conclusion) and merge watcher (`until gh pr view <N> --json state -q .state | grep -qE 'MERGED|CLOSED'; do sleep 120; done`, run_in_background).
- [ ] **Step 5:** Report PR URL + CI status. **User merges** (never self-merge).

---

## Self-Review (controller, before dispatching)

1. **Spec coverage:** Rungs 1–6 → Tasks 2–3; WE1–3 → Task 4; Code → Task 5; C/B/P/A → Task 6; Summary → Task 7; appendix → Task 8; bib + CI → Task 9. Framework choice (potential outcomes lead → DAG bridge) realized in rung order (PO in 1–3, DAG/do in 4–5). All spec sections mapped. ✅
2. **Placeholders:** none — every anchor is a NumPy-verified number; every proof has its argument. ✅
3. **Consistency:** the confounded DGP ($\beta=2,\gamma=3,\delta=1$, seed 19, unit variances) is used identically across Rung 3, WE1–3, and the Code Tie-in; bias $=1.5$, naive $=3.5$, adjusted/interventional $=2.0$ everywhere. ✅
