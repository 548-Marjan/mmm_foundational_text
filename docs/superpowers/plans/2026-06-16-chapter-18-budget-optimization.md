# Chapter 18 — Budget Optimization & Scenario Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (streamlined variant — one implementer subagent per task, controller reviews math inline) to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 18, which runs Part IV's optimizer on Part V's *fitted, uncertain* response curve — closing the Part IV allocation loop, proving the envelope theorem and the optimizer's curse, making scenario analysis rigorous via the shadow price, and ending Part V on the wound that motivates Part VI.

**Architecture:** A single Quarto `.qmd` chapter following the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus an appendix solutions block and one new bib entry. All numeric anchors are NumPy-verified (see "Verified Anchors" below). Verification is **build-based**: NumPy anchor checks + a single headless-runnable `{python}` cell + CI `quarto render` (HTML+PDF). There are no pytest unit tests — the chapter's empirical gates replace TDD, exactly as in Chapters 1–17.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), Python (`numpy`, `scipy.optimize`, `matplotlib`) for the Code Tie-in, BibTeX.

---

## Spec

Design spec: `docs/superpowers/specs/2026-06-16-chapter-18-budget-optimization-design.md` (committed). Read it for the full rung rationale. This plan operationalizes it task-by-task.

## House rules (enforce every task)

- **Never** name PyMC-Marketing or any MMM/PPL/sampler library. `numpy`/`scipy`/`matplotlib` only.
- **Never** use `\begin{psmallmatrix}` (breaks CI PDF/LuaLaTeX). Use `bmatrix`/`pmatrix`/`smallmatrix`.
- KaTeX: multi-line math uses `aligned` **inside** `$ … $` (never bare `\begin{align}`); `$$` delimiters on their own lines; even count of `$$` per file.
- "Key identities" in the Summary must be a **bulleted list**, never a run-on paragraph (sibling-chapter formatting rule).
- Cross-reference sibling chapters **by name and post-renumber number** (see Cross-refs below).
- Commit with identity `jlh530i` / `jlh530i@gmail.com`:
  `git -c user.name='jlh530i' -c user.email='jlh530i@gmail.com' commit -m "..."`.

## Cross-references (post-renumber numbering)

- **Ch. 12 — Convexity:** equal-marginal / KKT optimum; the $b^\star=(7.2,1.8)$, $\lambda=0.3727$ anchor.
- **Ch. 13 — Linear Programming:** shadow price as an LP dual multiplier.
- **Ch. 14 — Constrained Nonlinear Optimization & SLSQP:** the solver; multiplier-as-shadow-price.
- **Ch. 15 — The MMM Data-Generating Process:** the saturating, adstocked response curve.
- **Ch. 16 — Building & Fitting an MMM:** the posterior over the curve; the **Fisher-information ridge** (weak identification) that drives the wound.
- **Part VI (Ch. 19–22):** named as the heal — causal foundations, quasi-experimental designs, calibration.

## Files

- **Modify:** `parts/05-mmm-modeling/04-budget-optimization.qmd` — replace the stub body; keep H1 `# Budget Optimization & Scenario Analysis`.
- **Modify:** `appendix/solutions.qmd` — append a `## Chapter 18 — Budget Optimization & Scenario Analysis` block at end of file (line 2774), gated by `::: {.content-visible when-meta="show-solutions"}`.
- **Modify:** `references.bib` — add `@smith2006` (optimizer's curse).
- Read-only voice exemplars: `parts/05-mmm-modeling/03-dlm.qmd` (Ch. 17, immediate predecessor), `parts/05-mmm-modeling/02-building-fitting.qmd` (Ch. 16, source of the ridge), `parts/04-optimization/03-slsqp.qmd` (Ch. 14, the optimizer being reused).

---

## Verified Anchors (all NumPy-checked — use these exact numbers)

**Square-root response model (the spine):** two channels, $r_i(b)=a_i\sqrt{b_i}$, $a=(2,1)$, budget $B=9$, $S=\sum_i a_i^2 = 5$.

| Quantity | Value |
|---|---|
| Optimum $b^\star = B\,a^2/S$ | $(7.2,\ 1.8)$ |
| Shadow price $\lambda^\star = \sqrt S/(2\sqrt B) = 1/\sqrt{7.2}$ | $0.3727$ |
| Value $V^\star(9)=\sqrt{S B}=\sqrt{45}$ | $6.7082$ |
| Value function | $V^\star(B)=\sqrt{5B}$ |
| Envelope $dV^\star/dB = \sqrt5/(2\sqrt B)$ at $B=9$ | $0.3727 = \lambda^\star$ |

**Scenario table (WE2):**

| Scenario | Result |
|---|---|
| Budget $+20\%$, $B=10.8$ | exact $V^\star=\sqrt{54}=7.3485$; first-order $6.7082+0.3727\cdot1.8=7.3790$ (overshoots by concavity) |
| Budget $-20\%$, $B=7.2$ | $V^\star=\sqrt{36}=6.0000$ |
| Channel-1 cap $b_1\le5$ | $b=(5,4)$, $V=2\sqrt5+2=6.4721$ |
| Channel-2 dark $b_2=0$ | $b=(9,0)$, $V=2\sqrt9=6.0000$ |
| **Frontier** | $V^\star(B)=\sqrt{5B}$, marginal $\lambda^\star(B)=\sqrt5/(2\sqrt B)$ |
| **Inverse** (hit $V_0=8$) | $B=V_0^2/S=12.8$; $\lambda^\star=\sqrt5/(2\sqrt{12.8})=5/16=0.3125$ (exact) |

**WE3 — ridge posterior + optimizer's curse.** Posterior on $a$ centered $\bar a=(2,1)$ with covariance
`Cov = 0.45**2 * outer(u,u) + 0.08**2 * outer(v,v)`, where `u=(1,-1)/sqrt2` (attribution-contrast — the **poorly identified** ridge direction) and `v=(1,1)/sqrt2` (total effect — well identified). Sampler: `rng=np.random.default_rng(18)`, `N=200_000` draws, keep draws with both coords `>0` (≈99.9% kept). Optimum is closed-form per draw: $b^\star(a)=B\,a^2/S(a)$, $S(a)=a_1^2+a_2^2$ (no per-draw SLSQP needed).

| Quantity | Value (seed=18, N=200k) |
|---|---|
| $b_1^\star$ posterior median | $\approx 7.20$ |
| $b_1^\star$ 5th / 95th percentile | $\approx 4.35\ /\ 8.69$ |
| $\lambda^\star$ 5th / 95th percentile | $\approx 0.343\ /\ 0.431$ |
| $\max_b \mathbb E[V] = \sqrt{S B}$, $\bar S=5$ | $6.7082$ |
| $\mathbb E[\max_b V] = \sqrt B\,\mathbb E\lVert a\rVert$ | $\approx 6.8276$ |
| **Optimizer's-curse / EVPI gap** $\mathbb E[\max]-\max\mathbb E$ | $\approx 0.1194\ (\ge 0)$ |
| Jensen check: $\mathbb E\lVert a\rVert\ge\lVert\bar a\rVert$ | $2.276 \ge 2.236$ |

Treat the percentile/MC numbers as "$\approx$" in prose (they depend on seed+N); assert them in code with tolerances (e.g. `EVPI gap > 0.10`, `median b1 within 0.1 of 7.2`).

**WE3 — two equal-fit / different-optimum curves (the wound; exact, seed-free).** Put both curves on the same response-capacity circle $S=5.0625$ (so $V=\sqrt{9\cdot5.0625}=6.75$ exactly), at angles $\theta=15^\circ$ and $75^\circ$: $a=(\sqrt S\cos\theta,\sqrt S\sin\theta)$.

| Curve | $a$ | fit $V$ | optimum $b^\star$ |
|---|---|---|---|
| $\theta=15^\circ$ | $(2.1733,\ 0.5823)$ | $6.7500$ | $(8.3971,\ 0.6029)$ |
| $\theta=75^\circ$ | $(0.5823,\ 2.1733)$ | $6.7500$ | $(0.6029,\ 8.3971)$ |

Identical in-sample fit, mirror-image allocations — *good fit ≠ good decision*.

---

## Task 1 — Front matter + Motivation

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (replace stub `## Motivation`; keep H1).

- [ ] **Step 1: Write Motivation.** 3–5 short paragraphs establishing:
  - Parts IV and V each delivered half the budgeting answer: Part IV the *optimizer* (for a known curve), Part V the *curve* (fitted, uncertain). This chapter joins them.
  - The driving question: *"We can solve for the optimal budget. But should we trust the answer?"*
  - Preview the two new ideas beyond Part IV: (a) the **shadow price** $\lambda^\star=dV^\star/dB$ that powers scenario analysis, and (b) the **optimizer's curse** — optimizing an estimated curve overstates realized performance, by an amount equal to what an experiment would be worth.
  - Name the payoff arc: recover Part IV's $(7.2,1.8)$ on the fitted curve, then watch the Ch. 16 ridge turn that confident answer into a wide one — the **wound** Part VI heals.
  - MMM anchoring: this is the chapter a practitioner actually runs to set next quarter's spend.
- [ ] **Step 2: Lint** — even `$$` count, no `\begin{align}`, no `psmallmatrix`, H1 intact.
- [ ] **Step 3: Commit** — `feat(ch18): motivation`.

## Task 2 — Theory rungs 1–3 (problem setup; equal-marginal applied; envelope theorem P2)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (`## Theory & Proofs`).

- [ ] **Step 1: Rung 1 — the budgeting problem on a fitted curve.** State
  $\max_{b\ge0,\ \mathbf1^\top b=B} V(b;\theta)=\sum_i r_i(b_i;\theta_i)$ with $r_i$ the fitted saturating MMM response (Ch. 15) and $\theta$ carrying the Ch. 16 posterior. Emphasize the shift from Part IV: the objective is now *random* because $\theta$ is uncertain. Introduce the square-root surrogate $r_i(b)=a_i\sqrt b$ as the analytically tractable stand-in used throughout (concave, saturating, matches Part IV's anchor) — state plainly it is the toy that makes every number checkable by hand while preserving the saturating shape.
- [ ] **Step 2: Rung 2 — equal-marginal, applied (not re-proved).** State the equal-marginal / KKT condition $r_i'(b_i^\star)=\lambda^\star\ \forall i$ (cite **Ch. 12** by name; do not re-derive). For $r_i=a_i\sqrt b$ this gives $a_i/(2\sqrt{b_i^\star})=\lambda^\star$, hence $b_i^\star=B a_i^2/S$, recovering $b^\star=(7.2,1.8)$, $\lambda=0.3727$. New observation: the condition holds **per posterior draw**, so $p(\theta)$ induces $p(b^\star,\lambda^\star)$ — sets up Rungs 4–5.
- [ ] **Step 3: Rung 3 — Proof P2 (keystone): envelope theorem / shadow price.** State the theorem: at the optimum, $dV^\star/dB=\lambda^\star$. Prove via the envelope theorem on $L(b,\lambda)=V(b)+\lambda(B-\mathbf1^\top b)$: $dV^\star/dB=\partial L/\partial B|_\star=\lambda^\star$ (the $\partial/\partial b$ terms vanish by stationarity/KKT — cite **Ch. 12**; the active equality keeps feasibility). Use `aligned` inside `$ … $`. Then anchor analytically: $V^\star(B)=\sqrt{SB}$ (derive: substitute $b_i^\star=Ba_i^2/S$ into $\sum a_i\sqrt{b_i^\star}=\sqrt{B/S}\sum a_i^2=\sqrt{BS}$), so $dV^\star/dB=\sqrt S/(2\sqrt B)=\sqrt5/6=0.3727=\lambda^\star$ at $B=9$. State the interpretation: the multiplier we keep citing **is** the marginal value of a budget dollar — the engine of scenario analysis (Rung 6).
- [ ] **Step 4: Lint + commit** — `feat(ch18): theory rungs 1-3 (setup, equal-marginal, envelope theorem)`.

## Task 3 — Theory rungs 4–7 (propagation; optimizer's curse P1; scenario analysis; the wound)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (continue `## Theory & Proofs`).

- [ ] **Step 1: Rung 4 — propagating curve uncertainty.** Push $p(\theta\mid\text{data})$ (Ch. 16) through the optimizer: $b^\star(\theta),\lambda^\star(\theta),V^\star(\theta)$ each get a posterior. The Ch. 16 **Fisher-information ridge** (collinear adstock/saturation directions weakly identified) becomes a *wide* allocation posterior. State that point estimation collapses this to one confident allocation — the certainty-equivalent / plug-in plan — discarding real decision risk.
- [ ] **Step 2: Rung 5 — Proof P1 (keystone): the optimizer's curse.** State
  $\mathbb E_\theta[\max_b V(b;\theta)]\ge\max_b\mathbb E_\theta[V(b;\theta)]$. Prove: for any fixed $b'$, $\max_b V(b;\theta)\ge V(b';\theta)$ pointwise in $\theta$; take $\mathbb E_\theta$; then maximize the RHS over $b'$ (a max–expectation / Jensen-type argument; always true, no convexity-in-$\theta$ needed). Give the **three readings**: (1) **over-optimism** — optimizing on a point estimate and reporting $\max_b V(b;\hat\theta)$ overstates realized expected response; (2) **value of information** — the gap is the EVPI, exactly what resolving $\theta$ is worth; (3) **bridge to Part VI** — that gap is largest precisely on the weakly-identified ridge, so a confidently-wrong allocation is most costly where observational data is least able to resolve it. Cite `@smith2006` (optimizer's curse) and `@gelman2013` (Bayesian decision/EVPI).
- [ ] **Step 3: Rung 6 — scenario analysis via the shadow price.** What-if budgeting = re-optimize under perturbations. The envelope theorem gives the first-order answer *before* re-solving: $\Delta V\approx\lambda^\star\Delta B$; re-optimizing exposes where **curvature**/binding caps make it diverge. Then the two value-function tools: (i) **frontier** $V^\star(B)=\sqrt{SB}$ with decreasing marginal $\lambda^\star(B)=\sqrt S/(2\sqrt B)$ — the efficiency curve; (ii) **inverse** — to hit target $V_0$, solve $B=V_0^2/S$ (anchor: $V_0=8\Rightarrow B=12.8$, $\lambda^\star=5/16=0.3125$). Frontier and inverse are one identity read two ways.
- [ ] **Step 4: Rung 7 — the wound (closing rung).** Combine Rungs 4–6 under uncertainty: the *deterministic* scenario table looks crisp; the *posterior* version's rankings overlap and the shadow price carries a wide band — *the what-if table is more confident than the data warrants*. Two posterior-plausible curves give materially different optima (forward-reference WE3's mirror-image $b^\star$): **good fit ≠ good decision**. State plainly: observational MMM weakly identifies the curve, so the budget decision inherits that weak identification; the only way to shrink the EVPI gap is to *intervene and measure*. Hand to **Part VI** by name (causal foundations → quasi-experimental designs → calibration).
- [ ] **Step 5: Lint + commit** — `feat(ch18): theory rungs 4-7 (propagation, optimizer's curse, scenario analysis, the wound)`.

## Task 4 — Worked Examples (WE1, WE2, WE3)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (`## Worked Examples`).

- [ ] **Step 1: WE1 — closing the Part IV loop (deterministic).** Solve $r_i=a_i\sqrt b$, $a=(2,1)$, $B=9$ by the equal-marginal rule → $b^\star=(7.2,1.8)$, $\lambda=0.3727$, $V^\star=6.7082$ (show the arithmetic: $b_i^\star=Ba_i^2/S$, $5b_2=9$, etc.). Then verify the envelope theorem: $V^\star(B)=\sqrt{5B}$, $dV^\star/dB|_{B=9}=\sqrt5/6=0.3727=\lambda^\star$. State the payoff: the MMM allocation is *the same number* Part IV converged to — Part V supplied the curve Part IV assumed.
- [ ] **Step 2: WE2 — scenario analysis via the shadow price (deterministic).** Build the scenario table (use the Verified Anchors): $B=10.8$ (exact $7.3485$ vs first-order $7.3790$), $B=7.2$ ($6.0000$), cap $b_1\le5$ → $(5,4)$, $V=6.4721$, channel-2-dark → $(9,0)$, $V=6.0000$. Show the shadow price predicts the first column before any re-solve, and curvature explains the gap. Then the frontier ($V^\star(B)=\sqrt{5B}$) and the inverse ($V_0=8\Rightarrow B=12.8$, $\lambda^\star=0.3125$, verified by re-optimizing to $V^\star=8$).
- [ ] **Step 3: WE3 — the wound (under uncertainty).** Three parts: (a) **allocation posterior** — the ridge posterior on $a$ (describe the $u/v$ decomposition: total well-identified, attribution-contrast poorly identified) gives $b_1^\star$ median $\approx7.2$, 5/95 band $\approx[4.35,8.69]$. (b) **optimizer's curse numerically** — $\mathbb E[\max]\approx6.8276$ vs $\max\mathbb E=6.7082$, **gap $\approx0.1194$** = EVPI = what an experiment is worth (and the same number shows the in-sample optimum overstates the at-the-truth response). (c) **two plausible curves, two optima** — the exact $\theta=15^\circ/75^\circ$ pair: identical fit $V=6.75$, mirror-image $b^\star=(8.40,0.60)$ vs $(0.60,8.40)$. *Good fit ≠ good decision.*
- [ ] **Step 4: Lint + commit** — `feat(ch18): worked examples (Part IV loop, scenario analysis, the wound)`.

## Task 5 — Code Tie-in (controller writes directly)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (`## Code Tie-in`).

> Controller writes this task directly (streamlined pattern), then verifies headless. Single `{python}` cell, `numpy`+`scipy.optimize`+`matplotlib`, every claim asserted, figures end `plt.show()`.

- [ ] **Step 1: Write the cell** doing, in order:
  1. **Deterministic optimum (WE1):** `scipy.optimize.minimize(method="SLSQP")` on $-\sum a_i\sqrt{b_i}$, $a=(2,1)$, equality $\mathbf1^\top b=9$, bounds $b\ge0$. Assert `b*≈(7.2,1.8)` (atol 1e-2), response `≈6.7082` (1e-3), recovered multiplier `≈0.3727`. (Closes the Part IV loop on a live solver.)
  2. **Envelope check:** finite-difference $dV^\star/dB$ across a budget sweep using the analytic $V^\star(B)=\sqrt{5B}$; assert it tracks $\lambda^\star(B)=\sqrt5/(2\sqrt B)$ (atol 1e-3). Figure: $V^\star(B)$ with slope $=\lambda^\star$ annotated at $B=9$.
  3. **Scenario table + frontier + inverse (WE2):** compute the four scenario rows (closed-form), print first-order vs actual $\Delta V$; assert the anchors ($7.3485$, $6.0$, $6.4721$, $6.0$). Solve inverse $B=V_0^2/5$ for $V_0=8$; assert `B==12.8`, `lambda==0.3125`, and re-optimized `V*≈8`.
  4. **Uncertainty + optimizer's curse (WE3):** `rng=np.random.default_rng(18)`, ridge `Cov`, `N=200_000` draws, keep `>0`; closed-form `b*=B*a**2/S` per draw. Plot $b_1^\star$ posterior (hist) with the deterministic $7.2$ overlaid (the wound visible). Compute `Emax=mean(sqrt(S*B))`, `maxE=sqrt(5*B)`; assert `Emax-maxE>0.10`; print as EVPI. Construct the $\theta=15^\circ/75^\circ$ pair; assert both fit `≈6.75` and optima `≈(8.397,0.603)`/`(0.603,8.397)`.
- [ ] **Step 2: Verify headless** — `cd` to a temp dir, extract the cell, run `MPLBACKEND=Agg python3 cell.py`; confirm all asserts pass and it prints the EVPI gap.
- [ ] **Step 3: Lint + commit** — `feat(ch18): code tie-in — SLSQP optimum, envelope, scenarios, optimizer's curse`.

## Task 6 — Exercises (C / B / P / A — controller writes directly)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (`## Exercises`). Self-contained; **no inline solution links**.

- [ ] **Step 1: Write exercises.**
  - **C (Conceptual):** Why does optimizing on a point estimate overstate realized response? What does the EVPI gap mean operationally and why does an experiment shrink it? Why is good in-sample fit insufficient for a budgeting decision?
  - **B (By hand):** On $r_i=a_i\sqrt b$, $a=(2,1)$: derive $V^\star(B)=\sqrt{5B}$; confirm $dV^\star/dB=\lambda^\star$; compute budget-$+20\%$ first-order vs exact $\Delta V$; solve the inverse — budget to reach $V_0=8$ — and read off its shadow price ($B=12.8$, $\lambda^\star=5/16$).
  - **P (Prove it):** (i) Prove $\mathbb E[\max_b V]\ge\max_b\mathbb E[V]$ from the pointwise bound. (ii) Prove $dV^\star/dB=\lambda^\star$ for the equality-constrained problem (envelope theorem).
  - **A (Applied/code):** Extend the Code Tie-in — widen the ridge (raise `s_along`) and show the EVPI gap and allocation bands grow; add a minimum-spend floor $b_i\ge f$ and read its shadow price.
- [ ] **Step 2: Lint + commit** — `feat(ch18): exercises (C/B/P/A)`.

## Task 7 — Summary (auto-included; controller writes directly)

**Files:** Modify `parts/05-mmm-modeling/04-budget-optimization.qmd` (`## Summary` or the chapter's auto-included summary block — match Ch. 17's exact mechanism).

- [ ] **Step 1: Write Summary.** Bulleted **Key concepts** (the join of Part IV optimizer + Part V curve; shadow price; scenario analysis = shadow price at work; optimizer's curse; EVPI; the wound → Part VI). Bulleted **Key identities** (inline math, **bulleted** — not a paragraph):
  - equal-marginal: $r_i'(b_i^\star)=\lambda^\star$ for all funded $i$;
  - envelope / shadow price: $dV^\star/dB=\lambda^\star$;
  - value function: $V^\star(B)=\sqrt{SB}$ (with $S=\sum_i a_i^2$);
  - inverse scenario: $B=V_0^2/S$;
  - optimizer's curse: $\mathbb E[\max_b V]\ge\max_b\mathbb E[V]$, with the gap = EVPI.
- [ ] **Step 2: Lint + commit** — `feat(ch18): summary`.

## Task 8 — Appendix solutions block

**Files:** Modify `appendix/solutions.qmd` (append at end of file, after the Ch. 17 block; line ~2774).

- [ ] **Step 1: Append the block** `## Chapter 18 — Budget Optimization & Scenario Analysis`, wrapped in `::: {.content-visible when-meta="show-solutions"}` … `:::`, in chapter order (it is the last chapter). Provide full solutions:
  - **C:** prose answers (over-optimism = the in-sample optimum is selected for looking good on noisy estimates; EVPI = value of perfect information, shrinks because an experiment replaces the prior with data; fit measures description not decision-robustness).
  - **B:** full arithmetic — $V^\star(B)=\sqrt{5B}$ derivation; $dV^\star/dB=\sqrt5/(2\sqrt B)=0.3727$ at $B=9$; $+20\%$: first-order $7.3790$ vs exact $7.3485$; inverse $B=64/5=12.8$, $\lambda^\star=\sqrt5/(2\sqrt{12.8})=5/16=0.3125$.
  - **P:** the two proofs in full (optimizer's-curse pointwise→expectation→max argument; envelope theorem on the Lagrangian).
  - **A:** describe expected code behavior — wider ridge ⇒ larger EVPI gap and wider $b^\star$ band; a floor constraint's multiplier is its shadow price (value of relaxing the floor by a dollar).
- [ ] **Step 2: Lint** (even `$$` count in appendix; gated block closed) **+ commit** — `feat(ch18): appendix solutions`.

## Task 9 — Bibliography + final review + PR

**Files:** Modify `references.bib`; final lint of `parts/05-mmm-modeling/04-budget-optimization.qmd`.

- [ ] **Step 1: Add bib entry** to `references.bib`:
  ```bibtex
  @article{smith2006,
    title={The optimizer's curse: Skepticism and postdecision surprise in decision analysis},
    author={Smith, James E. and Winkler, Robert L.},
    journal={Management Science},
    volume={52}, number={3}, pages={311--322}, year={2006},
    publisher={INFORMS}
  }
  ```
  Confirm every `@key` cited in the chapter exists in `references.bib` (`@smith2006`, `@gelman2013`, optionally `@boyd2004`/`@nocedal2006`).
- [ ] **Step 2: Full chapter lint** — six template headings present and in order (Motivation, Theory & Proofs, Worked Examples, Code Tie-in, Exercises, Summary); H1 `# Budget Optimization & Scenario Analysis` intact; even `$$` count; no bare `\begin{align}`; no `\begin{psmallmatrix}`; no banned library names (word-boundary grep for `stan`, `pymc`, `orbit`, etc.); all citation keys valid.
- [ ] **Step 3: Re-run the code cell headless** one final time — confirm all asserts pass.
- [ ] **Step 4: Push branch + open PR** against `main` with `gh pr create`; summarize the chapter (two keystone proofs, the Part IV loop closure, scenario frontier+inverse, the wound → Part VI). Start a background CI-render watcher and a merge watcher:
  - render: poll `gh run list`/`gh run watch` for the PR's `quarto render` to a green conclusion.
  - merge: `until gh pr view <N> --json state -q .state | grep -qE 'MERGED|CLOSED'; do sleep 120; done` (run_in_background).
- [ ] **Step 5:** Report the PR URL and CI status. **User merges the PR** (never self-merge).

---

## Self-Review (controller, before dispatching)

1. **Spec coverage:** Rungs 1–7 → Tasks 2–3; WE1–3 → Task 4; Code Tie-in → Task 5; C/B/P/A → Task 6; Summary → Task 7; appendix → Task 8; `@smith2006` + CI → Task 9; scenario frontier+inverse → Rung 6 (Task 3) + WE2 (Task 4) + code step 3 (Task 5). All spec sections mapped. ✅
2. **Placeholders:** none — every anchor is a concrete NumPy-verified number; every proof has its argument sketched. ✅
3. **Consistency:** the square-root model $r_i=a_i\sqrt b$, $a=(2,1)$, $B=9$, $S=5$ is used identically across all tasks; ridge posterior params (`s_along=0.45`, `s_perp=0.08`, seed 18, N=200k) and the $\theta=15^\circ/75^\circ$ pair are stated once in Verified Anchors and referenced thereafter. ✅
