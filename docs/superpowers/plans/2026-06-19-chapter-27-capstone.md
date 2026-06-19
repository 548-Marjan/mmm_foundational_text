# Chapter 27 — Putting It All Together (Capstone) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Write Chapter 27, the capstone: assemble Ch.1–26 into one concrete end-to-end MMM pipeline (simulate → fit → check → calibrate → optimize → test) realized as a literate walkthrough, then close the decision loop. No new theorems.

**Architecture:** Replace the stub `parts/08-capstone/01-putting-it-together.qmd` with: Motivation → **The Assembled System** (six `### Stage N` subsections, each prose → small `{python}` cell → hand-off; cells share one Quarto kernel) → **Closing the Loop** → Exercises (C/A-weighted) → Summary + book sendoff. Append a gated appendix block. No new bib entries expected.

**Tech Stack:** Quarto `.qmd` (KaTeX), multiple `{python}` cells, `numpy` + `scipy.optimize` (SLSQP) + `matplotlib` + stdlib only.

---

## Conventions (enforced every task)

- **No new math.** This chapter cites results, never re-proves them. **No `$\blacksquare$`, no proof ladders.** Any displayed math uses Part IV spacing (`$$` on own lines, blank lines around displays, even `$$` count).
- **House rules.** `numpy`/`scipy`/`matplotlib`/stdlib only. **Never name** PyMC-Marketing or any MMM/PPL/sampler/causal library (PyMC, Stan, Orbit, numpyro, pyro, lightweight_mmm, causalimpact, networkx, hypothesis, pytest, sqlalchemy, kafka, pandas). Name methods, not authors/tools. No `\begin{psmallmatrix}`; `aligned` only inside `$$`.
- **Confidentiality guardrail (locked).** "Closing the Loop" is generic decision-theory principle only — "rank channels by shadow-price-weighted EVPI." NO firm-specific heuristic, production hyperparameters/priors/tuning, internal architecture/naming, or real data. Synthetic numbers only.
- **Locked anchors (NumPy-verified, seed 0; the whole pipeline is in `/tmp/ch27_pipeline.py`, verified headless):** corr(f₁,f₂)=0.842, κ(X)=32.2, κ(XᵀX)≈κ(X)²≈1037; posterior mean a=[2.367, 0.590] (contrast wrong = the wound), contrast-dir std 0.542 vs other-dir 0.158; EVPI 0.1432 (>0) → 0.0072 after one lift fold; SLSQP b\*=[7.259, 1.741], V\*=6.714, λ\*=0.373; EVPI-monotone invariant True; loop staircase 0.1432 → 0.0072 → 0.0060.
- **Determinism.** Single seed `np.random.default_rng(0)` for the DGP; EVPI Monte-Carlo uses its own fixed seed (3) inside the helper; ndraw=400_000; asserts use tolerances ~0.01 (sized to MC noise, Ch.26 Rung 4), never exact equality on MC quantities.
- **Git identity** jlh530i / jlh530i@gmail.com. Worktree `ch27-capstone` (off latest main; Ch.26 already in appendix, no rebase). Branch→PR→user merges; one complete commit-set.

---

## Task 1: Motivation

**Files:** Modify `parts/08-capstone/01-putting-it-together.qmd` (keep H1 `# Putting It All Together`; no anchors line).

- [ ] **Step 1:** Replace the stub callout + template with `## Motivation` (3–4 short paragraphs):
  - Parts I–VII built and grounded every component; the capstone assembles them into one concrete pipeline and runs a full turn of the decision loop.
  - This chapter proves nothing new — its job is composition and judgment. It deliberately breaks the six-section template: the body is a literate walkthrough.
  - State the driving question: *"You have all the pieces — a fitted model, a calibration store, an optimizer, a test harness. How do they compose into one system that turns data into a budget decision, and improves itself?"*
  - Preview the six stages and the closing loop in one sentence.
- [ ] **Step 2:** Verify: `grep -nE '^# |^## ' file` shows H1 + `## Motivation`, no stub callout.
- [ ] **Step 3:** Commit `feat(ch27): motivation`.

---

## Task 2: The Assembled System — Stages 1–3 (DGP, fit, find-the-wound)

**Files:** Modify the chapter (add `## The Assembled System` + Stages 1–3, each prose → `{python}` cell → hand-off).

The verified code for all stages lives in `/tmp/ch27_pipeline.py` (headless-clean). Split it across cells; **cells share one Quarto kernel**, so later cells use earlier variables (`X`, `a_hat`, `Ca`, `B`, `evpi_from_posterior`, …). Cell 1 carries the imports + `adstock`; the EVPI helper is defined in the Stage-3 cell (before first use).

- [ ] **Step 1: `## The Assembled System` intro** — one short paragraph: one synthetic two-channel MMM carried end to end; each stage names the chapter it stands on; the wiring is the point.

- [ ] **Step 2: Stage 1 — Simulate the DGP.** Prose: weekly spend for two channels, geometric **adstock** + sqrt **saturation**, noise → sales (Ch.14/15); the two channels' spends are collinear, so the spend/attribution *contrast* will be poorly identified — the ridge. Cell (from `/tmp/ch27_pipeline.py` Stage 1; imports here): build `x1,x2,f1,f2,X,y`; print `T`, `corr(f1,f2)=0.842`, `κ(X)=32.2`; `assert X.shape==(156,2)` and `abs(np.corrcoef(f1,f2)[0,1]-0.842)<0.01`. Hand-off: this is the ground truth the pipeline must recover.

- [ ] **Step 3: Stage 2 — Fit (conjugate).** Prose: transform features, closed-form **conjugate/ridge posterior** over (base, a₁, a₂) — deterministic, yields mean *and covariance* (Ch.5/16; samplers Ch.8–9 are the general alternative). Cell: build `Xd, cov, mean, a_hat, Ca`; print mean a=[2.367, 0.590], κ(XᵀX)≈κ(X)² (compute `cond(X)` and `cond(X.T@X)` on the **same 2-col `X`** so the squaring is exact: 32.2 and ≈1037), and the eigen-split of `Ca` (contrast-dir std 0.542 vs 0.158). `assert` κ(X.T@X) ≈ cond(X)**2 within 1%; `assert` contrast std > 3× the other. Hand-off: we have a posterior — but the point estimate's contrast is untrustworthy; is it good enough to *act on*?

- [ ] **Step 4: Stage 3 — Check + find the wound.** Prose: a residual check confirms fit, but the contrast is the high-κ ridge direction Ch.23 named ($\kappa(X^\top X)\approx\kappa(X)^2$); propagate the posterior through the B=9 allocation → wide budget posterior, **optimizer's-curse EVPI > 0** (Ch.18). Define the **Jensen-safe** EVPI helper (closed form `EVPI = sqrt(B)*(E||a|| - ||mu||) >= 0`, **no positivity filter**) in this cell. Print `EVPI=0.1432`; `assert 0.13 < evpi0 < 0.16`. Hand-off: the model fits but the decision is not yet trustworthy — the EVPI says an experiment is worth running.

- [ ] **Step 5:** Verify each cell extracts and runs headless (see Task 6 verification); even `$$` count.
- [ ] **Step 6:** Commit `feat(ch27): assembled system stages 1-3 (DGP, fit, wound)`.

---

## Task 3: The Assembled System — Stages 4–6 (calibrate, optimize, test)

**Files:** Modify the chapter (add Stages 4–6).

- [ ] **Step 1: Stage 4 — Calibrate from the prior store.** Prose: a **lift test** measures the poorly-identified contrast direction (Ch.20 estimand, Ch.21 mechanics); **fold its precision** into the posterior the way Ch.22's store does; EVPI falls. Cell (`fold_lift` helper + apply on `u=(1,-1)`, `s=0.10`): print `EVPI: 0.1432 -> 0.0072`; `assert evpi1 < evpi0` and `evpi1 < 0.02`. Hand-off: the posterior is now tight enough to act on — the payoff of the store.

- [ ] **Step 2: Stage 5 — Optimize (SLSQP).** Prose: run **SLSQP** allocation (Ch.13) on the now-tight posterior mean, recover the split and the budget **shadow price** λ\* as the common marginal return (envelope theorem, Ch.18). Cell (`neg_response`, `solve_budget`, apply to `a_hat1`, B=9): print b\*=[7.259, 1.741], V\*=6.714, λ\*=0.373; `assert np.allclose(b_star,[7.259,1.741],atol=0.02)`, `abs(V_star-6.714)<0.01`, `np.allclose(lam,lam[0],atol=1e-3)` and `abs(lam[0]-0.373)<0.01`. Hand-off: a defensible budget — and a shadow price per channel.

- [ ] **Step 3: Stage 6 — Test the loop.** Prose: wrap the pipeline in a **metamorphic test** (Ch.26) — assert the invariant the system must honor: EVPI is **monotone** under an information fold (Ch.22), standing on Ch.23/24/25 grounding. Cell: `assert evpi1 <= evpi0 + 1e-9`; print the invariant holds. Hand-off: the assembled system is verified to honor its own theorems.

- [ ] **Step 4:** Verify; even `$$` count.
- [ ] **Step 5:** Commit `feat(ch27): assembled system stages 4-6 (calibrate, optimize, test)`.

---

## Task 4: Closing the Loop (+ summary figure)

**Files:** Modify the chapter (add `## Closing the Loop`).

- [ ] **Step 1:** Prose (generic-principle, guardrail-compliant): shadow prices from Stage 5, weighted by each channel's remaining posterior uncertainty, rank which channel's uncertainty is **most costly to the decision** — that names the next experiment. Fold it back (Stage 4 again); by Ch.22 monotonicity the **EVPI falls monotonically** each turn. State plainly this is value-of-information decision theory, not an operational recipe; the closed loop data → posterior → calibration → allocation → measurement → store **is the book.**
- [ ] **Step 2:** Cell: run one more turn (fold `v=(1,1)`, `s=0.10`) → print staircase `0.1432 -> 0.0072 -> 0.0060`; `assert evpi2 <= evpi1 + 1e-9`. Then the **summary figure**: EVPI staircase across loop turns (step plot, turns 0/1/2 → 0.143/0.007/0.006), and/or before/after allocation. End the figure cell with `plt.show()`.
- [ ] **Step 3:** Verify headless run produces the figure with no exception.
- [ ] **Step 4:** Commit `feat(ch27): closing the loop + summary figure`.

---

## Task 5: Exercises + Summary

**Files:** Modify the chapter (add `## Exercises`, `## Summary`).

- [ ] **Step 1: `## Exercises`** (C/B/P/A headings; C/A-weighted; self-contained, no solution links):
  - **C (emphasized):** Walk the loop in words — what does each stage hand the next? The optimizer says reallocate 30% to search; name two things that would make you *not* trust that (and the stage that would catch each). Why does EVPI fall monotonically, and what does that guarantee?
  - **B:** Given a (pre, post) EVPI pair and a shadow-price vector, say which channel to experiment on next and the expected EVPI after one more fold.
  - **P (light):** Show the loop is consistent — a fold never increases EVPI (cite Ch.22's Loewner-monotonicity; application, no new proof).
  - **A (emphasized):** Add a third channel to the DGP and rerun the whole pipeline; vary the lift-test precision and trace allocation + EVPI; swap the metamorphic invariant for another from Ch.26 and confirm it still passes.
- [ ] **Step 2: `## Summary`** — bulleted **Key concepts** (the assembled pipeline; each stage's owner chapter; the closed decision loop; EVPI-driven experiment selection; tests as the loop's guard) + bulleted **Key takeaways** (the loop is contractive in EVPI; the theorems are the spec *and* the tests; a budget decision is only as trustworthy as its least-identified direction). End with a short **book sendoff**.
- [ ] **Step 3:** Verify section headings present; even `$$` count; **no `\blacksquare`** anywhere in the chapter (`grep -c blacksquare` → 0).
- [ ] **Step 4:** Commit `feat(ch27): exercises + summary + book sendoff`.

---

## Task 6: Appendix solutions

**Files:** Modify `appendix/solutions.qmd` (append `## Chapter 27 — Putting It All Together` after the Ch.26 block, inside the gated div).

- [ ] **Step 1:** Confirm insertion point: `grep -n '^## Chapter 26' appendix/solutions.qmd` exists; find the final `:::` closing the `content-visible` div; insert before it.
- [ ] **Step 2:** Append the Ch.27 block, full C/B/P/A (C and A carry the most weight):
  - **C:** the loop narrative — each stage's hand-off; two trust-breakers (e.g. a failed posterior-predictive check at Stage 3; a wide post-calibration contrast still giving EVPI ≫ 0 → don't trust the split yet) and the catching stage; EVPI monotone because each fold only *adds* precision (Ch.22 Loewner), guaranteeing the loop never makes the decision worse.
  - **B:** pick the channel with the largest shadow-price-weighted posterior variance; expected post-fold EVPI from the precision addition.
  - **P:** one paragraph applying Ch.22 monotonicity (no new proof, **no `$\blacksquare$`** — or a single cited application line).
  - **A:** sketch the three-channel extension and what moves.
- [ ] **Step 3:** Verify even `$$` count in `appendix/solutions.qmd`; block sits before the final `:::`.
- [ ] **Step 4:** Commit `feat(ch27): appendix solutions (C/B/P/A)`.

---

## Task 7: Final review + PR

**Files:** read-only review + push.

- [ ] **Step 1: Bib check** — `grep -oE '@[a-z0-9]+' chapter` ; confirm every cited key exists in `references.bib`. Add nothing unless a cited result lacks an entry (it should not — all cited results are from prior chapters).
- [ ] **Step 2: Lint:**
  ```bash
  grep -c '\$\$' parts/08-capstone/01-putting-it-together.qmd     # even
  grep -c '\$\$' appendix/solutions.qmd                            # even
  grep -c blacksquare parts/08-capstone/01-putting-it-together.qmd # 0 (no new proofs)
  grep -nE 'begin\{align\}|psmallmatrix' parts/08-capstone/01-putting-it-together.qmd  # none
  grep -niE 'pymc|stan|hypothesis|pytest|networkx|pandas|numpyro|pyro|sqlalchemy|kafka|lightweight' parts/08-capstone/01-putting-it-together.qmd  # none (real libs; ignore con-stan-t/in-stan-ce false positives)
  grep -nE '^## ' parts/08-capstone/01-putting-it-together.qmd     # Motivation, The Assembled System, Closing the Loop, Exercises, Summary
  ```
  Confirm the confidentiality guardrail by reading "Closing the Loop": principle-level only.
- [ ] **Step 3: Headless run** — extract ALL `{python}` cells in document order, concatenate, run `MPLBACKEND=Agg python3` (mimicking the shared kernel). Expect all asserts pass and the staircase/allocation/EVPI anchors print. (Verified reference: `/tmp/ch27_pipeline.py`.)
- [ ] **Step 4: Push + PR:**
  ```bash
  git push -u origin worktree-ch27-capstone
  gh pr create --title "Chapter 27 — Putting It All Together (Capstone)" --body "<summary + test plan; note: final chapter, completes the book>"
  ```
- [ ] **Step 5:** Start a background watcher on the render check + a merge-watcher; report PR number; watch render to green. Do **not** push follow-up commits after opening. User merges — on merge **the book is complete.**

---

## Self-Review

1. **Spec coverage:** Motivation (T1), 6 stages (T2–T3), Closing the Loop + figure (T4), Exercises C/A-weighted + Summary + sendoff (T5), appendix (T6), bib/lint/headless/PR (T7). ✓
2. **Placeholders:** every stage's prose + cell is specified against the verified `/tmp/ch27_pipeline.py`; anchors pinned. ✓
3. **Anchor consistency:** corr 0.842, κ(X) 32.2, κ² ≈1037, mean [2.367,0.590], EVPI 0.1432→0.0072, b\* [7.259,1.741], λ\* 0.373, staircase 0.1432→0.0072→0.0060 — identical across plan, cells, prose, appendix. ✓
4. **House rules + guardrail:** no banned libs, even `$$`, no `\blacksquare` (no new proofs), Closing-the-Loop principle-only, deterministic seed, branch→PR→user-merges. ✓
