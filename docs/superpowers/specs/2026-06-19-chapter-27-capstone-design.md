# Chapter 27 — Putting It All Together (Capstone): Design Spec

**Date:** 2026-06-19
**Part:** VIII — Capstone (the final chapter; closes the book)
**File:** `parts/08-capstone/01-putting-it-together.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This is the capstone. Parts I–VII built and grounded every component of a marketing-mix-modeling
system; this chapter assembles them into **one concrete, end-to-end pipeline** and runs a single full
turn of the decision loop the book has been building toward: simulate → fit → check → calibrate →
optimize → test, and then close the loop by letting the optimizer tell you which experiment to run next.

**It introduces no new mathematics.** Its job is integration and judgment, not new theorems. Every stage
*names and cites* the result it leans on (proved in an earlier chapter) and re-exercises it in context;
nothing is re-proved. This is a deliberate, one-time departure from the book's six-section chapter
template — the proof-heavy Part VII has just ended, and the last chapter should read like a practitioner
closing the loop, not a mathematician.

**Center of gravity (locked, via brainstorming):** one concrete pipeline, end to end, realized as a
**literate walkthrough** — prose explains each stage, a small code cell runs it, a hand-off sentence
carries the result into the next stage. The intellectual payoff is the **wiring**: how the keystones from
each part interlock into a closed decision loop.

**Driving question:** *"You have all the pieces — a fitted model, a calibration store, an optimizer, a
test harness. How do they compose into one system that turns data into a budget decision, and improves
itself?"*

**No anchors line** in this chapter (it is not anchored to a single canonical source); it synthesizes the
whole book.

**Styling (locked):** Part IV formula/spacing discipline still applies to any displayed math
(`$$` on their own lines, blank lines around displays). But there are **no proof ladders and no
`$\blacksquare$`** — the body is stage subsections, not theorem rungs.

---

## 2. Scope discipline & confidentiality guardrail

- **Re-uses only what the book already taught.** Every component is standard, published technique already
  developed in Ch.1–26 (conjugate Bayesian regression, adstock/saturation, SLSQP allocation, the envelope
  theorem / shadow prices, precision-pooled lift calibration, EVPI, metamorphic testing). The capstone
  introduces no new disclosure beyond the assembly.
- **Confidentiality guardrail (locked, user directive).** The "Closing the Loop" section stays at the
  **textbook-principle** level — decision-theory "rank channels by shadow-price-weighted EVPI." It must
  **NOT** include: any firm-specific experiment-selection heuristic; production hyperparameters, priors,
  or tuning; internal architecture, naming, or workflow that mirrors a real implementation; or anything
  but **synthetic** data and numbers. Principle, not playbook.
- **House rules (critical).** `numpy` + `scipy` + `matplotlib` + Python stdlib only. **Never name
  PyMC-Marketing or any MMM/PPL/sampler/causal library** (no PyMC, Stan, Orbit, numpyro, pyro,
  lightweight_mmm, causalimpact, networkx, hypothesis, pytest, sqlalchemy, kafka, pandas). Name methods,
  not authors/tools. NEVER `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their
  own lines; even `$$` count.

---

## 3. Chapter structure (template deliberately bent)

The rigid six-section template is replaced by:

1. **Motivation** (short) — the pieces exist; the question is composition and self-improvement.
2. **The Assembled System** — the literate walkthrough; this single section absorbs what were *Theory &
   Proofs + Worked Examples + Code Tie-in*. Six **stage subsections**, each `### Stage N — Title`, each a
   **prose → small `{python}` cell → hand-off** triple. Cells share one Quarto kernel (state flows
   forward); each cell prints and `assert`s its one key number.
3. **Closing the Loop** — the one genuinely new (systems) idea, generic-principle level: shadow prices
   rank which channel's uncertainty is most costly → that names the next experiment → it folds back into
   the store → EVPI falls monotonically (Ch.22). A small cell runs one extra turn and shows the EVPI
   staircase; the summary figure lives here.
4. **Exercises** — C/B/P/A headings retained for consistency, but **weighted to C (judgment) and A
   (extend the pipeline)**; minimal/no P (the proofs were earned in earlier chapters).
5. **Summary** — bulleted **Key concepts** + bulleted **Key takeaways** (not "Key identities" — this
   chapter has few new identities), closing with a short **"where to go next"** sendoff that ends the book.

---

## 4. The Assembled System — the six stages

One synthetic two-channel MMM carried end to end. Each stage cites its owner chapter(s).

**Stage 1 — Simulate the data-generating process.** Generate weekly spend for two channels, apply
**adstock** (geometric carryover) and **saturation** (a concave response), add noise, produce sales
(Ch.14 DGP, Ch.15). Print a few rows of the synthetic series; `assert` shapes and a known summary stat.
Hand-off: this is the ground truth the rest of the pipeline must recover.

**Stage 2 — Fit the model.** Transform the features (adstock+saturation), then compute a **closed-form
conjugate / ridge posterior** for the coefficients — deterministic, fast, and it yields a real posterior
mean *and covariance* to propagate (Ch.5 Bayesian inference, Ch.16 ridge; the sampler chapters Ch.8–9 are
cited as the general-purpose alternative). `assert` the posterior mean recovers the DGP coefficients
within tolerance. Hand-off: we have a posterior, but is it good enough to *act on*?

**Stage 3 — Check, and find the wound.** A posterior-predictive residual check (Ch.10) confirms fit, but
the spend/attribution **contrast** is poorly identified — exactly the high-condition-number ridge
direction Ch.23 named ($\kappa(X^\top X) \approx \kappa(X)^2$). Propagate the posterior through the
allocation to show the budget posterior is wide and the **optimizer's-curse EVPI is positive** (Ch.18).
`assert` EVPI $> 0$ (anchor $\approx 0.12$). Hand-off: the model fits but the *decision* is not yet
trustworthy — the EVPI says an experiment is worth running.

**Stage 4 — Calibrate from the prior store.** Estimate a **lift-test factor** (Ch.20 interventional
estimand, Ch.21 mechanics) and **fold its precision** into the posterior the way Ch.22's store does;
watch the EVPI fall (anchor $0.12 \to 0.03$). `assert` the post-fold EVPI $<$ pre-fold EVPI and matches
the anchor. Hand-off: the posterior is now tight enough that the allocation will be trustworthy — the
payoff of having built the store.

**Stage 5 — Optimize the budget.** Run **SLSQP** allocation on the now-tight posterior mean (Ch.13
solver, Ch.18 allocation), recover the optimal split and the **budget shadow price** $\lambda^\star$ as
the common marginal return (envelope theorem). `assert` the allocation and $\lambda^\star$ anchors.
Hand-off: we have a defensible budget — and a shadow price per channel.

**Stage 6 — Test the loop.** Wrap the pipeline in a **metamorphic test** (Ch.26): assert an invariant the
whole system must honor — e.g. permutation-invariance of the calibration fold, or monotonicity of EVPI
under added information — standing on the Ch.23 (conditioning), Ch.24 (acyclic module order), Ch.25
(idempotent/replayable store) grounding. `assert` the invariant holds. Hand-off: the assembled system is
verified to honor its own theorems.

---

## 5. Closing the Loop (generic-principle)

The shadow prices from Stage 5, weighted by each channel's remaining posterior uncertainty, rank which
channel's uncertainty is **most costly to the decision** — that is the channel to experiment on next.
Folding that next experiment into the store (Stage 4 again) tightens the posterior further, and by the
Ch.22 monotonicity result the **EVPI falls monotonically** with each turn. A short cell runs one
additional turn and prints the EVPI staircase (e.g. $0.12 \to 0.03 \to 0.018$); the **summary figure**
shows the EVPI staircase across loop turns (and/or the before/after allocation posterior).

State the principle plainly and stop there: this is decision theory (value of information driving
experiment choice), not any specific operational recipe. This closed loop — data → posterior →
calibration → allocation → measurement → back to the store — **is the book**.

---

## 6. Code realization & verification

- **Multiple small `{python}` cells** (one per stage + the loop cell), interleaved with prose, sharing one
  Quarto kernel so later cells use earlier variables. This intentionally breaks the book's prior
  one-cell-per-chapter convention — appropriate and expected for a literate capstone.
- Each cell **prints its key number(s) and `assert`s them**, so CI catches any drift.
- **Fit is closed-form conjugate/ridge** — fully deterministic, no seed-flakiness; any stochastic step
  (DGP noise, Monte-Carlo EVPI) is seeded `np.random.default_rng(27)` with tolerances sized to sampling
  noise (Ch.26 Rung 4).
- `numpy` + `scipy.optimize` (SLSQP, the Ch.13 routine) + `matplotlib` + stdlib only. One summary figure;
  cells/figure end appropriately for headless run (`MPLBACKEND=Agg python3`).
- **All numeric anchors NumPy-verified before writing** (DGP coeffs recovered; EVPI staircase
  $0.12 \to 0.03 \to \dots$; SLSQP allocation + $\lambda^\star$; metamorphic invariant holds). The plan
  pins exact anchor numbers.

---

## 7. Exercises (C / B / P / A — self-contained, no inline solution links; C/A-weighted)

- **C (judgment, emphasized):** Walk the loop in words — what does each stage hand the next? The optimizer
  says reallocate 30% to search; name two things that would make you *not* trust that recommendation
  (point to the stage that would catch each). Why does the EVPI fall monotonically, and what does that
  guarantee about the loop?
- **B:** Given a (pre, post) EVPI pair and a shadow-price vector, compute which channel to experiment on
  next and the expected EVPI after one more fold.
- **P (light):** Show the assembled loop is *consistent* — adding an experiment never increases EVPI
  (cite Ch.22's Loewner-monotonicity result; no new proof required, just the application).
- **A (extend the pipeline, emphasized):** Add a third channel to the DGP and rerun the whole pipeline;
  change the lift-test precision and trace how the allocation and EVPI move; swap the metamorphic invariant
  for a different one from Ch.26 and confirm it still passes.

---

## 8. Appendix solutions

Append `## Chapter 27 — Putting It All Together` to `appendix/solutions.qmd`, in chapter order (after the
Ch.26 block), inside the existing `content-visible when-meta="show-solutions"` div. Full C/B/P/A; the C
answers carry the most weight (judgment narrative), the A answer sketches the extended pipeline.

---

## 9. Summary (auto-included)

Bulleted **Key concepts** (the assembled pipeline; each stage's owner chapter; the closed decision loop;
EVPI-driven experiment selection; tests as the loop's guard) + bulleted **Key takeaways** (the loop is
contractive in EVPI; the theorems are the spec and the tests; a budget decision is only as trustworthy as
its least-identified direction). End with a short **book sendoff**.

---

## 10. Cross-references

Touches every part: I (linear algebra/calculus/probability underpin all stages), II (regression/Bayes —
the fit), III (sampling — cited as the general fit alternative; model checking), IV (convexity/LP/SLSQP —
the allocator), V (DGP/building/DLM/budget optimization), VI (causal foundations/quasi-experimental/
advanced calibration/prior store — the calibration stage and the loop), VII (CS foundations/architecture/
data engineering/testing — the grounding and the Stage-6 test).

---

## 11. Bibliography

No new references expected — the chapter cites results already in `references.bib`. The plan confirms every
`[@key]` used resolves; add nothing unless a cited result lacks an entry (it should not).

---

## 12. Verification (the real gate)

- Every numeric anchor verified in NumPy; the chained cells run top-to-bottom headless
  (`MPLBACKEND=Agg python3`) with all `assert`s passing and no third-party library.
- KaTeX/structure lint: even `$$` count, no `\begin{align}`, no `psmallmatrix`, valid citation keys, no
  banned library names, no `$\blacksquare$` (this chapter proves nothing new), Part IV display spacing.
- Confidentiality guardrail satisfied: "Closing the Loop" is principle-level only; no firm-specific
  heuristic/hyperparameters/architecture/real data.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges. This is the final chapter;
  on merge the book is complete.

---

## 13. Build note

Built in a fresh `ch27-capstone` worktree off the current `main` (which already contains the Ch.24/25/26
appendix blocks), so the Ch.27 solutions block appends in correct chapter order with no rebase needed.
