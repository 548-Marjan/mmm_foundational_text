# Chapter 18 — Budget Optimization & Scenario Analysis: Design Spec

**Date:** 2026-06-16
**Part:** V — MMM Modeling (closing chapter)
**File:** `parts/05-mmm-modeling/04-budget-optimization.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter **closes Part V on the wound**. Parts IV and V have, by now, each delivered half of a
budgeting answer and left the other half implicit:

- **Part IV (Ch. 12–14)** built the *optimization machinery* — the equal-marginal / KKT condition
  (Ch. 12), LP duality and shadow prices (Ch. 13), and the SLSQP solver (Ch. 14) — and proved all
  three converge on the same allocation $b^\star=(7.2,1.8)$, $\lambda=0.3727$, response $6.7082$ for a
  fixed, *known* response curve.
- **Part V (Ch. 15–17)** built the *response curve itself* — the saturating, adstocked
  data-generating process (Ch. 15), its Bayesian fit and the Fisher-information **ridge** that makes
  the curve only weakly identified from observational spend (Ch. 16), and its time-varying
  generalization (Ch. 17).

This chapter joins them: it runs Part IV's optimizer on Part V's *fitted, uncertain* curve. The new
content is **not** re-deriving optimization (we cite Ch. 12–14) — it is the **decision-under-uncertainty
layer** that fitting a curve makes both possible and necessary, and the honest accounting of what
that uncertainty does to the budgeting decision.

The chapter ends on a **wound**: optimizing spend over a response curve the observational time series
only weakly identifies produces an allocation that *looks* decisive but is over-confident, and an
in-sample optimum that *overstates* the response it will actually realize. The size of that gap is
exactly the value of running an experiment to resolve the uncertainty — which is what **Part VI**
(causal foundations, quasi-experimental designs, calibration) exists to do. The reorg note is
explicit: "Part V's budget-optimization chapter ends on a *wound* … Part VI exists to heal it."

**Through-line driving question:** *"We can solve for the optimal budget. But should we trust the
answer?"*

---

## 2. What this chapter does NOT do (scope discipline)

- It does **not** re-prove the equal-marginal / KKT condition, LP duality, or SLSQP convergence.
  Those are Part IV results; we **cite** Ch. 12–14 by name and **apply** them.
- It does **not** introduce a new solver. The Code Tie-in uses `scipy.optimize.minimize(method="SLSQP")`
  (citable; a general tool, not an MMM/PPL library) exactly as Ch. 14 did.
- It does **not** re-derive the saturation/adstock curve. That is Ch. 15's DGP; we reuse it.
- House rules apply: **never** name PyMC-Marketing or any MMM/PPL/sampler library; `numpy`/`scipy`/
  `matplotlib` only. **Never** use `\begin{psmallmatrix}` (PDF-breaking). KaTeX rules:
  `aligned` inside `$ … $` only, `$$` delimiters on their own lines, even `$$`-count per file.

---

## 3. Theory & Proofs — the rungs

The chapter has **two proof keystones** plus an applied bridge from Part IV. The numbering below is
the intended rung order; "(keystone)" marks the two load-bearing proofs.

**Rung 1 — The budgeting problem, restated on a fitted curve.**
Set up $\max_{b\ge 0,\ \mathbf 1^\top b = B}\ V(b;\theta)=\sum_i r_i(b_i;\theta_i)$ where $r_i$ is the
*fitted* MMM response (saturating, post-adstock; Ch. 15) and $\theta$ are its **estimated**
parameters carrying a posterior (Ch. 16). Contrast with Part IV, where $\theta$ was *given*. This is
the whole shift: the objective is now random because $\theta$ is uncertain.

**Rung 2 — Equal-marginal allocation, applied (not re-proved).**
State the equal-marginal / KKT optimum from Ch. 12 and note it recovers $b^\star=(7.2,1.8)$,
$\lambda=0.3727$ on Part IV's square-root anchor. New observation: the condition holds **per posterior
draw**, so a posterior over $\theta$ induces a posterior over $(b^\star,\lambda^\star)$. Light;
mostly a bridge that closes the Part IV loop and sets up Rung 4.

**Rung 3 — Proof P2 (keystone): the envelope theorem / shadow price.**
At the optimum, the marginal value of the budget equals the optimal multiplier:
$$\frac{dV^\star}{dB}=\lambda^\star.$$
Proof via the envelope theorem on the Lagrangian $L(b,\lambda)=V(b)+\lambda(B-\mathbf 1^\top b)$:
$dV^\star/dB=\partial L/\partial B$ evaluated at the optimum $=\lambda^\star$ (the primal-feasibility
derivatives vanish by stationarity / KKT). Rests on Ch. 12's KKT; the **new content is the
shadow-price reading** — the multiplier we keep citing *is* the answer to "what is the next budget
dollar worth?" **Anchor (analytic, exact):** for $r_i(b)=a_i\sqrt b$, $a=(2,1)$, the value function
is $V^\star(B)=\sqrt{S\,B}$ with $S=\sum_i a_i^2=5$, so $V^\star(9)=\sqrt{45}=6.7082$ and
$dV^\star/dB=\sqrt S/(2\sqrt B)=\sqrt5/6=0.3727=\lambda^\star$ — the shadow price *is* the budget
derivative, on the very anchor Part IV converged to. This rung is the mathematical engine of
**scenario analysis** (Rung 6).

**Rung 4 — Propagating curve uncertainty to allocation uncertainty.**
Push the Ch. 16 posterior over $\theta$ through the optimizer: $\theta\sim p(\theta\mid\text{data})
\Rightarrow b^\star(\theta),\ \lambda^\star(\theta),\ V^\star(\theta)$ each have a posterior. The
weakly-identified **Fisher ridge** of Ch. 16 (collinear adstock/saturation directions) becomes a
*wide* posterior on the allocation. State that point estimation collapses this to a single confident
allocation — the certainty-equivalent (plug-in) plan — and that doing so discards real decision risk.

**Rung 5 — Proof P1 (keystone): the optimizer's curse.**
The central result. For an objective $V(b;\theta)$ with uncertain $\theta$,
$$\mathbb{E}_\theta\!\left[\max_b V(b;\theta)\right]\ \ge\ \max_b \mathbb{E}_\theta\!\left[V(b;\theta)\right].$$
*Proof:* for any fixed $b'$, $\max_b V(b;\theta)\ge V(b';\theta)$ pointwise in $\theta$; take
$\mathbb E_\theta$ of both sides, then maximize the right side over $b'$. (A clean max–expectation /
Jensen-type argument; always true, no convexity-in-$\theta$ assumption.) **Three decision-theoretic
readings**, all flowing from this one inequality:
  1. **Over-optimism.** Optimizing on a *point estimate* and reporting that value
     ($\max_b V(b;\hat\theta)$) overstates the response the chosen plan will realize in expectation
     — the in-sample optimum is an optimistically biased estimate of out-of-sample performance.
  2. **The value of resolving uncertainty.** The gap is the expected value of perfect information
     (EVPI): LHS is the value if you could optimize *knowing* $\theta$; RHS is the best you can do
     *before* knowing it. The difference is exactly what an experiment that pins down $\theta$ is
     worth.
  3. **The bridge to Part VI.** That EVPI gap is large **precisely when the curve is weakly
     identified** (the Ch. 16 ridge). A confidently-wrong allocation is most costly exactly where the
     observational data is least able to resolve it — which is the case for the EVPI-paying
     experiments of Part VI.

**Rung 6 — Scenario analysis, made rigorous by the envelope theorem.**
"What-if" budgeting: re-optimize under perturbed budgets and constraints (budget $\pm\Delta B$, a
channel cap $b_i\le c$, a channel going dark $b_i=0$). The envelope theorem (Rung 3) gives the
**first-order** scenario answer *before* re-optimizing: $\Delta V\approx \lambda^\star\,\Delta B$.
Re-optimizing then verifies it and exposes where **curvature** (concavity / a binding cap) makes the
realized $\Delta V$ diverge from the linear prediction. Scenario analysis is thus not a detour: it is
the shadow price put to work. Two further scenario tools follow directly from the value function:
  - **The response-vs-budget frontier.** Sweeping $B$ traces $V^\star(B)=\sqrt{SB}$ and its slope
    $\lambda^\star(B)=\sqrt S/(2\sqrt B)$ — the *efficiency curve* of the spend plan, with the shadow
    price as its (monotonically decreasing) marginal. Every budget question is a point on this curve.
  - **The inverse scenario.** "What budget hits a target response $V_0$?" inverts the value function:
    $B=V_0^2/S$. At that budget the shadow price is read straight off $\lambda^\star(B)$. **Anchor:**
    to hit $V_0=8$ on the $a=(2,1)$ curve ($S=5$) requires $B=64/5=12.8$, where
    $\lambda^\star=\sqrt5/(2\sqrt{12.8})=5/16=0.3125$ — exact. The frontier and its inverse are the
    same identity read in two directions.

**Rung 7 — The wound (the chapter's closing rung).**
Combine Rungs 4–6 under uncertainty. The *deterministic* scenario table looks crisp and decisive; the
*posterior* version shows the scenario rankings overlap and the shadow price carries a wide band —
**the what-if table is more confident than the data warrants.** Two posterior-plausible curves give
materially different optima: *good fit ≠ good decision*. State plainly: observational MMM weakly
identifies the response curve, so the budgeting decision inherits that weak identification; the only
way to shrink the EVPI gap is to *intervene* and measure. Hand off to **Part VI** by name (causal
foundations → quasi-experimental designs → calibration), which heals the wound by supplying the
experimental identification this chapter shows is missing.

---

## 4. Worked Examples

**WE1 — Closing the Part IV loop (deterministic).**
Solve the *known-curve* budget problem and recover Part IV's answer exactly: $r_i(b)=a_i\sqrt b$,
$a=(2,1)$, $B=9$ ⇒ $b^\star=(7.2,1.8)$, $\lambda=0.3727$, $V^\star=6.7082$. Then **verify the envelope
theorem numerically**: $V^\star(B)=\sqrt{5B}$, and a finite-difference $dV^\star/dB$ at $B=9$ matches
$\lambda^\star=0.3727$. Payoff: the MMM allocation is *the same number* Part IV converged to — Part V
has now supplied the curve Part IV assumed. (By-hand-checkable; all anchors analytic.)

**WE2 — Scenario analysis via the shadow price (deterministic).**
On the same curve, run a scenario table: budget $+20\%$ ($B=10.8$), $-20\%$ ($B=7.2$), a channel-1
cap $b_1\le 5$, and channel-2-dark ($b_2=0$). For each, give the **first-order** envelope prediction
$\Delta V\approx\lambda^\star\Delta B$ and the **re-optimized actual**, and show the curvature gap.
Anchors (analytic, to verify in plan): $V^\star(10.8)=\sqrt{54}=7.3485$ (first-order predicts
$6.7082+0.3727\cdot1.8=7.3790$ — overshoots by concavity); $V^\star(7.2)=\sqrt{36}=6.0000$;
cap $b_1\le5\Rightarrow b=(5,4)$, $V=2\sqrt5+2=6.4721$; channel-2-dark $\Rightarrow b=(9,0)$,
$V=2\sqrt9=6.0000$. The shadow price predicts the table's first column before any re-solve.
Then two scenario tools from the value function (Rung 6): (i) the **frontier** $V^\star(B)=\sqrt{5B}$
with marginal $\lambda^\star(B)=\sqrt5/(2\sqrt B)$ swept over $B$; and (ii) the **inverse** — to hit
target $V_0=8$, solve $B=V_0^2/5=12.8$, where $\lambda^\star=5/16=0.3125$ (exact). Verify the inverse
by re-optimizing at $B=12.8$ and confirming $V^\star=8$.

**WE3 — The wound (under posterior uncertainty).**
Now the curve is *uncertain*. Place a posterior on the response coefficients reflecting the Ch. 16
**ridge** (two coefficients only weakly separated by the data). Draw from it and, per draw, re-optimize:
  - **Allocation posterior:** $b^\star$ and $\lambda^\star$ get wide bands — the confident
    deterministic optimum is one draw among many.
  - **Optimizer's curse, numerically:** estimate $\mathbb E_\theta[\max_b V]$ and
    $\max_b \mathbb E_\theta[V]$ and show the gap $\ge 0$ (Rung 5); report it as EVPI — what an
    experiment is worth.
  - **Two plausible curves, two optima:** pick two posterior-plausible parameter draws with nearly
    equal in-sample fit but materially different $b^\star$; tabulate "good fit, different decision."
Exact anchors pinned in the plan via NumPy (seeded posterior draws).

---

## 5. Code Tie-in

A single runnable `{python}` cell (NumPy + `scipy.optimize` + Matplotlib; verified headless under
`MPLBACKEND=Agg`). It does, in order:

1. **Deterministic optimum (WE1):** SLSQP on $\sum a_i\sqrt{b_i}$, $a=(2,1)$, $\mathbf1^\top b=9$;
   assert $b^\star\approx(7.2,1.8)$, response $\approx6.7082$, recovered multiplier
   $\approx0.3727$ — closing the Part IV loop on a live solver.
2. **Envelope check:** finite-difference $dV^\star/dB$ across a budget sweep; assert it tracks
   $\lambda^\star(B)=\sqrt5/(2\sqrt B)$. Figure: $V^\star(B)$ with slope $=\lambda^\star$ annotated.
3. **Scenario table + frontier + inverse (WE2):** re-optimize under $\pm20\%$, channel cap,
   channel-dark; print first-order vs actual $\Delta V$; assert the analytic anchors above. Sweep $B$
   to plot the **frontier** $V^\star(B)$ with slope $\lambda^\star(B)$ annotated; solve the **inverse**
   $B=V_0^2/5$ for $V_0=8$ and assert $B=12.8$, $\lambda^\star=0.3125$, re-optimized $V^\star=8$.
4. **Uncertainty propagation + optimizer's curse (WE3):** seed a posterior over the coefficients
   (Ch. 16 ridge), re-optimize per draw; plot the $b^\star$ posterior; compute
   $\mathbb E[\max]-\max\mathbb E \ge 0$ and print it as EVPI; show two equal-fit / different-optimum
   draws. Figure: allocation posterior + the deterministic point estimate overlaid (the wound made
   visible).

Every numeric claim asserted against NumPy. Figures end with `plt.show()`.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C (Conceptual):** Why does optimizing on a point estimate overstate realized response? What does
  the EVPI gap mean operationally, and why does it shrink when you run an experiment? Why is "good
  in-sample fit" insufficient for a budgeting decision?
- **B (By hand):** On $r_i=a_i\sqrt b$, $a=(2,1)$: derive $V^\star(B)=\sqrt{5B}$, confirm
  $dV^\star/dB=\lambda^\star$, compute the budget-$+20\%$ first-order vs exact $\Delta V$, and solve the
  inverse scenario — the budget $B$ needed to reach $V_0=8$ — reading off its shadow price ($B=12.8$,
  $\lambda^\star=5/16$).
- **P (Prove it):** Prove the optimizer's-curse inequality $\mathbb E[\max_b V]\ge\max_b\mathbb E[V]$
  from the pointwise bound; and prove the envelope identity $dV^\star/dB=\lambda^\star$ for the
  equality-constrained problem.
- **A (Applied / code):** Extend the Code Tie-in — widen the posterior (a flatter ridge) and show the
  EVPI gap and allocation bands grow; add a minimum-spend floor constraint and read its shadow price.

---

## 7. Appendix solutions

Add a `## Chapter 18 — Budget Optimization & Scenario Analysis` block to `appendix/solutions.qmd`,
**in chapter order** (after the Ch. 17 DLM block), gated by
`::: {.content-visible when-meta="show-solutions"}`. Full solutions for C/B/P/A, with the P-block
carrying the two proofs in full and the B-block the $\sqrt{5B}$ arithmetic.

---

## 8. Summary (auto-included)

Bulleted **Key concepts** and bulleted **Key identities** (inline math only — matches sibling
chapters; Key identities must be a bulleted list, never a run-on paragraph). Key identities to
include: equal-marginal $r_i'(b_i^\star)=\lambda^\star\ \forall i$; envelope $dV^\star/dB=\lambda^\star$;
value function $V^\star(B)=\sqrt{SB}$; optimizer's curse $\mathbb E[\max]\ge\max\mathbb E$; EVPI as the
gap.

---

## 9. Cross-references (by name, post-renumber numbering)

- **Ch. 12 (Convexity):** equal-marginal / KKT optimum, the $b^\star=(7.2,1.8)$ anchor.
- **Ch. 13 (Linear Programming):** shadow price as a dual / LP multiplier.
- **Ch. 14 (SLSQP):** the solver and the multiplier-as-shadow-price reading.
- **Ch. 15 (MMM DGP):** the saturating, adstocked response curve being optimized.
- **Ch. 16 (Building & Fitting):** the posterior over the curve and the **Fisher-information ridge**
  (weak identification) that drives the wound.
- **Part VI (Ch. 19–22):** named as the heal — causal foundations, quasi-experimental designs,
  calibration — supplying the experimental identification the EVPI gap shows is missing.

---

## 10. Bibliography

Reuse existing keys where possible (`@boyd2004` for convex optimization, `@gelman2013` for Bayesian
decision). Add if needed for the optimizer's curse / decision under uncertainty: a
decision-theory/value-of-information reference (verify key exists in `references.bib`; add a standard
entry if not — e.g. Smith & Winkler on the optimizer's curse, Raiffa/Schlaifer on EVPI). The plan
will confirm and add exact BibTeX.

---

## 11. Future work (deferred, not in this chapter)

This chapter treats scenario analysis through the value function and the shadow price: the scenario
table, the response-vs-budget **frontier**, and the **inverse** ("budget to hit a target"). A fuller,
first-class **scenario-planning** treatment is deliberately deferred — to be added to the book later
as either an expanded section here or a dedicated chapter: minimum-spend **floors**, **channel-to-
channel reallocation** under a fixed total, multi-period / pacing scenarios, and the **uncertainty
overlay** applied across all of them (posterior bands on every scenario in the table). Logged so it is
not lost.

---

## 12. Verification (the real gate)

- Every numeric anchor NumPy-verified before it is written.
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`).
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no
  `\begin{align}`, no `psmallmatrix`, valid citation keys.
- **CI `quarto render` (HTML + PDF) green on the PR** — Quarto is not installed locally; this is the
  gate. User merges the PR.
