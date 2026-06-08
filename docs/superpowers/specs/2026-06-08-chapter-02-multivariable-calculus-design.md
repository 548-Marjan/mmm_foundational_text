# Chapter 2 — Multivariable Calculus & the Optimization Toolkit: Design Spec

**Date:** 2026-06-08
**Branch:** `chapter-2-brainstorm`
**File to author:** `parts/01-foundations/02-multivariable-calculus.qmd`
**Canonical anchors:** Nocedal & Wright (1999/2006), Boyd & Vandenberghe (appendix A).

## Goal

Replace the scaffolded stub with a complete chapter that reactivates multivariable
calculus and turns it into the **optimization toolkit** the rest of the book runs
on. Scope is the **calculus *of* optimization** — gradient, Jacobian, Hessian,
multivariate Taylor, and the optimality/convexity conditions that rest on them —
**not** optimization algorithms (those are Part IV). Deliverable this cycle is a
**full chapter draft**: prose, proofs, worked examples, a runnable Python tie-in,
and full C/B/P/A exercises.

The chapter follows the book's fixed template:
Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises (C/B/P/A).

## Through-line (anchor example)

Every abstraction lands on one concrete object: the **SSE loss surface**
`L(β) = ‖y − Xβ‖²`, where `X` is the same media-mix design matrix from Chapter 1.
Chapter 1 closed by asserting that `XᵀX` is the curvature (Hessian) of the
regression loss; Chapter 2 *cashes that claim out with calculus*. Each tool earns
its place on `L`:

- **Gradient** — `∇L(β) = −2Xᵀ(y − Xβ)`; setting `∇L = 0` **re-derives the normal
  equations** `XᵀXβ̂ = Xᵀy`, this time by calculus rather than by projection. The
  book now has two independent derivations of its keystone, one per chapter.
- **Jacobian** — earns its place via the **residual map** `r(β) = y − Xβ`, whose
  Jacobian is `−X`; then `∇L = 2Jᵀr`. This is the **Gauss–Newton** shape, a
  named-and-deferred preview of how the *nonlinear* saturating MMM curves (Part V)
  are fit.
- **Hessian** — `∇²L(β) = 2XᵀX`, literally Chapter 1's symmetric PSD matrix, now
  revealed as curvature.
- **Convexity** — because `2XᵀX ⪰ 0`, `L` is convex, so its critical point is the
  *global* minimum. This is the guarantee the rest of the book leans on whenever it
  fits a least-squares or Gaussian-likelihood model.

The spine is: *how do you descend a loss surface, and when can you trust that the
bottom you reach is the only bottom?* Slope = gradient, "only bottom" = convexity.

## Section-by-section content

### Motivation
Open on the loss surface: *"Chapter 1 found `β̂` by projecting `y` onto the column
space. But how would you find it if you could only feel the slope under your feet —
and how do you know the bottom you reach is the only bottom?"* Slope is the
gradient; "only bottom" is convexity. That frames the whole toolkit as answering
*how to descend* and *when descent is trustworthy*.

### Theory & Proofs
A ladder, each rung a downstream payoff. Four results are **proved at full rigor**
(marked ★); mechanical results are stated and cited.

1. Partial derivatives, the **gradient**, the directional derivative
   `D_u f = ∇f · u`. **★ Proof: ∇f is the direction of steepest ascent** and is
   orthogonal to level sets (via Cauchy–Schwarz on `D_u f`).
2. The **Jacobian** of a vector-valued map; the **multivariate chain rule** (stated,
   with the residual-map application `r(β) = y − Xβ`, `J = −X`, worked through to
   `∇L = 2Jᵀr`).
3. The **Hessian**; **multivariate Taylor to second order** (stated; the integral /
   Lagrange remainder form is *cited* to Nocedal & Wright, not proved).
4. **★ First-order optimality:** an interior local minimizer has `∇f = 0` (proved).
5. **★ Second-order sufficient condition:** `∇f = 0` together with a positive
   definite Hessian ⇒ strict local minimum (proved, via the second-order Taylor
   expansion and Chapter 1's PD machinery).
6. **★ Convexity:** for `C²` `f`, convex ⇔ Hessian PSD everywhere (proved); applied
   to show `L(β)` is convex because `∇²L = 2XᵀX ⪰ 0`, so its critical point is
   global.
7. **Named-and-deferred to Part IV:** gradient descent, Lagrange multipliers, and
   the KKT conditions are defined *by name only*, each with a forward pointer. No
   algorithms in this chapter.

### Worked Examples
- (a) For a small 3-week, 2-channel `X`, compute `∇L(β)` and `∇²L` by hand, set
  `∇L = 0`, and recover the *same* `β̂` Chapter 1 solved by projection — the two
  derivations meeting on one number.
- (b) A 2-D scalar `f` with a **saddle point**: show `∇f = 0` holds but the Hessian
  is indefinite, illustrating why the first-order condition alone is insufficient
  and motivating the second-order test.

### Code Tie-in
NumPy / SciPy, runnable and self-contained:
- Build the synthetic `X`; form `∇L(β)` and `∇²L` in closed form.
- **Verify the gradient numerically** (finite differences vs. the analytic
  gradient) — a discipline that pays off for every model fit later in the book.
- Evaluate the Hessian's eigenvalues to confirm PSD / convexity, and tie the
  smallest eigenvalue back to Chapter 1's collinearity story: a near-zero
  eigenvalue is a flat valley in `L`, which is exactly what makes `β̂` unstable and
  descent slow.

### Exercises

Proofs taught in the body are **not** re-assigned as exercises. The P tier instead
poses *adjacent* proofs that extend or complement the four taught results, so the
reader does genuinely new reasoning on the same machinery.

- **C — Conceptual:** what the gradient and Hessian mean geometrically on a loss
  surface; why convexity is what lets you trust a single fitted `β̂`; what a saddle
  point is and why `∇f = 0` does not by itself certify a minimum.
- **B — By hand:** compute gradients and Hessians of small quadratic and
  non-quadratic functions; classify critical points (min / max / saddle) via the
  Hessian; hand-compute `∇L` for a tiny `X` and solve `∇L = 0`.
- **P — Prove it** (each extends, rather than repeats, a body proof):
  - **P1 — Second-order *necessary* condition** (companion to the taught
    *sufficient* one): prove that at an interior local minimum the Hessian is PSD
    (`∇²f ⪰ 0`). The gap between PSD-necessary and PD-sufficient is exactly where
    saddles and flat valleys live.
  - **P2 — First-order characterization of convexity:** prove that for `C¹` `f`,
    convexity ⇔ `f(z) ≥ f(β) + ∇f(β)ᵀ(z − β)` (the graph lies above every tangent
    plane) — a *different* characterization than the second-order Hessian test
    taught in the body, and the inequality the rest of the book actually uses.
  - **P3 — Convex ⇒ local min is global; strict convex ⇒ unique minimizer:** extend
    the taught convexity result to its payoff, then connect it to the through-line —
    `β̂` is the *unique* global minimizer of `L` exactly when `XᵀX` is PD (full
    column rank), tying back to Chapter 1's collinearity story.
  - **P4 — Convexity preserved under affine precomposition:** prove that if `g` is
    convex then `β ↦ g(Xβ + c)` is convex, and use it to show `L(β) = ‖y − Xβ‖²` is
    convex *without differentiating* — an alternative route to the body's
    Hessian-based argument.
- **A — Applied / code:** implement a finite-difference gradient check against the
  analytic `∇L`; hand-step a descent down the SSE bowl and watch convergence slow
  as the condition number of `XᵀX` worsens.

Exercise solutions go in the shared `appendix/solutions.qmd`, gated by the
`show-solutions` metadata flag (no inline links, per the preface).

## Rigor level

Chapter 2 continues Chapter 1's ramp: intuition and geometry first (Nocedal &
Wright's applied flavor), with four genuinely illuminating results proved at full
rigor (steepest ascent, first-order optimality, second-order sufficiency, convex ⇔
PSD Hessian) and mechanical results (Taylor's remainder) cited. It sits slightly
above Chapter 1 in rigor and foreshadows the heavier optimization treatment in
Part IV, matching the preface's "prove where illuminating, cite where mechanical."

## Conventions & constraints

- Quarto `.qmd`, KaTeX math (HTML) / LaTeX (PDF); number-sections on. Keep math to
  constructs KaTeX renders.
- Citations via `references.bib`: add/confirm keys for Nocedal & Wright and Boyd &
  Vandenberghe; reuse Chapter 1 keys (`@strang2016`) for the linear-algebra
  call-backs as needed.
- Code is original, minimal, runnable against pinned `requirements.txt`; no
  proprietary data or priors.
- Match the established voice of `index.qmd` and Chapter 1 (direct, no hand-waving),
  and explicitly thread the Chapter 1 → Chapter 2 continuity (`XᵀX` as Hessian).

## Out of scope (YAGNI)

Optimization *algorithms* (gradient descent, Newton, line search — Part IV);
constrained optimization, Lagrange multipliers, and KKT (Part IV); integration and
change-of-variables for densities (deferred to the probability/statistics chapter,
where it is actually needed); vector calculus (divergence, curl); and the
differential geometry of manifolds. The MMM response-curve through-line (saturating
transforms) and generic abstract scalar-field treatment are deliberately deferred to
later chapters where they earn their place.

## Success criteria

- Chapter renders cleanly (`quarto render`) in HTML and PDF with no math/build
  errors.
- Every theory subsection ties back to the SSE-loss-surface through-line, and the
  Chapter 1 → Chapter 2 continuity (`∇²L = 2XᵀX`) is explicit.
- The four named proofs are present and correct in the body, and are *not* duplicated
  as exercises.
- The four P exercises pose distinct, adjacent proofs (necessary condition,
  first-order convexity characterization, global/unique minimizer, affine-composition
  convexity) with matching solutions in the appendix.
- The code tie-in runs top-to-bottom and demonstrates the analytic-vs-numeric
  gradient check and the eigenvalue/conditioning tie-back to Chapter 1.
- All four exercise tiers are populated, with matching solutions in the appendix.
