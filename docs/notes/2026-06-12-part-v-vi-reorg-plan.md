# Reorganization Plan — Part V (Modeling) & Part VI (Causal Grounding & Calibration)

**Date:** 2026-06-12
**Context:** Working session that designed an advanced-material expansion of the book (dynamic
linear models, causal inference, quasi-experimental designs, advanced calibration, and a
"prior store" technique) and reorganized the chapter structure to accommodate it. This note is
the **source of truth** for the agreed structure so any session writing chapters follows it.
All four open structural decisions were resolved with the author and are **locked**.

---

## Summary of the change

Expand the book from 21 → **26 chapters**. The new material reshapes the V/VI boundary around a
single clean line:

- **Part V = build the best *observational* model you can** (now includes DLM).
- **Part VI = ground that model in *causal* evidence and make it compound** (causal foundations,
  quasi-experimental designs, advanced calibration, the prior store).

**Narrative spine:** build → make dynamic → optimize → hit the wall → understand why (causality)
→ get exogenous evidence (designs) → fold it in (calibration) → make it compound + close the loop
(prior store) → productionize → integrate.

The key structural insight: Part V's budget-optimization chapter ends on a *wound* — you are
optimizing spend over a response curve the observational time series only weakly identifies. Part
VI exists to heal that wound.

---

## Finalized table of contents (26 chapters)

Changed/new items are marked.

**Part I — Mathematical Foundations**
1. Linear Algebra
2. Multivariable Calculus
3. Probability & Statistics

**Part II — Regression & Bayesian Inference**
4. Linear Regression — **NEW closing section**: plant endogeneity ("a fitted coefficient is
   causal only under conditions we'll formalize later"). A teaser, not the full treatment.
5. Bayesian Inference
6. Hierarchical Regression

**Part III — Computation: Sampling**
7. Markov Chains
8. MCMC
9. HMC & NUTS
10. Model Checking & Diagnostics

**Part IV — Optimization**
11. Convexity
12. Linear Programming
13. Constrained Optimization (SLSQP)

**Part V — Marketing Mix Modeling: Modeling**
14. The MMM Data-Generating Process
15. Building & Fitting an MMM
16. **Dynamic Linear Models & State-Space MMM** — *MOVED here from the originally-proposed
    Part VI.* It is a modeling extension and belongs with its siblings.
17. Budget Optimization — *ends Part V on the wound (optimizing over a weakly-identified curve).*

**Part VI — Causal Grounding & Calibration**
18. Causal Inference Foundations — *opens VI. Decision: kept here, NOT moved into Part II.
    Answers "why can't the data pin the curve?" → spend is endogenous.*
19. Quasi-Experimental Design
20. Advanced Calibration
21. **The Prior Store & the Calibration–Optimization Loop** — *NEW synthesis chapter, the
    conceptual climax.*

**Part VII — Software Engineering & Computer Science**
22. CS Foundations
23. Software Architecture
24. Data Engineering — *the prior store is implemented here as a versioned data product.*
25. Testing & Reliability

**Part VIII — Capstone**
26. Putting It All Together — *spine: the MMM as a compounding causal-evidence system.*

---

## Resolved decisions (locked)

1. **Causal Inference Foundations stays at the head of Part VI** (Ch. 18), not relocated to the
   inference arc. Rationale: causality is most teachable once the reader has a concrete MMM with
   visibly endogenous spend, and keeping it adjacent to the designs that use it (Ch. 19) avoids a
   12-chapter gap between learning potential outcomes and applying them. Mitigation: a short
   endogeneity **teaser in Ch. 4**.
2. **DLM lives in Part V** (Ch. 16) as a modeling chapter, upstream of the two later chapters that
   reuse its machinery (BSTS in Ch. 19, the prior store in Ch. 21). **Required internal structure:
   classical DLM first, then Bayesian DLM** (see below).
3. **The prior store gets its own chapter** (Ch. 21), paired with the calibration–optimization
   loop. Carries real math (power-prior weighting, DLM-as-prior, hierarchical pooling). The
   capstone (Ch. 26) runs it end-to-end rather than re-deriving it.
4. **Part VI → VII handoff via the prior store as a data product.** VI defines it conceptually,
   VII (Data Engineering, Ch. 24) implements it with schema/versioning/lineage, VIII integrates.

---

## Per-chapter content notes

### Ch. 16 — Dynamic Linear Models & State-Space MMM

Present **classical DLM, then Bayesian DLM**:

- **Classical.** State-space form: observation `y_t = F_t θ_t + v_t`, `v_t ~ N(0, V_t)`; state
  `θ_t = G_t θ_{t-1} + w_t`, `w_t ~ N(0, W_t)`. Kalman filter (forward predict→update) gives
  closed-form Gaussian `p(θ_t | y_{1:t})`; RTS smoother gives `p(θ_t | y_{1:T})`. Variances
  `(V, W)` via MLE (prediction-error decomposition) or discount factors (West & Harrison).
- **Bayesian.** Priors on `(V, W)` and `θ_0`; full posterior over states *and* hyperparameters via
  **FFBS** (Forward-Filtering Backward-Sampling — present it as "the Kalman filter you just
  learned, run forward then sampled backward," not a new algorithm) + a Gibbs/MCMC step for the
  variances.

Three payoffs to surface explicitly:
- (a) It **rhymes with Ch. 4 → Ch. 5** (frequentist → Bayesian) at the time-series level.
- (b) It is the **cleanest conjugacy-vs-MCMC boundary lesson**: with Gaussian errors and known
  variances the Kalman recursions are exact (no sampling); make the variances unknown and you fall
  back to MCMC. Retroactively sharpens Part III.
- (c) It **splits the downstream dependencies cleanly**: classical Kalman → BSTS (Ch. 19);
  Bayesian DLM → prior store (Ch. 21).

### Ch. 19 — Quasi-Experimental Design

Organize the chapter around a **secant / tangent / validation taxonomy** — classify each design by
where its causal estimand lands on the response curve `S`:

| Calibration target | Constrains | Designs |
|---|---|---|
| **Secant** (difference / average slope) | level change over an interval | DiD, synthetic control, ITS, geo lift, **IV** |
| **Tangent** (point derivative / curvature) | slope (or kink) at one `x` | RDD, RKD |
| **Validation** | nothing in `θ` | any estimand that isn't a functional of `S` |

Full worked treatment of **DiD, synthetic control, ITS, IV, and BSTS/CausalImpact**. IV is the
endogeneity-of-spend lens (pairs with Ch. 18); BSTS builds its counterfactual in the classical
state-space machinery from Ch. 16. **No matching/propensity methods** (decided out of scope).
RDD/RKD, switchback, and staggered-adoption DiD are brief mentions only.

Key subtlety to make explicit: **IV is a secant in disguise.** The Wald estimator is an average
derivative over the complier spend range, which by the fundamental theorem of calculus equals a
secant slope — so it calibrates as a secant, not a true point-tangent. The genuine point-local
designs are RDD (level at the cutoff) and RKD (kink = curvature at the cutoff).

### Ch. 20 — Advanced Calibration

**Hard constraint: explain lift-test calibration as a *generic technique*. Do NOT name
PyMC-Marketing (or any library) in the book text.** It is a method, not a product.

Mechanism — add a likelihood factor tying the model-implied lift to the measured lift, keyed by
the taxonomy:
- **Secant:** `Δ̂ ~ Normal( S(x+δ; θ) − S(x; θ), s² )` (or a positivity-respecting Gamma).
- **Tangent:** match `S'(x; θ) · δ`.
- **Validation:** post-hoc reconciliation / posterior-predictive check only — no likelihood term.

Frame it as **Bayesian evidence synthesis**, not "experiment-as-prior-on-ROAS" (the
likelihood-factor view is cleaner and avoids double-counting).

Backbone metaphor (carry it through the chapter): the discrete lift is the **secant** of the
response curve; the marginal form is the **tangent** (`δ → 0`). Develop via the mean value
theorem; the sign of the secant–tangent gap under concavity (Taylor: `δ S'(x) + (δ²/2) S''(ξ)`,
`S'' < 0`); paired up/down tests bracketing `S'(x)`; and test-design ↔ identification (one tangent
under-identifies the shape — you need a spread of constraints, or a wide secant, to pin the
saturation parameters `κ`, `α`).

### Ch. 21 — The Prior Store & the Calibration–Optimization Loop

**Named contribution; conceptual climax.** A curated, versioned, *accumulating* repository of
priors on each channel's response curve, fed by the **business's existing quasi-experimental
backlog** (geo holdouts, ITS readouts, IV-style natural experiments, operational RDDs) that is
otherwise discarded after a one-off readout. Components:

1. **Schema = the secant/tangent/validation taxonomy.** Each study normalized to a functional of
   `θ` + operating point `x` + SE `s` + validity tier + provenance.
2. **Sequential update.** Posterior at period *t* becomes the prior at *t+1*.
3. **Power-prior down-weighting.** `p(θ | data, store) ∝ p₀(θ) · L_obs(θ) · ∏_k N(ŷ_k | g_k(θ),
   s_k²)^{w_k}` with `w_k ∈ [0,1]` = design credibility × temporal decay `ρ^{(t_now − t_k)}`.
4. **DLM tie-in.** Place store entries as timestamped observations feeding a Bayesian-DLM smoother
   over `θ_t` (reuses Ch. 16's Bayesian DLM). Static power-prior calibration is the random-walk-
   variance → 0 limit.
5. **Hierarchical pooling** across conflicting studies (reuses Ch. 6 random-effects meta-analysis).
6. **Validation-tier entries** (`w = 0`) act as posterior-predictive checks / conflict alarms;
   they never enter the likelihood.

**The calibration–optimization loop:** the optimizer (Part IV) flags high-leverage spend regions
→ aim the next experiment there (test-design ↔ identification) → estimate (Ch. 19) → update the
store (Ch. 20) → re-optimize. A closed active-learning loop.

Honest framing for the text: theoretically this is operationalized cumulative meta-analysis /
power priors / Bayesian evidence synthesis. The contribution is the **systematization** (schema,
weighting policy, decay model, governance), not novel theory — say so plainly.

---

## Dependency graph (new/changed edges)

```
Ch. 4 (endogeneity teaser) ─┐
                            ├─→ Ch. 18 Causal Foundations ─→ Ch. 19 QE Designs ─→ Ch. 20 Calibration ─→ Ch. 21 Prior Store
Ch. 6 Hierarchical ─────────┘                                      │                                          ↑
                                                                   │ (BSTS reuses classical Kalman)           │ (hierarchical pooling)
Ch. 16 classical Kalman ───────────────────────────────────────────┘                                          │
Ch. 16 Bayesian DLM (FFBS) ──────────────────────────────────────────────────────────────────────────────────┘
Ch. 13/IV Optimization ─→ Ch. 17 Budget Opt (the wound) ─→ Ch. 21 (the loop closes back to optimization)
Ch. 21 Prior Store ─→ Ch. 24 Data Engineering (implemented as a data product) ─→ Ch. 26 Capstone
```

---

## Sequencing for writing

Parts I–III are substantially drafted; IV–VIII (including all of the new material) are currently
stubs. Recommended order:

1. **Ch. 20 (Advanced Calibration) spec first** — it is the conceptual keystone and depends only
   on already-drafted chapters (5, 6, 15). Getting it right tells you how much causal scaffolding
   Ch. 18 must carry.
2. Then Ch. 16 (DLM), since it is the upstream dependency for both Ch. 19 (BSTS) and Ch. 21.
3. Then Ch. 18 → 19 → 21 in narrative order.

Per the project's chapter-authoring workflow: each chapter follows the fixed template (learning
objectives → theory → worked example → runnable code tie-in → C/B/P/A tiered exercises → summary),
with appendix solutions, verified by `quarto render` in CI.

---

## `_quarto.yml` impact (not yet applied)

The book's `_quarto.yml` `chapters:` list will need: Part V gains a `16-dlm.qmd`; a new Part VI
(`causal-foundations`, `quasi-experimental`, `advanced-calibration`, `prior-store`) inserted before
the current Software Engineering part; current Parts VI/VII renumber to VII/VIII. Do this as part
of scaffolding the new chapter files, not in this planning note.
