# Chapter 20 — Quasi-Experimental Design: Design Spec

**Date:** 2026-06-17
**Part:** VI — Causal Grounding & Calibration (second chapter)
**File:** `parts/06-causal-calibration/02-quasi-experimental.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

Chapter 19 proved that observational MMM cannot identify the causal response curve when spend is
endogenous and the confounder is unobserved. This chapter supplies the escape: **quasi-experimental
designs** that generate (or exploit) spend variation independent of demand, producing *interventional*
identification the observational series cannot.

The organizing idea — and the chapter's spine — is the **secant / tangent / validation taxonomy**.
Every quasi-experimental design produces an estimand that is a **functional of the response curve $S$**,
and the taxonomy classifies *which* functional:

| Calibration target | Constrains | Designs |
|---|---|---|
| **Secant** (difference / average slope) | a level change over a spend interval | DiD, synthetic control, ITS, geo lift, **IV** |
| **Tangent** (point derivative / curvature) | the slope $S'(x)$ or curvature $S''(x)$ at one $x$ | RDD, RKD |
| **Validation** | nothing in the curve parameters | any estimand that isn't a functional of $S$ |

This taxonomy is not just pedagogy: it is the **schema Chapter 21 ingests**. Calibration adds a
likelihood factor keyed by the taxonomy (secant → match $S(x{+}\delta)-S(x)$; tangent → match
$S'(x)\,\delta$; validation → posterior-predictive check only). So this chapter is the bridge from
"run an experiment" to "fold it into the model."

**Through-line driving question:** *"An experiment gives you one number — but a number about what part
of the curve?"* The taxonomy is the answer.

**Anchor curve (reused from Chapter 18):** $S(x) = 2\sqrt{x}$, a concave saturating response. Using a
*nonlinear* curve is essential — it is the only setting in which secant and tangent differ, which is
the whole content of the taxonomy.

---

## 2. Scope discipline

- **Full worked/code treatment:** DiD, synthetic control, ITS, IV, and BSTS (BSTS = the state-space
  way to build an ITS counterfactual, **explicitly reusing Chapter 17's Kalman filter**).
- **Two formal proofs only** (the keystones): DiD identification under parallel trends, and
  IV-is-a-secant. Synthetic control and ITS/BSTS get constructive worked treatment, not separate
  theorems.
- **Brief mention only:** RDD/RKD (the genuine tangent/curvature designs), switchback, and
  staggered-adoption DiD.
- **Out of scope (decided):** matching / propensity-score methods.
- House rules: **never** name PyMC-Marketing or any MMM/PPL/sampler library (BSTS/CausalImpact are
  described as *methods*, never as a product); `numpy`/`scipy`/`matplotlib` only. **Never**
  `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$ … $` only; `$$` on their own lines; even `$$`
  count.

---

## 3. Theory & Proofs — the rungs (6)

**Rung 1 — The secant / tangent / validation taxonomy.**
Define the response curve $S(x)$ and its functionals. A **secant** is the average slope
$[S(x_1)-S(x_0)]/(x_1-x_0)$ over an interval (equivalently, a level change $S(x_1)-S(x_0)$); a
**tangent** is a point derivative $S'(x)$ (or curvature $S''(x)$); a **validation** target is any
quantity that is not a functional of $S$ (e.g. a brand-awareness lift that does not constrain the
sales-response parameters). State the chapter's thesis: every design below identifies one of these,
and *which* one is what the experiment can legitimately tell you about the curve.

**Rung 2 — Proof P1 (KEYSTONE): difference-in-differences identifies the ATT under parallel trends.**
Potential-outcomes panel: regions $r\in\{T,C\}$, periods $\{\text{pre},\text{post}\}$; the treated
region $T$ receives a spend change in the post period. **Parallel-trends assumption:** the
counterfactual (untreated) trend of the treated region equals the observed trend of the control,
$\mathbb E[Y_{T,\text{post}}(0)-Y_{T,\text{pre}}(0)] = \mathbb E[Y_{C,\text{post}}(0)-Y_{C,\text{pre}}(0)]$.
**Theorem:** the DiD estimand
$\big(\mathbb E[Y_{T,\text{post}}]-\mathbb E[Y_{T,\text{pre}}]\big) - \big(\mathbb E[Y_{C,\text{post}}]-\mathbb E[Y_{C,\text{pre}}]\big)$
equals the ATT $\mathbb E[Y_{T,\text{post}}(1)-Y_{T,\text{post}}(0)]$. **Proof:** substitute observed =
potential outcomes, add and subtract $Y_{T,\text{post}}(0)$, and use parallel trends to cancel the
common trend; the region fixed effect (baseline demand — the Chapter 19 confounder) differences out.
**MMM reading:** in a geo experiment the ATT is the **lift** $S(x_1)-S(x_0)$ — a *secant* of the
response curve — and DiD recovers it while removing the baseline-demand confound that biases a naive
cross-section. **Anchor:** $S(x)=2\sqrt x$, treated spend $x_0=4\to x_1=9$, baseline demand
$\alpha_T=10,\alpha_C=8$: naive post comparison $=4$ (lift $2$ + baseline gap $2$), DiD $=2 = S(9)-S(4)$.

**Rung 3 — Proof P2 (KEYSTONE): IV / the Wald estimator is a secant in disguise.**
Instrument $W$ (e.g. randomized encouragement, an auction-price shock) shifts spend without affecting
sales except through spend. **Wald estimator:**
$\hat\beta_{\text{IV}} = \dfrac{\mathbb E[Y\mid W{=}1]-\mathbb E[Y\mid W{=}0]}{\mathbb E[X\mid W{=}1]-\mathbb E[X\mid W{=}0]}$.
Under instrument validity (relevance, exclusion, independence, monotonicity) it identifies the LATE
over compliers. **Theorem:** for a smooth response $Y_i(x)=S(x)+\alpha_i$, the complier's outcome
change is $S(x_1)-S(x_0)=\int_{x_0}^{x_1}S'(u)\,du$, so
$$
\hat\beta_{\text{IV}} = \frac{S(x_1)-S(x_0)}{x_1-x_0} = \frac{1}{x_1-x_0}\int_{x_0}^{x_1} S'(u)\,du,
$$
the **average derivative of $S$ over the complier spend interval** — a *secant slope* by the
fundamental theorem of calculus, equal to the tangent $S'(\xi)$ at some interior $\xi$ (mean value
theorem) but **not** the tangent at the operating point you care about. **Anchor:** $x_0=4,x_1=9$:
$\hat\beta_{\text{IV}} = 2/5 = 0.4 = S'(6.25)$, distinct from $S'(4)=0.5$ and $S'(9)=0.333$. Corrects
the common misconception that IV yields a point marginal-return; cashes Chapter 19's IV teaser.

**Rung 4 — Synthetic control: a constructed counterfactual (secant).**
When there is one treated unit and many controls, build a weighted combination of control regions —
convex weights $w\ge 0$, $\sum w=1$, chosen to match the treated region's **pre-period** outcomes — and
use it as the counterfactual; the post-period gap is the lift. State the identifying intuition (a
factor/latent-confounder model: matching the pre-period trajectory matches the latent demand factors)
and the Abadie convex-hull / pre-fit conditions, without a formal theorem. The estimand is again a
**secant** (level change). Anchor: controls' pre matched to treated pre $=14$; post gap $=2$.

**Rung 5 — Interrupted time series and BSTS: a state-space counterfactual (secant).**
ITS fits the treated unit's **pre-intervention** trajectory and extrapolates it as the counterfactual;
the post-intervention gap is the lift. The modern construction (BSTS / structural-time-series
counterfactual) builds that extrapolation with the **Chapter 17 local-level Kalman filter** — make the
reuse explicit: the filtered/forecast level *is* the counterfactual, its forecast variance gives the
uncertainty band. The estimand is a **secant**. (Describe BSTS/CausalImpact as a method; never as a
product.) Anchor: pre-period level $\approx 14$, forecast counterfactual $14$, observed post $16$,
gap $2$.

**Rung 6 — The genuine tangent designs (brief), and the handoff.**
RDD identifies the response **level** at a cutoff and RKD the **kink** (curvature $S''$) — the only
designs whose estimand is a true point-local **tangent**, obtained as $\delta\to 0$ rather than over an
interval. Note switchback and staggered-adoption DiD as variants, briefly. Close: each design hands
Chapter 21 a **calibrated functional of $S$** tagged by the taxonomy (secant / tangent / validation),
which is exactly the schema the calibration likelihood and the prior store will consume —
intervene → identify → **calibrate** → re-optimize.

---

## 4. Worked Examples

**WE1 — Geo-lift difference-in-differences (secant).**
The 2×2 panel with the anchor numbers: treated region spend $4\to9$, control held at $4$, baseline
demand $\alpha_T=10,\alpha_C=8$. Tabulate $Y$ (pre/post × T/C): $Y_{T}=14\to16$, $Y_C=12\to12$. Show the
**naive** post cross-section $16-12=4$ is biased (lift $2$ + baseline gap $2$); the **DiD**
$(16-14)-(12-12)=2$ recovers the lift $S(9)-S(4)=2$, differencing out the confound. Verify parallel
trends holds by construction, and identify the lift as a **secant** (slope $2/5=0.4$).

**WE2 — Instrumental variables / Wald (secant in disguise).**
An instrument shifts complier spend $4\to9$. Compute $\hat\beta_{\text{IV}} = \Delta Y/\Delta X = 2/5
= 0.4$. Show it equals the secant slope and, by MVT, the tangent at $\xi=6.25$ — but **not** $S'(4)=0.5$
or $S'(9)=0.333$. Contrast with Chapter 19's confounded naive slope (which a non-experimental
regression would return). The lesson: IV gives an interval-average slope; to act at a specific
operating point you must know where on the curve that average sits.

**WE3 — Synthetic control (constructed counterfactual).**
One treated region, several controls with different baselines. Choose convex weights matching the
treated region's pre-period level ($14$); the synthetic control stays $\approx14$ post (controls
unchanged) while the treated region rises to $16$, so the gap $=2=$ the lift. Note how the pre-period
fit is the identification lever (good pre-fit ⇒ credible counterfactual), and that the estimand is a
secant. (Exact weights pinned in the plan via constrained least squares.)

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + `scipy` + `matplotlib`; verified headless). On one
simulated geo panel built from $S(x)=2\sqrt x$ + region fixed effects + noise (seeded), it:
1. **DiD:** computes the 2×2 estimate; asserts $\approx 2$ and that the naive cross-section $\approx 4$
   is biased.
2. **IV / Wald:** instrument shifting spend $4\to9$; asserts $\hat\beta_{\text{IV}}\approx 0.4$ = secant,
   prints $S'(4),S'(9)$ to show it is neither.
3. **Synthetic control:** convex weights via constrained least squares (`scipy.optimize.nnls` or
   `minimize` on the simplex) matching the treated pre-period; asserts post gap $\approx 2$.
4. **ITS / BSTS counterfactual:** the Chapter 17 local-level Kalman filter on the treated pre-period,
   forecast extrapolated as the counterfactual; asserts post gap $\approx 2$ with a forecast band.
5. Figure: the four designs' recovered lifts side by side (all $\approx 2$), and a panel showing the
   secant of $S$ on $[4,9]$ against the tangents at $4,9$ — the taxonomy, drawn.
Every numeric claim asserted against the analytic anchors.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why is every quasi-experimental estimand a functional of $S$, and why does the taxonomy
  matter for using the result? Why is IV's number a secant, not a point marginal return? What does
  parallel trends assume, and when is it violated in a geo experiment?
- **B:** On $S(x)=2\sqrt x$: compute the lift and secant slope for $x_0=4,x_1=9$; find the MVT point
  $\xi$ where the tangent equals the secant; show the DiD arithmetic on the 2×2 panel and that the
  naive comparison is biased by the baseline gap.
- **P:** (i) Prove DiD identifies the ATT under parallel trends. (ii) Prove the Wald estimand equals
  the average derivative $\frac{1}{x_1-x_0}\int_{x_0}^{x_1}S'(u)\,du$ (the secant), stating the
  instrument-validity conditions used.
- **A:** Extend the Code Tie-in — break parallel trends (give the control a different trend) and show
  DiD bias; widen the instrument's complier interval and show the Wald secant drifts; add a second
  control region and show synthetic-control pre-fit improves the counterfactual.

---

## 7. Appendix solutions

Append `## Chapter 20 — Quasi-Experimental Design` to `appendix/solutions.qmd`, **in chapter order**
(after the Chapter 19 block), inside the existing `content-visible` gated div. Full C/B/P/A; the
P-block carries both proofs (DiD parallel-trends; Wald = average-derivative secant via FTC).

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities: secant
$[S(x_1)-S(x_0)]/(x_1-x_0)$; DiD $=(\bar Y_{T,\text{post}}-\bar Y_{T,\text{pre}})-(\bar Y_{C,\text{post}}-\bar Y_{C,\text{pre}})$;
parallel trends; Wald $=\Delta Y/\Delta X = \frac{1}{x_1-x_0}\int S'$; tangent $S'(x)$ vs curvature
$S''(x)$. Close tying forward to Chapter 21 (calibration ingests the tagged functional) and back to
Chapter 18 (the curve being identified is the one the optimizer needs).

---

## 9. Cross-references (post-renumber numbering)

- **Chapter 17 (Dynamic Linear Models):** the local-level Kalman filter reused for the ITS/BSTS
  counterfactual.
- **Chapter 18 (Budget Optimization):** the response curve $S$ these designs identify is the one the
  optimizer consumes.
- **Chapter 19 (Causal Inference Foundations):** the confounding these designs defeat; the IV teaser
  cashed here.
- **Chapter 21 (Advanced Calibration):** ingests each design's estimand as a taxonomy-tagged
  likelihood factor.

---

## 10. Bibliography

Reuse `@angrist2009` (already in `references.bib` from Chapter 19) for DiD/IV. Add:
- `@abadie2010` — Abadie, Diamond & Hainmueller, synthetic control (JASA 2010).
- `@brodersen2015` — Brodersen et al., Bayesian structural time-series / CausalImpact (Annals of
  Applied Statistics 2015).
The plan confirms exact BibTeX and that all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor NumPy-verified (done: lift $2$, secant $0.4$, $\xi=6.25$, naive $4$, DiD $2$,
  Wald $\approx0.4$, SC gap $2$).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`).
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no
  `\begin{align}`, no `psmallmatrix`, valid citation keys, no banned library names.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
