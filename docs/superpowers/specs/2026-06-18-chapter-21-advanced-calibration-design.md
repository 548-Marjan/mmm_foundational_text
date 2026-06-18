# Chapter 21 — Advanced Calibration: Design Spec

**Date:** 2026-06-18
**Part:** VI — Causal Grounding & Calibration (third chapter; the conceptual keystone)
**File:** `parts/06-causal-calibration/03-advanced-calibration.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This is the chapter where Part VI pays off. Chapter 19 proved observational MMM cannot identify the
response curve under unobserved confounding (the structural wound); Chapter 20 produced
*interventional* estimands — taxonomy-tagged functionals of the curve (secant / tangent / validation).
This chapter **folds those estimands into the model**: it adds a **calibration likelihood factor**
that ties the model-implied functional to the measured one, compressing the posterior in exactly the
direction the observational data left flat — **healing the Chapter 16 ridge and the Chapter 18 wound**.
Then it hands the accumulating, weighted collection of these constraints to Chapter 22 (the prior
store).

The chapter closes the book's central loop numerically: **intervene → identify → calibrate →
re-optimize.** The EVPI gap Chapter 18 could only *name* is here made to *shrink*.

**Through-line driving question:** *"You ran the experiment and got a number. How do you put it into the
model — and how much does it help?"*

**Framing (locked):** Bayesian **evidence synthesis** — the experiment enters as a *likelihood factor*
on a functional of $\theta$, not as a hand-tuned "prior on ROAS." The likelihood-factor view composes
cleanly with the observational likelihood and avoids double-counting.

**Center of gravity (locked):** posterior compression ("calibration heals the ridge") is the central
keystone; the secant–tangent matching bias is the second keystone; power-prior down-weighting is kept
**brief** and deferred to Chapter 22.

**Anchor curve:** $S(x;\theta)=\theta\sqrt{x}$ (the Chapter 18/20 surrogate with an explicit scale
parameter $\theta$; truth $\theta=2$), plus a two-parameter coefficient model $(a_1,a_2)$ for the
ridge (reused from Chapter 18's WE3).

---

## 2. Scope discipline

- Does **not** re-derive the secant/tangent/validation taxonomy (Chapter 20 owns it); it *uses* it as
  the key for the calibration factor.
- Does **not** build the prior-store data product (schema, versioning, governance) — Chapter 22.
  Power-prior credibility/decay weighting is introduced briefly (Rung 5) and handed off.
- **House rule (critical, from the reorg note):** explain lift-test calibration as a **generic
  technique**. **Never** name PyMC-Marketing or any MMM/PPL/sampler library. It is a method, not a
  product. `numpy`/`scipy`/`matplotlib` only.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$ … $` only; `$$` on their own lines; even
  `$$` count.

---

## 3. Theory & Proofs — the rungs (5)

**Rung 1 — The calibration likelihood factor (taxonomy-keyed).**
Bayesian evidence synthesis: with observational data $D$ and an experimental estimate $\hat g$ of a
functional $g(\theta)$ of the response curve,
$$
p(\theta \mid D, \hat g) \;\propto\; p(\theta)\,\underbrace{L_{\text{obs}}(\theta\mid D)}_{\text{the MMM fit}}\,\underbrace{L_{\text{exp}}(\hat g \mid \theta)}_{\text{calibration factor}} .
$$
The calibration factor is keyed by Chapter 20's taxonomy:
- **Secant:** $\hat\Delta \sim \mathcal N\!\big(S(x{+}\delta;\theta) - S(x;\theta),\, s^2\big)$ (a
  positivity-respecting Gamma is the practical variant).
- **Tangent:** match the marginal, $\widehat{S'}\cdot\delta \sim \mathcal N\!\big(S'(x;\theta)\,\delta,\, s^2\big)$.
- **Validation:** *no* likelihood term — a posterior-predictive check / conflict alarm only.
State plainly: the experiment is *information about a functional of $\theta$*, entering through the
likelihood; this is cleaner than re-expressing it as a prior and never double-counts the same data.

**Rung 2 — Proof P1 (KEYSTONE): a calibration factor compresses the posterior in the constrained
direction (the ridge heals).**
Linear-Gaussian model: prior/observational posterior $\theta \sim \mathcal N(m_{\text{obs}},
\Lambda_{\text{obs}}^{-1})$; a calibration measurement of a linear functional $c^\top\theta$ with
$\hat g \sim \mathcal N(c^\top\theta, s^2)$. **Theorem:** the updated posterior precision adds,
$$
\Lambda_{\text{post}} = \Lambda_{\text{obs}} + \tfrac{1}{s^2}\,c\,c^\top,
$$
so the posterior variance **along the measured direction $c$ strictly decreases**, and most where
$\Lambda_{\text{obs}}$ is smallest — the **ridge**. **Proof:** complete the square in the Gaussian
log-posterior; the calibration factor contributes $\tfrac{1}{s^2}cc^\top$ to the precision (a rank-one
update), and by the matrix-determinant/Sherman–Morrison view the variance along $c$ drops from
$1/\Lambda_{\text{obs},c}$ to $1/(\Lambda_{\text{obs},c}+1/s^2)$. **MMM reading:** the Chapter 16
Fisher-information ridge is the near-flat eigendirection of $\Lambda_{\text{obs}}$; a Chapter 20
experiment aimed at that direction supplies the missing precision, and the Chapter 18 allocation
posterior — wide because of the ridge — narrows. **Anchor (exact):** ridge precision
$\Lambda_{\text{obs}}=1$ (variance $1$, std $1$) along the contrast direction $u=(1,-1)/\sqrt2$, a
calibration precision $1/s^2=3$ → $\Lambda_{\text{post}}=4$, variance $0.25$, std $0.5$: the ridge std
**halves**, variance drops $4\times$ (the precision ratio).

**Rung 3 — Proof P2 (KEYSTONE): the secant–tangent matching relationship.**
A finite lift test measures a **secant**; the optimizer wants a **tangent**. **Theorem (Taylor with
remainder):** for $S\in C^2$,
$$
\frac{S(x{+}\delta)-S(x)}{\delta} = S'(x) + \frac{\delta}{2}\,S''(\xi), \qquad \xi\in(x,x{+}\delta).
$$
Under concavity ($S''<0$) the secant **underestimates** the tangent $S'(x)$ at the lower endpoint by
$\tfrac{\delta}{2}|S''(\xi)|$ — an $O(\delta)$ bias. **Corollary (paired bracketing):** averaging a
matched up-test and down-test about $x$ cancels the $O(\delta)$ term, recovering $S'(x)$ to
$O(\delta^2)$:
$\tfrac12\big(\tfrac{S(x+\delta)-S(x)}{\delta}+\tfrac{S(x)-S(x-\delta)}{\delta}\big)=S'(x)+O(\delta^2)$.
This is the rigorous cost of calibrating a tangent with a secant, and the design fix. **Anchors:**
$S=2\sqrt x$, $x=4$, $\delta=5$: secant $0.4$, tangent $S'(4)=0.5$, gap $-0.1$, with $S''(\xi)=-0.04$
at $\xi\approx5.39$. Paired about $x=9$, $\delta=4$: up-secant $0.303$, down-secant $0.382$, average
$0.342$ vs $S'(9)=0.333$ — the average is second-order accurate where each secant is off by $\approx0.04$.

**Rung 4 — Identification: test-design ↔ curve shape.**
A saturating curve has shape parameters (a scale $\theta$, and saturation parameters $\kappa,\alpha$ in
the general Hill form). **One** secant or tangent pins **one** functional — it under-identifies the
shape: many curves pass through a single measured lift. State the identification principle: to pin the
*shape* you need a **spread of constraints** — lifts at several spend levels, or a wide secant plus a
tangent — so the experiments jointly determine $(\theta,\kappa,\alpha)$. This is the **active-learning
hook**: the Chapter 18 optimizer flags the high-leverage spend region, and the next experiment is
aimed there (closing the loop). Connect the geometry: each constraint is a surface in parameter space;
identification is their intersection collapsing to a point.

**Rung 5 — Power-prior down-weighting (brief) and the handoff to Chapter 22.**
Not every experiment deserves full weight: a stale or low-credibility study should enter
down-weighted. Introduce the **power-prior** factor $L_{\text{exp}}(\hat g\mid\theta)^{w}$ with
$w\in[0,1]$ = (design credibility) × (temporal decay $\rho^{\,t_{\text{now}}-t_k}$); $w=1$ is full
trust, $w=0$ is a validation-tier check that never enters the likelihood. State that a *collection* of
such weighted factors, accumulated and versioned, is the **prior store** — developed as a data product
in Chapter 22 (schema = the secant/tangent/validation taxonomy; sequential update; hierarchical pooling
of conflicting studies). Keep this brief: the mechanism here, the system there. Close on the full loop:
**intervene (Ch. 19–20) → identify the functional → calibrate (here) → re-optimize (Ch. 18)** — and the
EVPI gap, named in Chapter 18, shrinks.

---

## 4. Worked Examples

**WE1 — Secant calibration tightens the curve (numeric).**
Scale model $S(x;\theta)=\theta\sqrt x$ (truth $\theta=2$). Observational posterior wide:
$\theta\sim\mathcal N(2,0.5^2)$ (precision $4$). Chapter 20's geo lift measured $\hat\Delta=2$ over
$[4,9]$, and $S(9;\theta)-S(4;\theta)=\theta(\sqrt9-\sqrt4)=\theta$, so the secant factor is
$\hat\Delta\sim\mathcal N(\theta, s^2)$ with $s=0.2$ (precision $25$). Posterior precision
$4+25=29$, variance $1/29=0.0345$, std $0.186$, mean $(2/0.25+2/0.04)/29=2.0$: the lift pins the scale,
tightening $\theta$ from std $0.5$ to $0.186$ (a $2.7\times$ reduction). The fitted curve's marginal
returns — what the optimizer differentiates — tighten with it.

**WE2 — The secant–tangent gap, and paired bracketing.**
On $S=2\sqrt x$: a lift test moving spend $4\to9$ returns the secant $0.4$, but the marginal return at
$x=4$ is $S'(4)=0.5$ — calibrating $S'(4)$ with this secant biases it low by $0.1$, the
$\tfrac{\delta}{2}S''(\xi)$ term ($\xi\approx5.39$). Show the paired fix about $x=9$, $\delta=4$:
up-secant $0.303$, down-secant $0.382$, average $0.342$ vs $S'(9)=0.333$ — the average brackets the
tangent and is second-order accurate. Lesson: match the estimand class to the parameter you are
calibrating, or pay a known-sign, $O(\delta)$ bias.

**WE3 — Healing the ridge (the wound closes).**
Reuse Chapter 18's ridge posterior on $(a_1,a_2)$ centered $(2,1)$: the contrast direction
$u=(1,-1)/\sqrt2$ has std $0.45$ (the ridge), giving the Chapter 18 EVPI $\approx0.12$. A geo
experiment supplies a calibration factor on the contrast (it measures the attribution split the
observational data could not), **halving** the ridge std to $0.225$. Propagate through the Chapter 18
optimizer: the allocation posterior narrows and the optimizer's-curse / **EVPI gap shrinks from
$\approx0.12$ to $\approx0.03$** (a $\sim4\times$ reduction, tracking the variance). State plainly: this
is the wound of Chapter 18 closing — the calibration bought back the identification the observational
series lacked.

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + `matplotlib`; verified headless). It:
1. **Precision addition (P1):** build $\Lambda_{\text{obs}}$ with a flat ridge direction; add a
   rank-one calibration precision $\tfrac{1}{s^2}cc^\top$; print the ridge variance before/after and
   assert the $4\times$ drop; figure: posterior covariance ellipse before (elongated ridge) vs after
   (compressed).
2. **Secant calibration (WE1):** conjugate Gaussian update of $\theta$ from the lift; assert posterior
   mean $\approx2$, std $\approx0.186$.
3. **Secant–tangent (WE2):** compute the secant, the tangent, the gap, and the paired-test average;
   assert the average is closer to $S'(9)$ than either single secant.
4. **EVPI heals (WE3):** Monte-Carlo EVPI (the Chapter 18 computation) at ridge std $0.45$ then
   $0.225$; assert it drops from $\approx0.12$ to $\approx0.03$; figure: EVPI vs ridge std with the
   before/after points marked.
Every numeric claim asserted; seed `np.random.default_rng(21)`.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why does an experiment enter as a likelihood factor rather than a prior on ROAS, and what does
  the likelihood-factor view avoid? Why does a calibration constraint help most in the ridge direction?
  Why does one lift under-identify the curve shape?
- **B:** On $S=2\sqrt x$: compute the secant–tangent gap for $x=4,\delta=5$ and the paired-test average
  about $x=9,\delta=4$; and do the conjugate Gaussian update of WE1 (precision addition $4+25=29$).
- **P:** (i) Prove the posterior-precision addition $\Lambda_{\text{post}}=\Lambda_{\text{obs}}+\tfrac1{s^2}cc^\top$
  and that the variance along $c$ strictly decreases. (ii) Prove the secant–tangent Taylor relation and
  the second-order accuracy of the paired average.
- **A:** Extend the Code Tie-in — vary the calibration precision $1/s^2$ and plot the residual ridge
  variance (and the EVPI) against it; add a *second* calibration constraint at a different spend level
  and show the curve *shape* (two parameters) becomes identified where one constraint did not.

---

## 7. Appendix solutions

Append `## Chapter 21 — Advanced Calibration` to `appendix/solutions.qmd`, **in chapter order** (after
the Chapter 20 block), inside the existing `content-visible` gated div. Full C/B/P/A; the P-block
carries both proofs (precision-addition / variance decrease; secant–tangent Taylor + paired bracketing).

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities: the
evidence-synthesis posterior $p(\theta\mid D,\hat g)\propto p(\theta)L_{\text{obs}}L_{\text{exp}}$;
secant factor $\hat\Delta\sim\mathcal N(S(x{+}\delta;\theta)-S(x;\theta),s^2)$; precision addition
$\Lambda_{\text{post}}=\Lambda_{\text{obs}}+\tfrac1{s^2}cc^\top$; secant–tangent
$\frac{S(x+\delta)-S(x)}{\delta}=S'(x)+\frac{\delta}{2}S''(\xi)$; power-prior weight $w\in[0,1]$. Close
tying back to Chapter 18 (the EVPI gap shrinks) and forward to Chapter 22 (the prior store).

---

## 9. Cross-references (post-renumber numbering)

- **Chapter 16 (Building & Fitting):** the Fisher-information ridge this chapter fills in.
- **Chapter 18 (Budget Optimization):** the wound / EVPI gap that shrinks; the allocation posterior
  that narrows.
- **Chapter 19 (Causal Foundations):** the confounding that made calibration necessary.
- **Chapter 20 (Quasi-Experimental Design):** the taxonomy-tagged estimands ingested here.
- **Chapter 22 (The Prior Store):** the accumulating, weighted, versioned collection of these
  constraints as a data product.

---

## 10. Bibliography

Reuse `@gelman2013` (Bayesian decision / evidence synthesis / informative priors — already in
`references.bib`). Add for the power-prior brief:
- `@ibrahim2000` — Ibrahim & Chen, "Power Prior Distributions for Regression Models" (Statistical
  Science 2000).
The plan confirms exact BibTeX and that all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor NumPy-verified (done: precision addition variance $1\to0.25$; secant $0.4$ vs
  tangent $0.5$; paired average $0.342$ vs $0.333$; WE1 std $0.5\to0.186$; EVPI $0.12\to0.03$).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`).
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no `\begin{align}`,
  no `psmallmatrix`, valid citation keys, no banned library names.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
