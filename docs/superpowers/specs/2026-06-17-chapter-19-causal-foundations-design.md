# Chapter 19 — Causal Inference Foundations: Design Spec

**Date:** 2026-06-17
**Part:** VI — Causal Grounding & Calibration (opening chapter)
**File:** `parts/06-causal-calibration/01-causal-foundations.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter **opens Part VI** and answers the question Chapter 18 left open as a wound: *why can't
observational MMM pin the response curve?* The answer is **endogeneity** — spend is not randomly
assigned. Managers spend more when they already expect high sales, so spend is correlated with
demand drivers (seasonality, promotions, competitive actions, pricing) that *also* move sales. That
correlation makes the fitted response coefficient biased **as a property of the data-generating
process, not the sample size** — no amount of observational data removes it.

Chapter 18 framed the wound as *weak identification* (a wide posterior ridge). This chapter sharpens
it: a meaningful part of that ridge is **structural confounding bias**, and the only escape is to
change how the data is generated — to *intervene* (Chapter 20's quasi-experimental designs) and fold
the result back into the model (Chapter 21's calibration). Causal foundations is the conceptual hinge
that makes Chapters 20–21 *necessary* rather than optional.

**Through-line driving question:** *"A fitted coefficient is causal only under conditions we must
name."* This chapter names them.

**Framework choice (locked):** lead with **potential outcomes** (Neyman–Rubin) — the bias is most
concrete as omitted-variable bias in the fitted MMM regression — then bridge to **structural causal
models / DAGs and the do-operator** (Pearl) as the graphical language that makes the identification
condition visual and general, proving the two are equivalent (ignorability ⇔ backdoor criterion).

---

## 2. Scope discipline (what this chapter does NOT do)

- It does **not** develop the quasi-experimental designs (DiD, synthetic control, ITS, IV, BSTS,
  RDD/RKD) — those are **Chapter 20**. IV is *teased* here only as "the natural escape when the
  confounder is unobserved," not worked.
- It does **not** build the calibration likelihood or the prior store — **Chapters 21–22**.
- **No matching / propensity-score methods** (decided out of scope for the book).
- House rules: **never** name PyMC-Marketing or any MMM/PPL/sampler library; `numpy`/`scipy`/
  `matplotlib` only. **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$ … $` only; `$$` on
  their own lines; even `$$` count.

---

## 3. Theory & Proofs — the rungs

**Rung 1 — The causal question in MMM (potential outcomes).**
Define potential outcomes: $Y_i(x)$ is period $i$'s sales *if* spend were set to $x$. The
dose–response (response curve) is $\mu(x) = \mathbb E[Y_i(x)]$ — exactly the curve the MMM wants and
the optimizer of Chapter 18 needs. The **fundamental problem of causal inference**: we observe only
$Y_i(X_i)$ at the spend $X_i$ that actually occurred, never the counterfactual. The fitted MMM
regresses observed sales on observed spend; whether its slope equals $d\mu/dx$ depends on *how spend
was chosen*. Spend is not randomized — it is an endogenous management decision.

**Rung 2 — Ignorability and the observational/interventional gap.**
Define **unconfoundedness / ignorability**: $Y_i(x) \perp X_i \mid Z_i$ for a covariate set $Z$.
Interpretation: among periods with the same demand state $Z$, the level of spend is "as good as
random." Under ignorability (plus positivity/overlap), the causal dose–response is recovered by
adjustment: $\mu(x) = \mathbb E_Z[\,\mathbb E[Y \mid X=x, Z]\,]$. Without it, the observational
regression $\mathbb E[Y\mid X=x]$ measures something else. State the estimand the chapter targets:
the response slope $d\mu/dx$ (in MMM terms, the marginal return the optimizer consumes).

**Rung 3 — Proof P1 (KEYSTONE): the confounding / omitted-variable-bias theorem.**
Linear structural MMM: $Y = \beta X + \gamma Z + \varepsilon$ (true response slope $\beta$),
with endogenous spend $X = \delta Z + \nu$ (spend tracks demand $Z$). The naive regression of $Y$ on
$X$ alone has
$$
\operatorname{plim}\hat\beta_{\text{naive}} = \frac{\operatorname{cov}(X,Y)}{\operatorname{var}(X)}
 = \beta + \gamma\,\frac{\operatorname{cov}(X,Z)}{\operatorname{var}(X)}.
$$
**Proof:** $\operatorname{cov}(X,Y) = \beta\operatorname{var}(X) + \gamma\operatorname{cov}(X,Z) +
\operatorname{cov}(X,\varepsilon)$; with $\varepsilon$ orthogonal to $X$, divide by
$\operatorname{var}(X)$. The bias term $\gamma\,\operatorname{cov}(X,Z)/\operatorname{var}(X)$ is zero
**iff** $\operatorname{cov}(X,Z)=0$ (spend exogenous) or $\gamma=0$ ($Z$ does not affect sales). The
MMM reading: when spend is endogenous and the demand driver moves sales, the fitted response slope is
biased, and **the bias does not shrink with more data** — it is a property of the joint distribution,
the rigorous statement of Chapter 18's wound. **Anchor:** $\beta=2$, $\gamma=3$, $\delta=1$,
$\operatorname{var}(Z)=\operatorname{var}(\nu)=1 \Rightarrow$ bias $=3\cdot\frac{1}{2}=1.5$, so
$\operatorname{plim}\hat\beta_{\text{naive}} = 3.5$ — a $75\%$ overstatement of the true $\beta=2$.

**Rung 4 — Structural causal models, DAGs, and the do-operator (the bridge).**
Introduce the SCM: each variable a function of its parents plus noise; the DAG draws the parent
arrows. The confounded MMM is the **fork** $X \leftarrow Z \rightarrow Y$ with $X \rightarrow Y$:
$Z$ opens a **backdoor path** $X \leftarrow Z \rightarrow Y$ that carries non-causal association.
Define the **do-operator**: $P(Y \mid \mathrm{do}(X=x))$ is the distribution in the *mutilated* graph
where $X$'s incoming arrows are cut and $X$ is set to $x$ (the **truncated factorization**). Contrast
$P(Y\mid \mathrm{do}(x))$ (intervention — cut $Z\to X$) with $P(Y\mid x)$ (observation — $Z$ still
covaries with $X$). In the linear-Gaussian MMM the interventional slope is $\beta$ and the
observational slope is $\beta + \text{bias}$ — the same gap as Rung 3, now graphical.

**Rung 5 — Proof P2: backdoor adjustment and identification.**
**Theorem (backdoor adjustment).** If $Z$ satisfies the backdoor criterion relative to $(X,Y)$ —
blocks every backdoor path and contains no descendant of $X$ — then
$$
P(Y \mid \mathrm{do}(X=x)) = \sum_z P(Y \mid X=x, Z=z)\,P(Z=z).
$$
**Proof** via truncated factorization: the interventional joint deletes $X$'s mechanism, so
$P(\cdot\mid \mathrm{do}(x)) = \big[\prod_{V\neq X} P(V\mid \mathrm{pa}_V)\big]_{X=x}$; marginalize out
$Z$ and use that $Z$'s distribution is unchanged by the intervention. **Bridge:** show the graphical
backdoor criterion is exactly the potential-outcomes **ignorability** $Y(x)\perp X\mid Z$ — two
languages, one condition. Conclusion: observational MMM recovers the causal curve **iff** a sufficient
adjustment set is observed and conditioned on.

**Rung 6 — The wound, restated causally, and the handoff.**
In MMM the sufficient adjustment set is usually **not observed**: you cannot measure every demand
shock, competitive action, or pricing decision that drove both spend and sales. When $Z$ is
unobserved the backdoor cannot be closed, ignorability fails, and **no observational estimator
recovers $\beta$** — the bias is irreducible. This is the causal content of Chapter 16's ridge and
Chapter 18's wound. The escape is to change the data-generating process: **Chapter 20** supplies
*interventional* identification (geo experiments, holdouts, IV-style natural experiments — IV is the
unobserved-confounder escape, teased here), and **Chapter 21** folds those estimates back into the
posterior as calibration. Name the loop: intervene → identify → calibrate → re-optimize.

---

## 4. Worked Examples

**WE1 — A confounded MMM, and adjustment in action (numeric).**
Concrete DGP: demand $Z\sim N(0,1)$; spend $X = Z + \nu$ ($\nu\sim N(0,1)$ — managers spend more in
high-demand periods); sales $Y = 2X + 3Z + \varepsilon$ (true response slope $\beta=2$). Show:
(i) the naive OLS slope of $Y$ on $X$ is $\approx 3.50$ — biased high by $1.5$, matching the analytic
$\gamma\,\operatorname{cov}(X,Z)/\operatorname{var}(X) = 3\cdot\tfrac12$; (ii) regressing $Y$ on
$X$ **and** the observed confounder $Z$ recovers $\hat\beta\approx 2.00$ (the backdoor adjustment of
Rung 5, made arithmetic). The fitted curve is wrong by $75\%$ until the confounder is controlled.

**WE2 — The unobserved confounder (the wound).**
Same DGP, but now $Z$ is *unobserved* (a latent demand index the analyst never measured). The only
feasible regression is the naive one, stuck at $\hat\beta\approx 3.50$; the adjustment of WE1 is
unavailable. Quantify the irreducible bias and state plainly: every observational estimator built
from $(X,Y)$ alone inherits it — this is Chapter 18's wound in causal language, and Chapter 16's
ridge is its posterior shadow. Tease: an instrument (a driver of $X$ that does not affect $Y$ except
through $X$) would escape, but that is Chapter 20.

**WE3 — Observation vs intervention via the do-operator.**
Compute both slopes on the same DGP: the **observational** slope $dE[Y\mid X=x]/dx \approx 3.5$ versus
the **interventional** slope $dE[Y\mid \mathrm{do}(X=x)]/dx = 2$. Show the intervention numerically by
*severing* $Z\to X$ — regenerate spend independent of demand ($X=\nu$, randomized assignment), keep
the same structural $Y=2X+3Z+\varepsilon$, and recover $\hat\beta\approx 2.0$ from the naive
regression. The act of randomizing spend closes the backdoor by construction — the mathematical
content of "run an experiment."

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + `numpy.linalg.lstsq` + `matplotlib`; verified headless).
Seeded `np.random.default_rng(19)`, $N$ large. It:
1. Simulates the confounded DGP ($\beta=2,\gamma=3,\delta=1$).
2. **Naive** OLS ($Y\sim X$): recovers $\approx 3.50$; asserts the analytic OVB $=1.5$.
3. **Adjusted** OLS ($Y\sim X+Z$): recovers $\approx 2.00$ (backdoor adjustment).
4. **Hidden confounder:** show adjustment is unavailable; naive stuck at $3.50$.
5. **Intervention (do):** regenerate $X$ independent of $Z$, re-fit naive: recovers $\approx 2.00$.
6. Figure: scatter of $(X,Y)$ with the steep observational fit ($3.5$) and the true/interventional
   slope ($2.0$) overlaid; optionally colored by $Z$ to show confounding visually.
Every numeric claim asserted against the analytic values.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why is a fitted MMM coefficient not automatically causal? State the fundamental problem of
  causal inference in MMM terms. Why does confounding bias *not* vanish with more data? What does
  "ignorability" mean operationally for spend, and why is it usually indefensible in observational
  MMM?
- **B:** Derive $\operatorname{plim}\hat\beta_{\text{naive}} = \beta + \gamma\operatorname{cov}(X,Z)/
  \operatorname{var}(X)$ for the linear MMM; evaluate the bias for the WE1 parameters; show it is zero
  when $\operatorname{cov}(X,Z)=0$.
- **P:** (i) Prove the OVB theorem (P1). (ii) Prove the backdoor adjustment formula (P2) via
  truncated factorization, and show ignorability $Y(x)\perp X\mid Z$ is equivalent to $Z$ blocking
  the backdoor path in the fork $X\leftarrow Z\rightarrow Y$.
- **A:** Extend the Code Tie-in — sweep the confounding strength $\delta$ and plot naive bias vs
  $\delta$ (linear in $\operatorname{cov}(X,Z)$); add a second observed confounder and show the
  adjustment set must be *complete*; simulate an instrument and recover $\beta$ (a Chapter 20 teaser).

---

## 7. Appendix solutions

Append a `## Chapter 19 — Causal Inference Foundations` block to `appendix/solutions.qmd`, **in
chapter order** (after the Chapter 18 block), inside the existing `content-visible` gated div. Full
C/B/P/A solutions; the P-block carries both proofs in full (OVB plim derivation; truncated-
factorization backdoor proof + the ignorability⇔backdoor equivalence).

---

## 8. Summary (auto-included)

Bulleted **Key concepts** and bulleted **Key identities** (inline math; identities bulleted, never a
run-on paragraph). Identities: potential-outcome dose–response $\mu(x)=\mathbb E[Y(x)]$; ignorability
$Y(x)\perp X\mid Z$; OVB $\operatorname{plim}\hat\beta=\beta+\gamma\operatorname{cov}(X,Z)/
\operatorname{var}(X)$; backdoor adjustment $P(Y\mid\mathrm{do}(x))=\sum_z P(Y\mid x,z)P(z)$;
observation-vs-intervention slope gap. Close by tying forward: Chapter 20 (interventional
identification) and Chapter 21 (calibration) heal the wound this chapter diagnoses.

---

## 9. Cross-references (post-renumber numbering)

- **Chapter 4 (Linear Regression):** the endogeneity teaser planted there ("a fitted coefficient is
  causal only under conditions we'll formalize later") is cashed out here.
- **Chapter 16 (Building & Fitting):** the Fisher-information ridge — its causal reading is
  confounding when the unidentified direction is a demand confounder.
- **Chapter 18 (Budget Optimization):** the wound (optimizing over a weakly/biasedly identified
  curve) is what this chapter explains structurally.
- **Chapter 20 (Quasi-Experimental Designs):** supplies the interventional identification (IV teased
  here as the unobserved-confounder escape).
- **Chapter 21 (Advanced Calibration):** folds interventional estimates back into the posterior.

---

## 10. Bibliography

Add to `references.bib` (confirm keys do not already exist):
- `@imbens2015` — Imbens & Rubin, *Causal Inference for Statistics, Social, and Biomedical Sciences*
  (2015) — potential outcomes.
- `@angrist2009` — Angrist & Pischke, *Mostly Harmless Econometrics* (2009) — OVB, regression
  causality, IV.
- `@pearl2009` — Pearl, *Causality* (2nd ed., 2009) — SCMs, do-calculus, backdoor.
Reuse `@gelman2013` where a Bayesian-decision reference is apt. The plan confirms exact BibTeX.

---

## 11. Verification (the real gate)

- Every numeric anchor NumPy-verified before written (done: naive $3.50$, OVB $1.5$, adjusted $2.00$,
  do-slope $2.0$ at seed 19).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`).
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no
  `\begin{align}`, no `psmallmatrix`, valid citation keys, no banned library names.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
