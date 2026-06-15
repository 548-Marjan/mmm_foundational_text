# Chapter 14 — The MMM Data-Generating Process: Design Spec

**Date:** 2026-06-14
**Status:** Approved (brainstorm complete). Next: implementation plan via `writing-plans`.
**Part:** V — Marketing Mix Modeling: Modeling (opening chapter).
**File:** `parts/05-mmm-modeling/01-mmm-dgp.qmd` (currently a 33-line template stub).
**Canonical anchor:** Jin et al. 2017, *Bayesian Methods for Media Mix Modeling with
Carryover and Shape Effects* (`@jin2017`). **Do NOT cite `@pymcmarketing` or name any
library** — locked house rule.

---

## 1. Role in the book

Chapter 14 opens Part V. It is the chapter where the marketing response curve stops
being *given* and becomes *defined*: the generative story of how spend produces sales,
specified **before** any fitting (estimation/inference is Chapter 15). Two structural
payoffs it must deliver:

1. **Supply Part IV's missing floor.** Chapters 11–13 optimized over a Hill saturation
   curve that was *used before it was defined*. Chapter 14 is its definitional source,
   using the **same parameterization** Part IV used so the cross-reference is exact.
2. **End on the wound.** The DGP is constructed so that observational, collinear,
   seasonally-confounded spend **cannot** identify it. This is the door Part VI
   (Ch. 18–21) is built to walk through. Ch. 14 *opens* the wound rigorously; it does
   **not** resolve it (no IV, no experiments, no priors-as-fix here).

The chapter follows the fixed template: Motivation → Theory & Proofs → Worked Examples →
Code Tie-in → Exercises (C/B/P/A) → Summary (auto-included), with appendix solutions
gated by `show-solutions`.

## 2. Scope decisions (locked in brainstorming)

- **Identification depth:** a *compact but real* non-identifiability theorem — state it,
  prove it by construction, and stop, deferring the cure explicitly to Part VI. Not a
  sprawling identification analysis (that would poach Part VI).
- **Non-media components** (baseline, trend, seasonality, controls): first-class in the
  structural equation but **lightly developed** — defined and motivated in a sentence or
  two each, no proofs. They exist to (a) write the full generative model and (b) power
  the "spend correlated with seasonal demand" confounding in the identification theorem.
  Full theory-and-proof weight goes to the two media transforms.
- **Media transforms get the weight:** adstock (carryover) and saturation (shape).

## 3. Theory rungs

1. **The structural MMM equation.**
   $$y_t = \text{baseline}_t + \text{trend}_t + \text{seasonality}_t
        + \sum_c \beta_c\, g\big(a_\lambda(x_{c,t})\big) + \gamma^\top z_t + \varepsilon_t.$$
   The generative model as a whole; each term named and motivated. Non-media terms light.
2. **Adstock / carryover.** Geometric adstock $a_t = x_t + \lambda a_{t-1}$,
   $\lambda \in [0,1)$; normalized form (weights sum to 1); the **delayed-peak** variant
   of Jin et al. (carryover that peaks a few periods out) described qualitatively.
3. **Saturation / shape.** The Hill curve $g(x;K,n) = \dfrac{x^n}{K^n + x^n}$ defined
   here, with the **same symbols** Part IV uses ($K$ = half-saturation point, $n$ = shape
   exponent). $T$ denotes series length, so $n$ is free for the Hill exponent.
4. **Assembling the DGP** + the **order-of-operations** decision (adstock-then-saturate,
   per Jin et al.) — set up for proof P3.
5. **The identification wound.** Why observational spend (collinear across channels,
   confounded with seasonality, often at a narrow operating range) cannot recover the
   curve — set up for proof P4, hand off to Part VI.

## 4. Proof slate (four proofs)

| # | Theorem | Rung | Statement |
|---|---|---|---|
| **P1** | Adstock geometric-series identities | 2 | Normalized weights $(1-\lambda)\lambda^k$ sum to 1 (a constant input passes through at its level); unnormalized sustained-spend steady state $x/(1-\lambda)$; half-life $\ln 0.5/\ln\lambda$; mean lag $\lambda/(1-\lambda)$. Pure geometric-series argument. |
| **P2** | Hill shape theorem | 3 | $g$ is strictly increasing and bounded in $(0,1)$; for $n\le 1$ it is concave on $(0,\infty)$ (diminishing returns from the first dollar); for $n>1$ it is S-shaped with a single inflection at $x^\star = K\big(\tfrac{n-1}{n+1}\big)^{1/n}$ and is concave on $(x^\star,\infty)$. **This is the diminishing-returns floor Part IV assumed.** |
| **P3** | Non-commutativity of adstock and saturation | 4 | $g\circ a_\lambda \neq a_\lambda\circ g$ in general, shown by an explicit numeric counterexample — so the transform ordering is a genuine modeling choice, not a notational convention. |
| **P4** | Non-identifiability (compact) | 5 | Under spend confined to a single operating point (or two channels in a fixed ratio), exhibit **distinct** parameter sets that produce **identical** $\mathbb{E}[y]$ on the observed support — e.g. a $\beta$–$K$ trade-off, or a contribution split between collinear channels. Prove by construction. Close with an explicit statement that breaking this requires the exogenous variation supplied in Part VI. |

## 5. Worked examples + anchor numbers

All numerics to be NumPy-verified in the plan.

- **WE1 — geometric adstock by hand.** $\lambda = 0.5$. Impulse spend
  $[100, 0, 0, \dots] \to [100, 50, 25, 12.5, \dots]$. Half-life $\ln 0.5/\ln 0.5 = 1$
  period; mean lag $\lambda/(1-\lambda) = 1$; sustained spend $100$/period $\to$ steady
  state $100/(1-0.5) = 200$.
- **WE2 — Hill curve.** $K = 3$, $n = 2 \Rightarrow$ inflection
  $x^\star = 3\,(1/3)^{1/2} = \sqrt 3 \approx 1.732$. **This is exactly the curve Part IV
  optimized over** (Ch. 11's inflection anchor was $\sqrt 3$) — closing the Part IV loop.
  Compute $g$ and marginal $g'$ at a couple of operating points.
- **WE3 — identification.** A concrete construction: two $(\beta, K)$ pairs (or two
  contribution splits between fixed-ratio channels) that fit the observed support
  identically, making P4 tangible.

## 6. Code tie-in

A single runnable `{python}` cell that **simulates the full DGP** over $T$ weeks:
baseline + trend + seasonality + two channels (adstock → saturate) + Gaussian noise → a
synthetic weekly sales series. It plots (a) the geometric adstock decay, (b) the Hill
saturation curves, and (c) the assembled sales series. The simulated series is the
**dataset Chapter 15 will fit** — an intentional hand-off, not a throwaway. Runs
top-to-bottom under `MPLBACKEND=Agg python3`; every pinned number asserted against NumPy.

## 7. Exercises, Summary, housekeeping

- **Exercises** C/B/P/A (self-contained, no inline solution links); appendix solutions
  appended to `appendix/solutions.qmd` as `## Chapter 14 — The MMM Data-Generating
  Process`, gated by `::: {.content-visible when-meta="show-solutions"}`, before the final
  closing `:::`.
- **Summary** auto-included; **Key concepts** and **Key identities** both as **bulleted
  lists** (inline math only) — per the formatting-consistency rule. No run-on paragraph.
- **Anchor line:** replace the stub's `*Canonical anchors: @jin2017; @pymcmarketing.*`
  with `*Canonical anchors: @jin2017.*` (scrub the library citation).
- **Bundled fix:** correct the **two** Chapter 11 cross-references that attribute the Hill
  curve to "Chapter 5" — re-point both to Chapter 14, now its real source. They are at
  `parts/04-optimization/01-convexity.qmd:19` ("the diminishing-returns saturation curve
  introduced in Chapter 5") and `:69` ("the **Hill saturation curve** introduced in
  Chapter 5"). **Leave `:200` untouched** — that "Chapter 5" reference points at the
  ridge/MAP Gaussian model, which genuinely is Chapter 5 and is correct.

## 8. Conventions (enforced per task)

- KaTeX: `aligned` inside `$$ … $$` (never bare `\begin{align}`); `$$` delimiters on their
  own lines; even `$`-count per line. **Never** `\begin{psmallmatrix}` (PDF/LuaLaTeX
  undefined — passes KaTeX/HTML, breaks the PDF build); use `bmatrix`/`pmatrix`/`smallmatrix`.
- Build-based verification: the chapter's code cell runs headless; the real gate is CI
  `quarto render` (HTML + PDF) on the PR — watch **both** jobs to green before reporting done.
- Library-agnostic: name no MMM library in the text. `numpy`/`scipy`/`matplotlib` are fine.

## 9. Critical files

- `parts/05-mmm-modeling/01-mmm-dgp.qmd` — the chapter (replace stub body; keep H1
  `# The MMM Data-Generating Process`; fix the anchor line).
- `appendix/solutions.qmd` — append the Chapter 14 solutions block before the final `:::`.
- `parts/04-optimization/01-convexity.qmd` — the two Hill cross-reference fixes (lines 19
  and 69); leave line 200 (ridge/MAP = Chapter 5, correct) untouched.
- Voice/structure exemplars (read-only): `parts/04-optimization/03-slsqp.qmd` (immediate
  predecessor), `parts/04-optimization/01-convexity.qmd` (defines the Part IV Hill usage).
