# Chapter 14 — The MMM Data-Generating Process Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 14, the opener of Part V, which defines the MMM data-generating process (adstock carryover + Hill saturation atop lightly-developed baseline/trend/seasonality), proves four properties of it, simulates it into the dataset Chapter 15 will fit, and ends on the identification "wound" that Part VI heals.

**Architecture:** Replace the template stub `parts/05-mmm-modeling/01-mmm-dgp.qmd` with a full chapter following the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary). Append a gated solutions block to `appendix/solutions.qmd`. Two bundled housekeeping edits: scrub the `@pymcmarketing` anchor; fix the two Chapter 11 Hill cross-references. The chapter is the definitional source of the Hill curve Part IV already optimizes.

**Tech Stack:** Quarto (`.qmd`), KaTeX math, one `{python}` cell (numpy + matplotlib, run headless under `MPLBACKEND=Agg`). Real verification gate is CI `quarto render` (HTML + PDF).

---

## Conventions (enforce in EVERY task)

- **KaTeX:** display math uses `$$ … $$` with delimiters **on their own lines**; multi-line uses `aligned` **inside** `$$ … $$` — never bare `\begin{align}`. Every line has an **even** count of `$`. Inline math uses single `$`.
- **PDF safety:** **never** use `\begin{psmallmatrix}` (undefined in the CI LuaLaTeX/PDF build — it renders in KaTeX so HTML passes but the PDF stage fails). Use `bmatrix`, `pmatrix`, or `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]`. (No matrices are expected in this chapter, but hold the rule.)
- **Library-agnostic:** name **no** MMM library anywhere in the prose. `numpy`/`scipy`/`matplotlib` are fine to name. The only canonical anchor is `@jin2017`.
- **Notation (fixed for the whole chapter):** $K$ = half-saturation point, $n$ = Hill shape exponent, $\lambda$ = adstock decay, $T$ = series length (so $n$ is free for the exponent), $\beta_c$ = channel effect, $g$ = saturation, $a_\lambda$ = adstock.
- **Verification:** the single code cell must run top-to-bottom under `MPLBACKEND=Agg python3` with all asserts passing; end figures with `plt.show()`. Pinned numbers below are NumPy-verified — keep them exact.
- **Commits:** one commit per task, message prefix `feat(ch14): …`. Use `git -C <worktree>` or ensure CWD is the worktree root. Identity is already set (`jlh530i` / `jlh530i@gmail.com`).

## Verified anchor numbers (do not change)

- **Adstock** $\lambda=0.5$: impulse $[100,0,0,\dots]\to[100,50,25,12.5,\dots]$; half-life $\ln 0.5/\ln 0.5 = 1$; mean lag $\lambda/(1-\lambda)=1$; sustained-spend steady state $100/(1-0.5)=200$.
- **Hill** $K=3,\ n=2$: inflection $x^\star=3(1/3)^{1/2}=\sqrt3\approx1.7320508$; $g(K)=g(3)=0.5$; $g(\sqrt3)=0.25$; $g'(3)=1/6\approx0.1666667$; $g'(\sqrt3)=\sqrt3/8\approx0.2165064$.
- **Non-commutativity** (spend $(3,0)$, $\lambda=0.5$, $K=3,n=2$): adstock-then-saturate period-1 value $g(1.5)=0.2$; saturate-then-adstock period-1 value $0.5\cdot g(3)=0.25$. $0.2\neq0.25$.
- **Non-identifiability** (single operating point $x_0=3$): $(\beta,K,n)=(1,3,2)$ gives contribution $1\cdot g(3)=0.5$; $(\beta,K,n)=(2,\,3\sqrt3,\,2)$ gives $2\cdot g(3;3\sqrt3,2)=2\cdot0.25=0.5$. Identical on the observed support; $3\sqrt3\approx5.196152$.

## File structure

- **`parts/05-mmm-modeling/01-mmm-dgp.qmd`** — the chapter. Replace the stub body; keep H1 `# The MMM Data-Generating Process`; fix the anchor line.
- **`appendix/solutions.qmd`** — append `## Chapter 14 — The MMM Data-Generating Process` solutions block (gated) before the file's final closing `:::`.
- **`parts/04-optimization/01-convexity.qmd`** — fix the two Hill cross-references (lines ~19 and ~69 say "Chapter 5"; re-point to Chapter 14). Leave the line ~200 ridge/MAP "Chapter 5" reference untouched (it is correct).
- Read-only exemplars: `parts/04-optimization/03-slsqp.qmd` (predecessor voice/structure), `parts/04-optimization/01-convexity.qmd` (how Part IV uses the Hill curve).

---

### Task 1: Front matter, anchor scrub, and Motivation

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (replace stub header block + write Motivation)

- [ ] **Step 1: Fix the H1 and anchor line.** Keep the H1 `# The MMM Data-Generating Process`. Replace the stub anchor line `*Canonical anchors: @jin2017; @pymcmarketing.*` with exactly:

```
*Canonical anchor: @jin2017.*
```

Delete the `::: {.callout-note … Stub …}` block entirely.

- [ ] **Step 2: Write the `## Motivation` section.** 3–5 paragraphs. Content beats to hit:
  - Part IV took the response curve $r_i(\cdot)$ as *given* and optimized over it. This chapter writes down where that curve comes from — the generative model of how weekly spend becomes weekly sales.
  - The MMM is a structural regression: sales = baseline demand the brand would have had anyway + the incremental effect of each channel's spend, passed through two nonlinear transforms (carryover in time, diminishing returns in amount) + noise.
  - Name the two transforms as the chapter's stars: **adstock** (today's spend keeps working for weeks) and **saturation** (the tenth million dollars buys less than the first). Foreshadow that everything else — trend, seasonality, controls — is scaffolding.
  - The driving tension to plant for the payoff: once we can *simulate* this process we can also ask the hard question — given only observed spend and sales, can we *recover* the curve? The answer (no, not from observational data alone) is the wound the chapter ends on and Part VI heals.
  - One sentence connecting forward: the synthetic series built here is the dataset Chapter 15 fits.

- [ ] **Step 3: Verify and commit.**

Run: `grep -c '\$' parts/05-mmm-modeling/01-mmm-dgp.qmd` → expect an even number. Confirm `pymcmarketing` is gone: `grep -c pymcmarketing parts/05-mmm-modeling/01-mmm-dgp.qmd` → expect `0`.

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): front matter, anchor scrub, and Motivation"
```

---

### Task 2: Theory rungs 1–2 (structural equation; adstock) + Proof P1

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (add `## Theory & Proofs` and first two rungs)

- [ ] **Step 1: Open `## Theory & Proofs` and write Rung 1 — the structural MMM equation.** Present the generative model:

$$
y_t = \underbrace{\alpha + \tau_t + \sigma_t}_{\text{baseline}} \;+\; \sum_{c=1}^{C} \beta_c\, g\!\big(a_\lambda(x_{c,t})\big) \;+\; \gamma^\top z_t \;+\; \varepsilon_t,
\qquad \varepsilon_t \sim \mathcal{N}(0,\sigma_\varepsilon^2).
$$

Name each term in one or two sentences: $\alpha$ intercept, $\tau_t$ slow trend, $\sigma_t$ seasonality (e.g. yearly), $z_t$ exogenous controls (price, holidays), $x_{c,t}$ channel-$c$ spend, $\beta_c$ its effect size, $g\circ a_\lambda$ the saturation-of-adstock transform. State plainly that the non-media terms are kept deliberately simple here — they matter mostly because they share variance with spend (the confounding seen later).

- [ ] **Step 2: Write Rung 2 — adstock / carryover.** Define geometric adstock two equivalent ways:

$$
a_\lambda(x)_t \;=\; x_t + \lambda\, a_\lambda(x)_{t-1} \;=\; \sum_{k=0}^{\infty} \lambda^k x_{t-k}, \qquad \lambda \in [0,1).
$$

Define the **normalized** form with weights $w_k = (1-\lambda)\lambda^k$ (so $\sum_k w_k = 1$, a weighted moving average that preserves the level of a constant input). Describe the **delayed-peak** variant of @jin2017 qualitatively in 2–3 sentences (carryover that rises to a peak a few weeks out before decaying, capturing ad campaigns that build before they fade) — no proof for the delayed form.

- [ ] **Step 3: Write Proof P1 — adstock geometric-series identities.** State as a proposition, then prove. Provide this skeleton in full prose:
  - **Normalization.** $\sum_{k=0}^\infty w_k = (1-\lambda)\sum_{k=0}^\infty \lambda^k = (1-\lambda)\cdot\frac{1}{1-\lambda}=1$ (geometric series, valid since $|\lambda|<1$). Hence for constant input $x_t\equiv x$, normalized adstock $=x$.
  - **Sustained-spend steady state.** Unnormalized, $x_t\equiv x \Rightarrow a = x\sum_{k\ge0}\lambda^k = x/(1-\lambda)$.
  - **Impulse response.** $x_0=1$, else $0 \Rightarrow a_k=\lambda^k$; its total mass is $1/(1-\lambda)$.
  - **Half-life.** Smallest lag with $\lambda^t \le \tfrac12$ solves $t = \ln(1/2)/\ln\lambda = \log_\lambda \tfrac12$.
  - **Mean lag.** $\sum_{k\ge0} k\,w_k = (1-\lambda)\sum_{k\ge0} k\lambda^k = (1-\lambda)\cdot\frac{\lambda}{(1-\lambda)^2} = \frac{\lambda}{1-\lambda}$, using $\sum_{k\ge0} k\lambda^k = \lambda/(1-\lambda)^2$.
  - Close by evaluating at $\lambda=0.5$: half-life $1$, mean lag $1$, steady state $2x$ — the WE1 numbers.

- [ ] **Step 4: Verify and commit.** Even `$`-count; no `\begin{align}` (`grep -c 'begin{align}' …` → 0); no `psmallmatrix`.

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): structural equation, adstock, and the carryover identities (P1)"
```

---

### Task 3: Theory rungs 3–4 (saturation; assembling + ordering) + Proofs P2, P3

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd`

- [ ] **Step 1: Write Rung 3 — saturation / shape.** Define the Hill curve exactly as Part IV uses it:

$$
g(x; K, n) = \frac{x^{n}}{K^{n} + x^{n}}, \qquad x \ge 0,\; K>0,\; n>0.
$$

State the reading of each parameter: $K$ is the half-saturation point ($g(K)=\tfrac12$), $n$ controls steepness/shape. One sentence noting this is the very curve Chapters 11–13 optimized over — defined here at last.

- [ ] **Step 2: Write Proof P2 — Hill shape theorem.** State: $g$ is strictly increasing and bounded in $(0,1)$; for $n\le1$ it is concave on $(0,\infty)$; for $n>1$ it has a single inflection at $x^\star=K\big(\tfrac{n-1}{n+1}\big)^{1/n}$, convex on $(0,x^\star)$ and concave on $(x^\star,\infty)$. Proof skeleton in full prose:
  - **Monotone/bounded.** $g'(x)=\dfrac{nK^n x^{n-1}}{(K^n+x^n)^2}>0$; limits $g(0^+)=0$, $g(\infty)=1$; $g(K)=\tfrac12$.
  - **Curvature.** Differentiate to
  $$
  g''(x) = \frac{nK^n x^{n-2}}{(K^n+x^n)^3}\Big[(n-1)K^n-(n+1)x^n\Big].
  $$
  The prefactor is positive for $x>0$, so $\operatorname{sign} g'' = \operatorname{sign}\big[(n-1)K^n-(n+1)x^n\big]$.
  - **Case $n\le1$:** bracket $<0$ for all $x>0$ ⇒ $g''<0$ ⇒ concave (diminishing returns from the first dollar).
  - **Case $n>1$:** bracket $=0$ at $x^n=K^n\frac{n-1}{n+1}$, i.e. $x^\star=K\big(\frac{n-1}{n+1}\big)^{1/n}$; positive below (convex toe), negative above (concave). Single inflection.
  - Evaluate at $K=3,n=2$: $x^\star=\sqrt3$, $g(\sqrt3)=0.25$, $g'(3)=1/6$, $g'(\sqrt3)=\sqrt3/8$.
  - One sentence: this is the concavity Part IV *assumed* when it invoked equal-marginal allocation; for $n>1$ the convex toe is exactly why Ch. 11–13 had to worry about non-convexity.

- [ ] **Step 3: Write Rung 4 — assembling the DGP and the order of operations.** Spend → adstock → saturation → scale by $\beta_c$ → sum → add baseline/controls/noise. State that @jin2017 apply **adstock first, then saturation**, and that this ordering is a real choice, motivating P3.

- [ ] **Step 4: Write Proof P3 — non-commutativity.** State: $g\circ a_\lambda \neq a_\lambda\circ g$ in general. Prove by explicit counterexample (full arithmetic shown): spend $(x_0,x_1)=(3,0)$, $\lambda=0.5$, $K=3$, $n=2$.
  - Adstock-then-saturate: $a=(3,\,1.5)$, output $\big(g(3),g(1.5)\big)=(0.5,\,0.2)$ since $g(1.5)=\frac{2.25}{9+2.25}=0.2$.
  - Saturate-then-adstock: $g$ first gives $(0.5,\,0)$, then adstock gives $(0.5,\,0.25)$.
  - Period 1: $0.2 \ne 0.25$. Conclude the transform ordering changes the model, so it is a modeling decision, not a convention. (Note they agree at the impulse period — the gap appears in carryover.)

- [ ] **Step 5: Verify and commit.** Even `$`-count; no `\begin{align}`/`psmallmatrix`.

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): Hill saturation (P2) and transform non-commutativity (P3)"
```

---

### Task 4: Theory rung 5 (identification wound) + Proof P4

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd`

- [ ] **Step 1: Write Rung 5 — the identification wound.** Motivate in prose: real spend is (a) often confined to a narrow operating range, (b) collinear across channels (budgets move together), and (c) confounded with the baseline (spend ramps up in high-demand seasons). Each makes recovering the curve from observed $(x,y)$ alone ill-posed. Frame the cleanest case formally for P4.

- [ ] **Step 2: Write Proof P4 — non-identifiability at a single operating point.** State: if observed spend is confined to a single point $x_0$ (or, more generally, the regressor lies on a lower-dimensional set than the parameter count), the parameters $(\beta,K,n)$ are not identified — the data constrain only the scalar $\beta\,g(x_0;K,n)$. Prove by construction:
  - Fix $x_0=3$. The contribution is $\beta\,g(3;K,n)$, a single number.
  - $(\beta,K,n)=(1,3,2)$ gives $1\cdot g(3;3,2)=0.5$.
  - $(\beta,K,n)=(2,\,3\sqrt3,\,2)$ gives $2\cdot g(3;3\sqrt3,2)=2\cdot\frac{9}{27+9}=2\cdot0.25=0.5$.
  - Two distinct parameter vectors, identical expected sales on the observed support ⇒ the likelihood is flat along a ridge; no amount of data at $x_0$ separates them. Generalize in a sentence: identifying the *shape* $(K,n)$ requires spend that *varies* across the curve.
  - **Close on the wound (do not resolve):** state explicitly that breaking this degeneracy needs exogenous variation in spend — randomized or quasi-experimental — which is exactly what Part VI (Chapters 18–21) supplies. One sentence; no IV/experiment/prior mechanics here.

- [ ] **Step 3: Verify and commit.**

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): the identification wound (P4), deferred to Part VI"
```

---

### Task 5: Worked Examples

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (add `## Worked Examples`)

- [ ] **Step 1: WE1 — geometric adstock by hand.** $\lambda=0.5$, impulse spend $[100,0,0,0]$. Show the recursion producing $[100,50,25,12.5]$; compute half-life $=1$ week, mean lag $=1$ week; show sustained spend $100$/week → steady state $200$. Walk the arithmetic so a reader can reproduce it without code.

- [ ] **Step 2: WE2 — the Hill curve at $K=3,n=2$.** Compute $g(3)=0.5$ (half-saturation), the inflection $x^\star=\sqrt3\approx1.732$ with height $g(\sqrt3)=0.25$, and the marginals $g'(3)=1/6\approx0.167$ and $g'(\sqrt3)=\sqrt3/8\approx0.217$. State explicitly: this is the exact curve whose optimum Part IV computed (Ch. 11's $\sqrt3$ inflection) — the loop closes.

- [ ] **Step 3: WE3 — non-identifiability made concrete.** Take a single observed operating point $x_0=3$ with measured contribution $0.5$. Show $(\beta,K)=(1,3)$ and $(\beta,K)=(2,3\sqrt3)$ both reproduce it (holding $n=2$). State the takeaway: one operating point pins one point on the curve, never its shape — a preview of why Chapter 20 needs experiments spread across spend levels.

- [ ] **Step 4: Verify and commit.** Even `$`-count.

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): worked examples (adstock, Hill, identification)"
```

---

### Task 6: Code Tie-in — simulate the DGP

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (add `## Code Tie-in` with one `{python}` cell)

- [ ] **Step 1: Write the section intro** (2–3 sentences): we now turn the equations into a generator and produce the synthetic weekly series Chapter 15 will fit. Pin a seed for reproducibility.

- [ ] **Step 2: Write the single `{python}` cell.** It must (a) define `geometric_adstock`, `hill`, and `hill_prime`; (b) assert every WE/proof anchor; (c) simulate a full 2-year DGP; (d) plot three panels. Use exactly this cell (verified headless):

```{python}
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(14)

# --- transforms ---
def geometric_adstock(x, lam):
    out = np.empty_like(x, dtype=float)
    prev = 0.0
    for t, xt in enumerate(x):
        prev = xt + lam * prev
        out[t] = prev
    return out

def hill(x, K, n):
    x = np.asarray(x, dtype=float)
    return x**n / (K**n + x**n)

def hill_prime(x, K, n):
    x = np.asarray(x, dtype=float)
    return n * K**n * x**(n - 1) / (K**n + x**n)**2

# --- P1 / WE1: adstock identities at lambda = 0.5 ---
imp = geometric_adstock(np.array([100, 0, 0, 0]), 0.5)
assert np.allclose(imp, [100, 50, 25, 12.5])
assert np.isclose(np.log(0.5) / np.log(0.5), 1.0)          # half-life
assert np.isclose(0.5 / (1 - 0.5), 1.0)                    # mean lag
assert np.isclose(100 / (1 - 0.5), 200.0)                  # steady state

# --- P2 / WE2: Hill at K=3, n=2 ---
xstar = 3 * (1/3)**0.5
assert np.isclose(xstar, np.sqrt(3))
assert np.isclose(hill(3, 3, 2), 0.5)
assert np.isclose(hill(xstar, 3, 2), 0.25)
assert np.isclose(hill_prime(3, 3, 2), 1/6)
assert np.isclose(hill_prime(xstar, 3, 2), np.sqrt(3) / 8)

# --- P3: non-commutativity, spend (3, 0), lambda = 0.5 ---
spend2 = np.array([3.0, 0.0])
A = hill(geometric_adstock(spend2, 0.5), 3, 2)            # adstock then saturate
B = geometric_adstock(hill(spend2, 3, 2), 0.5)           # saturate then adstock
assert np.isclose(A[1], 0.2) and np.isclose(B[1], 0.25)

# --- P4: non-identifiability at x0 = 3 ---
assert np.isclose(1 * hill(3, 3, 2), 0.5)
assert np.isclose(2 * hill(3, 3 * np.sqrt(3), 2), 0.5)

# --- full DGP: 2 years of weekly sales ---
T = 104
t = np.arange(T)
baseline = 10.0 + 0.02 * t + 2.0 * np.sin(2 * np.pi * t / 52)     # intercept + trend + yearly seasonality
x1 = np.clip(3.0 + 1.5 * np.sin(2 * np.pi * t / 52) + rng.normal(0, 0.5, T), 0, None)  # seasonal channel
x2 = np.clip(rng.gamma(2.0, 1.0, T), 0, None)                    # bursty channel
beta = (4.0, 3.0)
contrib1 = beta[0] * hill(geometric_adstock(x1, 0.5), 3, 2)
contrib2 = beta[1] * hill(geometric_adstock(x2, 0.3), 2, 1.5)
noise = rng.normal(0, 0.3, T)
sales = baseline + contrib1 + contrib2 + noise
assert sales.shape == (T,)
assert np.all(np.diff(hill(np.linspace(0.01, 20, 200), 3, 2)) > 0)   # saturation is monotone

# --- figures ---
fig, ax = plt.subplots(1, 3, figsize=(15, 4))

ax[0].stem(range(6), geometric_adstock(np.array([100, 0, 0, 0, 0, 0]), 0.5))
ax[0].set_title("Geometric adstock impulse (λ=0.5)")
ax[0].set_xlabel("weeks since spend"); ax[0].set_ylabel("carryover")

xs = np.linspace(0, 12, 200)
ax[1].plot(xs, hill(xs, 3, 2), label="K=3, n=2 (S-shaped)")
ax[1].plot(xs, hill(xs, 2, 1.5), label="K=2, n=1.5")
ax[1].axvline(np.sqrt(3), ls="--", color="grey")
ax[1].set_title("Hill saturation"); ax[1].set_xlabel("adstocked spend"); ax[1].legend()

ax[2].plot(t, sales, label="sales")
ax[2].plot(t, baseline, ls="--", label="baseline")
ax[2].set_title("Simulated DGP (104 weeks)"); ax[2].set_xlabel("week"); ax[2].legend()

plt.tight_layout()
plt.show()
```

- [ ] **Step 3: Run it headless.** Extract the cell to `/tmp/ch14_cell.py` and run `MPLBACKEND=Agg python3 /tmp/ch14_cell.py`. Expected: no output, no assertion error, exit 0.

- [ ] **Step 4: Commit.**

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): code tie-in simulating the DGP into Ch15's dataset"
```

---

### Task 7: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (add `## Exercises` with the four subsections)

> **Controller note:** per the established workflow, the controller writes Exercises directly (not via subagent). Self-contained — no inline solution links. Provide 2–3 per tier. Use only concepts defined in this chapter.

- [ ] **Step 1: Write `### C — Conceptual / Reading Comprehension`.** e.g.: why does normalized adstock preserve a constant input's level while unnormalized adstock inflates it? In words, what does $n>1$ buy you (and cost you) versus $n\le1$? Why does a single operating point identify only one number?

- [ ] **Step 2: Write `### B — By Hand`.** e.g.: compute the $\lambda=0.7$ adstock half-life and mean lag; find the Hill inflection for $K=5,n=3$; given spend $(2,0)$, $\lambda=0.5$, $K=2,n=2$, compute both transform orderings and show they differ.

- [ ] **Step 3: Write `### P — Prove It`.** e.g.: prove $g(K)=\tfrac12$ for all $n$; prove the mean-lag identity $\lambda/(1-\lambda)$ from $\sum k\lambda^k$; prove that for $n\le1$ the Hill curve is concave on $(0,\infty)$ (no inflection).

- [ ] **Step 4: Write `### A — Applied / Code`.** e.g.: extend the simulation to three channels and plot contributions; empirically locate the Hill inflection by finite differences and match $x^\star$; demonstrate the identification ridge by fitting $(\beta,K)$ to single-operating-point data and showing the flat objective.

- [ ] **Step 5: Commit.**

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): exercises (C/B/P/A)"
```

---

### Task 8: Summary (auto-included)

**Files:**
- Modify: `parts/05-mmm-modeling/01-mmm-dgp.qmd` (add `## Summary`)

> **Controller note:** controller writes the Summary directly. **Key concepts** and **Key identities** are BOTH bulleted lists (per the formatting-consistency rule — never a run-on paragraph). Key identities are inline math only.

- [ ] **Step 1: Write `## Summary`** with a 1–2 sentence lead, then:
  - `**Key concepts**` — a bulleted list (one bullet each): the structural MMM equation; adstock/carryover as a geometric filter; Hill saturation and the $n\le1$ vs $n>1$ shape split; transform ordering is a modeling choice; the identification wound (observational spend can't pin the curve) and its deferral to Part VI; the synthetic DGP is Chapter 15's data.
  - `**Key identities**` — a bulleted list (short bold label + inline math each): adstock recursion $a_t=x_t+\lambda a_{t-1}$ and steady state $x/(1-\lambda)$; half-life $\ln 0.5/\ln\lambda$; mean lag $\lambda/(1-\lambda)$; Hill $g(x;K,n)=x^n/(K^n+x^n)$ with $g(K)=\tfrac12$; inflection $x^\star=K\big(\tfrac{n-1}{n+1}\big)^{1/n}$ for $n>1$; curvature sign $\propto (n-1)K^n-(n+1)x^n$; non-identification ridge $\beta\,g(x_0;K,n)=\text{const}$.
  - One closing sentence handing off to Chapter 15 (fitting this DGP).

- [ ] **Step 2: Verify** the Summary uses bulleted lists (mirror `parts/04-optimization/03-slsqp.qmd`'s now-bulleted Summary) and even `$`-count. Commit.

```bash
git add parts/05-mmm-modeling/01-mmm-dgp.qmd
git commit -m "feat(ch14): summary"
```

---

### Task 9: Appendix solutions block

**Files:**
- Modify: `appendix/solutions.qmd` (append a gated Chapter 14 block before the final closing `:::`)

- [ ] **Step 1: Find the insertion point.** Open `appendix/solutions.qmd`; locate the final closing `:::` of the last chapter's `content-visible` block. Insert the new block immediately before the file's terminal `:::` so it stays inside the solutions gating, matching the Chapter 13 block's structure.

- [ ] **Step 2: Write the block.** Header `## Chapter 14 — The MMM Data-Generating Process`, wrapped in `::: {.content-visible when-meta="show-solutions"}` … `:::`. Provide full solutions for every C/B/P/A exercise authored in Task 7. For the B/P numeric ones, show worked arithmetic; for A, give a short verified `{python}` cell (run it headless before committing). Pin: $\lambda=0.7$ half-life $=\ln0.5/\ln0.7\approx1.9434$, mean lag $0.7/0.3\approx2.3333$; $K=5,n=3$ inflection $5(2/4)^{1/3}=5\cdot0.7937\approx3.9685$.

- [ ] **Step 3: Verify.** Confirm the gating opens and closes exactly once for the new block (balanced `:::`); if the A solution has a code cell, run it headless (exit 0). Commit.

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch14): appendix solutions"
```

---

### Task 10: Bundled housekeeping — Chapter 11 cross-reference fixes

**Files:**
- Modify: `parts/04-optimization/01-convexity.qmd` (two edits; leave the ridge reference alone)

- [ ] **Step 1: Fix line ~19.** Change "the diminishing-returns saturation curve introduced in Chapter 5" → "the diminishing-returns saturation curve introduced in Chapter 14". (Confirm context: it is the budget-allocation paragraph about concave $r_i$.)

- [ ] **Step 2: Fix line ~69.** Change "The canonical example is the **Hill saturation curve** introduced in Chapter 5," → "… introduced in Chapter 14,". (Confirm context: the concavity definition paragraph.)

- [ ] **Step 3: Confirm the ridge reference is untouched.** `grep -n "Chapter 5" parts/04-optimization/01-convexity.qmd` should still show the line ~200 reference (the ridge/MAP Gaussian model — genuinely Chapter 5). Verify only the two Hill references changed.

- [ ] **Step 4: Commit.**

```bash
git add parts/04-optimization/01-convexity.qmd
git commit -m "fix(ch11): re-point Hill cross-references to Chapter 14 (its real source)"
```

---

### Task 11: Final review, headless re-run, PR

**Files:**
- All of the above.

- [ ] **Step 1: Structure/KaTeX lint.** Confirm in `parts/05-mmm-modeling/01-mmm-dgp.qmd`: H1 present; the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises`, `## Summary`); even `$`-count for the file; zero `begin{align}`; zero `psmallmatrix`; zero `pymcmarketing`. Run:

```bash
f=parts/05-mmm-modeling/01-mmm-dgp.qmd
grep -c '\$' "$f"; grep -c 'begin{align}' "$f"; grep -c 'psmallmatrix' "$f"; grep -c pymcmarketing "$f"
grep -nE '^#( |$)|^## ' "$f"
```
Expected: even number; `0`; `0`; `0`; H1 then the six `##` headings in order.

- [ ] **Step 2: Re-run the code cell headless.** Re-extract and run under `MPLBACKEND=Agg python3` — exit 0, all asserts pass.

- [ ] **Step 3: Push and open PR.**

```bash
git push -u origin worktree-ch14-mmm-dgp
gh pr create --base main --head worktree-ch14-mmm-dgp \
  --title "feat: Chapter 14 — The MMM Data-Generating Process" \
  --body "Part V opener. Defines the MMM DGP (adstock + Hill saturation + light baseline/trend/seasonality); proves the adstock geometric-series identities, the Hill shape theorem, transform non-commutativity, and a compact non-identifiability result (deferred to Part VI); simulates the DGP into Ch15's dataset. Bundles: scrub the pymc-marketing anchor; fix two Ch11 Hill cross-references. CI quarto render (HTML+PDF) is the gate."
```

- [ ] **Step 4: Watch CI to green.** Watch BOTH the HTML and PDF render jobs of `render.yml` on the PR to a green conclusion before reporting done. The `deploy` job is skipped on PR branches (runs on merge to main) — do not wait on it. If the PDF stage fails on a LaTeX/KaTeX gap (e.g. an environment undefined in LuaLaTeX), fix and re-push; HTML-side lint cannot catch PDF-only gaps. Then hand off to the user to merge.

---

## Self-Review (controller, before dispatching Task 1)

**Spec coverage:** Motivation (T1) ✓; structural equation + non-media scaffolding (T2) ✓; adstock + P1 (T2) ✓; Hill + P2 (T3) ✓; ordering + P3 (T3) ✓; identification wound + P4, cure deferred (T4) ✓; WE1–3 with the Part IV loop-closure (T5) ✓; full-DGP code as Ch15 hand-off (T6) ✓; C/B/P/A exercises (T7) ✓; bulleted Summary (T8) ✓; gated appendix solutions (T9) ✓; anchor scrub (T1) ✓; two Ch11 cross-ref fixes, ridge left alone (T10) ✓; CI HTML+PDF gate (T11) ✓.

**Placeholder scan:** no TBD/TODO; all anchor numbers NumPy-verified and inlined; the code cell is complete and was run headless.

**Consistency:** notation fixed once ($K,n,\lambda,T,\beta,g,a_\lambda$) and used identically across tasks; `geometric_adstock`/`hill`/`hill_prime` defined once in T6 and reused; proof statements match the worked-example numbers.
