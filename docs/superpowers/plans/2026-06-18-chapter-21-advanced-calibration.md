# Chapter 21 — Advanced Calibration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author Chapter 21 (Advanced Calibration) — Part VI's keystone — folding Chapter 20's taxonomy-tagged experimental estimands into the MMM posterior as a calibration likelihood factor (Bayesian evidence synthesis), proving the posterior-compression and secant–tangent results, and numerically closing the Chapter 18 EVPI gap.

**Architecture:** Replace the stub body of `parts/06-causal-calibration/03-advanced-calibration.qmd` with the six-heading template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary). Append a Chapter 21 solutions block to `appendix/solutions.qmd`. Add one BibTeX entry. Verify via headless code run + CI `quarto render`, ship as a PR.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), `numpy`/`scipy`/`matplotlib` only, `references.bib`.

---

## Conventions (enforced every task)

- **HOUSE RULES (critical):** NEVER name PyMC-Marketing or any MMM/PPL/sampler library (no PyMC, Stan, Orbit, numpyro, pyro, lightweight_mmm, causalimpact). Calibration/lift-testing is a **method**, not a product. `numpy`/`scipy`/`matplotlib` only.
- **KaTeX:** use `aligned` inside `$$ … $$` only — never bare `\begin{align}`. `$$` delimiters on their own lines. **Even** `$$` count per file. NEVER `\begin{psmallmatrix}`.
- Proofs end with `$\blacksquare$`.
- "Key identities" in the Summary must be a **bulleted list**, not a run-on paragraph.
- Keep H1 `# Advanced Calibration` and the anchors line `*Canonical anchors: Bayesian evidence synthesis; informative priors (Gelman et al., BDA3).*`. Remove the stub callout.
- Citation keys: `@gelman2013` (exists), `@ibrahim2000` (added in T9).
- Commit identity `jlh530i` / `jlh530i@gmail.com`. Branch → PR → user merges; never self-merge; never commit to main.

## NumPy-verified anchors (all confirmed before planning)

- **P1 precision addition:** $\Lambda_{\text{obs}}=1$, $1/s^2=3 \Rightarrow \Lambda_{\text{post}}=4$; variance $1\to0.25$, std $1\to0.5$ (halves, 4× variance drop).
- **WE1 secant calibration:** prior precision $4$ ($\theta\sim\mathcal N(2,0.5^2)$), lift precision $25$ ($s=0.2$); posterior precision $29$, variance $0.03448$, std $0.18570$, mean $2.0$.
- **P2 / WE2 secant–tangent:** $S=2\sqrt x$. $x=4,\delta=5$: secant $0.4$, tangent $S'(4)=0.5$, gap $-0.1$; $S''(\xi)=-0.04 \Rightarrow \xi\approx5.386$. Paired about $x=9,\delta=4$: up $0.30278$, down $0.38197$, average $0.34237$ vs $S'(9)=0.33333$.
- **WE3 / EVPI heal:** reusing the Chapter 18 Monte-Carlo (abar $(2,1)$, $B=9$, contrast $u=(1,-1)/\sqrt2$, well-identified $v$ std $0.08$, seed `default_rng(21)`): ridge std $0.45 \Rightarrow$ EVPI $0.1197$; ridge std $0.225 \Rightarrow$ EVPI $0.0315$ (~3.8× reduction).

## Critical files

- `parts/06-causal-calibration/03-advanced-calibration.qmd` — the chapter (currently a stub; replace body).
- `appendix/solutions.qmd` — append Chapter 21 block before the final `:::` (currently last line, file is 2880 lines; one big `content-visible` gated div).
- `references.bib` — add `@ibrahim2000`.
- Voice/structure exemplars (read-only): `parts/06-causal-calibration/02-quasi-experimental.qmd` (Ch20, immediate predecessor; owns the secant/tangent/validation taxonomy), `parts/05-mmm-modeling/04-budget-optimization.qmd` (Ch18; the EVPI/ridge wound and the exact Monte-Carlo to reuse in WE3).
- Authoritative design: `docs/superpowers/specs/2026-06-18-chapter-21-advanced-calibration-design.md`.

---

### Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

- [ ] **Step 1: Strip the stub.** Keep lines 1–3 (H1 + blank + anchors line). Delete the stub callout (lines 5–7) and the placeholder italic lines under each heading. Keep the six template headings (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises` with its four `### C/B/P/A` subheads, `## Summary`) in order.

- [ ] **Step 2: Write Motivation (3–4 paragraphs).** The driving question: *"You ran the experiment and got a number. How do you put it into the model — and how much does it help?"* Cover, in prose matching Ch20's register:
  - Chapter 19 proved observational MMM cannot identify the response curve under unobserved confounding (the structural wound); Chapter 20 produced *interventional* estimands — taxonomy-tagged functionals of the curve (secant / tangent / validation).
  - The naive move — hand-tune a "prior on ROAS" until the model agrees with the experiment — double-counts data and is undisciplined. The disciplined move is **Bayesian evidence synthesis**: the experiment enters as a **likelihood factor** on a *functional* of $\theta$, composing cleanly with the observational likelihood.
  - Preview the payoff: this factor compresses the posterior in exactly the direction the observational data left flat — healing the Chapter 16 ridge and shrinking the Chapter 18 EVPI gap (which Chapter 18 could only *name*). Close the loop: **intervene → identify → calibrate → re-optimize.**
  - Name the anchor: $S(x;\theta)=\theta\sqrt x$ (truth $\theta=2$), plus the two-parameter $(a_1,a_2)$ ridge from Chapter 18's WE3.
  - Do NOT name any library. Calibration is a method.

- [ ] **Step 3: Verify** even `$$` count in the file so far (`grep -c '\$\$'` is even), headings present and ordered. Commit: `feat(ch21): motivation`.

---

### Task 2: Theory & Proofs — Rungs 1–3 (the two keystones)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

Write three numbered rungs under `## Theory & Proofs`. Match Ch20's rung prose density and proof rigor.

- [ ] **Step 1: Rung 1 — The calibration likelihood factor (taxonomy-keyed).** Bayesian evidence synthesis. State
$$
p(\theta \mid D, \hat g) \;\propto\; p(\theta)\,L_{\text{obs}}(\theta\mid D)\,L_{\text{exp}}(\hat g \mid \theta),
$$
with $L_{\text{obs}}$ = the MMM fit, $L_{\text{exp}}$ = the calibration factor on a functional $g(\theta)$ of the curve. Key the factor by Chapter 20's taxonomy:
  - **Secant:** $\hat\Delta \sim \mathcal N\!\big(S(x{+}\delta;\theta) - S(x;\theta),\, s^2\big)$ (note a positivity-respecting Gamma as the practical variant, one sentence).
  - **Tangent:** match the marginal, $\widehat{S'}\cdot\delta \sim \mathcal N\!\big(S'(x;\theta)\,\delta,\, s^2\big)$.
  - **Validation:** *no* likelihood term — a posterior-predictive check / conflict alarm only.
  State plainly: the experiment is information about a *functional* of $\theta$, entering through the likelihood; cleaner than re-expressing it as a prior, and never double-counts the same data.

- [ ] **Step 2: Rung 2 — Proof P1 (KEYSTONE): a calibration factor compresses the posterior in the constrained direction.** Linear-Gaussian setup: observational posterior $\theta \sim \mathcal N(m_{\text{obs}}, \Lambda_{\text{obs}}^{-1})$; calibration of a linear functional $c^\top\theta$ with $\hat g \sim \mathcal N(c^\top\theta, s^2)$. **Theorem:**
$$
\Lambda_{\text{post}} = \Lambda_{\text{obs}} + \tfrac{1}{s^2}\,c\,c^\top.
$$
  **Proof:** complete the square in the Gaussian log-posterior; the calibration factor contributes the rank-one $\tfrac{1}{s^2}cc^\top$ to the precision; by the Sherman–Morrison / scalar projection view the variance along $c$ drops from $1/\Lambda_{\text{obs},c}$ to $1/(\Lambda_{\text{obs},c}+1/s^2)$ — a strict decrease, largest where $\Lambda_{\text{obs}}$ is smallest (the ridge). End `$\blacksquare$`. **MMM reading:** the Chapter 16 Fisher-information ridge is the near-flat eigendirection of $\Lambda_{\text{obs}}$; a Chapter 20 experiment aimed at it supplies the missing precision, and the wide Chapter 18 allocation posterior narrows. **Anchor (state exactly):** $\Lambda_{\text{obs}}=1$ along $u=(1,-1)/\sqrt2$, $1/s^2=3 \Rightarrow \Lambda_{\text{post}}=4$, variance $1\to0.25$, std halves, variance drops $4\times$.

- [ ] **Step 3: Rung 3 — Proof P2 (KEYSTONE): the secant–tangent matching relationship.** A finite lift test measures a **secant**; the optimizer wants a **tangent**. **Theorem (Taylor with remainder),** $S\in C^2$:
$$
\frac{S(x{+}\delta)-S(x)}{\delta} = S'(x) + \frac{\delta}{2}\,S''(\xi), \qquad \xi\in(x,x{+}\delta).
$$
  Under concavity ($S''<0$) the secant **underestimates** the lower-endpoint tangent by $\tfrac{\delta}{2}|S''(\xi)|$ — an $O(\delta)$ bias. **Corollary (paired bracketing):** averaging a matched up- and down-test about $x$ cancels the $O(\delta)$ term:
$$
\tfrac12\!\left(\frac{S(x{+}\delta)-S(x)}{\delta}+\frac{S(x)-S(x{-}\delta)}{\delta}\right)=S'(x)+O(\delta^2).
$$
  Prove both (FTC + MVT for the first; second-order Taylor expansions of the two secants for the cancellation). End `$\blacksquare$`. **Anchors:** $S=2\sqrt x$, $x=4,\delta=5$: secant $0.4$, tangent $0.5$, gap $-0.1$, $S''(\xi)=-0.04$ at $\xi\approx5.39$; paired about $x=9,\delta=4$: up $0.303$, down $0.382$, average $0.342$ vs $S'(9)=0.333$.

- [ ] **Step 4: Verify** even `$$` count, no bare `\begin{align}`, no `psmallmatrix`, both proofs end `$\blacksquare$`. Commit: `feat(ch21): theory rungs 1-3 (calibration factor + two keystone proofs)`.

---

### Task 3: Theory & Proofs — Rungs 4–5 (identification + power-prior handoff)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

- [ ] **Step 1: Rung 4 — Identification: test-design ↔ curve shape.** A saturating curve has shape parameters (scale $\theta$; saturation $\kappa,\alpha$ in the general Hill form). **One** secant or tangent pins **one** functional — it under-identifies the shape (many curves pass through a single measured lift). State the identification principle: to pin the *shape* you need a **spread of constraints** — lifts at several spend levels, or a wide secant plus a tangent — so the experiments jointly determine $(\theta,\kappa,\alpha)$. This is the **active-learning hook**: the Chapter 18 optimizer flags the high-leverage spend region; the next experiment is aimed there (closing the loop). Geometry: each constraint is a surface in parameter space; identification is their intersection collapsing to a point.

- [ ] **Step 2: Rung 5 — Power-prior down-weighting (brief) + handoff to Chapter 22.** Not every experiment deserves full weight. Introduce the **power-prior** factor $L_{\text{exp}}(\hat g\mid\theta)^{w}$ with $w\in[0,1]$ = (design credibility) × (temporal decay $\rho^{\,t_{\text{now}}-t_k}$); $w=1$ full trust, $w=0$ a validation-tier check that never enters the likelihood [@ibrahim2000]. State that a *collection* of such weighted factors, accumulated and versioned, is the **prior store** — developed as a data product in Chapter 22 (schema = the secant/tangent/validation taxonomy; sequential update; hierarchical pooling of conflicting studies). Keep brief: the mechanism here, the system there. Close on the full loop — **intervene (Ch. 19–20) → identify the functional → calibrate (here) → re-optimize (Ch. 18)** — and the Chapter 18 EVPI gap shrinks.

- [ ] **Step 3: Verify** even `$$` count, citation `@ibrahim2000` used (resolved in T9). Commit: `feat(ch21): theory rungs 4-5 (identification + power-prior handoff)`.

---

### Task 4: Worked Examples (WE1–WE3)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

Write three worked examples under `## Worked Examples`, each carrying its result into arithmetic. Match Ch20 WE prose.

- [ ] **Step 1: WE1 — Secant calibration tightens the curve.** Scale model $S(x;\theta)=\theta\sqrt x$ (truth $\theta=2$). Observational posterior $\theta\sim\mathcal N(2,0.5^2)$ (precision $4$). Chapter 20's geo lift measured $\hat\Delta=2$ over $[4,9]$; since $S(9;\theta)-S(4;\theta)=\theta(\sqrt9-\sqrt4)=\theta$, the secant factor is $\hat\Delta\sim\mathcal N(\theta,s^2)$ with $s=0.2$ (precision $25$). Posterior precision $4+25=29$, variance $0.0345$, std $0.186$, mean $(2/0.25+2/0.04)/29=2.0$. The lift pins the scale: std $0.5\to0.186$ (a $2.7\times$ reduction). State the marginal returns the optimizer differentiates tighten with it.

- [ ] **Step 2: WE2 — The secant–tangent gap, and paired bracketing.** On $S=2\sqrt x$: a lift moving spend $4\to9$ returns secant $0.4$, but the marginal return at $x=4$ is $S'(4)=0.5$ — calibrating $S'(4)$ with this secant biases it low by $0.1$, the $\tfrac{\delta}{2}S''(\xi)$ term ($\xi\approx5.39$). Show the paired fix about $x=9,\delta=4$: up-secant $0.303$, down-secant $0.382$, average $0.342$ vs $S'(9)=0.333$ — second-order accurate. Lesson: match the estimand class to the parameter you calibrate, or pay a known-sign $O(\delta)$ bias.

- [ ] **Step 3: WE3 — Healing the ridge (the wound closes).** Reuse Chapter 18's ridge posterior on $(a_1,a_2)$ centered $(2,1)$: contrast $u=(1,-1)/\sqrt2$ has std $0.45$ (the ridge), giving the Chapter 18 EVPI $\approx0.12$. A geo experiment supplies a calibration factor on the contrast (measuring the attribution split the observational data could not), **halving** the ridge std to $0.225$. Propagate through the Chapter 18 optimizer: the allocation posterior narrows and the **EVPI gap shrinks from $\approx0.12$ to $\approx0.03$** (~4×, tracking variance). State plainly: this is the Chapter 18 wound closing — calibration bought back the identification the observational series lacked.

- [ ] **Step 4: Verify** even `$$` count; the three numbers match the anchors. Commit: `feat(ch21): worked examples WE1-WE3`.

---

### Task 5: Code Tie-in (single runnable cell)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

- [ ] **Step 1: Write a single ```{python}``` cell** under `## Code Tie-in`, prefaced by a short prose paragraph (Ch20-style) describing what it does. `numpy` + `matplotlib` only; seed `np.random.default_rng(21)`; figures end `plt.show()`. Four numbered blocks, every numeric claim `assert`ed:
  1. **Precision addition (P1):** build $\Lambda_{\text{obs}}$ with a flat ridge direction $u=(1,-1)/\sqrt2$ (eigenvalue $1$) and a stiff direction $v$; add rank-one $\tfrac1{s^2}cc^\top$ with $c=u$, $1/s^2=3$; print ridge variance before ($1.0$) / after ($0.25$); `assert` the $4\times$ drop. Figure: posterior covariance ellipse before (elongated) vs after (compressed).
  2. **Secant calibration (WE1):** conjugate Gaussian update of $\theta$ (precision $4+25=29$); `assert` posterior mean $\approx2$ and std $\approx0.186$.
  3. **Secant–tangent (WE2):** compute secant, tangent $S'(9)$, gap, and paired-test average; `assert` the paired average ($0.342$) is closer to $S'(9)=0.333$ than either single secant ($0.303$, $0.382$).
  4. **EVPI heals (WE3):** the Chapter 18 Monte-Carlo EVPI at ridge std $0.45$ then $0.225$ (abar $(2,1)$, $B=9$, $v$-std $0.08$, keep positive draws); `assert` it drops from $\approx0.12$ to $\approx0.03$ and that the second is at least ~3× smaller. Figure: EVPI vs ridge std with before/after points marked.

- [ ] **Step 2: Extract and run headless.** Copy the cell body to a temp `.py`, run `MPLBACKEND=Agg python3 /tmp/ch21_code.py`. Expected: all asserts pass; prints show variance $1.0\to0.25$, WE1 std $\approx0.186$, paired avg $\approx0.342$, EVPI $\approx0.12\to\approx0.03$. Fix until clean.

- [ ] **Step 3: Verify** even `$$` count; no banned libraries (`grep -iE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight'` returns nothing). Commit: `feat(ch21): code tie-in (precision addition, secant calibration, EVPI heal)`.

---

### Task 6: Exercises (C / B / P / A — self-contained)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

- [ ] **Step 1: Write the four exercise blocks** under the existing `### C/B/P/A` subheads. No inline solution links (solutions live in the appendix, T8).
  - **C (Conceptual):** (C1) Why does an experiment enter as a likelihood factor rather than a prior on ROAS, and what does the likelihood-factor view avoid? (C2) Why does a calibration constraint help most in the ridge direction? (C3) Why does one lift under-identify the curve shape?
  - **B (By hand):** (B1) On $S=2\sqrt x$ compute the secant–tangent gap for $x=4,\delta=5$. (B2) Compute the paired-test average about $x=9,\delta=4$ and compare to $S'(9)$. (B3) Do the conjugate Gaussian update of WE1 (precision $4+25=29$; report posterior std and mean).
  - **P (Prove it):** (P1) Prove the posterior-precision addition $\Lambda_{\text{post}}=\Lambda_{\text{obs}}+\tfrac1{s^2}cc^\top$ and that the variance along $c$ strictly decreases. (P2) Prove the secant–tangent Taylor relation and the second-order accuracy of the paired average.
  - **A (Applied):** (A1) Extend the Code Tie-in — vary the calibration precision $1/s^2$ and plot residual ridge variance (and EVPI) against it. (A2) Add a *second* calibration constraint at a different spend level and show the curve *shape* (two parameters) becomes identified where one constraint did not.

- [ ] **Step 2: Verify** even `$$` count. Commit: `feat(ch21): exercises C/B/P/A`.

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/06-causal-calibration/03-advanced-calibration.qmd`

- [ ] **Step 1: Write `## Summary`** with a one-paragraph wrap, then **bulleted** "Key concepts" and **bulleted** "Key identities" (inline math). Identities (each a bullet):
  - evidence-synthesis posterior $p(\theta\mid D,\hat g)\propto p(\theta)\,L_{\text{obs}}\,L_{\text{exp}}$
  - secant factor $\hat\Delta\sim\mathcal N(S(x{+}\delta;\theta)-S(x;\theta),s^2)$
  - precision addition $\Lambda_{\text{post}}=\Lambda_{\text{obs}}+\tfrac1{s^2}cc^\top$
  - secant–tangent $\frac{S(x+\delta)-S(x)}{\delta}=S'(x)+\frac{\delta}{2}S''(\xi)$
  - power-prior weight $w\in[0,1]$
  Close tying back to Chapter 18 (the EVPI gap shrinks) and forward to Chapter 22 (the prior store).

- [ ] **Step 2: Verify** Key identities is a bulleted list (not a paragraph), even `$$` count. Commit: `feat(ch21): summary`.

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd`

- [ ] **Step 1: Append `## Chapter 21 — Advanced Calibration`** immediately before the final `:::` (currently the last line; the file is one big `content-visible` gated div). Provide full C/B/P/A solutions:
  - **C1–C3:** likelihood-factor avoids double-counting + composes with the observational likelihood (not a hand-tuned prior); calibration helps most where $\Lambda_{\text{obs}}$ is smallest (the ridge) because the rank-one update's relative variance reduction is largest there; one lift fixes one functional, leaving the multi-parameter shape under-determined.
  - **B1:** secant $0.4$, tangent $S'(4)=0.5$, gap $-0.1 = \tfrac{\delta}{2}S''(\xi)$, $\xi\approx5.39$. **B2:** up $0.30278$, down $0.38197$, average $0.34237$ vs $S'(9)=0.33333$ (second-order accurate). **B3:** precision $4+25=29$, variance $0.03448$, std $0.18570$, mean $(8+50)/29=2.0$.
  - **P1:** complete-the-square proof of $\Lambda_{\text{post}}=\Lambda_{\text{obs}}+\tfrac1{s^2}cc^\top$ + scalar projection / Sherman–Morrison strict variance decrease along $c$. End `$\blacksquare$`. **P2:** FTC+MVT for $\frac{S(x+\delta)-S(x)}{\delta}=S'(\xi)$ form and Taylor-remainder for $S'(x)+\tfrac{\delta}{2}S''(\xi)$; second-order Taylor of both paired secants to show the $O(\delta)$ terms cancel. End `$\blacksquare$`.
  - **A1, A2:** describe the expected plots (residual ridge variance $\propto 1/(\Lambda_{\text{obs}}+1/s^2)$ falling as $1/s^2$ grows; EVPI falling with it; a second constraint at another spend level intersects the first surface so the two-parameter shape is pinned where one constraint left a curve of solutions).

- [ ] **Step 2: Verify** even `$$` count in `appendix/solutions.qmd`; the block sits before the final `:::`; proofs end `$\blacksquare$`. Commit: `feat(ch21): appendix solutions (C/B/P/A)`.

---

### Task 9: Bibliography + final review + PR

**Files:**
- Modify: `references.bib`
- Verify: `parts/06-causal-calibration/03-advanced-calibration.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: Add the BibTeX entry** to `references.bib`:

```bibtex
@article{ibrahim2000,
  author  = {Ibrahim, Joseph G. and Chen, Ming-Hui},
  title   = {Power Prior Distributions for Regression Models},
  journal = {Statistical Science},
  year    = {2000},
  volume  = {15},
  number  = {1},
  pages   = {46--60}
}
```

- [ ] **Step 2: Lint pass on both files.** Confirm: even `$$` count in the chapter and in the appendix; no bare `\begin{align}`; no `\begin{psmallmatrix}`; six template headings present and in order; H1 `# Advanced Calibration` and anchors line intact; stub callout gone; all proofs end `$\blacksquare$`; Key identities bulleted; only `@gelman2013`/`@ibrahim2000` cited and both resolve in `references.bib`; no banned library names anywhere (`grep -iE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight|marketing'`).

- [ ] **Step 3: Re-run the extracted code cell headless** one final time (`MPLBACKEND=Agg python3`); all asserts pass.

- [ ] **Step 4: Verify local CI hooks active** (`git config --get core.hooksPath` points at the tracked hooks dir; if not, run the repo installer). Then commit: `feat(ch21): power-prior reference (Ibrahim & Chen 2000)`.

- [ ] **Step 5: Push the branch and open a PR** against `main` (title `feat: Chapter 21 — Advanced Calibration`; body summarizing the two keystones + the EVPI heal). Then start the background merge-watcher:
  `until gh pr view <N> --json state -q .state | grep -qE 'MERGED|CLOSED'; do sleep 120; done` (run via Bash `run_in_background: true`).

- [ ] **Step 6: Watch CI `quarto render` (HTML + PDF) to a green conclusion.** This is the real gate. If it fails, read the log, fix, push, re-watch. Report the green PR to the user; the user merges.

---

## Self-Review (done at plan time)

- **Spec coverage:** Rungs 1–5 → Tasks 2–3; WE1–3 → Task 4; code 4 blocks → Task 5; exercises C/B/P/A → Task 6; summary → Task 7; appendix → Task 8; bib `@ibrahim2000` → Task 9. All spec sections mapped.
- **Anchor consistency:** P1 (1→0.25, std halves), WE1 (29 / 0.186 / 2.0), P2/WE2 (0.4 vs 0.5; 0.342 vs 0.333), WE3/EVPI (0.1197→0.0315) are identical across Theory, Worked Examples, Code, and Appendix tasks — all NumPy-verified above.
- **No placeholders:** every task carries the actual math/numbers/commands.
- **House rules:** library ban, KaTeX rules, even `$$`, bulleted Key identities, `$\blacksquare$`, branch→PR→user-merges all enforced per task.
