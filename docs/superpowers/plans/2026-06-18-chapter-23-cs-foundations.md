# Chapter 23 — CS Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author Chapter 23 — Part VII's opening chapter — grounding the book's computation in floating-point arithmetic and algorithmic cost, with the keystone $\kappa(X^\top X)=\kappa(X)^2$ tying the Chapter 16 identifiability ridge to numerical ill-conditioning.

**Architecture:** Replace the stub body of `parts/07-swe-cs/01-cs-foundations.qmd` with the six-heading template. Append a Chapter 23 solutions block to `appendix/solutions.qmd`. No bib additions. Verify via headless code run + CI `quarto render`, ship as a PR.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), `numpy`/`scipy`/`matplotlib` only, `references.bib`.

---

## STYLING CONVENTIONS — match Part IV (`parts/04-optimization/01-convexity.qmd`) — USER REQUIREMENT

The user explicitly asked that formula styling and especially **proof spacing** match Part IV. Follow these exactly:

- **Rungs are `### Rung N — Title` H3 headings** (not inline `**Rung N — ...**` bold). This is the Part IV convention and gives natural vertical separation. (Chapters 19–22 used inline bold; Chapter 23 follows Part IV per the user.)
- **Theorem / Lemma / Proposition** lead a paragraph as `**Theorem (name).**` The statement may end with a colon introducing a centered display.
- **Proofs are AIRY, step-by-step.** Open with `**Proof.**` then a short framing line. Each logical step is its **own short paragraph** that ends in a colon or transition, followed by a **blank line**, a `$$` display **with each `$$` on its own line**, another **blank line**, then the next step's paragraph. NEVER cram a multi-step derivation into one dense paragraph. Model the rhythm on Part IV Rung 3/Rung 4 proofs.
- **Multi-line math** uses `\begin{aligned} … \end{aligned}` INSIDE `$$ … $$` (e.g. stacked gradient/derivative lines, paired inequalities). NEVER a bare `\begin{align}`.
- Use `\quad \text{for all } …` style annotations inside displays and `\bigl( … \bigr)` sizing where it aids readability, as Part IV does.
- Every display has a blank line before the opening `$$` and after the closing `$$`.
- Proofs end with `$\blacksquare$`.

## Other conventions (enforced every task)

- **HOUSE RULES (critical):** NEVER name PyMC-Marketing or any MMM/PPL/sampler library (no PyMC, Stan, Orbit, numpyro, pyro, lightweight_mmm, causalimpact). `numpy`/`scipy`/`matplotlib` only.
- **KaTeX:** `aligned` inside `$$ … $$` only. `$$` delimiters on their own lines. **Even** `$$` count per file. NEVER `\begin{psmallmatrix}`.
- "Key identities" in the Summary must be a **bulleted list**.
- Keep H1 `# CS Foundations` and the anchors line `*Canonical anchors: CLRS; Goldberg (floating-point).*`. Remove the stub callout.
- Citation keys (all already in `references.bib`): `@clrs2009`, `@goldberg1991`, `@trefethen1997`.
- **Wall-clock timing is illustrative only — never hard-asserted** (threaded BLAS makes it noisy).
- Commit identity `jlh530i` / `jlh530i@gmail.com`. Branch → PR → user merges; never self-merge; never commit to main.

## NumPy-verified anchors (all confirmed before planning)

- **Machine epsilon:** $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22\times10^{-16}$; unit roundoff $u = 2^{-53} \approx 1.11\times10^{-16}$. $0.1+0.2 \ne 0.3$; $1+u = 1$; $1+\varepsilon_{\text{mach}} \ne 1$.
- **WE1 cancellation:** $x = 10^8 + \{1,2,3\}$, population variance $2/3 \approx 0.667$. Naive one-pass $\overline{x^2}-\bar x^2 = \mathbf{0.0}$ (total cancellation); two-pass $\overline{(x-\bar x)^2} = 0.667$.
- **WE2 keystone:** near-collinear $X$ (50×2, second column = first + $10^{-6}$ noise): $\kappa(X)\approx1.03\times10^6$, $\kappa(X^\top X)\approx1.05\times10^{12}$, ratio $\kappa(X^\top X)/\kappa(X)^2 = 1.0000$. Known $\beta=(2,-1)$: normal-equations solve error $\approx7.9\times10^{-5}$, QR/lstsq error $\approx7.2\times10^{-11}$.
- **WE3 complexity:** Cholesky $\approx n^3/3$ flops; analytical flop ratio per doubling $= 8\times$ (wall-clock illustrative only).

## Critical files

- `parts/07-swe-cs/01-cs-foundations.qmd` — the chapter (currently a stub; replace body).
- `appendix/solutions.qmd` — append Chapter 23 block before the final `:::` (one `content-visible` gated div).
- Styling exemplar (read-only, AUTHORITATIVE for spacing): `parts/04-optimization/01-convexity.qmd` (Part IV — Rung headings + airy proofs).
- Content/voice exemplars (read-only): `parts/06-causal-calibration/04-prior-store.qmd` (Ch22), `parts/05-mmm-modeling/04-budget-optimization.qmd` (Ch18, the ridge), `parts/01-foundations/01-linear-algebra.qmd` (Ch1, SVD).
- Authoritative design: `docs/superpowers/specs/2026-06-18-chapter-23-cs-foundations-design.md`.

---

### Task 1: Front matter + Motivation

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`

- [ ] **Step 1: Strip the stub.** Keep lines 1–3 (H1 + blank + anchors line). Delete the stub callout and placeholder italics. Keep the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises` with its four `### C/B/P/A` subheads, `## Summary`).

- [ ] **Step 2: Write Motivation (3–4 paragraphs).** Driving question: *"Your model fit, sampler, and optimizer all run on finite-precision arithmetic at finite cost. When does that silently corrupt the answer — and what does it cost to scale?"* Cover:
  - Parts I–VI built models and ran them on a substrate taken for granted: finite-precision arithmetic at finite cost. This chapter makes both precise — opening Part VII, which engineers that substrate.
  - Two pillars: floating-point (when arithmetic silently corrupts an answer) and algorithmic cost (what it takes to scale). The center of gravity is floating-point and conditioning.
  - Preview the unifying payoff: the Chapter 16 identifiability ridge is *also* a numerical conditioning problem — near-collinear spend columns make $X$ nearly singular, and a naive solve loses precision in exactly the direction the data fail to identify.
  - Anchor everything to the MMM stack already built (regression normal equations, Gaussian posteriors, samplers, optimizers, Kalman/prior-store recursions). Not a general CS survey.
  - Do NOT name any library.

- [ ] **Step 3: Verify** even `$$` count, headings present/ordered. Commit: `feat(ch23): motivation`.

---

### Task 2: Theory & Proofs — Rungs 1–3 (cost model, IEEE-754, cancellation)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`. **Use `### Rung N — Title` H3 headings and Part IV airy proof spacing.** Open `## Theory & Proofs` with one short orienting paragraph listing the five rungs (as Part IV / Ch22 do).

- [ ] **Step 1: `### Rung 1 — The cost model and asymptotic complexity`.** RAM model; big-O / $\Theta$ / $\Omega$. Costs of the MMM primitives: matrix–vector $O(n^2)$, matrix–matrix $O(n^3)$, Cholesky $O(n^3/3)$, triangular solve $O(n^2)$, per-sample MCMC cost. **`**Proposition (Cholesky flop count).**`** then airy **Proof**: the $k$-th column costs $\approx k$ multiply–adds over $\approx k$ rows; summing,

  $$
  \sum_{k=1}^{n} k^2 = \frac{n(n+1)(2n+1)}{6} \approx \frac{n^3}{3},
  $$

  end `$\blacksquare$`. One sentence placing sorting / data structures as broader CLRS context, not developed. MMM reading: a $d$-channel Gaussian posterior solve is $O(d^3)$, times draws when sampled.

- [ ] **Step 2: `### Rung 2 — The IEEE-754 model and the rounding bound`.** Double = $\pm m\cdot 2^e$, 52-bit mantissa; $\varepsilon_{\text{mach}} = 2^{-52}$, unit roundoff $u = 2^{-53}$. **`**Theorem (floating-point axiom).**`** $\mathrm{fl}(x \mathbin{op} y) = (x \mathbin{op} y)(1+\delta)$, $|\delta|\le u$. Airy **Proof**: rounding to nearest gives absolute error $\le$ half the spacing at that magnitude, so relative error $\le u$ (display the bound). End `$\blacksquare$`. Note errors compose over many operations.

- [ ] **Step 3: `### Rung 3 — Catastrophic cancellation`.** Subtracting near-equal numbers explodes *relative* error. Airy derivation: with $\hat x = x(1+\delta_x)$, $\hat y = y(1+\delta_y)$, step to the bound

  $$
  \frac{|(\hat x - \hat y) - (x-y)|}{|x-y|} \le \frac{|x|+|y|}{|x-y|}\,u,
  $$

  which blows up as $x\to y$. Anchor: naive one-pass variance $\overline{x^2}-\bar x^2$ vs two-pass $\overline{(x-\bar x)^2}$.

- [ ] **Step 4: Verify** even `$$` count, no bare `\begin{align}`/`psmallmatrix`, Rung 1 & 2 proofs end `$\blacksquare$`, proof spacing is airy (blank lines around every `$$`). Commit: `feat(ch23): theory rungs 1-3 (cost model, IEEE-754, cancellation)`.

---

### Task 3: Theory & Proofs — Rungs 4–5 (keystone conditioning + backward stability)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`. Part IV styling.

- [ ] **Step 1: `### Rung 4 — The condition number and $\kappa(X^\top X)=\kappa(X)^2$ (keystone)`.** Define $\kappa(A) = \lVert A\rVert\,\lVert A^{-1}\rVert = \sigma_{\max}/\sigma_{\min}$. State the forward-error bound (relative solution error $\le \kappa(A)$ × relative perturbation; $\log_{10}\kappa$ digits lost). **`**Theorem.**`** $\kappa(X^\top X) = \kappa(X)^2$. Airy **Proof** via the SVD: with $X = U\Sigma V^\top$,

  $$
  X^\top X = V\Sigma^\top U^\top U \Sigma V^\top = V \Sigma^2 V^\top,
  $$

  so the singular values of $X^\top X$ are $\sigma_i(X)^2$ and

  $$
  \kappa(X^\top X) = \frac{\sigma_{\max}^2}{\sigma_{\min}^2} = \kappa(X)^2 .
  $$

  End `$\blacksquare$`. Conclusion: forming/solving the normal equations doubles the digits lost; use QR/SVD ($\kappa(X)$) not $X^\top X$ ($\kappa(X)^2$). **MMM tie (unifying observation):** the Chapter 16 ridge is near-collinear columns, $\sigma_{\min}(X)\approx0$, $\kappa(X)$ enormous; the identifiability ridge *is* numerical ill-conditioning, and the direction calibration (Chapter 21) compresses is the direction a naive solve loses.

- [ ] **Step 2: `### Rung 5 — Backward stability and the practical upshot`.** Define backward stability (computed solution is exact for a slightly perturbed problem). Cholesky/QR solves are backward stable, so achievable accuracy is $\sim\kappa(A)\,u$ — cannot beat conditioning, but a stable algorithm does not worsen it [@trefethen1997]. Practical rules: never invert (solve); never form $X^\top X$ (QR); sum stably; two-pass variance. Close forward to Part VII (architecture/data-eng/testing protect cost and correctness at scale).

- [ ] **Step 3: Verify** even `$$` count, `@trefethen1997` used, Rung 4 proof ends `$\blacksquare$`, airy spacing. Commit: `feat(ch23): theory rungs 4-5 (keystone conditioning + backward stability)`.

---

### Task 4: Worked Examples (WE1–WE3)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`. Use `### WE1/WE2/WE3 — Title` subheadings and Part IV airy display spacing. Use the EXACT anchor numbers.

- [ ] **Step 1: WE1 — Catastrophic cancellation in the variance.** $x = 10^8 + \{1,2,3\}$, population variance $2/3 \approx 0.667$; naive one-pass returns $\mathbf{0.0}$ (total cancellation — variance below the resolution near $10^{16}$), two-pass returns $0.667$. Lesson: algebraic equivalence is not numerical equivalence.

- [ ] **Step 2: WE2 — $\kappa(X^\top X)=\kappa(X)^2$ and the digits lost (keystone).** Near-collinear $X$: $\kappa(X)\approx1.0\times10^6$, $\kappa(X^\top X)\approx1.05\times10^{12}$, ratio $1.0000$. Known $\beta=(2,-1)$: normal-equations error $\approx7.9\times10^{-5}$, QR error $\approx7.2\times10^{-11}$ — a $\approx10^6$ factor, exactly the extra $\kappa(X)$. The ridge as a numerical object.

- [ ] **Step 3: WE3 — The $O(n^3)$ cost of the Gaussian solve.** Cholesky $\approx n^3/3$ flops; tabulate analytical flop counts at growing $n$, confirm $8\times$ per doubling; mention wall-clock is illustrative/noisy. MMM reading: doubling channels octuples per-solve cost, times draws when sampled.

- [ ] **Step 4: Verify** even `$$` count; numbers match anchors. Commit: `feat(ch23): worked examples WE1-WE3`.

---

### Task 5: Code Tie-in (single runnable cell)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`

- [ ] **Step 1: Write a single ```{python}``` cell** under `## Code Tie-in` with a short prose preface. `numpy` + `matplotlib`; seed `np.random.default_rng(23)`; figures end `plt.show()`. Four blocks:
  1. **Machine epsilon:** print $\varepsilon_{\text{mach}}$, $u$; `assert (0.1+0.2) != 0.3`, `assert 1.0 + 2.0**-53 == 1.0`, `assert 1.0 + 2.0**-52 != 1.0`.
  2. **Cancellation (WE1):** `x = 1e8 + np.array([1.,2.,3.])`; naive `(x**2).mean()-x.mean()**2`, two-pass `((x-x.mean())**2).mean()`; `assert naive == 0.0`; `assert abs(twopass - 2/3) < 1e-9`.
  3. **Conditioning (WE2):** build near-collinear `X` (`t = linspace`, second col = `t + 1e-6*randn`); `kX = np.linalg.cond(X)`, `kXtX = np.linalg.cond(X.T@X)`; `assert abs(kXtX/kX**2 - 1) < 0.05`; solve known `beta=(2,-1)` via `np.linalg.solve(X.T@X, X.T@b)` and via `np.linalg.lstsq(X,b)`; `assert err_qr < err_ne` and `assert err_ne/err_qr > 1e3`. Figure: solve error vs $\kappa(X)$ sweep (normal-eqns $\sim\kappa^2 u$ vs QR $\sim\kappa u$), log–log.
  4. **Complexity (WE3):** analytical flop counts `n**3/3` at `n in [...]`; `assert` doubling ratio $\approx 8$ (`abs(flops(2n)/flops(n) - 8) < 1e-6`); time `np.linalg.cholesky` and plot as illustration (NO timing assert). Figure: cost vs $n$ log–log with slope-3 reference.

- [ ] **Step 2: Extract and run headless.** Copy to `/tmp/ch23_code.py`, run `MPLBACKEND=Agg python3`. Expected: all asserts pass; prints show $u\approx1.11e{-16}$, naive variance $0.0$ / two-pass $0.667$, $\kappa$ ratio $\approx1$, QR error $\ll$ normal-eqns error, flop ratio $8$. Fix until clean.

- [ ] **Step 3: Verify** even `$$` count; no banned libraries (`grep -inwE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm'`). Commit: `feat(ch23): code tie-in (epsilon, cancellation, conditioning, complexity)`.

---

### Task 6: Exercises (C / B / P / A — self-contained)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`

- [ ] **Step 1: Write the four exercise blocks** under the existing `### C/B/P/A` subheads.
  - **C:** (C1) Why never invert a matrix or form $X^\top X$ to solve least squares? (C2) What is machine epsilon, and why is catastrophic cancellation about *relative* error? (C3) In what sense is the Chapter 16 ridge simultaneously a statistical and a numerical object?
  - **B:** (B1) Compute the unit roundoff for double precision and state how many decimal digits it corresponds to. (B2) For a small $2\times2$ near-collinear $X$, compute $\kappa(X)$ and $\kappa(X^\top X)$ from the singular values and the digits lost by each route. (B3) Evaluate naive vs two-pass variance on a small large-mean dataset by hand.
  - **P:** (P1) Prove $\kappa(X^\top X)=\kappa(X)^2$ via the SVD. (P2) Prove the rounding bound $\mathrm{fl}(x\mathbin{op}y)=(x\mathbin{op}y)(1+\delta)$, $|\delta|\le u$, and derive the cancellation relative-error bound $\frac{|x|+|y|}{|x-y|}u$.
  - **A:** (A1) Sweep the conditioning of $X$ and plot realized normal-equations and QR solve error vs $\kappa(X)$, confirming the $\kappa^2 u$ and $\kappa u$ slopes. (A2) Add a calibration-style rank-one precision increment (Chapter 21) and show it lowers $\kappa$ and the realized error together — the ridge healed numerically as well as statistically.

- [ ] **Step 2: Verify** even `$$` count. Commit: `feat(ch23): exercises C/B/P/A`.

---

### Task 7: Summary (auto-included)

**Files:** Modify `parts/07-swe-cs/01-cs-foundations.qmd`

- [ ] **Step 1: Write `## Summary`** with a one-paragraph wrap, then **bulleted** "Key concepts" and **bulleted** "Key identities" (inline math). Identities (each a bullet):
  - unit roundoff $u = 2^{-53} \approx 1.11\times10^{-16}$
  - floating-point axiom $\mathrm{fl}(x\mathbin{op}y) = (x\mathbin{op}y)(1+\delta)$, $|\delta|\le u$
  - cancellation bound $\frac{|x|+|y|}{|x-y|}\,u$
  - condition number $\kappa(A) = \sigma_{\max}/\sigma_{\min}$ and forward-error bound (relative error $\le \kappa\cdot$ perturbation)
  - keystone $\kappa(X^\top X) = \kappa(X)^2$
  - Cholesky cost $\approx n^3/3$
  Close tying back to Chapter 16 (ridge as ill-conditioning) and forward to Part VII.

- [ ] **Step 2: Verify** Key identities bulleted; even `$$` count. Commit: `feat(ch23): summary`.

---

### Task 8: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Append `## Chapter 23 — CS Foundations`** immediately before the final `:::`. Use Part IV airy proof spacing in the P-block. Full C/B/P/A:
  - **C1–C3:** forming/inverting $X^\top X$ squares conditioning ($\kappa^2$) and doubles digits lost, so solve via QR; machine epsilon is the gap above $1$ ($2^{-52}$), cancellation is relative because the absolute operand errors are fixed while the *difference* shrinks; the ridge is statistical (unidentified contrast direction) and numerical (near-singular $X$, huge $\kappa$) — the same direction.
  - **B1:** $u = 2^{-53} \approx 1.11\times10^{-16}$, $\approx 15$–16 decimal digits. **B2:** for $X$ with singular values $\sigma_{\max}, \sigma_{\min}$: $\kappa(X)=\sigma_{\max}/\sigma_{\min}$, $\kappa(X^\top X)=(\sigma_{\max}/\sigma_{\min})^2$; digits lost $\log_{10}\kappa$ vs $2\log_{10}\kappa$. **B3:** large-mean naive vs two-pass worked on small numbers.
  - **P1:** SVD proof of $\kappa(X^\top X)=\kappa(X)^2$ (airy). **P2:** rounding bound (round-to-nearest, half-spacing) and cancellation bound derivation (airy). Both end `$\blacksquare$`.
  - **A1, A2:** describe the expected log–log slopes (normal-eqns slope $\propto\kappa^2 u$, QR $\propto\kappa u$); the rank-one increment raises $\sigma_{\min}$, lowering $\kappa$ and realized error together.

- [ ] **Step 2: Verify** even `$$` count in `appendix/solutions.qmd`; block before final `:::`; proofs end `$\blacksquare$`. Commit: `feat(ch23): appendix solutions (C/B/P/A)`.

---

### Task 9: Final review + PR (no bib additions needed)

**Files:** verify `parts/07-swe-cs/01-cs-foundations.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: Confirm citations resolve.** `@clrs2009`, `@goldberg1991`, `@trefethen1997` all already in `references.bib` (verified at plan time) — no additions.

- [ ] **Step 2: Lint pass on both files.** Confirm: even `$$` count in chapter and appendix; no bare `\begin{align}`; no `\begin{psmallmatrix}`; six template headings present and in order; H1 and anchors line intact; stub callout gone; all proofs end `$\blacksquare$`; Key identities bulleted; only `@clrs2009`/`@goldberg1991`/`@trefethen1997` cited and all resolve; no banned library names (`grep -inwE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm'`). **Styling check:** rungs are `### Rung N` H3 headings, proofs are airy (blank lines around every `$$`), multi-line math uses `\begin{aligned}` — matching Part IV.

- [ ] **Step 3: Re-run the extracted code cell headless** (`MPLBACKEND=Agg python3`); all (non-timing) asserts pass.

- [ ] **Step 4: Commit** any final lint fixes: `chore(ch23): final lint`.

- [ ] **Step 5: Push the branch and open a PR** against `main` (title `feat: Chapter 23 — CS Foundations`; body summarizing the keystone $\kappa(X^\top X)=\kappa(X)^2$ + the ridge-as-ill-conditioning tie + Part VII opening). Start the background merge-watcher.

- [ ] **Step 6: Watch CI `quarto render` (HTML + PDF) to a green conclusion** — the real gate. If it fails, read the log, fix, push, re-watch. Report the green PR; the user merges.

---

## Self-Review (done at plan time)

- **Spec coverage:** Rungs 1–5 → Tasks 2–3; WE1–3 → Task 4; code 4 blocks → Task 5; exercises C/B/P/A → Task 6; summary → Task 7; appendix → Task 8. All spec sections mapped.
- **Anchor consistency:** machine eps ($u=1.11e{-16}$), WE1 (naive $0.0$ / two-pass $0.667$), WE2 ($\kappa(X)\approx10^6$, $\kappa(X^\top X)\approx10^{12}$, errors $7.9e{-5}$ vs $7.2e{-11}$), WE3 (flop ratio $8\times$) identical across Theory, Worked Examples, Code, Appendix — all NumPy-verified above.
- **Styling:** Part IV `### Rung` headings + airy proof spacing + `\begin{aligned}` baked into Tasks 2, 3, 8 and the final-review styling check (Task 9 Step 2), per the user's explicit request.
- **No placeholders; house rules** (library ban, KaTeX, even `$$`, bulleted Key identities, `$\blacksquare$`, no timing asserts, branch→PR→user-merges) enforced per task.
