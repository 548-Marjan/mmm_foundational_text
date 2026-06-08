# Chapter 1 — Linear Algebra Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the stub `parts/01-foundations/01-linear-algebra.qmd` with a complete, MMM-targeted linear algebra chapter (theory + three proofs, worked examples, a runnable NumPy/SciPy tie-in, and C/B/P/A exercises with appendix solutions).

**Architecture:** A single Quarto `.qmd` file authored section by section against the book's fixed template, with the media-mix **design matrix `X`** as the through-line. Exercise solutions live in the shared `appendix/solutions.qmd`, gated by the `show-solutions` metadata flag. Verification is build-based, not unit-test-based.

**Tech Stack:** Quarto (book project), KaTeX math, Python 3.11 with NumPy / SciPy (pinned in `requirements.txt`), `references.bib` (keys `@strang2016`, `@axler2015`, `@trefethen1997`).

---

## Verification model (read first)

This is prose + math + one code cell, so there is no pytest suite. Each task is
verified by one or more of:

- **PY** — extract the chapter's Python and run it standalone:
  `python3 /tmp/ch01_codecell.py` → expect clean exit, printed output as described.
- **MATH** — eyeball that every `$...$` / `$$...$$` is balanced and uses only
  KaTeX-supported constructs (no `\begin{align}` without `aligned`, no custom
  macros). Quick check: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
  returns an **even** number.
- **RENDER** — full `quarto render` (HTML + PDF). Quarto is **not installed
  locally**; this gate runs in CI (`.github/workflows/render.yml`) on the PR. If a
  local check is wanted, install Quarto first (`brew install quarto`).

Commit after every task. Keep the working tree clean between tasks.

## File structure

- **Modify:** `parts/01-foundations/01-linear-algebra.qmd` — the whole chapter
  (currently a stub; replace section-by-section, preserving the H2 headings already
  present so the template is honored).
- **Modify:** `appendix/solutions.qmd` — append a "Chapter 1 — Linear Algebra"
  section holding worked solutions, gated by `show-solutions`.
- **Reference (do not edit):** `references.bib`, `index.qmd`, `_quarto.yml`.

Heading skeleton to preserve in the chapter file (from the stub):
`# Linear Algebra` → `## Motivation` → `## Theory & Proofs` → `## Worked Examples`
→ `## Code Tie-in` → `## Exercises` with `### C…`, `### B…`, `### P…`, `### A…`.

## Notation conventions (use consistently everywhere)

- Design matrix `$X \in \mathbb{R}^{n \times p}$` (n weeks, p columns), response
  `$y \in \mathbb{R}^n$`, coefficients `$\beta \in \mathbb{R}^p$`, fit
  `$\hat y = X\hat\beta$`.
- Transpose `$X^\top$`; Gram matrix `$X^\top X$`; column space `$\mathcal{C}(X)$`;
  null space `$\mathcal{N}(X)$`; rank `$r$`.
- Eigenpairs `$(\lambda_i, v_i)$`; SVD `$X = U\Sigma V^\top$`; Cholesky
  `$A = LL^\top$` with `$L$` lower-triangular, positive diagonal.
- Inner product `$\langle u, v\rangle = u^\top v$`, norm `$\lVert v\rVert$`.

---

### Task 1: Motivation

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (replace the `## Motivation`
  stub body; remove the `.callout-note` "Stub" block at the top).

- [ ] **Step 1: Remove the stub callout and write the Motivation section**

Delete lines 5–7 (the `::: {.callout-note ...} Stub. ... :::` block). Keep
`# Linear Algebra` and the canonical-anchors line. Under `## Motivation`, write
~250–350 words that:
- Introduce `X` concretely: rows = weeks, columns = `[intercept, TV spend, search
  spend, price, seasonality]`, `y` = weekly sales.
- Pose the driving question: *"You ran TV and search together every week. How much
  did each drive sales — and when can you even tell them apart?"*
- Name the three things the chapter delivers against that question: (1) the column
  space = "sales patterns the media can explain," (2) least squares = projecting
  `y` onto it, (3) the eigenvalues of `$X^\top X$` = whether the answer is stable.
- Close with a one-line map of the theory ladder and a forward pointer to Part II
  (regression) and Part III (sampling). Cite `@strang2016` for the column picture.

- [ ] **Step 2: MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: an even number (0 is fine if Motivation uses only inline math).

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): write Motivation section, drop stub callout"
```

---

### Task 2: Theory rungs 1–3 (vectors → subspaces/rank → projection & normal equations)

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (under `## Theory & Proofs`).

- [ ] **Step 1: Write rung 1 — vectors, matrices, matrix–vector product as a linear combination of columns**

Define vectors/matrices on the `X` object. State the key identity explicitly:
`$$X\beta = \sum_{j=1}^{p} \beta_j\, x_j,\qquad x_j = \text{the } j\text{-th column of } X.$$`
Interpret: a fitted sales curve is a weighted blend of channel columns. ~200 words.

- [ ] **Step 2: Write rung 2 — column space, null space, rank, and collinearity**

Define `$\mathcal{C}(X)$`, `$\mathcal{N}(X)$`, and rank `$r \le \min(n,p)$`. State
the equivalence to surface for MMM:
> Two channels are perfectly collinear ⟺ their columns are linearly dependent ⟺
> `$\operatorname{rank}(X) < p$` ⟺ `$X^\top X$` is singular ⟺ `$\beta$` is not
> identified.
Give the near-collinear (vs exactly collinear) caveat: rank stays `$p$` but
`$X^\top X$` is ill-conditioned. ~250 words.

- [ ] **Step 3: Write rung 3 — orthogonality, projection, and the normal equations (PROOF 1)**

Define orthogonal projection of `$y$` onto `$\mathcal{C}(X)$`. State and **prove**
the keystone:

> **Theorem (Normal equations).** If `$X$` has full column rank, the unique
> minimizer of `$\lVert y - X\beta\rVert^2$` is `$\hat\beta = (X^\top X)^{-1}X^\top
> y$`, equivalently the solution of `$X^\top X\hat\beta = X^\top y$`, and the
> residual `$y - X\hat\beta$` is orthogonal to every column of `$X$`.

Proof to include (geometry + calculus, either is fine; use the orthogonality one):
expand `$f(\beta)=\lVert y-X\beta\rVert^2$`, set `$\nabla f = -2X^\top(y-X\beta)=0$`,
giving `$X^\top X\beta = X^\top y$`; full column rank ⟹ `$X^\top X$` invertible (forward-
ref to Task 4's PD result) ⟹ unique solution; convexity (Hessian `$2X^\top X \succeq
0$`) ⟹ it is the global min. Note residual orthogonality `$X^\top(y-X\hat\beta)=0$`.
~300 words incl. proof.

- [ ] **Step 4: MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: even number.

- [ ] **Step 5: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): theory rungs 1-3 incl. normal-equations proof"
```

---

### Task 3: Theory rungs 4–5 (spectral theorem → positive definiteness)

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (under `## Theory & Proofs`).

- [ ] **Step 1: Write rung 4 — eigenvalues/eigenvectors and the spectral theorem (PROOF 2)**

Define `$Av=\lambda v$`. State and prove (at the real-symmetric level):

> **Theorem (Spectral theorem, symmetric case).** A real symmetric matrix
> `$A=A^\top$` has only real eigenvalues, and eigenvectors for distinct eigenvalues
> are orthogonal; hence `$A = Q\Lambda Q^\top$` for an orthogonal `$Q$`.

Proof to include for the provable core: (i) **real eigenvalues** — for `$Av=\lambda
v$`, `$\bar v^\top A v = \lambda \lVert v\rVert^2$` and equals its own conjugate
since `$A$` real symmetric ⟹ `$\lambda\in\mathbb{R}$`. (ii) **orthogonality** — if
`$Av_1=\lambda_1 v_1$`, `$Av_2=\lambda_2 v_2$`, `$\lambda_1\ne\lambda_2$`, then
`$\lambda_1 v_1^\top v_2 = (Av_1)^\top v_2 = v_1^\top A v_2 = \lambda_2 v_1^\top
v_2$` ⟹ `$v_1^\top v_2 = 0$`. State the full orthonormal-basis / `$Q\Lambda Q^\top$`
form and **cite `@axler2015`** for the complete diagonalizability argument (induction
on invariant subspaces). Tie back: `$X^\top X$` is symmetric, so it always
diagonalizes this way. ~320 words incl. proof.

- [ ] **Step 2: Write rung 5 — positive (semi)definiteness and the quadratic form (PROOF 3)**

Define `$A\succeq 0$` (PSD) and `$A\succ 0$` (PD) via `$z^\top A z \ge 0$` /
`$>0$ for $z\ne 0$`. State and prove:

> **Proposition.** A symmetric `$A$` is positive definite ⟺ all its eigenvalues are
> positive. Consequently `$A\succ0 \Rightarrow A$` is invertible.

Proof: use `$A=Q\Lambda Q^\top$` from rung 4; substitute `$w=Q^\top z$` so `$z^\top
A z = w^\top \Lambda w = \sum_i \lambda_i w_i^2$`. This is `$>0$` for all `$w\ne0$`
⟺ every `$\lambda_i>0$`. Positive eigenvalues ⟹ `$\det A=\prod\lambda_i>0$` ⟹
invertible. Tie back: `$z^\top X^\top X z = \lVert Xz\rVert^2 \ge 0$` always (so
`$X^\top X\succeq0$`), and `$>0$` exactly when `$X$` has full column rank — closing
the loop with Task 2's identification claim and Task 2's forward reference.
~300 words incl. proof.

- [ ] **Step 3: MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: even number.

- [ ] **Step 4: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): theory rungs 4-5, spectral theorem + PD proofs"
```

---

### Task 4: Theory rungs 6–7 (SVD, Cholesky)

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (under `## Theory & Proofs`).

- [ ] **Step 1: Write rung 6 — the SVD**

State `$X = U\Sigma V^\top$` (`$U\in\mathbb{R}^{n\times n}$`,
`$V\in\mathbb{R}^{p\times p}$` orthogonal, `$\Sigma$` diagonal with singular values
`$\sigma_1\ge\cdots\ge\sigma_r>0$`). Give the geometric reading (rotate–scale–rotate)
and the two key links: `$X^\top X = V\Sigma^2 V^\top$` so `$\sigma_i^2$` are the
eigenvalues of the Gram matrix, and the **condition number** `$\kappa(X) =
\sigma_1/\sigma_r$` measures collinearity / least-squares sensitivity. **Cite
`@trefethen1997`** for existence (do not prove). ~250 words.

- [ ] **Step 2: Write rung 7 — Cholesky factorization**

State: every SPD `$A$` has a unique `$A=LL^\top$`, `$L$` lower-triangular with
positive diagonal. Give the one-paragraph "why it matters" pointer: solving
`$X^\top X\hat\beta = X^\top y$` stably, and **generating correlated draws** — if
`$\Sigma = LL^\top$` and `$z\sim\mathcal N(0,I)$`, then `$Lz\sim\mathcal N(0,\Sigma)$`
— forward-referencing Part II (Bayesian posteriors) and Part III (HMC mass
matrices). ~200 words.

- [ ] **Step 3: MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: even number.

- [ ] **Step 4: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): theory rungs 6-7, SVD and Cholesky"
```

---

### Task 5: Worked Examples

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (replace `## Worked Examples`
  stub body).

- [ ] **Step 1: Write Example (a) — normal equations by hand**

Use a tiny 3-week, 2-column matrix (intercept + one channel) so arithmetic is
clean, e.g.
`$$X=\begin{bmatrix}1&1\\1&2\\1&3\end{bmatrix},\quad y=\begin{bmatrix}1\\2\\2\end{bmatrix}.$$`
Compute `$X^\top X=\begin{bmatrix}3&6\\6&14\end{bmatrix}$`,
`$X^\top y=\begin{bmatrix}5\\11\end{bmatrix}$`, solve `$\hat\beta$` by hand
(`$\det=6$`, invert), and report `$\hat\beta=(2/3,\,1/2)$`. Show the residual is
orthogonal to both columns as a sanity check. (Verify these numbers in Task 6's code
cell.)

- [ ] **Step 2: Write Example (b) — near-collinearity blows up the estimate**

Replace the channel column with two nearly-parallel channels (e.g. columns
`$[1,2,3]^\top$` and `$[1,2,3.01]^\top$`). Show `$X^\top X$` is nearly singular:
compute its eigenvalues (one near zero), report the large condition number, and
contrast `$\hat\beta$` before/after a tiny perturbation of `$y$` to show the
estimate is unstable. Keep numbers explicit so they match the code cell.

- [ ] **Step 3: PY check (pre-verify the example arithmetic)**

Write `/tmp/ch01_examples_check.py` reproducing both examples with NumPy and
`print` the computed `$\hat\beta$`, eigenvalues, and condition numbers.
Run: `python3 /tmp/ch01_examples_check.py`
Expected: Example (a) prints `beta = [0.667 0.5]` (≈), residual·columns ≈ 0;
Example (b) prints a small eigenvalue (~1e-2 or less) and large `cond` (≫ 1e2).
Fix the prose numbers to match the printout if they differ.

- [ ] **Step 4: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): worked examples (normal equations + collinearity)"
```

---

### Task 6: Code Tie-in (runnable NumPy/SciPy cell)

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (replace `## Code Tie-in`
  stub body with one ```` ```{python} ```` cell).

- [ ] **Step 1: Write the code cell**

Insert a single executable Quarto Python cell. Exact content (adjust comments to
match prose voice, keep the logic):

````markdown
```{python}
import numpy as np
from scipy.linalg import cholesky

rng = np.random.default_rng(0)
n = 200

def make_X(corr):
    """Design matrix: intercept, TV, and a search column 'corr'-correlated with TV."""
    tv = rng.normal(size=n)
    search = corr * tv + np.sqrt(1 - corr**2) * rng.normal(size=n)
    return np.column_stack([np.ones(n), tv, search])

beta_true = np.array([1.0, 2.0, -0.5])

def solve_three_ways(X, y):
    normal = np.linalg.solve(X.T @ X, X.T @ y)         # normal equations
    lstsq  = np.linalg.lstsq(X, y, rcond=None)[0]       # QR-based
    U, s, Vt = np.linalg.svd(X, full_matrices=False)    # SVD
    svd = Vt.T @ ((U.T @ y) / s)
    return normal, lstsq, svd

for corr in [0.3, 0.999]:
    X = make_X(corr)
    y = X @ beta_true + 0.1 * rng.normal(size=n)
    normal, lstsq, svd = solve_three_ways(X, y)
    G = X.T @ X
    eig = np.linalg.eigvalsh(G)
    print(f"corr={corr}: cond(X)={np.linalg.cond(X):.1f}, "
          f"min eig(XtX)={eig.min():.3g}")
    print(f"  normal vs lstsq max|diff| = {np.max(np.abs(normal - lstsq)):.2e}")
    print(f"  lstsq  vs svd   max|diff| = {np.max(np.abs(lstsq  - svd)):.2e}")

# Cholesky of an SPD Gram matrix: the factorization behind correlated sampling.
X = make_X(0.3)
L = cholesky(X.T @ X, lower=True)
print("Cholesky reconstructs XtX:", np.allclose(L @ L.T, X.T @ X))
```
````

Add 2–3 sentences after the cell reading the output: at `corr=0.3` all three solvers
agree and `cond` is small; at `corr=0.999` `cond` explodes and the normal-equations
route loses the most precision, motivating QR/SVD — the numerical payoff of the
theory.

- [ ] **Step 2: PY check — extract and run the cell standalone**

Copy the cell body (between the ```` ```{python} ```` fences) into
`/tmp/ch01_codecell.py`.
Run: `python3 /tmp/ch01_codecell.py`
Expected: clean exit (code 0); prints two `corr=` blocks where `corr=0.3` shows
`cond(X)` ~ single/low-double digits and tiny diffs, `corr=0.999` shows large
`cond(X)` (≥ 1e2) and a larger normal-vs-lstsq diff; final line
`Cholesky reconstructs XtX: True`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): runnable NumPy/SciPy code tie-in"
```

---

### Task 7: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/01-foundations/01-linear-algebra.qmd` (fill the four `###` exercise
  subsections).

- [ ] **Step 1: Write the four tiers**

Populate each tier with 2–3 numbered exercises. Concretely:

- **C — Conceptual:** (C1) In words, what is `$\mathcal{C}(X)$` for a media matrix
  and what does it mean for `$y\notin\mathcal{C}(X)$`? (C2) Why does perfect
  collinearity make `$\beta$` unidentifiable — phrase it via columns *and* via
  `$X^\top X$`. (C3) What does `$\kappa(X)$` being large imply for an estimated ROI?
- **B — By hand:** (B1) For `$X=[[1,0],[1,1],[1,2]]$`, `$y=[1,0,2]^\top$`, solve the
  normal equations. (B2) Find eigenvalues/eigenvectors of
  `$\begin{bmatrix}2&1\\1&2\end{bmatrix}$` and confirm orthogonality. (B3) Show
  `$\begin{bmatrix}2&1\\1&2\end{bmatrix}\succ0$` via the eigenvalue test.
- **P — Prove it:** (P1) Prove `$\hat\beta=(X^\top X)^{-1}X^\top y$` minimizes
  `$\lVert y-X\beta\rVert^2$` (full column rank). (P2) Prove eigenvectors of a
  symmetric matrix for distinct eigenvalues are orthogonal. (P3) Prove
  `$A\succ0\Rightarrow A$` invertible.
- **A — Applied / code:** (A1) Extend the code cell: sweep `corr` over
  `[0, 0.5, 0.9, 0.99, 0.999]` and plot `$\kappa(X)$` vs the empirical variance of
  `$\hat\beta$` over 100 resamples. (A2) Draw `$10^4$` samples from
  `$\mathcal N(0,\Sigma)$` with `$\Sigma=\begin{bmatrix}1&0.8\\0.8&1\end{bmatrix}$`
  using its Cholesky factor; confirm the empirical covariance matches.

Each exercise is self-contained; **no inline links to solutions** (per the preface).

- [ ] **Step 2: MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: even number.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/01-linear-algebra.qmd
git commit -m "feat(ch01): exercises across C/B/P/A tiers"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (append a Chapter 1 section, gated by
  `show-solutions`).

- [ ] **Step 1: Inspect the existing solutions file structure**

Run: `sed -n '1,40p' appendix/solutions.qmd`
Match its existing gating idiom (e.g. a `{.content-visible when-meta="show-solutions"}`
block or `eval` on the metadata flag). Reuse whatever pattern is already there; do
not invent a new one.

- [ ] **Step 2: Write solutions for B and P tiers (and A sanity answers)**

Add `## Chapter 1 — Linear Algebra` with worked solutions for B1–B3, P1–P3, and the
expected numeric/plot outcomes for A1–A2. For B1 give `$\hat\beta$` explicitly; for
P1–P3 give the full proofs (they mirror the chapter's three proofs — restate them
self-contained). Wrap the section in the same `show-solutions` gate the file already
uses so reader builds omit it.

- [ ] **Step 3: PY check — verify B1 and A2 numbers**

Add B1 and A2 computations to `/tmp/ch01_sol_check.py` and run:
`python3 /tmp/ch01_sol_check.py`
Expected: B1 `$\hat\beta$` matches the prose; A2 empirical covariance ≈
`$\begin{bmatrix}1&0.8\\0.8&1\end{bmatrix}$` within Monte-Carlo error.

- [ ] **Step 4: Commit**

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch01): appendix solutions for B/P/A exercises"
```

---

### Task 9: Full-chapter integration check

**Files:**
- Possibly modify: `parts/01-foundations/01-linear-algebra.qmd`,
  `appendix/solutions.qmd` (only fixes surfaced by the checks below).

- [ ] **Step 1: Cross-reference and notation sweep**

Read the whole chapter once. Confirm: notation matches the conventions table
everywhere; the three forward/back references between Tasks 2, 3, and 5 actually
resolve (identification ↔ PD ↔ full column rank); every `###` exercise tier is
populated; no leftover stub text ("_Why this matters…_", "Stub", etc.).
Run: `grep -ni 'stub\|_why this matters\|TODO\|TBD' parts/01-foundations/01-linear-algebra.qmd`
Expected: no matches.

- [ ] **Step 2: Final MATH check**

Run: `grep -c '\$\$' parts/01-foundations/01-linear-algebra.qmd`
Expected: even number. Spot-check that no `$$` block uses unsupported KaTeX.

- [ ] **Step 3: RENDER gate**

If Quarto is installed locally: `quarto render parts/01-foundations/01-linear-algebra.qmd`
→ expect clean HTML build, code cell executes, no math errors. Otherwise this gate
runs in CI on the PR (`.github/workflows/render.yml`) — note in the PR that local
Quarto was unavailable and CI is the render authority.

- [ ] **Step 4: Commit any fixes**

```bash
git add -A
git commit -m "fix(ch01): integration sweep — notation, refs, stub cleanup"
```

---

## Self-Review (completed by plan author)

**Spec coverage:** Through-line → Tasks 1–6. Theory ladder rungs 1–7 → Tasks 2–4.
Three named proofs (normal equations, spectral, PD) → Tasks 2, 3. Worked examples
(a)/(b) → Task 5. NumPy/SciPy tie-in with three solvers + Cholesky → Task 6. C/B/P/A
exercises → Task 7. Appendix solutions gated by `show-solutions` → Task 8. Render +
notation success criteria → Task 9. SVD existence cited not proved → Task 4. Out-of-
scope items are simply never introduced. All spec sections map to a task.

**Placeholder scan:** No "TBD/TODO/handle edge cases" left in steps; the code cell is
given in full; proof contents are stated explicitly rather than deferred.

**Type/notation consistency:** `make_X`, `solve_three_ways`, `beta_true` names are
consistent across Task 6 and referenced by Tasks 5/7/8; `$X^\top X$`, `$\hat\beta$`,
`$\kappa(X)$`, `$A=LL^\top$` used uniformly per the conventions table.
