# Chapter 2 — Multivariable Calculus & the Optimization Toolkit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/01-foundations/02-multivariable-calculus.qmd` with a complete chapter — Motivation, Theory & Proofs (four full proofs), Worked Examples, a runnable Python code tie-in, and C/B/P/A exercises — plus matching appendix solutions, all rendering cleanly under Quarto.

**Architecture:** The chapter is authored section-by-section into one Quarto `.qmd`, each section a separate commit. The anchor object throughout is the SSE loss surface `L(β) = ‖y − Xβ‖²` from Chapter 1's design matrix `X`. Math is KaTeX-safe; the single code cell is original NumPy/SciPy verified to run standalone; exercise solutions live in the shared `appendix/solutions.qmd` gated by `show-solutions`.

**Tech Stack:** Quarto (book project), KaTeX (HTML math) / LaTeX (PDF), Python 3 with NumPy/SciPy (pinned in `requirements.txt`), BibTeX (`references.bib`, keys `@nocedal2006`, `@boyd2004`, `@strang2016` already present).

**Design spec:** `docs/superpowers/specs/2026-06-08-chapter-02-multivariable-calculus-design.md`

---

## Environment prerequisites

- **Quarto must be on PATH.** It is not installed in the authoring sandbox; install from <https://quarto.org> (or `brew install quarto`) before running render verifications. Confirm with `quarto --version`.
- Python deps: `python3 -m pip install -r requirements.txt` (only NumPy/SciPy are needed for this chapter's code).
- All paths below are relative to the repo root `/Users/jameshenson/Documents/mmm_foundational_text`.

## Notation conventions (use consistently in every section)

- Loss: `L(β) = ‖y − Xβ‖²` (no ½ factor, so `∇L = −2Xᵀ(y − Xβ)`, `∇²L = 2XᵀX`). Keep this convention everywhere — do not silently switch to a ½ version.
- `X` is `n×p` (n weeks, p channels+controls), `y ∈ ℝⁿ`, `β ∈ ℝᵖ`.
- Residual map `r(β) = y − Xβ`, Jacobian `J = ∂r/∂β = −X`.
- Math delimiters: inline `$...$`, display `$$...$$`. Use only KaTeX-supported macros (`\nabla`, `\succeq`, `\top`, `\mathbf`, `\langle\rangle`, `\|\cdot\|`). Avoid `\xrightarrow`, custom `\newcommand` in-line, and AMS environments KaTeX can't parse; use `\begin{aligned}...\end{aligned}` for multi-line.

## File Structure

- **Modify (replace stub):** `parts/01-foundations/02-multivariable-calculus.qmd` — the whole chapter. One file, one responsibility (the chapter).
- **Modify:** `appendix/solutions.qmd` — append a "Chapter 2" solutions block inside the existing `::: {.content-visible when-meta="show-solutions"}` region.
- **No change:** `references.bib` (`@nocedal2006`, `@boyd2004` already present), `_quarto.yml` (chapter already registered).

---

### Task 1: Motivation section

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Replace the stub body, keeping the H1 title, with the Motivation section**

Open the file. Keep the existing first line `# Multivariable Calculus & the Optimization Toolkit` and the `*Canonical anchors...*` line. Delete the `::: {.callout-note ...} **Stub.** ... :::` block and every `_italic placeholder_` line under each `##` heading, leaving the `##` headings in place as we fill them. Replace the `## Motivation` body with:

```markdown
## Motivation

Chapter 1 found the least-squares coefficients $\hat\beta$ by projecting the
sales vector $y$ onto the column space of the design matrix $X$ — a geometric
move. But suppose you could not see the whole geometry at once. Suppose you
could only stand at some guess $\beta$ and *feel the slope of the error under
your feet*. Two questions decide whether you can still find $\hat\beta$:

1. **Which way is downhill, and how steep?** That is the **gradient**.
2. **If you reach a flat spot, is it the bottom — and is it the *only*
   bottom?** That is **convexity**.

This chapter builds the multivariable-calculus answer to both, on one recurring
object: the sum-of-squares loss surface
$$
L(\beta) = \lVert y - X\beta \rVert^2 .
$$
By the end we will have re-derived the normal equations $X^\top X\hat\beta =
X^\top y$ a second time — by calculus rather than projection — and shown that
the matrix $X^\top X$ that diagnosed collinearity in Chapter 1 is, exactly, the
*curvature* of $L$. The tools (gradient, Jacobian, Hessian, Taylor expansion,
the optimality and convexity conditions) are the same ones every model in this
book is fit with [@nocedal2006; @boyd2004].
```

- [ ] **Step 2: Verify the section renders (Quarto available)**

Run: `quarto render parts/01-foundations/02-multivariable-calculus.qmd --to html`
Expected: exits 0, produces HTML, no KaTeX "ParseError" in output.

If Quarto is not yet installed, defer this check to Task 10 and instead do a math-delimiter sanity check:
Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0, 'unbalanced \$\$'; print('delimiters balanced')"`
Expected: prints `delimiters balanced`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: motivation — the loss surface, slope and convexity"
```

---

### Task 2: Theory — gradient & steepest-ascent proof

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Write the gradient subsection under `## Theory & Proofs`**

Replace the `## Theory & Proofs` placeholder with the heading plus this first subsection:

```markdown
## Theory & Proofs

We climb a ladder of tools, each rung paying off on $L$. Four results are
proved in full; mechanical results are stated and cited.

### Gradient and directional derivative

For $f:\mathbb{R}^p\to\mathbb{R}$ differentiable at $\beta$, the **gradient**
$\nabla f(\beta)$ is the vector of partial derivatives
$\big(\partial f/\partial\beta_1,\dots,\partial f/\partial\beta_p\big)^\top$.
The **directional derivative** in a unit direction $u$ ($\lVert u\rVert=1$) is
the rate of change of $f$ along $u$:
$$
D_u f(\beta) = \lim_{h\to 0}\frac{f(\beta+hu)-f(\beta)}{h}
             = \nabla f(\beta)^\top u .
$$

::: {.callout-note appearance="simple"}
**Theorem (steepest ascent).** Among all unit directions $u$, the directional
derivative $D_u f(\beta)=\nabla f(\beta)^\top u$ is largest when $u$ points
along $\nabla f(\beta)$, and the level set of $f$ through $\beta$ is orthogonal
to $\nabla f(\beta)$.

**Proof.** By Cauchy–Schwarz, $\nabla f(\beta)^\top u \le
\lVert\nabla f(\beta)\rVert\,\lVert u\rVert = \lVert\nabla f(\beta)\rVert$, with
equality iff $u$ is a nonnegative multiple of $\nabla f(\beta)$; the unit
maximizer is therefore $u^\star=\nabla f(\beta)/\lVert\nabla f(\beta)\rVert$
(assuming $\nabla f(\beta)\neq 0$). For the orthogonality claim, let
$\gamma(t)$ be any curve lying in the level set $\{x:f(x)=f(\beta)\}$ with
$\gamma(0)=\beta$. Then $f(\gamma(t))$ is constant, so differentiating and
using the chain rule, $0=\tfrac{d}{dt}f(\gamma(t))\big|_{0}
=\nabla f(\beta)^\top\gamma'(0)$. Hence $\nabla f(\beta)$ is orthogonal to every
tangent direction of the level set. $\qquad\blacksquare$
:::

So the negative gradient is the locally steepest *downhill* direction — the
fact every gradient-based fitting method in this book exploits.
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: theory — gradient, directional derivative, steepest-ascent proof"
```

---

### Task 3: Theory — Jacobian, chain rule, residual map (Gauss–Newton preview)

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Append the Jacobian subsection after the gradient subsection**

```markdown
### Jacobian and the chain rule

A vector-valued map $r:\mathbb{R}^p\to\mathbb{R}^n$ has **Jacobian**
$J(\beta)\in\mathbb{R}^{n\times p}$ with entries $J_{ij}=\partial r_i/\partial
\beta_j$. The **multivariate chain rule** says that for a composition
$g\circ r$ with $g:\mathbb{R}^n\to\mathbb{R}$,
$$
\nabla_\beta\,(g\circ r)(\beta) = J(\beta)^\top\,\nabla g\big(r(\beta)\big).
$$

Apply it to least squares. The **residual map** is $r(\beta)=y-X\beta$, a linear
(hence affine) map with constant Jacobian
$$
J=\frac{\partial r}{\partial\beta} = -X,
$$
and $L(\beta)=g(r(\beta))$ with $g(r)=\lVert r\rVert^2$, $\nabla g(r)=2r$.
The chain rule gives the gradient of the loss directly:
$$
\nabla L(\beta) = J^\top\nabla g(r) = (-X)^\top\,2\,(y-X\beta)
               = -2X^\top(y-X\beta).
$$
The form $\nabla L = 2J^\top r$ is the **Gauss–Newton** shape. Here $J=-X$ is
constant because the model is linear; when we later fit the *saturating* media
curves of Part V, $J$ will depend on $\beta$, but this same expression is what
those nonlinear least-squares fits descend.
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: theory — Jacobian, chain rule, residual map / Gauss-Newton gradient"
```

---

### Task 4: Theory — Hessian, Taylor, first-order & second-order optimality proofs

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Append the Hessian / optimality subsection**

```markdown
### Hessian, Taylor expansion, and optimality

The **Hessian** $\nabla^2 f(\beta)\in\mathbb{R}^{p\times p}$ collects the second
partials $\big[\nabla^2 f\big]_{ij}=\partial^2 f/\partial\beta_i\partial\beta_j$;
for $C^2$ functions it is symmetric (Clairaut). The **second-order Taylor
expansion** about $\beta$ is
$$
f(\beta+\delta) = f(\beta) + \nabla f(\beta)^\top\delta
  + \tfrac12\,\delta^\top\nabla^2 f(\beta)\,\delta + o(\lVert\delta\rVert^2),
$$
with the exact (Lagrange/integral) remainder cited to [@nocedal2006, Thm. 2.1].
For the loss, differentiating $\nabla L(\beta)=-2X^\top(y-X\beta)$ once more,
$$
\nabla^2 L(\beta) = 2X^\top X,
$$
constant in $\beta$ — and *exactly* the symmetric PSD matrix Chapter 1 used to
diagnose collinearity. Its eigenvalues are the curvatures of the bowl $L$.

::: {.callout-note appearance="simple"}
**Theorem (first-order necessary condition).** If $\beta^\star$ is an interior
local minimizer of a differentiable $f$, then $\nabla f(\beta^\star)=0$.

**Proof.** Fix any direction $u$. The scalar $\phi(t)=f(\beta^\star+tu)$ has a
local minimum at $t=0$, so $\phi'(0)=\nabla f(\beta^\star)^\top u=0$. As $u$ was
arbitrary, every component of $\nabla f(\beta^\star)$ vanishes. $\qquad\blacksquare$
:::

::: {.callout-note appearance="simple"}
**Theorem (second-order sufficient condition).** If $\nabla f(\beta^\star)=0$
and $\nabla^2 f(\beta^\star)$ is positive definite, then $\beta^\star$ is a
strict local minimizer.

**Proof.** Let $\lambda_{\min}>0$ be the smallest eigenvalue of
$H=\nabla^2 f(\beta^\star)$ (real and positive by the spectral theorem and PD,
from Chapter 1). With $\nabla f(\beta^\star)=0$, Taylor gives
$$
f(\beta^\star+\delta)-f(\beta^\star)
 = \tfrac12\,\delta^\top H\,\delta + o(\lVert\delta\rVert^2)
 \ge \tfrac12\lambda_{\min}\lVert\delta\rVert^2 + o(\lVert\delta\rVert^2).
$$
For $\lVert\delta\rVert$ small enough the right side is positive, so
$f(\beta^\star+\delta)>f(\beta^\star)$ for all small $\delta\neq 0$.
$\qquad\blacksquare$
:::
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: theory — Hessian, Taylor, first- and second-order optimality proofs"
```

---

### Task 5: Theory — convexity ⇔ PSD Hessian proof + deferred pointers

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Append the convexity subsection and the deferral note**

```markdown
### Convexity and the global guarantee

A function $f$ is **convex** if for all $x,z$ and $t\in[0,1]$,
$f(tx+(1-t)z)\le t f(x)+(1-t)f(z)$.

::: {.callout-note appearance="simple"}
**Theorem (second-order characterization).** A $C^2$ function $f$ on an open
convex domain is convex if and only if $\nabla^2 f(x)\succeq 0$ (positive
semidefinite) for every $x$.

**Proof.** ($\Leftarrow$) For any $x,z$ let $\delta=z-x$ and
$\psi(t)=f(x+t\delta)$ on $[0,1]$. Then $\psi''(t)=\delta^\top\nabla^2
f(x+t\delta)\,\delta\ge 0$, so $\psi$ is a convex scalar function, which is
exactly the convexity inequality for $f$ along the segment.
($\Rightarrow$) If some $\nabla^2 f(x_0)$ had a negative eigenvalue with
eigenvector $v$, then $\psi(t)=f(x_0+tv)$ would satisfy $\psi''(0)=v^\top
\nabla^2 f(x_0)v<0$, making $\psi$ strictly concave near $0$ and violating
convexity along that segment. $\qquad\blacksquare$
:::

Apply this to the loss. Since $\nabla^2 L(\beta)=2X^\top X$ and
$v^\top X^\top X v=\lVert Xv\rVert^2\ge 0$ for all $v$, the Hessian is PSD
everywhere, so **$L$ is convex**. Therefore any $\beta$ with $\nabla L(\beta)=0$
— any solution of the normal equations — is a *global* minimizer, not merely a
local one. When $X$ has full column rank, $X^\top X$ is positive definite, $L$
is strictly convex, and $\hat\beta$ is the *unique* global minimizer; rank
deficiency (collinearity) is precisely when that guarantee weakens.

### What this chapter does not cover

Descending $L$ in practice — **gradient descent**, Newton's method, line
search — and **constrained** optimization via **Lagrange multipliers** and the
**KKT conditions** are the subject of Part IV [@nocedal2006; @boyd2004]; we name
them here only to mark where this calculus leads. Integration and
change-of-variables for probability densities arrive in the
probability/statistics chapter, where they are needed.
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: theory — convexity <=> PSD Hessian proof, global-min guarantee, deferrals"
```

---

### Task 6: Worked Examples

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Replace the `## Worked Examples` placeholder**

```markdown
## Worked Examples

**(a) Gradient, Hessian, and recovering $\hat\beta$ by calculus.** Take three
weeks, two channels (no intercept, to keep arithmetic light):
$$
X=\begin{bmatrix}1&0\\1&1\\0&1\end{bmatrix},\qquad
y=\begin{bmatrix}1\\3\\2\end{bmatrix}.
$$
Then
$$
X^\top X=\begin{bmatrix}2&1\\1&2\end{bmatrix},\qquad
X^\top y=\begin{bmatrix}4\\5\end{bmatrix}.
$$
The gradient is $\nabla L(\beta)=-2X^\top(y-X\beta)=2\big(X^\top X\beta-X^\top
y\big)$. Setting $\nabla L=0$ is the normal equations $X^\top X\beta=X^\top y$:
$$
\begin{aligned}
2\beta_1+\beta_2 &= 4\\
\beta_1+2\beta_2 &= 5
\end{aligned}
\quad\Rightarrow\quad
\hat\beta=\begin{bmatrix}1\\2\end{bmatrix}.
$$
The Hessian $\nabla^2 L=2X^\top X=\begin{bmatrix}4&2\\2&4\end{bmatrix}$ has
eigenvalues $2$ and $6$ (both positive), so by the second-order sufficient
condition $\hat\beta$ is a strict — and, by convexity, global and unique —
minimizer. This is the same $\hat\beta$ the projection of Chapter 1 would
return.

**(b) A saddle: why $\nabla f=0$ is not enough.** Let
$f(\beta)=\beta_1^2-\beta_2^2$. Then $\nabla f=(2\beta_1,-2\beta_2)^\top$, which
vanishes at the origin, yet
$\nabla^2 f=\begin{bmatrix}2&0\\0&-2\end{bmatrix}$ is **indefinite**
(eigenvalues $+2,-2$). Along $\beta_2=0$ the origin looks like a minimum; along
$\beta_1=0$ it looks like a maximum. The first-order condition flags it as a
candidate; only the Hessian reveals it is neither — exactly the gap the
second-order test closes.
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Verify the worked arithmetic is correct**

Run:
```bash
python3 -c "
import numpy as np
X=np.array([[1,0],[1,1],[0,1]]); y=np.array([1,3,2])
b=np.linalg.solve(X.T@X, X.T@y)
assert np.allclose(b,[1,2]), b
assert sorted(np.linalg.eigvalsh(2*X.T@X).round(6))==[2.0,6.0]
print('worked example (a) checks out:', b)
"
```
Expected: prints `worked example (a) checks out: [1. 2.]`.

- [ ] **Step 4: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: worked examples — calculus recovers beta-hat; a saddle point"
```

---

### Task 7: Code Tie-in (runnable, verified standalone)

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Replace the `## Code Tie-in` placeholder with prose and one Python cell**

```markdown
## Code Tie-in

We form the analytic gradient and Hessian of $L$, check the gradient against
finite differences (the habit that catches a hundred modeling bugs later), and
read the Hessian's eigenvalues as curvature — tying the smallest one back to
Chapter 1's collinearity story.

```{python}
import numpy as np

rng = np.random.default_rng(0)

# Synthetic media design matrix: 40 weeks, intercept + 2 channels.
n = 40
intercept = np.ones(n)
tv = rng.normal(size=n)
search = rng.normal(size=n)
X = np.column_stack([intercept, tv, search])
beta_true = np.array([2.0, 1.5, -0.5])
y = X @ beta_true + rng.normal(scale=0.5, size=n)

def loss(beta):
    r = y - X @ beta
    return r @ r

def grad_analytic(beta):
    return -2.0 * X.T @ (y - X @ beta)

def grad_numeric(beta, eps=1e-6):
    g = np.zeros_like(beta)
    for j in range(len(beta)):
        e = np.zeros_like(beta); e[j] = eps
        g[j] = (loss(beta + e) - loss(beta - e)) / (2 * eps)
    return g

beta0 = np.zeros(3)
print("max |analytic - numeric| gradient:",
      np.max(np.abs(grad_analytic(beta0) - grad_numeric(beta0))))

# Hessian is constant: 2 X^T X. Its eigenvalues are the curvatures of the bowl.
H = 2.0 * X.T @ X
eig = np.linalg.eigvalsh(H)
print("Hessian eigenvalues:", np.round(eig, 3))
print("all positive (=> strictly convex, unique minimizer):", bool(np.all(eig > 0)))

# Solve the normal equations (the calculus critical point) and confirm.
beta_hat = np.linalg.solve(X.T @ X, X.T @ y)
print("||gradient at beta_hat||:", np.linalg.norm(grad_analytic(beta_hat)))
print("condition number of X^T X:", np.linalg.cond(X.T @ X))
```

The numeric and analytic gradients agree to roughly $10^{-5}$; every Hessian
eigenvalue is positive, so $L$ is a strictly convex bowl with a unique bottom;
and $\nabla L(\hat\beta)\approx 0$ confirms the normal-equation solution is that
bottom. If the two channels were made nearly collinear, the smallest eigenvalue
— and the condition number — would blow up: a flat valley in $L$, which is what
makes $\hat\beta$ unstable, exactly as in Chapter 1.
```

Note: the ```` ```{python} ```` fence above is the executable Quarto cell. Keep the closing ```` ``` ```` of that cell, then the trailing prose, then the closing ```` ``` ```` of the outer markdown block in this plan is not part of the file — only the cell fence goes in the `.qmd`.

- [ ] **Step 2: Verify the code cell runs standalone**

Extract the cell body to a temp file and run it (this does not require Quarto):
```bash
python3 - <<'PY'
import re, subprocess, sys, tempfile, os
src = open('parts/01-foundations/02-multivariable-calculus.qmd').read()
m = re.search(r'```\{python\}\n(.*?)\n```', src, re.S)
assert m, "no python cell found"
code = m.group(1)
fd, path = tempfile.mkstemp(suffix='.py'); os.write(fd, code.encode()); os.close(fd)
r = subprocess.run([sys.executable, path], capture_output=True, text=True)
print(r.stdout)
if r.returncode != 0:
    print(r.stderr); sys.exit(1)
assert "all positive (=> strictly convex, unique minimizer): True" in r.stdout
print("CODE CELL OK")
PY
```
Expected: prints the cell's output and `CODE CELL OK`. The `max |analytic - numeric|` line should be a small number (< 1e-3).

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: code tie-in — analytic vs numeric gradient, Hessian eigenvalues, conditioning"
```

---

### Task 8: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/01-foundations/02-multivariable-calculus.qmd`

- [ ] **Step 1: Replace the `## Exercises` block (and its four tier sub-headings) with full exercises**

```markdown
## Exercises

### C — Conceptual / Reading Comprehension

1. In one sentence each, say what $\nabla L(\beta)$ and $\nabla^2 L(\beta)$ tell
   you geometrically about the loss surface.
2. Why does convexity of $L$ let you trust a *single* fitted $\hat\beta$, rather
   than worrying you have found one of many local minima?
3. What is a saddle point, and why is $\nabla f=0$ alone not enough to certify a
   minimum? Refer to Worked Example (b).

### B — By Hand

1. For $f(\beta)=3\beta_1^2+\beta_2^2-\beta_1\beta_2$, compute $\nabla f$ and
   $\nabla^2 f$, find the unique critical point, and classify it via the Hessian.
2. For the design matrix $X=\begin{bmatrix}1&2\\1&0\\1&1\end{bmatrix}$ and
   $y=(2,1,2)^\top$, write $\nabla L(\beta)=0$ explicitly and solve for
   $\hat\beta$.
3. Give a $2\times2$ symmetric matrix that is PSD but not PD, and exhibit the
   direction $v$ with $v^\top H v=0$. Relate this to a perfectly collinear pair
   of channels.

### P — Prove It

These extend the chapter's proofs rather than repeat them.

1. **(Second-order necessary condition.)** Prove that if $\beta^\star$ is an
   interior local minimizer of a $C^2$ function $f$, then $\nabla^2
   f(\beta^\star)\succeq 0$. (Hint: apply the scalar restriction $\phi(t)=
   f(\beta^\star+tu)$ and use $\phi''(0)\ge 0$.) Explain why this gives only
   PSD, while the sufficient condition in the text needed PD.
2. **(First-order characterization of convexity.)** Prove that a $C^1$ function
   $f$ is convex if and only if $f(z)\ge f(\beta)+\nabla f(\beta)^\top(z-\beta)$
   for all $\beta,z$ — i.e. the graph lies above each tangent plane.
3. **(Global, unique minimizer.)** Prove that for a convex $f$ every local
   minimizer is global, and that if $f$ is *strictly* convex the minimizer is
   unique. Then state precisely the condition on $X$ under which $\hat\beta$ is
   the unique global minimizer of $L$, and justify it.
4. **(Convexity under affine precomposition.)** Prove that if $g:\mathbb{R}^n\to
   \mathbb{R}$ is convex and $A,c$ are a fixed matrix and vector, then
   $\beta\mapsto g(A\beta+c)$ is convex. Use this — together with the convexity
   of $r\mapsto\lVert r\rVert^2$ — to show $L(\beta)=\lVert y-X\beta\rVert^2$ is
   convex *without computing any derivatives*.

### A — Applied / Code

1. Write a finite-difference gradient checker `grad_numeric(f, beta)` and use it
   to verify the analytic gradient of $L$ at three random points; report the
   maximum discrepancy.
2. Build an $X$ whose two channels are increasingly correlated (e.g.
   `search = rho*tv + sqrt(1-rho**2)*noise` for $\rho\in\{0,0.9,0.99,0.999\}$).
   For each $\rho$, print the smallest eigenvalue of $X^\top X$ and the condition
   number, and take ten fixed-step gradient-descent steps on $L$ from
   $\beta=0$; report how the loss after ten steps degrades as $\rho\to1$. Relate
   the slow-down to the flat valley the small eigenvalue describes.
```

- [ ] **Step 2: Math-delimiter sanity check**

Run: `python3 -c "t=open('parts/01-foundations/02-multivariable-calculus.qmd').read(); assert t.count('$$')%2==0; print('ok')"`
Expected: prints `ok`.

- [ ] **Step 3: Commit**

```bash
git add parts/01-foundations/02-multivariable-calculus.qmd
git commit -m "ch2: exercises — C/B/P/A; P-tier poses adjacent proofs"
```

---

### Task 9: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd`

- [ ] **Step 1: Insert a Chapter 2 solutions section inside the `content-visible` block**

In `appendix/solutions.qmd`, locate the line `### How solutions are organized`
and insert the following block **immediately before** it (so it stays inside the
`::: {.content-visible when-meta="show-solutions"}` region):

```markdown
## Chapter 2 — Multivariable Calculus & the Optimization Toolkit {.unnumbered}

**C.1** $\nabla L$ points in the steepest-uphill direction and is zero only at
the minimizer; $\nabla^2 L=2X^\top X$ is the (constant) curvature of the bowl —
its eigenvalues say how sharply $L$ rises in each principal direction.
**C.2** Convexity makes every critical point a global minimum, so a single
$\nabla L=0$ solution is *the* answer; there are no other valleys to miss.
**C.3** A saddle has $\nabla f=0$ but the Hessian is indefinite, so the point is
a minimum along some directions and a maximum along others (Worked Example (b));
$\nabla f=0$ only marks candidates.

**B.1** $\nabla f=(6\beta_1-\beta_2,\;2\beta_2-\beta_1)^\top$,
$\nabla^2 f=\begin{bmatrix}6&-1\\-1&2\end{bmatrix}$. Setting $\nabla f=0$ gives
$\beta=(0,0)$; the Hessian has positive eigenvalues (trace $8>0$, determinant
$11>0$), so the origin is a strict local — and, since $f$ is a PD quadratic,
global — minimum.
**B.2** $X^\top X=\begin{bmatrix}3&3\\3&5\end{bmatrix}$,
$X^\top y=\begin{bmatrix}5&6\end{bmatrix}^\top$; solving
$3\beta_1+3\beta_2=5,\;3\beta_1+5\beta_2=6$ gives
$\hat\beta=(\tfrac{7}{6},\tfrac12)^\top$.
**B.3** $H=\begin{bmatrix}1&1\\1&1\end{bmatrix}$ is PSD (eigenvalues $0,2$) but
not PD; $v=(1,-1)^\top$ gives $v^\top H v=0$. This is two channels that are
identical columns — a perfectly collinear pair, along which $L$ is flat.

**P.1** Restrict to $\phi(t)=f(\beta^\star+tu)$. As $t=0$ is a local min,
$\phi''(0)=u^\top\nabla^2 f(\beta^\star)u\ge 0$ for every $u$, i.e.
$\nabla^2 f(\beta^\star)\succeq 0$. We get only "$\ge 0$": a flat direction
($u^\top Hu=0$) is consistent with a minimum, so PSD is necessary but PD is what
*guarantees* one.
**P.2** ($\Rightarrow$) Convexity gives, for $t\in(0,1]$, $f(\beta+t(z-\beta))
\le f(\beta)+t\,[f(z)-f(\beta)]$; rearranging and letting $t\to0^+$ yields
$\nabla f(\beta)^\top(z-\beta)\le f(z)-f(\beta)$. ($\Leftarrow$) Given the
tangent inequality, take any $x,z$ and $w=tx+(1-t)z$; apply the inequality at
$\beta=w$ to both $x$ and $z$, multiply by $t$ and $1-t$ and add; the gradient
terms cancel and convexity of $f$ follows.
**P.3** If $\beta^\star$ is a local but not global min, some $z$ has $f(z)<
f(\beta^\star)$; convexity along the segment gives $f(\beta^\star+t(z-\beta^\star))
\le f(\beta^\star)+t[f(z)-f(\beta^\star)]<f(\beta^\star)$ for small $t>0$,
contradicting locality. For uniqueness under strict convexity, two distinct
minimizers $a\neq b$ with equal value give $f(\tfrac{a+b}2)<\tfrac12 f(a)+\tfrac12
f(b)=f(a)$, a contradiction. For $L$: $\hat\beta$ is the unique global minimizer
iff $X^\top X$ is PD, i.e. iff $X$ has full column rank (no collinearity), since
then $L$ is strictly convex.
**P.4** For $t\in[0,1]$, $A(t\beta_1+(1-t)\beta_2)+c=t(A\beta_1+c)+(1-t)(A\beta_2
+c)$, so $g(A\cdot+c)$ evaluated at the convex combination equals $g$ of the
convex combination of $A\beta_i+c$; convexity of $g$ then bounds it by the convex
combination of values, proving $\beta\mapsto g(A\beta+c)$ convex. Taking
$g(r)=\lVert r\rVert^2$ (convex), $A=-X$, $c=y$ gives $L(\beta)=g(-X\beta+y)$
convex with no derivatives.

**A.1** A central-difference checker matches the analytic gradient
$-2X^\top(y-X\beta)$ to $\approx10^{-5}$–$10^{-6}$ at random points (see the
chapter code cell for the pattern).
**A.2** As $\rho\to1$ the smallest eigenvalue of $X^\top X\to0$ and the
condition number $\to\infty$; fixed-step gradient descent makes almost no
progress along the flat valley, so the loss after ten steps stays far above its
minimum — the numerical face of collinearity from Chapter 1.
```

- [ ] **Step 2: Confirm the block is inside the show-solutions region**

Run:
```bash
python3 -c "
t=open('appendix/solutions.qmd').read()
i=t.index('content-visible when-meta=\"show-solutions\"')
j=t.index('Chapter 2 — Multivariable Calculus')
k=t.rindex(':::')
assert i < j < k, 'Chapter 2 block is outside the content-visible region'
print('solutions block correctly placed')
"
```
Expected: prints `solutions block correctly placed`.

- [ ] **Step 3: Commit**

```bash
git add appendix/solutions.qmd
git commit -m "ch2: appendix solutions for C/B/P/A exercises"
```

---

### Task 10: Full render verification (HTML + PDF)

**Files:**
- No edits expected; this task verifies the whole chapter builds.

- [ ] **Step 1: Confirm Quarto is installed**

Run: `quarto --version`
Expected: prints a version (e.g. `1.5.x`). If "command not found", install Quarto first (see Environment prerequisites) — this task cannot pass without it.

- [ ] **Step 2: Render the reader edition (HTML)**

Run: `quarto render --to html`
Expected: exits 0; no KaTeX `ParseError`, no `WARN` about undefined citation keys for `nocedal2006`/`boyd2004`, and `_book/parts/01-foundations/02-multivariable-calculus.html` exists.

- [ ] **Step 3: Render the instructor edition with solutions visible (HTML)**

Run: `quarto render --to html --metadata show-solutions:true`
Expected: exits 0; the rendered solutions page now contains the "Chapter 2 — Multivariable Calculus" heading. Verify:
`grep -q "Chapter 2 — Multivariable Calculus" _book/appendix/solutions.html && echo "solutions visible"`
Expected: prints `solutions visible`.

- [ ] **Step 4: Render PDF**

Run: `quarto render --to pdf`
Expected: exits 0; produces the book PDF with no LaTeX math errors. (Requires a TeX install; if unavailable in CI, document the skip and run locally before merge.)

- [ ] **Step 5: Final self-review against the spec**

Confirm each spec success criterion holds: through-line present in every theory subsection; `∇²L = 2XᵀX` continuity stated; four proofs in the body and NOT duplicated as exercises; four P exercises are the adjacent set (necessary condition, first-order convexity, global/unique, affine composition); code cell ran (Task 7); all four exercise tiers populated with matching solutions.

- [ ] **Step 6: Commit any freeze/render artifacts that should be tracked**

```bash
git add -A
git commit -m "ch2: render verification — HTML reader + instructor + PDF build clean" || echo "nothing to commit"
```

---

## Self-Review (performed during plan authoring)

- **Spec coverage:** Motivation (Task 1), four-rung theory ladder incl. all four proofs (Tasks 2–5), worked examples a/b (Task 6), code tie-in with numeric gradient check + eigenvalue/conditioning tie-back (Task 7), C/B/P/A with the adjacent P set (Task 8), appendix solutions (Task 9), clean HTML+PDF render (Task 10). All spec sections mapped.
- **Deferred item resolved:** `@nocedal2006` and `@boyd2004` already exist in `references.bib`; no bib task needed.
- **Placeholders:** none — every content step contains the actual prose, math, code, or exercises.
- **Consistency:** the no-½ loss convention (`∇L=−2Xᵀ(y−Xβ)`, `∇²L=2XᵀX`) is used identically across Tasks 1–9; `r(β)=y−Xβ`, `J=−X` consistent in Tasks 3–9; the P-tier proofs in Task 8 match their solutions in Task 9 one-for-one.
- **Environment caveat:** Quarto is absent from the authoring sandbox; render checks (Task 10) require installing it, while every intermediate task is verifiable without Quarto (delimiter checks + standalone Python execution).
