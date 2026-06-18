# Chapter 23 — CS Foundations: Design Spec

**Date:** 2026-06-18
**Part:** VII — Software Engineering & Computer Science (first chapter; the computational grounding)
**File:** `parts/07-swe-cs/01-cs-foundations.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter opens Part VII and grounds every computation the book has run on its actual substrate:
finite-precision arithmetic executed at finite cost. The matrix solves behind regression (Chapter 4)
and Gaussian posteriors (Chapters 5–6), the per-sample work of the samplers (Chapters 8–10), the
optimizers (Chapters 12–14), and the Kalman / prior-store recursions (Chapters 17, 22) all rest on
floating-point and on algorithmic cost. This chapter makes both precise.

**Center of gravity (locked, via AskUserQuestion):** floating-point and conditioning. The keystone is
**$\kappa(X^\top X) = \kappa(X)^2$** — forming and solving the normal equations squares the
conditioning, doubling the digits lost — with algorithmic complexity (big-O, the cost of the linear-
algebra and sampling primitives) as the supporting pillar.

**Unifying payoff:** the **Chapter 16 identifiability ridge is also a numerical conditioning problem.**
Near-collinear design columns make $\sigma_{\min}(X) \approx 0$, hence $\kappa(X)$ huge; the statistical
identifiability ridge *is* numerical ill-conditioning, and the very direction that calibration
(Chapter 21) compresses is the direction in which floating-point silently loses precision.

**Driving question (locked):** *"Your model fit, sampler, and optimizer all run on finite-precision
arithmetic at finite cost. When does that silently corrupt the answer — and what does it cost to
scale?"*

**Scope note (locked):** stays tightly anchored to the MMM computational stack — linear-algebra solves,
samplers, the ridge. NOT a general CLRS survey: sorting / data-structure theory gets at most a passing
mention in Rung 1.

---

## 2. Scope discipline

- Anchors the computation already done in Parts I–VI; does not introduce new MMM modeling.
- **House rule (critical):** describe everything as generic technique. **Never** name PyMC-Marketing or
  any MMM/PPL/sampler library. `numpy`/`scipy`/`matplotlib` only.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their own lines;
  even `$$` count.
- Wall-clock timing is **illustrative only** — never hard-asserted (BLAS threading makes it noisy);
  assert the analytical flop-count ratio instead.

---

## 3. Theory & Proofs — the rungs (5; keystone at Rung 4)

**Rung 1 — The cost model and asymptotic complexity (supporting pillar).**
The RAM cost model and the big-O / $\Theta$ / $\Omega$ hierarchy. The cost of the primitives the MMM
stack runs: matrix–vector $O(n^2)$, matrix–matrix $O(n^3)$, **Cholesky factorization $O(n^3/3)$**,
triangular solve $O(n^2)$, and the per-sample cost of MCMC (each leapfrog step a gradient evaluation,
so a run is $O(\text{draws} \times \text{leapfrog} \times \text{grad cost})$). **Proof (Cholesky flop
count):** the $k$-th column of the factorization costs $\approx k$ multiply–adds across $\approx k$
rows, so the total is $\sum_{k=1}^{n} k^2 \approx n^3/3$. State the MMM reading: a Gaussian posterior
solve over $d$ channels is $O(d^3)$, multiplied by the number of draws when sampled. Sorting / data-
structure theory mentioned in one sentence as the broader CLRS context, not developed.

**Rung 2 — The IEEE-754 model and the rounding bound.**
A double is $\pm\,m \cdot 2^{e}$ with a 52-bit mantissa; the gap between $1$ and the next
representable number is machine epsilon $\varepsilon_{\text{mach}} = 2^{-52} \approx 2.22\times10^{-16}$,
and the **unit roundoff** is $u = 2^{-53} \approx 1.11\times10^{-16}$. **The fundamental axiom of
floating-point:** for any basic operation, $\mathrm{fl}(x \mathbin{op} y) = (x \mathbin{op} y)(1+\delta)$
with $|\delta| \le u$. **Proof:** rounding to the nearest representable number incurs an absolute error
at most half the spacing at that magnitude, so the relative error is at most $u$. Consequence: a single
operation is accurate to relative error $u$, but errors compose over many operations.

**Rung 3 — Catastrophic cancellation.**
Subtracting two nearly-equal numbers explodes the *relative* error: small absolute rounding errors in
the operands, negligible relative to the operands, become large relative to the small difference.
**Derivation:** if $\hat x = x(1+\delta_x)$ and $\hat y = y(1+\delta_y)$ with $|\delta_\cdot| \le u$,
then the relative error of $\hat x - \hat y$ is bounded by $\frac{|x|+|y|}{|x-y|}\,u$, which blows up as
$x \to y$. **Anchor:** the naive one-pass variance $\overline{x^2} - \bar x^2$ versus the two-pass
$\overline{(x-\bar x)^2}$ — catastrophic when the mean is large relative to the spread.

**Rung 4 — Proof P (KEYSTONE): the condition number and $\kappa(X^\top X) = \kappa(X)^2$.**
Define the (2-norm) condition number $\kappa(A) = \lVert A\rVert\,\lVert A^{-1}\rVert =
\sigma_{\max}/\sigma_{\min}$. State the forward-error bound: the relative error in solving $Ax=b$ is
bounded by $\kappa(A)$ times the relative perturbation, so $\log_{10}\kappa$ digits are lost.
**Theorem:** $\kappa(X^\top X) = \kappa(X)^2$ in the 2-norm. **Proof:** with the SVD $X = U\Sigma V^\top$,
$X^\top X = V\Sigma^2 V^\top$, so the singular values of $X^\top X$ are $\sigma_i(X)^2$ and
$\kappa(X^\top X) = \sigma_{\max}^2/\sigma_{\min}^2 = \kappa(X)^2$. Therefore forming and solving the
normal equations squares the conditioning and doubles the digits lost; solve least squares via QR or
SVD (conditioning $\kappa(X)$) rather than $X^\top X$ (conditioning $\kappa(X)^2$). $\blacksquare$
**MMM tie (the unifying observation):** the Chapter 16 ridge is near-collinear spend columns, i.e.
$\sigma_{\min}(X) \approx 0$ and $\kappa(X)$ enormous; the statistical identifiability ridge *is*
numerical ill-conditioning. The eigendirection the calibration factor (Chapter 21) compresses is
exactly the direction in which a naive normal-equations solve loses all its digits.

**Rung 5 — Backward stability and the practical upshot (synthesis).**
Define backward stability: a computed solution is the exact solution of a slightly perturbed problem.
Cholesky and QR solves are backward stable, so the achievable accuracy is $\sim \kappa(A)\,u$ — you
cannot beat the conditioning, but a stable algorithm does not make it worse. State the practical rules
the rest of the stack relies on: never invert a matrix (solve instead), never form $X^\top X$ (use QR),
sum stably (Kahan / pairwise), prefer the two-pass variance. Close forward to Part VII: software
architecture, data engineering, and testing are the engineering that protect these numerical and cost
properties at scale.

---

## 4. Worked Examples

**WE1 — Catastrophic cancellation in the variance (numeric).**
Data $x = 10^8 + \{1, 2, 3\}$, whose population variance is $\tfrac{2}{3} \approx 0.667$. The naive
one-pass formula $\overline{x^2} - \bar x^2$ returns **exactly $0.0$** in double precision — total
cancellation, because the variance lies below the resolution of numbers near $10^{16}$. The two-pass
formula $\overline{(x-\bar x)^2}$ returns $0.667$ exactly. Lesson: algebraic equivalence is not
numerical equivalence; subtract small quantities before, not after, squaring.

**WE2 — $\kappa(X^\top X) = \kappa(X)^2$ and the digits lost (the keystone).**
A near-collinear design $X$ (two almost-parallel columns) has $\kappa(X) \approx 1.0\times10^{6}$ and
$\kappa(X^\top X) \approx 1.05\times10^{12}$ — the ratio is $1.0000$, the squaring exactly. Solving for
a known coefficient $\beta = (2, -1)$: the normal-equations solve has error $\approx 7.9\times10^{-5}$
while the QR/least-squares solve has error $\approx 7.2\times10^{-11}$ — a factor of $\approx 10^{6}$,
exactly the extra $\kappa(X)$ that forming $X^\top X$ introduced. This is the ridge as a numerical
object: the worse the identification, the more digits a naive solve discards.

**WE3 — The $O(n^3)$ cost of the Gaussian solve (supporting pillar).**
The cost of factorizing a $d \times d$ posterior precision grows as $n^3/3$ flops. Tabulate the
analytical flop counts at growing $n$ and confirm the $8\times$ growth per doubling; show wall-clock
timing alongside as illustration (noting it is noisy at small $n$ because of threaded BLAS). MMM
reading: doubling the channel count octuples the per-solve cost, and sampling multiplies that by the
draw count — why the dimension of the model is a first-order engineering constraint.

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + `matplotlib`; verified headless). It:
1. **Machine epsilon (Rung 2):** print $\varepsilon_{\text{mach}}$ and $u$; assert $0.1 + 0.2 \ne 0.3$,
   $1 + u = 1$, and $1 + \varepsilon_{\text{mach}} \ne 1$.
2. **Catastrophic cancellation (WE1):** naive vs two-pass variance on $x = 10^8 + \{1,2,3\}$; assert the
   naive result is $0.0$ and the two-pass is within $10^{-9}$ of $2/3$.
3. **Conditioning keystone (WE2):** build the near-collinear $X$; compute $\kappa(X)$, $\kappa(X^\top
   X)$; assert their ratio $\approx 1$ (squaring); solve via normal equations and via QR against a known
   $\beta$ and assert the QR error is orders of magnitude smaller. Figure: solve error vs $\kappa(X)$
   over a sweep, normal-equations ($\sim\kappa^2 u$) vs QR ($\sim\kappa u$).
4. **Complexity (WE3):** analytical $n^3/3$ flop counts at growing $n$; assert the doubling ratio
   $\approx 8$; show wall-clock timing as an illustrative figure (no hard timing assert). Figure:
   cost vs $n$ on log–log with a slope-3 reference line.
Every non-timing numeric claim asserted; seed `np.random.default_rng(23)`.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why should you never invert a matrix or form $X^\top X$ to solve least squares? What is machine
  epsilon, and why is catastrophic cancellation about *relative* error? In what sense is the Chapter 16
  ridge simultaneously a statistical and a numerical object?
- **B:** Compute the unit roundoff for double precision; for a small $2\times2$ near-collinear $X$,
  compute $\kappa(X)$ and $\kappa(X^\top X)$ from the singular values and the digits lost by each;
  evaluate naive vs two-pass variance on a small large-mean dataset by hand.
- **P:** (i) Prove $\kappa(X^\top X) = \kappa(X)^2$ via the SVD. (ii) Prove the floating-point rounding
  bound $\mathrm{fl}(x\mathbin{op}y) = (x\mathbin{op}y)(1+\delta)$, $|\delta|\le u$, and derive the
  catastrophic-cancellation relative-error bound $\frac{|x|+|y|}{|x-y|}u$.
- **A:** Extend the Code Tie-in — sweep the conditioning of $X$ and plot the realized solve error of the
  normal-equations and QR routes against $\kappa(X)$, confirming the $\kappa^2 u$ and $\kappa u$ slopes;
  then add a calibration-style rank-one precision increment (Chapter 21) and show it lowers $\kappa$ and
  the realized error together — the ridge healed numerically as well as statistically.

---

## 7. Appendix solutions

Append `## Chapter 23 — CS Foundations` to `appendix/solutions.qmd`, **in chapter order** (after the
Chapter 22 block), inside the existing `content-visible` gated div. Full C/B/P/A; the P-block carries
both proofs (the SVD proof of $\kappa(X^\top X)=\kappa(X)^2$; the rounding bound and the cancellation
relative-error bound).

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities: unit
roundoff $u = 2^{-53}$; the floating-point axiom $\mathrm{fl}(x\mathbin{op}y) = (x\mathbin{op}y)(1+\delta)$,
$|\delta|\le u$; cancellation bound $\frac{|x|+|y|}{|x-y|}u$; condition number
$\kappa(A)=\sigma_{\max}/\sigma_{\min}$ and the forward-error bound; the keystone
$\kappa(X^\top X)=\kappa(X)^2$; Cholesky cost $\approx n^3/3$. Close tying back to Chapter 16 (the ridge
as ill-conditioning) and forward to Part VII (architecture / data engineering / testing protect these
properties at scale).

---

## 9. Cross-references

- **Chapter 1 (Linear Algebra):** SVD and singular values reused for the conditioning proof.
- **Chapter 4 (Linear Regression):** the normal equations whose conditioning is squared.
- **Chapter 16 (Building & Fitting):** the Fisher-information ridge re-read as ill-conditioning.
- **Chapter 21 (Advanced Calibration):** the rank-one precision increment that lowers $\kappa$.
- **Part VII (Chapters 24–26):** architecture, data engineering, testing as protecting cost and
  numerical correctness at scale.

---

## 10. Bibliography

All keys already present in `references.bib` — no additions:
- `@clrs2009` — Cormen, Leiserson, Rivest & Stein, *Introduction to Algorithms* (complexity).
- `@goldberg1991` — Goldberg, "What Every Computer Scientist Should Know About Floating-Point
  Arithmetic" (the rounding model, cancellation).
- `@trefethen1997` — Trefethen & Bau, *Numerical Linear Algebra* (conditioning, backward stability).

---

## 11. Verification (the real gate)

- Every numeric anchor NumPy-verified (done: $u = 1.11\times10^{-16}$; naive variance $=0.0$ vs two-pass
  $0.667$; $\kappa(X)\approx10^6$, $\kappa(X^\top X)\approx10^{12}$, ratio $1.0000$; normal-eqns error
  $7.9\times10^{-5}$ vs QR $7.2\times10^{-11}$; Cholesky $n^3/3$ flop ratio $8\times$ per doubling).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`); no wall-clock asserts.
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no `\begin{align}`,
  no `psmallmatrix`, valid citation keys, no banned library names.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
