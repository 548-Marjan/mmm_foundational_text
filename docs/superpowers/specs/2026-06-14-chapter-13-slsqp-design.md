# Chapter 13 — Constrained Nonlinear Optimization & SLSQP: Design Spec

## Placement & Role

- **File:** `parts/04-optimization/03-slsqp.qmd`
- **Part:** IV — Optimization (the **final** chapter and capstone: Ch11 Convexity → Ch12 Linear Programming → **Ch13 SLSQP**)
- **Anchor:** Nocedal & Wright (`@nocedal2006`, primary); Boyd & Vandenberghe (`@boyd2004`).
- **Builds on:** Ch11 (KKT conditions, Lagrangian, complementary slackness, the local=global keystone, the equal-marginal-returns capstone, the non-convexity/multi-start motivation), Ch12 (the smooth/non-concave boundary where LP stops; KKT as the thing solvers target), Ch2 (gradients, Hessians, Taylor — Newton named but its convergence never proved), Ch1 (linear systems, Cholesky, least squares).
- **Hands off to:** Part V, Ch17 (Budget Optimization) — which applies this solver to a *fitted* MMM posterior; and more broadly the production budget optimizer of an MMM stack (described generically — **never name PyMC-Marketing or any library**, per the book's house rule).
- **What Ch11/Ch12 promised (must deliver):** Ch11 — "SLSQP attacks the KKT system iteratively... forms a quadratic model of the Lagrangian and a linearization of the constraints, solves that QP for a search direction... a Newton iteration on the KKT equations"; and "multi-start SLSQP" for non-convex S-curves. Ch12 — "the smooth, possibly non-concave world... is the domain of Chapter 13's SLSQP."

## Driving Question (the spine)

*How do you actually solve the budget problem when the response is smooth and possibly non-concave — the case LP could only approximate?* The answer is **Sequential Quadratic Programming**: treat the KKT conditions as a nonlinear system and solve it by **Newton's method**, which turns out to be exactly the repeated solution of a **quadratic-programming subproblem** (a quadratic model of the Lagrangian under linearized constraints). **SLSQP** is the production implementation — Kraft's algorithm — that recasts each QP as a linear least-squares problem, manages inequality constraints by an active set, and globalizes with a merit-function line search. Run on the real Hill-response budget problem, it recovers Ch11's exact equal-marginal optimum on concave response (no piecewise-linear approximation), and on S-shaped response it is launched from multiple starts to find the global optimum.

## Scope — Comprehensive (Part IV capstone; full Kraft internals + Newton convergence proof)

### Theory rungs

1. **The constrained NLP and the KKT conditions.** General smooth program $\min f(x)$ s.t. $c_i(x)=0\ (i\in\mathcal{E})$, $c_i(x)\le0\ (i\in\mathcal{I})$. The Lagrangian $\mathcal{L}(x,\lambda)=f(x)+\sum_i\lambda_i c_i(x)$. **First-order necessary (KKT) conditions** under a constraint qualification (**LICQ** — active constraint gradients linearly independent): stationarity $\nabla_x\mathcal{L}=0$, primal feasibility, dual feasibility $\lambda_{\mathcal{I}}\ge0$, complementary slackness. **Proof** of the equality-constrained first-order condition $\nabla f=-\sum\lambda_i\nabla c_i$ (Lagrange, via the tangent-space/implicit-function argument); the inequality case stated under LICQ with reference. Second-order sufficient condition (reduced Hessian of $\mathcal{L}$ positive definite on the critical cone) — stated. Tie to Ch11: there we proved KKT *sufficient* for convex programs; here are the *necessary* conditions for smooth (possibly non-convex) ones.

2. **The root-finding view.** A KKT point is a root of a nonlinear system. For the equality-constrained problem, $F(x,\lambda)=\begin{psmallmatrix}\nabla f(x)+\nabla c(x)^\top\lambda\\ c(x)\end{psmallmatrix}=0$ (NOTE: write with `\begin{bmatrix}` or `\begin{pmatrix}`, **never** `psmallmatrix` — PDF build). Solving the NLP $\equiv$ finding a root of $F$. Inequalities are handled by identifying the *active set* (next rungs). This reframing is the bridge from "optimize" to "solve equations," which is what makes Newton applicable.

3. **Newton's method and its quadratic convergence (proved).** Newton for a nonlinear system $F(z)=0$: $z_{k+1}=z_k-F'(z_k)^{-1}F(z_k)$. **Theorem (local quadratic convergence) + full proof:** if $F$ is $C^2$, $F(z^\star)=0$, and $F'(z^\star)$ is nonsingular, then for $z_0$ near $z^\star$, $\|z_{k+1}-z^\star\|\le C\|z_k-z^\star\|^2$ (Taylor expansion of $F$ about $z_k$, bound the remainder, use $F'$ invertibility). This is the engine behind SQP's speed — the digit-doubling per iteration — and it is *not proved anywhere else in the book*, so this rung closes that gap and underwrites the whole chapter. Quasi-Newton **superlinear** convergence stated for later.

4. **KEYSTONE — Sequential Quadratic Programming.** Apply Newton's method to the KKT system. **Theorem:** the Newton step for the KKT system is identical to the solution $(p,\lambda_{k+1})$ of the **QP subproblem**
$$
\min_p\ \nabla f(x_k)^\top p + \tfrac12 p^\top \nabla^2_{xx}\mathcal{L}(x_k,\lambda_k)\,p \quad\text{s.t.}\quad \nabla c_i(x_k)^\top p + c_i(x_k) = 0\ (i\in\mathcal{E}),\ \le 0\ (i\in\mathcal{I}).
$$
**Full proof** (write the KKT/Newton block system $\begin{bmatrix}\nabla^2_{xx}\mathcal{L} & \nabla c^\top\\ \nabla c & 0\end{bmatrix}\begin{bmatrix}p\\\lambda_{k+1}\end{bmatrix}=-\begin{bmatrix}\nabla f\\ c\end{bmatrix}$ and recognize it as the KKT system *of the QP*). This is exactly what Ch11 and Ch12 promised — "a quadratic model of the Lagrangian and a linearization of the constraints," "a Newton iteration on the KKT equations." The chapter's intellectual center.

5. **The Lagrangian Hessian and quasi-Newton.** $\nabla^2_{xx}\mathcal{L}$ is expensive and may be indefinite (bad for the QP). Replace it with a **BFGS** approximation $B_k\succ0$, updated from gradient differences with **damping** (Powell) to preserve positive-definiteness, so each QP subproblem is a well-posed **convex** QP with a unique solution. **Superlinear convergence** stated (`@nocedal2006`). Why $B_k\succ0$ also guarantees the QP direction is a descent direction for a merit function (sets up globalization).

6. **SLSQP internals (Kraft) — full treatment.** The production algorithm.
   - **(a) The least-squares recast (the "LS").** With $B_k\succ0$, factor $B_k=L_kL_k^\top$ (Cholesky, Ch1); the QP objective $\tfrac12 p^\top B_k p + \nabla f^\top p$ becomes $\tfrac12\|L_k^\top p + L_k^{-1}\nabla f\|^2 + \text{const}$, so the equality-constrained QP is a **linear least-squares problem with linear constraints** (LSEI/NNLS), solved stably without forming $B_k^{-1}$. This is the "Sequential **Least Squares**" in SLSQP.
   - **(b) Active-set strategy.** Inequalities are handled by maintaining a working set of constraints treated as active (equalities); the QP is solved over that set, multipliers checked for sign (dual feasibility), constraints added/dropped until KKT-consistent.
   - **(c) Merit-function line search (globalization).** From a poor start the full Newton step may not reduce infeasibility+cost; an $\ell_1$ **merit function** $\phi(x;\mu)=f(x)+\mu\sum_i|c_i(x)|_+$ measures progress, and a backtracking line search along the QP direction $p$ guarantees sufficient decrease, giving **global convergence** to a KKT point from arbitrary starts (state the descent property; mention the Maratos effect / second-order correction in one sentence). The result: local quadratic/superlinear speed *and* global reliability.

7. **The non-convexity reality.** SLSQP converges to a **local** KKT point — by Ch11's keystone, that is the global optimum *only if the problem is convex*. For a non-convex (S-shaped) budget problem, several local optima exist; SLSQP from one start finds one of them, and the answer depends on initialization. The remedy is **multi-start**: run from many initial allocations, keep the best KKT point. Honest statement of what the method does and does not guarantee (no global certificate off the convex world — the very boundary Ch11/Ch12 drew).

8. **MMM capstone — the real budget optimization.** The problem LP could only approximate, now solved exactly on the *smooth* response.
   - **(a) Concave case (vindication).** Smooth response $r_i(b)=a_i\sqrt b$ (or a concave Hill), budget $\mathbf{1}^\top b=B$, bounds $b\ge0$ (+ caps). SLSQP converges to the unique global optimum, and the returned **Lagrange multiplier on the budget equals Ch11's $\lambda$** — recovering $(7.2,1.8)$, $\lambda=0.3727$, response $6.7082$ **exactly**, with no piecewise-linear approximation. The equal-marginal-returns capstone, solved by the production algorithm.
   - **(b) S-shaped case (why we needed SLSQP).** Hill response with $n>1$ (convex toe), budget $B$. Multi-start SLSQP finds multiple local optima (e.g. interior vs corner), and the best is reported — the case LP and a single descent both mishandle.
   - Bridge to Part V: a *fitted* MMM posterior gives the response curves; Ch17 runs this solver (often over posterior draws) to allocate budget under uncertainty. The whole of Part IV — convex theory, LP duality, SQP — converges here on the algorithm a real MMM budget optimizer runs.

### Proof inventory (`\blacksquare`)

- **P1** — equality-constrained first-order (Lagrange) necessary conditions.
- **P3** — Newton's local quadratic convergence (the gap-closing proof; Q-decision B).
- **P4** — SQP step = Newton-on-KKT = QP subproblem (keystone).

Stated precisely with citation (not fully proved): the inequality KKT conditions under LICQ, second-order sufficient conditions, BFGS superlinear convergence, the merit-function global-convergence theorem (descent property shown, full proof cited).

## Running Example — One SQP step, the concave budget solve, and Newton digit-doubling

- **Worked Example 1 — one SQP/Newton-KKT step by hand.** A small equality-constrained NLP; form the QP subproblem (gradient, Lagrangian-Hessian or $B_0$, linearized constraint), solve the KKT block system, take the step. Use a case where the structure is transparent (e.g. a quadratic objective with a linear constraint, where SQP converges in **one** step — illustrating that the QP subproblem *is* the problem when it is already a QP — plus a nonlinear-constraint case taking a couple of steps).
- **Worked Example 2 — the concave budget optimum, exactly.** $r_i(b)=a_i\sqrt b$, $a=(2,1)$, $B=9$: write the KKT/equal-marginal condition $r_i'(b_i)=\lambda$, solve $b^\star=(7.2,1.8)$, $\lambda=0.3727$, response $6.7082$ — and note SLSQP returns precisely this, multiplier and all (the Ch11/Ch12 number, now from the real solver).
- **Worked Example 3 — Newton's quadratic convergence numerically.** A scalar root (e.g. solving $g(x)=0$) showing the error squaring each step (digit-doubling), making the Rung-3 theorem concrete.
- **Code Tie-in (single `{python}` cell, NumPy/Matplotlib/SciPy).** Uses `scipy.optimize.minimize(method="SLSQP")` — a general scientific solver, **not** a marketing library:
  - **(a)** Solve the concave budget problem; assert SLSQP recovers $(7.2,1.8)$ and that the budget multiplier (`res` KKT / finite-difference) $\approx0.3727$ — matching Ch11/Ch12 exactly, *no PWL*. Print iteration count.
  - **(b)** Solve the S-shaped (Hill $n=2$) budget problem from many starts; show several local optima and that multi-start finds the best (assert best $\ge$ any single start).
  - **(c)** A from-scratch Newton iteration on a scalar/2-D root showing the error squaring each step (the digit-doubling table), tying to the Rung-3 proof.
  - Figures (each `plt.show()`): convergence-error log plot (Newton vs linear); the S-curve objective over the budget split with the multi-start landing points. All numerics asserted; pin into prose.

## Template Conformance

Six H2 in order: **Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary** (Summary inline math only). Exercises C/B/P/A, self-contained. Appendix solutions under `## Chapter 13 — Constrained Nonlinear Optimization & SLSQP` (em-dash), gated by `show-solutions`, before the file's final `:::`.

## KaTeX / Build Conventions

- `aligned` inside `$$…$$`; never bare `\begin{align}`; `$$` on own lines; even `$` count.
- **No `\begin{psmallmatrix}`** (mathtools env undefined in the CI LuaLaTeX/PDF build — it broke Ch11's PDF; matrices appear often in this chapter, so use `\begin{bmatrix}`/`\begin{pmatrix}` throughout, or `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]`).
- Single `{python}` cell runs under `MPLBACKEND=Agg python3`; SciPy `scipy.optimize.minimize(method="SLSQP")` available.
- Real gate: CI `quarto render` (HTML **and PDF**) green on the PR — watch BOTH (matrix-heavy chapter ⇒ PDF risk).

## Citations

`@nocedal2006` (primary — KKT, Newton, SQP, BFGS, line search, merit functions, the SLSQP/SQP theory), `@boyd2004` (convexity/KKT cross-reference). No new bib keys. (Kraft's SLSQP report may be cited in prose by name without a bib entry, or note it as described in `@nocedal2006`.)

## Out of Scope (deliberately)

- Interior-point/barrier methods for NLP — one-sentence mention (the alternative globalization), not developed; Ch12 already placed interior-point for LP.
- Trust-region SQP variants — one-sentence mention vs the line-search globalization used here.
- Derivative-free / global optimizers (genetic, Bayesian optimization) — out of scope; multi-start is the global strategy taught.
- Automatic differentiation — note gradients are required and how they're obtained (analytic/finite-difference), not a full treatment.
