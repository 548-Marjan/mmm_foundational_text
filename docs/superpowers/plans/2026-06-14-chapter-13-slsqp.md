# Chapter 13 — Constrained Nonlinear Optimization & SLSQP Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 13 — the Part IV capstone: the general smooth constrained NLP and its KKT conditions, the root-finding view, Newton's method (with a quadratic-convergence proof), the keystone equivalence SQP-step = Newton-on-KKT = QP-subproblem, damped BFGS, the full SLSQP/Kraft internals (least-squares recast, active set, merit line search), the non-convex/multi-start reality, and the real Hill-response budget optimization that recovers Ch11/Ch12's exact equal-marginal optimum. Shipped as a PR whose CI `quarto render` (HTML **and PDF**) passes.

**Architecture:** Single `{python}` cell (NumPy/Matplotlib/SciPy `optimize.minimize(method="SLSQP")`). Three full proofs (Lagrange necessity, Newton quadratic convergence, the SQP=Newton-KKT=QP keystone). Fixed chapter template.

**Tech Stack:** Quarto `.qmd`, KaTeX, Python 3 (NumPy, Matplotlib, SciPy), `references.bib`.

---

## Conventions (enforce every task)

- **Worktree:** `~/.config/superpowers/worktrees/mmm_foundational_text/ch13-slsqp`. **Use explicit `git -C <worktree>` for every git command** (a bare `cd` may commit to main, where a hook blocks it). Identity set (`jlh530i`).
- **KaTeX:** `aligned` inside `$$…$$`; **never** bare `\begin{align}`; `$$` on own lines; even `$` count (`python3 -c "print(open(F).read().count(chr(36))%2)"` → 0). Inline math only in Summary.
- **PDF-SAFE LaTeX — CRITICAL this chapter:** it is matrix-heavy (KKT block systems). **NEVER use `\begin{psmallmatrix}`** (mathtools env undefined in CI LuaLaTeX/PDF — it silently passes HTML and breaks PDF, as in Ch11). Use `\begin{bmatrix}` or `\begin{pmatrix}` for ALL matrices/vectors; for an inline small matrix use `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]`. Check `grep -c psmallmatrix <file>` → 0 every task.
- **Chapter file:** `parts/04-optimization/03-slsqp.qmd`. Keep H1 `# Constrained Nonlinear Optimization & SLSQP` and the anchors italic line. Replace the stub body.
- **Citations:** `@nocedal2006` (primary), `@boyd2004`. No new bib keys. (Kraft's 1988 SLSQP report may be named in prose; attribute its theory to `@nocedal2006`.)
- **House rule:** describe the production budget optimizer generically; **never name PyMC-Marketing or any library.** `scipy.optimize.minimize(method="SLSQP")` is a general scientific tool and is fine to use/name.
- **Verification:** single `{python}` cell runs under `MPLBACKEND=Agg python3`; figures end `plt.show()`. Real gate: CI `quarto render` HTML+PDF green (watch BOTH — matrix-heavy ⇒ PDF risk).
- **Exemplars (read-only):** `parts/04-optimization/01-convexity.qmd` (Ch11 — KKT, Lagrangian, complementary slackness, equal-marginal capstone), `parts/04-optimization/02-linear-programming.qmd` (Ch12 — the boundary handing off to here, the (7.2,1.8)/λ number).

## Verified numeric anchors (NumPy/SciPy-checked)

- **WE1 — one SQP step.** $\min x_1^2 + x_2^2$ s.t. $x_1 + x_2 = 2$. Lagrangian $\mathcal{L}=x_1^2+x_2^2+\lambda(x_1+x_2-2)$, $\nabla^2_{xx}\mathcal{L}=2I$, $\nabla c=(1,1)$. From $x_0=(0,0),\lambda_0=0$: the Newton-KKT block system
$\begin{bmatrix}2&0&1\\0&2&1\\1&1&0\end{bmatrix}\begin{bmatrix}p_1\\p_2\\\lambda_1\end{bmatrix}=\begin{bmatrix}0\\0\\2\end{bmatrix}$
gives $p=(1,1)$, $\lambda_1=-2$, so $x_1=(1,1)$ — the exact optimum in **one step** (objective quadratic + constraint linear ⇒ QP subproblem is the problem). Optimal $f=2$.
- **WE3 — Newton digit-doubling.** $g(x)=x^2-2$, $x_{k+1}=\tfrac12(x_k+2/x_k)$, $x_0=1$: errors $|x_k-\sqrt2|=$ $4.14\text{e}{-1},\ 8.58\text{e}{-2},\ 2.45\text{e}{-3},\ 2.12\text{e}{-6},\ 1.59\text{e}{-12}$ — each $\approx C(\text{prev})^2$ with $C=1/(2\sqrt2)\approx0.354$.
- **Code (a) — concave budget, SLSQP.** $r_i(b)=a_i\sqrt b$, $a=(2,1)$, $B=9$; `minimize(method="SLSQP")` → $b^\star=(7.2,1.8)$ (to 4dp $7.2004,1.7996$), marginals $a_i/(2\sqrt{b_i})=(0.3727,0.3727)=\lambda$, response $6.7082$, ~7 iterations. **Matches Ch11/Ch12 exactly, no PWL.** ($\lambda=1/\sqrt{7.2}=0.3727$.)
- **Code (b) — S-curve multi-start.** Hill $r(b)=b^2/(3^2+b^2)$, two channels, $B=6$; multi-start SLSQP over the budget split: central starts → interior optimum $b=(3,3)$, total $1.0$; extreme starts ($\approx0,\approx6$) → corner, total $0.8$. Best $=1.0$; distinct optima $\{0.8,1.0\}$. (hill(3)+hill(3)=0.5+0.5=1.0; hill(6)=0.8.)

---

## Task 1: Front matter + Motivation

**Files:** Modify `parts/04-optimization/03-slsqp.qmd` (keep H1 + anchors line; write `## Motivation`).

- [ ] **Step 1: Write Motivation (~420–540 words).**
  - The Part IV capstone. Ch11 built the theory (KKT, local=global, equal-marginal capstone); Ch12 solved the *linear/piecewise-linear-concave* case exactly but hit a boundary: real response is **smooth and possibly non-concave** (S-shaped). This chapter solves *that* problem.
  - Driving question: *how do you actually compute the optimal budget when the response is a smooth curve?* Answer: treat the KKT conditions as a system of equations and solve them by **Newton's method** — which turns out to be the repeated solution of a **quadratic-programming subproblem** (Sequential Quadratic Programming). **SLSQP** is the production implementation.
  - Deliver on the promises: Ch11 said SLSQP is "a Newton iteration on the KKT equations"; Ch12 said it owns "the smooth, possibly non-concave world... launched from multiple starting points." This chapter proves both.
  - The payoff, stated: on concave response SLSQP recovers Ch11/Ch12's *exact* equal-marginal optimum (the same $(7.2,1.8)$, $\lambda=0.3727$) with no piecewise-linear approximation — the Lagrange multiplier it returns *is* $\lambda$; on S-shaped response it is run from many starts to find the global optimum.
  - Honesty: SLSQP finds a *local* KKT point — global only when the problem is convex (Ch11's keystone). Multi-start is the global strategy, not a guarantee.
  - Roadmap: KKT for the smooth NLP → the root-finding view → Newton & its quadratic convergence → the SQP keystone → BFGS → the SLSQP internals → the non-convex reality → the real budget optimization.
  - Position as the bridge to Part V: a fitted MMM gives the curves; Ch17 runs this solver.
- [ ] **Step 2: Lint + commit.** Even `$`, no bare align, no psmallmatrix. `git -C <WT> add … && git -C <WT> commit -m "feat(ch13): Motivation — solving the smooth constrained budget problem"`.

---

## Task 2: Theory rung 1 — The constrained NLP and KKT (Lagrange necessity proof)

**Files:** Modify chapter (begin `## Theory & Proofs`; rung 1).

- [ ] **Step 1: Open `## Theory & Proofs`** with a short ladder paragraph (NLP+KKT → root-finding → Newton+convergence → SQP keystone → BFGS → SLSQP internals → non-convex reality → budget capstone).

- [ ] **Step 2: Rung 1 — the smooth constrained program and its KKT conditions.**
  - The program: $\min_x f(x)$ s.t. $c_i(x)=0\ (i\in\mathcal{E})$, $c_i(x)\le0\ (i\in\mathcal{I})$, with $f,c_i$ smooth ($C^2$). The **Lagrangian** $\mathcal{L}(x,\lambda)=f(x)+\sum_i\lambda_i c_i(x)$.
  - **Constraint qualification (LICQ):** at a feasible $x^\star$, the gradients $\{\nabla c_i(x^\star):i\in\mathcal{E}\cup\mathcal{A}(x^\star)\}$ of active constraints are linearly independent (define active set $\mathcal{A}$).
  - **Theorem (first-order necessary / KKT).** If $x^\star$ is a local minimizer at which LICQ holds, there exist multipliers $\lambda^\star$ with: stationarity $\nabla_x\mathcal{L}(x^\star,\lambda^\star)=\nabla f(x^\star)+\sum_i\lambda_i^\star\nabla c_i(x^\star)=0$; primal feasibility; dual feasibility $\lambda_i^\star\ge0\ (i\in\mathcal{I})$; complementary slackness $\lambda_i^\star c_i(x^\star)=0$.
  - **Proof (∎) of the equality-constrained case** (the clean, illuminating core; inequality case stated under LICQ, cite `@nocedal2006`): at a local min subject to $c(x)=0$, consider any tangent direction $d$ with $\nabla c_i(x^\star)^\top d=0$ for all $i\in\mathcal{E}$. By LICQ and the implicit function theorem there is a feasible curve $x(t)$ with $x(0)=x^\star$, $x'(0)=d$; since $t=0$ minimizes $f(x(t))$, $\tfrac{d}{dt}f(x(t))|_0=\nabla f(x^\star)^\top d=0$. So $\nabla f(x^\star)\perp$ every tangent direction, i.e. $\nabla f(x^\star)\in\text{span}\{\nabla c_i(x^\star)\}$, giving $\nabla f(x^\star)=-\sum_i\lambda_i^\star\nabla c_i(x^\star)$ for some $\lambda^\star$. $\blacksquare$
  - **Second-order sufficient condition** (stated): if additionally the reduced Hessian $d^\top\nabla^2_{xx}\mathcal{L}\,d>0$ for all nonzero $d$ in the critical cone, $x^\star$ is a strict local minimizer.
  - **Tie to Ch11:** there KKT was proved *sufficient* for convex programs; here are the *necessary* conditions for smooth, possibly non-convex ones — the equations every solver targets.
- [ ] **Step 3: Lint + commit.** One `\blacksquare`. `-m "feat(ch13): Theory rung 1 — constrained NLP and KKT (Lagrange proof)"`.

---

## Task 3: Theory rungs 2–3 — Root-finding view & Newton's quadratic convergence (proof)

**Files:** Modify chapter (rungs 2–3).

- [ ] **Step 1: Rung 2 — the root-finding view.**
  - A KKT point is a **root of a nonlinear system**. For the equality-constrained problem, stack stationarity and feasibility:
$$
F(x,\lambda)=\begin{bmatrix}\nabla f(x)+\nabla c(x)^\top\lambda\\ c(x)\end{bmatrix}=0,
$$
where $\nabla c(x)$ is the Jacobian of the constraints. (Use `bmatrix`, never psmallmatrix.) Solving the NLP $\equiv$ finding a root of $F$. Inequalities are folded in by identifying the **active set** and treating active constraints as equalities (developed in the SLSQP rung).
  - This reframing — from "minimize" to "solve $F=0$" — is what makes **Newton's method** the natural engine; the rest of the chapter is consequences of this one move.

- [ ] **Step 2: Rung 3 — Newton's method and local quadratic convergence (full proof).**
  - Newton for $F(z)=0$ (here $z=(x,\lambda)$): linearize $F(z_k+\Delta)\approx F(z_k)+F'(z_k)\Delta$, set to $0$, step $z_{k+1}=z_k-F'(z_k)^{-1}F(z_k)$. ($F'$ is the Jacobian; for the KKT system it is the **KKT matrix** of the next rung.)
  - **Theorem (local quadratic convergence).** If $F\in C^2$, $F(z^\star)=0$, and $F'(z^\star)$ is nonsingular, then there is a neighborhood of $z^\star$ in which Newton's iterates satisfy $\|z_{k+1}-z^\star\|\le C\|z_k-z^\star\|^2$ for a constant $C$.
  - **Proof (∎).** Use `aligned`. Let $e_k=z_k-z^\star$. From $z_{k+1}-z^\star=e_k-F'(z_k)^{-1}F(z_k)=F'(z_k)^{-1}\big[F'(z_k)e_k-F(z_k)\big]$. Taylor-expand $F(z^\star)=0$ about $z_k$: $0=F(z_k)-F'(z_k)e_k+R_k$ with $\|R_k\|\le \tfrac12 L\|e_k\|^2$ ($L$ a Lipschitz bound on $F'$ near $z^\star$). Hence $F'(z_k)e_k-F(z_k)=R_k$, so $e_{k+1}=F'(z_k)^{-1}R_k$ and $\|e_{k+1}\|\le\|F'(z_k)^{-1}\|\cdot\tfrac12 L\|e_k\|^2$. By continuity $\|F'(z_k)^{-1}\|\le 2\|F'(z^\star)^{-1}\|$ near $z^\star$, giving $\|e_{k+1}\|\le C\|e_k\|^2$ with $C=L\|F'(z^\star)^{-1}\|$. $\blacksquare$
  - Significance: **digit-doubling** per iteration — the practical reason SQP/SLSQP is fast (forward-reference Worked Example 3). State that **quasi-Newton** (BFGS, Rung 5) trades the exact Jacobian for an approximation and achieves **superlinear** (not quite quadratic) convergence (`@nocedal2006`). Newton's convergence is assumed nowhere earlier in the book; this rung supplies it, and it is exactly the engine the chapter rides on.
- [ ] **Step 3: Lint + commit.** New `\blacksquare` (Newton). `-m "feat(ch13): Theory rungs 2-3 — root-finding view and Newton convergence (proof)"`.

---

## Task 4: Theory rung 4 — KEYSTONE: SQP = Newton-on-KKT = QP subproblem (proof)

**Files:** Modify chapter (rung 4, keystone).

- [ ] **Step 1: Rung 4 — Sequential Quadratic Programming.**
  - Apply Newton (Rung 3) to the KKT system $F(x,\lambda)=0$ (Rung 2). The Jacobian is the **KKT matrix**
$$
F'(x,\lambda)=\begin{bmatrix}\nabla^2_{xx}\mathcal{L}(x,\lambda) & \nabla c(x)^\top\\ \nabla c(x) & 0\end{bmatrix},
$$
so the Newton step $(p,\,\delta\lambda)$ solves, writing $\lambda_{k+1}=\lambda_k+\delta\lambda$,
$$
\begin{bmatrix}\nabla^2_{xx}\mathcal{L}_k & \nabla c_k^\top\\ \nabla c_k & 0\end{bmatrix}\begin{bmatrix}p\\ \lambda_{k+1}\end{bmatrix}=-\begin{bmatrix}\nabla f_k\\ c_k\end{bmatrix}.
$$
(Show the algebra that turns $\delta\lambda$ into $\lambda_{k+1}$: the top block is $\nabla^2_{xx}\mathcal{L}_k\,p+\nabla c_k^\top\lambda_{k+1}=-\nabla f_k$.)
  - **Theorem (SQP = Newton-on-KKT).** The Newton step above is exactly the primal–dual solution $(p,\lambda_{k+1})$ of the **quadratic-programming subproblem**
$$
\min_p\ \nabla f_k^\top p + \tfrac12 p^\top\nabla^2_{xx}\mathcal{L}_k\,p \quad\text{s.t.}\quad \nabla c_k\,p + c_k = 0\ (\mathcal{E}),\quad \nabla c_k\,p + c_k\le 0\ (\mathcal{I}).
$$
  - **Proof (∎).** Write the KKT conditions *of the QP*: its Lagrangian is $\nabla f_k^\top p+\tfrac12 p^\top\nabla^2_{xx}\mathcal{L}_k p+\mu^\top(\nabla c_k p+c_k)$; stationarity in $p$ gives $\nabla^2_{xx}\mathcal{L}_k\,p+\nabla f_k+\nabla c_k^\top\mu=0$, and the (linearized) constraints give $\nabla c_k p+c_k=0$. Identifying $\mu=\lambda_{k+1}$, this is *exactly* the block system above. So solving the QP subproblem $\equiv$ taking the Newton step on the KKT system. $\blacksquare$
  - This is the chapter's center and the redemption of Ch11/Ch12's promise: **SQP iterates by repeatedly minimizing a quadratic model of the Lagrangian subject to linearized constraints**, and that is a Newton iteration on the KKT equations — inheriting the quadratic convergence of Rung 3. Note the QP linearizes constraints and quadratically models the *Lagrangian* (not just $f$) — the constraint curvature enters through $\nabla^2_{xx}\mathcal{L}$.
- [ ] **Step 2: Lint + commit.** New `\blacksquare`. `-m "feat(ch13): Theory rung 4 — SQP = Newton-on-KKT = QP subproblem (keystone)"`.

---

## Task 5: Theory rungs 5–6a — Lagrangian Hessian/BFGS and the least-squares recast

**Files:** Modify chapter (rung 5 + rung 6 part (a)).

- [ ] **Step 1: Rung 5 — the Lagrangian Hessian and quasi-Newton (BFGS).**
  - Problem: $\nabla^2_{xx}\mathcal{L}_k$ is expensive (needs second derivatives of $f$ and every $c_i$) and may be **indefinite** away from a solution, which makes the QP subproblem nonconvex/ill-posed.
  - Fix: replace it with a **BFGS** approximation $B_k\succ0$, built from successive gradient differences $s_k=x_{k+1}-x_k$, $y_k=\nabla_x\mathcal{L}_{k+1}-\nabla_x\mathcal{L}_k$ via the BFGS update; **Powell's damping** modifies $y_k$ when $s_k^\top y_k$ is too small, guaranteeing $B_{k+1}\succ0$. Then every QP subproblem is a **convex** QP with a unique solution. State the BFGS update formula. **Superlinear convergence** stated (`@nocedal2006`).
  - Because $B_k\succ0$, the QP direction $p$ is a **descent direction** for a merit function (used by the line search in Rung 6c) — note this now.

- [ ] **Step 2: Rung 6(a) — the least-squares recast (the "LS" in SLSQP).**
  - With $B_k\succ0$, take the Cholesky factor $B_k=L_kL_k^\top$ (Ch1). The QP objective
$$
\tfrac12 p^\top B_k p+\nabla f_k^\top p=\tfrac12\big\|L_k^\top p + L_k^{-1}\nabla f_k\big\|^2 + \text{const},
$$
so minimizing the quadratic is minimizing a **Euclidean norm** — a linear **least-squares** objective in $p$. With the linearized constraints, the QP subproblem becomes a **linearly-constrained linear least-squares problem** (an LSEI/least-squares-with-equality-and-inequality problem, itself reducible to NNLS), solved by stable orthogonal factorizations **without ever forming $B_k^{-1}$ or the indefinite KKT matrix**. This numerically stable recast is the "Sequential **Least Squares**" of SLSQP (Kraft); attribute the algorithm to Kraft 1988, its theory to `@nocedal2006`. Show the completing-the-square algebra for the recast.
- [ ] **Step 3: Lint + commit.** `-m "feat(ch13): Theory rungs 5-6a — BFGS and the least-squares recast"`.

---

## Task 6: Theory rung 6b–c — Active set and merit-function globalization

**Files:** Modify chapter (rung 6 parts (b)–(c)).

- [ ] **Step 1: Rung 6(b) — the active-set strategy for inequalities.**
  - The QP subproblem has linearized **inequalities**; SLSQP solves it by an **active-set** method: maintain a **working set** of inequalities provisionally treated as equalities, solve the resulting equality-constrained least-squares QP, then check the QP multipliers — a negative multiplier on a working-set inequality means it should be **dropped** (dual infeasible), a violated inactive constraint should be **added**. Iterate until the working set is consistent with KKT (all multipliers $\ge0$, all constraints satisfied). This is complementary slackness (Ch11/Ch12) enforced combinatorially, per QP solve.

- [ ] **Step 2: Rung 6(c) — merit-function line search (global convergence).**
  - The full Newton/SQP step $p$ has quadratic convergence only **near** a solution; from a poor start it may increase infeasibility or cost. **Globalize** with a line search on a **merit function** that trades objective against constraint violation, e.g. the $\ell_1$ merit
$$
\phi(x;\mu)=f(x)+\mu\sum_{i\in\mathcal{E}}|c_i(x)|+\mu\sum_{i\in\mathcal{I}}\max\{0,c_i(x)\}.
$$
  - Take a step $x_{k+1}=x_k+\alpha p$ with $\alpha\in(0,1]$ chosen by backtracking to enforce sufficient decrease in $\phi$. **State the descent property:** for $\mu$ large enough (exceeding the largest multiplier magnitude), the SQP direction $p$ is a descent direction for $\phi$ (because $B_k\succ0$ and $p$ solves the linearized-constraint QP), so a step length exists that decreases $\phi$ — giving **global convergence to a KKT point from arbitrary starts**. Full proof cited (`@nocedal2006`); one sentence on the **Maratos effect** (the merit function can reject good steps near the solution, fixed by a second-order correction).
  - Synthesis: SLSQP = damped-BFGS SQP (quadratic model of the Lagrangian, Rung 4) + least-squares QP solve (Rung 6a) + active set (6b) + merit line search (6c) = **local quadratic/superlinear speed with global reliability**. This is the algorithm `scipy.optimize.minimize(method="SLSQP")` runs.
- [ ] **Step 3: Lint + commit.** `-m "feat(ch13): Theory rung 6b-c — active set and merit-function globalization"`.

---

## Task 7: Theory rungs 7–8 — The non-convex reality and the MMM budget capstone

**Files:** Modify chapter (rungs 7–8, completing Theory).

- [ ] **Step 1: Rung 7 — what SLSQP does and does not guarantee.**
  - SLSQP converges (globally, from any start) to a point satisfying the **KKT conditions** — a *local* optimum. By Ch11's keystone, a local KKT point is the **global** optimum **only if the problem is convex**. For a non-convex (S-shaped response) budget problem, several local KKT points exist; SLSQP from one start finds one, and *which* one depends on initialization.
  - The remedy is **multi-start**: run SLSQP from many initial allocations and keep the best KKT point. This is a heuristic for global optimality, not a certificate — the honest boundary Ch11/Ch12 drew. (One sentence: genuinely global methods — branch-and-bound, Bayesian optimization — exist but are out of scope; multi-start is the practical standard.)

- [ ] **Step 2: Rung 8 — the real budget optimization (capstone).**
  - The problem LP could only approximate, solved exactly on the **smooth** response.
  - **(a) Concave case (vindication).** $\max_b\sum_i r_i(b_i)$ s.t. $\mathbf{1}^\top b=B$, $b\ge0$ (+ caps), with $r_i$ smooth concave. SLSQP converges to the unique global optimum; its KKT stationarity is exactly Ch11's equal-marginal condition $r_i'(b_i^\star)=\lambda$, and the **Lagrange multiplier SLSQP returns on the budget constraint is $\lambda$**. For $r_i=a_i\sqrt b$, $a=(2,1)$, $B=9$: $b^\star=(7.2,1.8)$, $\lambda=0.3727$, response $6.7082$ — **identical** to Ch11's hand derivation and Ch12's refined-LP limit, now from the production solver with **no piecewise-linear approximation**. State this is the convergence of the whole Part IV arc: convex theory (Ch11) → LP duality (Ch12) → the SQP solver (Ch13), all yielding the same $\lambda$.
  - **(b) S-shaped case (why SLSQP).** Hill response with $n>1$ (convex toe), budget $B$. The objective is non-concave; multi-start SLSQP finds multiple local optima (interior vs corner) and reports the best — the case LP's "fill steepest first" and a single descent both get wrong. Concrete: $n=2$, $K=3$, two channels, $B=6$ → interior $(3,3)$ total $1.0$ beats corner total $0.8$; some starts land on the corner (forward-reference Code Tie-in).
  - **Bridge to Part V.** A *fitted* MMM posterior supplies the response curves $r_i$ (with uncertainty); Chapter 17 runs this very solver — often once per posterior draw — to allocate budget under uncertainty. End the Theory section: Part IV built the mathematics of optimization from convex foundations to the production algorithm; Part V turns it loose on a real, fitted model.
- [ ] **Step 3: Lint + commit.** `-m "feat(ch13): Theory rungs 7-8 — non-convex reality and the budget capstone"`.

---

## Task 8: Worked Examples

**Files:** Modify chapter (`## Worked Examples`). Reproduce verified anchors exactly. **No psmallmatrix — use bmatrix.**

- [ ] **Step 1: Example 1 — one SQP step by hand.** $\min x_1^2+x_2^2$ s.t. $x_1+x_2=2$. Show $\nabla f=(2x_1,2x_2)$, $\nabla^2_{xx}\mathcal{L}=2I$, $\nabla c=(1,1)$, $c=x_1+x_2-2$. From $x_0=(0,0),\lambda_0=0$ write the block system $\begin{bmatrix}2&0&1\\0&2&1\\1&1&0\end{bmatrix}\begin{bmatrix}p_1\\p_2\\\lambda_1\end{bmatrix}=\begin{bmatrix}0\\0\\2\end{bmatrix}$ (RHS $=-[\nabla f;c]=-[0,0,-2]$), solve $p=(1,1)$, $\lambda_1=-2$, so $x_1=(1,1)$ — the exact optimum in one step. Explain *why* one step suffices: the objective is quadratic and the constraint linear, so the QP subproblem *is* the original problem; SQP reduces to solving it once. Contrast: with a nonlinear constraint the linearization changes each iteration, needing several steps (shown in code).
- [ ] **Step 2: Example 2 — the concave budget optimum, exactly.** $r_i(b)=a_i\sqrt b$, $a=(2,1)$, $B=9$. KKT stationarity $r_i'(b_i)=\tfrac{a_i}{2\sqrt{b_i}}=\lambda$ ⇒ $b_i\propto a_i^2$ ⇒ $b^\star=(7.2,1.8)$, common $\lambda=\tfrac{2}{2\sqrt{7.2}}=0.3727$, response $2\sqrt{7.2}+\sqrt{1.8}=6.7082$. State that SLSQP returns precisely this allocation *and* this multiplier — the same number Ch11 derived by hand and Ch12 reached as a refined-LP dual price, now the exact smooth solve.
- [ ] **Step 3: Example 3 — Newton's quadratic convergence, concretely.** $g(x)=x^2-2$, Newton $x_{k+1}=\tfrac12(x_k+2/x_k)$, $x_0=1$. Tabulate $x_k$ and the error $|x_k-\sqrt2|$: $0.414,\ 0.0858,\ 0.00245,\ 2.12\text{e}{-6},\ 1.59\text{e}{-12}$ — the error roughly squares each step (digit-doubling), the concrete face of the Rung 3 theorem, and the reason SLSQP needs only a handful of iterations near the optimum.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch13): Worked Examples — one SQP step, concave budget, Newton convergence"`.

---

## Task 9: Code Tie-in

**Files:** Modify chapter (`## Code Tie-in`, single `{python}` cell, NumPy/Matplotlib/SciPy).

- [ ] **Step 1: Build the cell.** In order:
  1. Imports (`numpy`, `matplotlib.pyplot`, `from scipy.optimize import minimize`).
  2. **(a) Concave budget via SLSQP.** $r_i=a_i\sqrt b$, $a=(2,1)$, $B=9$; minimize $-(a\cdot\sqrt b)$ with `method="SLSQP"`, equality `sum(b)=B`, bounds `b>=0`. Print $b^\star$ (expect $(7.2,1.8)$), response (expect $6.7082$), `res.nit`. Compute the budget multiplier as the common marginal $a_i/(2\sqrt{b_i^\star})$ (expect $0.3727$) and confirm both channels' marginals agree (= $\lambda$). `assert np.allclose(b_star,[7.2,1.8],atol=1e-2)`, `assert abs(lam-0.3727)<1e-3`. Print a line: "SLSQP recovers Ch11/Ch12's optimum exactly, no PWL."
  3. **(b) S-curve multi-start.** Hill $r(b)=b^2/(9+b^2)$, two channels, $B=6$; objective over split $b_1\in[0,6]$, `minimize(neg, x0=[s], method="SLSQP", bounds=[(0,6)])` for ~12 starts in $[0.05,5.95]$; collect totals; print best and the distinct optima found ($\{0.8,1.0\}$). `assert` best $=$ max over starts and that $>1$ distinct optimum value appears (multi-start matters).
  4. **(c) Newton digit-doubling (from scratch).** Iterate $x_{k+1}=\tfrac12(x_k+2/x_k)$ from $x_0=1$; print a table of $x_k$ and error $|x_k-\sqrt2|$; `assert` the error at step 4 is $<10^{-10}$ (quadratic). Optionally a 2-D Newton on a small system to mirror the KKT view.
  5. **Figures (each `plt.show()`):** (a) Newton error vs iteration on a log scale (the quadratic plummet) optionally vs a linear-rate reference; (b) the S-curve total response over the budget split $b_1\in[0,6]$ with the multi-start landing points marked (interior max vs corners).
- [ ] **Step 2: Verify headless.** Extract to `/tmp/ch13_cell.py`; `cd <WT> && MPLBACKEND=Agg python3 /tmp/ch13_cell.py`; iterate until asserts pass and prints match. **Pin printed numbers into surrounding prose.**
- [ ] **Step 3: Lint + commit.** Even `$`; one ```{python} cell; no psmallmatrix in prose. `-m "feat(ch13): Code Tie-in — SLSQP budget solve, multi-start, Newton convergence"`.

---

## Task 10: Exercises (written directly by the controller)

**Files:** Modify chapter (`## Exercises`). Four tiers, self-contained. **No psmallmatrix.**

- [ ] **Step 1: Write C/B/P/A.**
  - **C (3):** why a KKT point found by SLSQP is global only for convex problems (Ch11 keystone) and what multi-start does about it; why SQP linearizes constraints but uses a quadratic model of the *Lagrangian* (not just $f$); what the merit-function line search buys over the raw Newton step.
  - **B (3):** take one SQP/Newton-KKT step on a small equality-constrained NLP by hand (fresh numbers; build and solve the block system); solve a 2-channel concave budget by the equal-marginal/KKT condition and identify $\lambda$; run two Newton iterations of a scalar root by hand and show the error squaring. **NumPy-verify each.**
  - **P (2–3):** prove the equality-constrained first-order (Lagrange) condition via the tangent-space argument; prove the SQP-step = QP-subproblem equivalence (write the QP's KKT and match the block system); prove Newton's local quadratic convergence (the Taylor-remainder bound).
  - **A (2–3):** solve a constrained budget problem with `scipy.optimize.minimize(method="SLSQP")` and verify the returned solution satisfies KKT (stationarity residual ≈ 0); demonstrate multi-start on an S-shaped problem and count distinct optima; implement a from-scratch Newton iteration and plot the quadratic error decay.
- [ ] **Step 2: Lint + commit.** `-m "feat(ch13): Exercises — C/B/P/A tiers"`.

---

## Task 11: Summary (written directly by the controller)

**Files:** Modify chapter (`## Summary`, auto-included; inline math only — no `$$`, no psmallmatrix).

- [ ] **Step 1: Key concepts + Key identities.**
  - **Concepts:** the smooth constrained NLP and KKT (necessary conditions, LICQ); the root-finding view; Newton's method and quadratic convergence; the SQP keystone (Newton-on-KKT = QP subproblem); the Lagrangian Hessian and damped BFGS (superlinear); the SLSQP internals (least-squares recast, active set, merit line search); local vs global and multi-start; the budget capstone where SLSQP recovers Ch11/Ch12's exact $\lambda$.
  - **Identities (inline):** Lagrangian $\mathcal{L}=f+\sum_i\lambda_i c_i$; KKT stationarity $\nabla f+\sum_i\lambda_i\nabla c_i=0$ with $\lambda_{\mathcal I}\ge0$, $\lambda_i c_i=0$; Newton step $z_{k+1}=z_k-F'(z_k)^{-1}F(z_k)$; quadratic convergence $\|e_{k+1}\|\le C\|e_k\|^2$; the SQP block system (write with inline `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]` or describe in words — keep inline); the QP subproblem $\min_p \nabla f^\top p+\tfrac12 p^\top\nabla^2_{xx}\mathcal{L}\,p$ s.t. linearized constraints; equal-marginal $r_i'(b_i^\star)=\lambda$. Note SLSQP returns the budget multiplier $=\lambda$.
- [ ] **Step 2: Lint + commit.** Inline math only. `-m "feat(ch13): Summary — key concepts and identities"`.

---

## Task 12: Appendix solutions block

**Files:** Modify `appendix/solutions.qmd` (insert before the final `:::`).

- [ ] **Step 1: Append `## Chapter 13 — Constrained Nonlinear Optimization & SLSQP`** with full worked solutions to every C/B/P/A exercise (NumPy/SciPy-verified; code cells use `minimize(method="SLSQP")` and from-scratch Newton), inside the `show-solutions` block, before the file's last `:::`. **No psmallmatrix** (use bmatrix). No inline links back.
- [ ] **Step 2: Lint + commit.** Even `$`; closes with single `:::`; Ch13 heading present (em-dash); `grep -c psmallmatrix` = 0; appendix code cells run headless. `-m "feat(ch13): appendix solutions for SLSQP"`.

---

## Task 13: Final review + PR + CI

- [ ] **Step 1: Chapter lint.** H1 `# Constrained Nonlinear Optimization & SLSQP`; six template H2 in order; no bare `\begin{align}`; even `$`; **`grep -c psmallmatrix` = 0**; `\blacksquare` count ≥3 (Lagrange, Newton, SQP keystone); citations resolve (`@nocedal2006`/`@boyd2004`).
- [ ] **Step 2: Appendix lint.** Even `$`; closes with single `:::`; Ch13 heading present; **no psmallmatrix**.
- [ ] **Step 3: Re-run chapter + appendix code cells headless** — asserts pass, prints match pinned numbers.
- [ ] **Step 4: Through-line check.** Ch11 KKT/keystone + Ch12 boundary paid off; SQP keystone delivered; concave budget recovers exact $\lambda=0.3727$; S-curve multi-start; bridge to Part V/Ch17.
- [ ] **Step 5: Push + PR.** `git -C <WT> push -u origin ch13-slsqp`; `gh pr create` (title "Chapter 13 — Constrained Nonlinear Optimization & SLSQP").
- [ ] **Step 6: Watch CI `quarto render` (HTML **and PDF**) to green** — matrix-heavy chapter, so the PDF stage is the real risk (psmallmatrix-class bugs surface only there). Confirm BOTH jobs green. Report PR link + status. Then the user merges.

---

## Self-Review (controller, before dispatching Task 1)

- **Spec coverage:** NLP+KKT ✓(T2) · Lagrange proof ✓(T2) · root-finding ✓(T3) · Newton+convergence proof ✓(T3) · SQP keystone proof ✓(T4) · BFGS ✓(T5) · LS recast ✓(T5) · active set ✓(T6) · merit line search ✓(T6) · non-convex/multi-start ✓(T7) · budget capstone ✓(T8) · 3 proofs ✓.
- **Placeholder scan:** none (PR body at T13).
- **Anchor consistency:** WE1 one-step $(1,1)/\lambda{=}-2$; WE3/Code(c) Newton errors to $1.6\text{e}{-12}$; Code(a) $(7.2,1.8)/\lambda{=}0.3727/6.7082$; Code(b) multi-start $\{0.8,1.0\}$ — all NumPy/SciPy-verified, used identically in Worked Examples (T8) and Code (T9); Exercises (T10) use fresh, separately-verified numbers.
- **PDF safety:** psmallmatrix banned and rechecked every task (this chapter has many matrices); CI watch covers HTML **and** PDF.
