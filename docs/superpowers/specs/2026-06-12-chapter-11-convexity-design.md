# Chapter 11 — Convexity Theory: Design Spec

## Placement & Role

- **File:** `parts/04-optimization/01-convexity.qmd`
- **Part:** IV — Optimization (this is the **opening** chapter of Part IV)
- **Anchor:** Boyd & Vandenberghe (`@boyd2004`); supporting `@nocedal2006`.
- **Predecessors it builds on:** Ch1 (linear algebra — Hessian, PSD matrices, eigenvalues), Ch2 (multivariable calculus — gradients, Hessians, Taylor), Ch5 (the Hill saturation curve $f(x;K,n)=x^n/(K^n+x^n)$, already canonical in the book).
- **Hands off to:** Ch12 (Linear Programming — applies KKT to polyhedra, develops LP duality) and Ch13 (SLSQP — the constrained algorithm that the MMM budget optimizer actually runs; arrives with its optimality theory already proven here). Part V's budget-optimization chapter is the ultimate payoff.
- **Why it leads Part IV:** convexity is the property that separates tractable optimization (a local optimum is the answer) from intractable optimization (a landscape of traps). The whole of Part IV — and the trustworthiness of any budget recommendation — rests on whether the problem is convex.

## Driving Question (the spine)

*Why is one optimization problem solvable and another a minefield?* The chapter answers: **convexity**. It builds convex analysis from sets to functions to constrained optimality, proves the two results that make optimization work — **a local minimum of a convex function is global** (keystone) and **at the optimal budget, marginal returns are equalized across channels** (capstone) — and then confronts the uncomfortable truth that real MMM response curves are *not* always convex, which is exactly why Part IV needs the machinery of Ch12–13.

## Scope — Comprehensive (Ch11 is the theory hub of Part IV)

### Theory rungs

1. **Convex sets.** Definition (the segment between any two points stays inside); examples that recur — hyperplanes, halfspaces, polyhedra, the norm ball, and the **budget feasible set** $\{b \in \mathbb{R}^n : b \ge 0,\ \mathbf{1}^\top b \le B\}$ (a simplex/polyhedron). Intersection of convex sets is convex (proof). The feasible set of the MMM budget problem is convex — this is what the later optimality theory will need.

2. **Convex functions.** Definition via the chord-above-graph inequality $f(\lambda x + (1-\lambda)y) \le \lambda f(x) + (1-\lambda)f(y)$; convex domain requirement; the **epigraph** characterization ($f$ convex $\iff$ $\text{epi}\,f$ is a convex set). **Concave** $= -$convex. Motivate with the saturation curve: a diminishing-returns response is *concave*, so maximizing it is a convex problem.

3. **First-order condition.** For differentiable $f$: convex $\iff f(y) \ge f(x) + \nabla f(x)^\top (y - x)$ for all $x, y$ — the graph lies above every tangent. **Full proof** (1-D via the secant-slope monotonicity, then the general case by restricting to the line through $x, y$). Geometric reading: the supporting hyperplane.

4. **Second-order condition.** For twice-differentiable $f$: convex $\iff \nabla^2 f(x) \succeq 0$ (Hessian positive semidefinite) on the domain. **Full proof** (Taylor with integral remainder / the 1-D restriction $g(t) = f(x + tv)$ and $g''(t) = v^\top \nabla^2 f\, v \ge 0$, reusing Ch1's PSD machinery). **Applied to the Hill curve:** compute $f''(x)$ for $f(x;K,n)$, show $f'' \le 0$ (concave) for $n=1$ on $x \ge 0$, and locate the inflection for $n>1$ (sets up the non-convexity rung).

5. **Operations that preserve convexity.** Nonnegative weighted sums, composition with an affine map, pointwise maximum/supremum, and the scalar-composition rules. **Proofs** for the load-bearing ones (nonneg sum, affine composition). Why this matters: it lets us certify the convexity of a *built-up* objective (a sum of per-channel responses, a regularized loss) without re-deriving from scratch — the negative log-posterior of Ch5's ridge model is convex by these rules.

6. **Jensen's inequality.** $f(\mathbb{E}[X]) \le \mathbb{E}[f(X)]$ for convex $f$ (reverse for concave). **Proof** (supporting-hyperplane / finite case by induction). MMM reading: because response is concave, **splitting spend across flights (averaging) beats lumping** — and the gap is the curvature, a quantitative statement of diminishing returns.

7. **KEYSTONE — local = global.** For a convex $f$ over a convex set, **every local minimum is a global minimum**, and for differentiable $f$ a stationary point $\nabla f(x^\star) = 0$ is sufficient for global optimality. **Full proof** (convex combination contradiction). This is the result that makes optimization *work*: it is why a downhill algorithm cannot be fooled.

8. **Constrained optimality — Lagrange & KKT.** The convex program $\min f(x)$ s.t. $g_i(x) \le 0,\ h_j(x) = 0$. The **Lagrangian**; **Lagrange multipliers** for equality constraints; the **KKT conditions** (stationarity, primal/dual feasibility, complementary slackness) and the theorem that, under a constraint qualification (Slater's condition, stated), KKT is **necessary and sufficient** for a global optimum of a convex problem. **Proof** of sufficiency for the convex case (the clean direction); necessity stated with Slater. This is the theory hub the rest of Part IV draws on.

9. **CAPSTONE — equal marginal returns.** The budget-allocation problem: $\max_{b \ge 0} \sum_i r_i(b_i)$ s.t. $\mathbf{1}^\top b = B$, with each response $r_i$ concave and increasing. **Derive via KKT** that at the optimum every funded channel satisfies $r_i'(b_i^\star) = \lambda$ (a common marginal return / shadow price of budget), and unfunded channels have $r_i'(0) \le \lambda$ — the **water-filling** solution. **Full proof** from Rung 8. This is the mathematical heart of MMM budget optimization and the exact condition Ch13's SLSQP converges to.

10. **The honest rung — when the MMM problem is *not* convex.** The Hill curve with exponent $n > 1$ is **S-shaped**: convex on a "toe" below the inflection point, concave above it. So a real response curve is not globally concave, the summed objective is **not concave**, and the budget problem is **non-convex** — the keystone's guarantee evaporates, multiple local optima appear, and a downhill method can land on a bad one. Explain what is lost (no global certificate), how practitioners cope (multi-start, good initialization, sometimes convex relaxations), and why this is precisely why Ch13 reaches for **SLSQP with multiple starts** rather than trusting a single descent. Bridge explicitly to Ch12 (LP, the convex special case where the geometry is fully understood) and Ch13.

### Proof inventory (the rigor the chapter is graded on)

Full theorem→proof (`\blacksquare`):
- **P-first-order** (Rung 3): convexity $\iff$ tangent-below-graph.
- **P-second-order** (Rung 4): convexity $\iff$ Hessian PSD.
- **P-Jensen** (Rung 6).
- **P-keystone** (Rung 7): local min $=$ global min.
- **P-KKT-sufficiency** (Rung 8): KKT $\Rightarrow$ global optimum for a convex program.
- **P-capstone** (Rung 9): equal-marginal-returns from KKT.

Slater/necessity of KKT and the full operations list are **stated precisely** with references, not all proved (standard, and proving every preservation rule would bloat the chapter).

## Running Example — Hill Saturation & the Two-Channel Budget

Continuity with the book's established saturation transform:

- **Response curve:** the Hill function $f(x;K,n) = \dfrac{x^n}{K^n + x^n}$ (already canonical from Ch5), scaled by a per-channel effectiveness. Concave on $x \ge 0$ for $n = 1$; S-shaped (convex toe) for $n > 1$.
- **Worked examples (hand-computable):**
  - (a) Convexity certificate via the Hessian: verify a simple built objective (e.g. a quadratic spend cost, or $-\log$ utility) is convex by sign of $f''$ / PSD Hessian, and locate the Hill inflection $x^\star = K\left(\frac{n-1}{n+1}\right)^{1/n}$ for $n > 1$.
  - (b) Equal-marginal-returns allocation by hand for **two channels** with concave responses (use tractable concave forms, e.g. $r_i(b) = a_i\sqrt{b}$ or $a_i\log(1 + b/k_i)$, so $r_i'(b) = \lambda$ solves in closed form); compute the split, the common marginal $\lambda$, and confirm it beats an equal/naive split.
  - (c) Show an $n>1$ Hill curve has an inflection (compute it) and argue the two-channel S-curve budget problem has multiple KKT points.
- **Code tie-in (single runnable `{python}` cell, NumPy/Matplotlib + SciPy):**
  - Plot concave ($n\le1$) vs S-shaped ($n>1$) Hill curves; mark the inflection.
  - Numerically verify convexity (Hessian eigenvalues $\ge 0$ / chord test) on a concave objective.
  - Solve the **concave** two-channel (or $n$-channel) budget allocation under $\mathbf{1}^\top b = B$ and confirm the optimizer equalizes marginal returns $r_i'(b_i^\star)$ across funded channels (the capstone, numerically).
  - Demonstrate **non-convexity**: on the S-shaped ($n>1$) problem, run a local optimizer from several starts and show it converges to *different* optima — the multi-start motivation for Ch13. All numeric claims NumPy/SciPy-verified before authoring; figures end `plt.show()`.

## Template Conformance

Six H2 sections in order: **Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary** (Summary auto-included, inline math only). Exercises C/B/P/A, self-contained (no inline solution links). Appendix solutions appended to `appendix/solutions.qmd` under `## Chapter 11 — Convexity Theory` (em-dash), gated by `::: {.content-visible when-meta="show-solutions"}`, inserted before the file's final `:::`.

## KaTeX / Build Conventions

- `aligned` inside `$$…$$`; never bare `\begin{align}`. `$$` delimiters on their own lines. Even `$` count.
- Single `{python}` cell runs top-to-bottom under `MPLBACKEND=Agg python3`; SciPy is available (in `requirements.txt`).
- Real gate: CI `quarto render` (HTML + PDF) green on the PR. Quarto not installed locally.

## Citations

- Existing and used here: `@boyd2004` (convex optimization — the anchor), `@nocedal2006` (numerical optimization — KKT, algorithms), `@gelman2013` where the convex-posterior link is made. No new bib keys required (Boyd & Vandenberghe and Nocedal & Wright cover sets, functions, conditions, Lagrange/KKT, duality).

## Out of Scope (deliberately)

- LP duality and the simplex method — that is Ch12.
- The SLSQP algorithm internals (sequential quadratic programming, line search) — that is Ch13.
- Subgradients / non-smooth convex analysis beyond a one-sentence mention (the book's objectives are smooth).
- Conic/semidefinite programming — out of the book's scope.
