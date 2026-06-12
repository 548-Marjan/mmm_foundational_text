# Chapter 12 — Linear Programming: Design Spec

## Placement & Role

- **File:** `parts/04-optimization/02-linear-programming.qmd`
- **Part:** IV — Optimization (the **middle** chapter: Ch11 Convexity → **Ch12 Linear Programming** → Ch13 Constrained Optimization / SLSQP)
- **Anchors:** Bertsimas & Tsitsiklis (`@bertsimas1997`, primary); Boyd & Vandenberghe (`@boyd2004`); Nocedal & Wright (`@nocedal2006`).
- **Builds on:** Ch11 (convex sets/polyhedra, the Fundamental theorem's reliance on convexity, KKT and complementary slackness, the equal-marginal-returns capstone, Lagrangian duality, Slater/strong duality), Ch1 (linear algebra — bases, rank, solving $Ax=b$).
- **Hands off to:** Ch13 (SLSQP) — the smooth nonlinear solver for when the response is genuinely curved/non-concave and the LP/PWL approximation no longer suffices.
- **Ch11's set-up to honor (verbatim intent):** "Chapter 12 takes $f$ and every $g_i$ to be **affine**: the feasible set becomes a polyhedron, complementary slackness becomes the active-set logic of vertices, and the KKT system collapses into the elegant symmetry of **linear-programming duality** ... the global optimum sits at a **vertex** and can be found with certainty."

## Driving Question (the spine)

*What happens to the optimization theory of Ch11 when everything is linear?* The answer is the most complete and tractable corner of optimization: the feasible set is a polyhedron, an optimum sits at a **vertex**, the **simplex** method walks vertex to vertex to find it, and Ch11's Lagrangian duality sharpens into **LP duality** — where the dual variables are **shadow prices**, and the shadow price on the budget is exactly the equal-marginal-return $\lambda$ of Ch11's capstone, now produced by an algorithm instead of by hand.

**Pedagogical arc (locked):** teach LP *on its own terms first* (standard form → polyhedral geometry → Fundamental Theorem → simplex → duality), **then** reveal the MMM payoff — the budget problem, under a **piecewise-linear concave** approximation of the saturation curve, *is* an LP — so $\lambda$ reappears as a dual price. Do NOT lead with the PWL trick.

## Scope — Comprehensive

### Theory rungs

1. **LP standard form.** $\min c^\top x$ s.t. $Ax = b,\ x \ge 0$ (and the inequality/canonical form $\min c^\top x$ s.t. $Ax \ge b,\ x\ge 0$); converting via **slack/surplus variables**; the feasible set is a polyhedron (cite Ch11). Why "linear program," and that maximization is $-\min(-c)$.

2. **Polyhedral geometry.** Vertices / **extreme points**; **basic feasible solutions (BFS)** — choose $m$ linearly independent columns of $A$ (a basis), solve for the basic variables, set the rest to zero, require $\ge 0$. **Theorem + proof:** $x$ is a vertex of the feasible polyhedron $\iff$ $x$ is a basic feasible solution. (The bridge from geometry to algebra that makes simplex possible.)

3. **Fundamental Theorem of LP.** **Theorem:** if an LP in standard form has an optimal solution, then it has an optimal solution at a vertex (BFS); if it has a feasible solution, it has a basic feasible one. **Full proof** (a non-vertex optimum lies in the relative interior of a face / on a segment between feasible points; moving along the edge direction in which the objective does not increase reaches a vertex with no worse value — finitely many vertices ⇒ an optimal vertex exists). State the trichotomy: an LP is infeasible, unbounded, or has a finite optimum attained at a vertex.

4. **The simplex method.** The algorithm as the constructive engine of the Fundamental Theorem: start at a BFS (vertex), compute **reduced costs**, if all are $\ge 0$ the vertex is optimal, else pick an **entering variable** (negative reduced cost) and a **leaving variable** via the **min-ratio test**, pivot to an adjacent BFS with no worse cost; repeat. Geometric reading: walking along edges of the polyhedron, downhill, to a corner. State termination (finiteness of vertices) and the degeneracy/cycling caveat (Bland's rule, stated). One **worked pivot** inline (full tableau in Worked Examples). Note complexity reality: exponential worst case, excellent in practice; mention interior-point as the polynomial alternative in one sentence.

5. **LP duality — construction & weak duality.** Every primal LP has a **dual** LP (give the standard primal–dual pair: primal $\min c^\top x$ s.t. $Ax \ge b, x\ge 0$ ↔ dual $\max b^\top y$ s.t. $A^\top y \le c, y \ge 0$). Derive it as the Lagrangian dual of Ch11 (the dual variables are the constraint multipliers). **Theorem (weak duality) + proof:** any dual-feasible $y$ lower-bounds any primal-feasible $x$: $b^\top y \le c^\top x$. Corollary: equal values ⇒ both optimal; primal unbounded ⇒ dual infeasible.

6. **Strong duality & complementary slackness.** **Theorem (strong duality):** if the primal has a finite optimum, so does the dual, and their optimal values are equal (zero gap). Prove it as the LP specialization of Ch11's KKT/Slater strong-duality result (or via simplex termination producing a dual-feasible $y$ from the optimal basis with matching value — give the cleaner of the two; cite `@bertsimas1997` for the full Farkas-based proof if abbreviated). **Theorem (complementary slackness):** primal-feasible $x$ and dual-feasible $y$ are both optimal $\iff$ for every $i$, $y_i(a_i^\top x - b_i) = 0$ and $x_j(c_j - (A^\top y)_j)=0$ — the LP face of Ch11's $\lambda_i g_i = 0$. **Proof.** This is the rung that closes the loop with Ch11.

7. **Shadow prices — the economic reading.** The optimal dual variable $y_i^\star$ is the **shadow price** of constraint $i$: the marginal change in the optimal objective per unit relaxation of $b_i$ (state the sensitivity result; tie to strong duality $\partial p^\star/\partial b_i = y_i^\star$ within a basis). A binding constraint has a positive shadow price; a slack one has zero (complementary slackness). For the budget constraint, the shadow price is the marginal response per extra dollar — **exactly Ch11's $\lambda$**.

8. **The MMM payoff — budget allocation as an LP (PWL bridge).** Approximate each channel's concave saturation response $r_i$ by a **piecewise-linear concave** function: split spend $b_i = \sum_k s_{ik}$ into segments $s_{ik}\in[0, w_{ik}]$ with **decreasing** per-segment slopes $m_{i1} > m_{i2} > \cdots$ (concavity). Then $\max \sum_{i,k} m_{ik} s_{ik}$ s.t. $\sum_{i,k} s_{ik} = B$, segment caps $0\le s_{ik}\le w_{ik}$, plus **planning constraints** (per-channel min/max spend) is a genuine **LP**. Because slopes decrease, simplex (or the dual) **fills the steepest segments first** — the diminishing-returns logic — and the **dual price on the budget equals the common marginal slope at the optimum: Ch11's equal-marginal-returns $\lambda$, rediscovered by LP duality.** Layered caps/minimums contribute their own shadow prices (the marginal value of relaxing a cap). This is the chapter's keystone application.

9. **The honest rung — where LP stops.** PWL-concave-as-LP is an *approximation*: refining the segmentation converges to the smooth optimum (state the error shrinks with segment width). But it rests on **concavity** — the trick works only because decreasing slopes let "fill steepest first" be optimal. For an **S-shaped** response (Ch11's $n>1$ Hill, convex toe), a PWL fit has *increasing-then-decreasing* slopes, "fill steepest first" is wrong, and the problem is **not** an LP — it becomes an integer/combinatorial or genuinely nonlinear program. That is the boundary: LP owns the linear/concave-separable world exactly; the smooth, possibly non-concave world is **Ch13's SLSQP**. Bridge explicitly.

### Proof inventory (`\blacksquare`)

- **P2** — vertex $\iff$ basic feasible solution.
- **P3** — Fundamental Theorem of LP (optimum at a vertex).
- **P5** — weak duality.
- **P6** — strong duality (LP specialization of Ch11 KKT) **and** complementary slackness.

Sensitivity/shadow-price result (Rung 7) stated precisely with a short derivation, not a full theorem; simplex termination/Bland's rule stated, not proved.

## Running Example — A Small LP, its Dual, and the PWL Budget Allocation

- **Worked Example 1 — a 2-variable LP by hand.** A small media LP (two channels, a budget and one cap) solved graphically: draw the feasible polygon, evaluate the objective at each vertex, identify the optimal vertex — the Fundamental Theorem made concrete. Then **one full simplex pivot sequence** on its standard form (slack variables, reduced costs, entering/leaving, pivot to the optimal BFS), confirming the graphical answer.
- **Worked Example 2 — the dual and complementary slackness.** Write the dual of Example 1's LP; solve it; verify **strong duality** (primal value = dual value) and read **complementary slackness** (which constraints bind, which dual prices are zero). Interpret the budget shadow price economically.
- **Worked Example 3 — PWL-concave budget allocation.** Two channels, each concave response approximated by 2–3 segments with decreasing slopes; set up the LP; solve (by the steepest-first logic, by hand); show the optimal allocation, that the **budget dual price = the common marginal slope = Ch11's $\lambda$**, and compare to a naive split. All numbers NumPy/SciPy-verified before authoring.
- **Code Tie-in (single `{python}` cell, NumPy/Matplotlib/SciPy):**
  - Build and solve a media LP with `scipy.optimize.linprog` (HiGHS); print the optimal vertex and objective.
  - Read the **dual values** (`result.ineqlin.marginals` / `eqlin.marginals`) and confirm they are the shadow prices; verify strong duality (primal = dual objective) and complementary slackness numerically.
  - Solve the **PWL-concave budget allocation** as an LP; show the budget dual price equals the common marginal slope; **refine the PWL** (more segments) and show the LP optimum converges to the smooth concave optimum from Ch11 (reusing the Hill/sqrt response), with the dual price → the smooth $\lambda$.
  - (Optional) show the S-curve case: a PWL fit of an $n>1$ Hill has non-monotone slopes, "fill steepest first" gives the wrong answer, motivating Ch13. Assertions pin the key equalities; figures end `plt.show()`.

## Template Conformance

Six H2 sections in order: **Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary** (Summary auto-included, inline math only). Exercises C/B/P/A, self-contained. Appendix solutions under `## Chapter 12 — Linear Programming` (em-dash), gated by `::: {.content-visible when-meta="show-solutions"}`, before the file's final `:::`.

## KaTeX / Build Conventions

- `aligned` inside `$$…$$`; never bare `\begin{align}`; `$$` on their own lines; even `$` count.
- **No `psmallmatrix`** (mathtools env absent in the PDF/LuaLaTeX build — it passed HTML but broke Ch11's PDF). Use `pmatrix`/`bmatrix`, or `smallmatrix` wrapped in `\bigl[ ... \bigr]`.
- Single `{python}` cell runs under `MPLBACKEND=Agg python3`; SciPy available (`scipy.optimize.linprog`, HiGHS).
- Real gate: CI `quarto render` (HTML **and PDF**) green on the PR. Quarto not installed locally — the PDF stage only runs in CI, so watch it.

## Citations

`@bertsimas1997` (LP theory — standard form, simplex, duality, the Farkas/strong-duality proofs), `@boyd2004` (duality as convex-program specialization), `@nocedal2006` (algorithms). No new bib keys.

## Out of Scope (deliberately)

- Integer / mixed-integer programming and branch-and-bound — one-sentence mention at the non-convex boundary (Rung 9), not developed.
- Network-flow / transportation LPs as a separate topic — the media LP is the only worked structure.
- Interior-point/barrier method internals — one-sentence mention (Rung 4) as the polynomial alternative; the smooth-KKT view is Ch13's.
- Sensitivity-analysis ranging beyond the basic shadow-price reading.
