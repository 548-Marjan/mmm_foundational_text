# Chapter 12 — Linear Programming Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 12 — Linear Programming, the middle chapter of Part IV: LP standard form, polyhedral geometry, the Fundamental Theorem (vertex optimality), the simplex method, and LP duality (weak/strong + complementary slackness) as the specialization of Ch11's KKT — then the MMM payoff where a piecewise-linear-concave budget problem is an LP whose budget dual price is Ch11's equal-marginal-return λ. Shipped as a PR whose CI `quarto render` (HTML **and PDF**) passes.

**Architecture:** Single `{python}` cell (NumPy/Matplotlib/SciPy `linprog`/HiGHS). Four full proofs (vertex⇔BFS, Fundamental Theorem, weak duality, strong duality + complementary slackness). Fixed chapter template.

**Tech Stack:** Quarto `.qmd`, KaTeX, Python 3 (NumPy, Matplotlib, SciPy), `references.bib`.

---

## Conventions (enforce every task)

- **Worktree:** `~/.config/superpowers/worktrees/mmm_foundational_text/ch12-linear-programming`. **Use explicit `git -C <worktree>` for every git command** (a bare `cd` may commit to main, where a hook blocks it). Identity set (`jlh530i`).
- **KaTeX:** `aligned` inside `$$…$$`; **never** bare `\begin{align}`; `$$` on their own lines; even `$` count (`python3 -c "print(open(F).read().count(chr(36))%2)"` → 0). Inline math only in Summary.
- **PDF-safe LaTeX:** **do NOT use `\begin{psmallmatrix}`** — that mathtools env is undefined in the CI LuaLaTeX build (it broke Ch11's PDF despite passing HTML). Use `pmatrix`/`bmatrix`, or `\bigl[\begin{smallmatrix}…\end{smallmatrix}\bigr]`.
- **Chapter file:** `parts/04-optimization/02-linear-programming.qmd`. Keep H1 `# Linear Programming` and the anchors italic line. Replace the stub body.
- **Citations:** `@bertsimas1997` (primary), `@boyd2004`, `@nocedal2006`. No new bib keys.
- **Verification:** single `{python}` cell runs under `MPLBACKEND=Agg python3`; figures end `plt.show()`. SciPy `linprog(method="highs")` available. Real gate: CI `quarto render` HTML+PDF green (PDF stage runs only in CI — watch it).
- **Exemplars (read-only):** `parts/04-optimization/01-convexity.qmd` (Ch11 — duality/KKT set-up, rung/proof voice), the Ch11 capstone (equal-marginal-returns λ) which this chapter rediscovers.

## Verified numeric anchors (SciPy-checked)

- **WE1 LP:** $\max 3x_1+2x_2$ s.t. $x_1+x_2\le4,\ x_1\le2,\ x\ge0$. Optimal vertex $x^\star=(2,2)$, value $10$. Vertices: $(0,0)\!\to\!0,(2,0)\!\to\!6,(2,2)\!\to\!10,(0,4)\!\to\!8$. **Dual:** $\min 4y_1+2y_2$ s.t. $y_1+y_2\ge3,\ y_1\ge2,\ y\ge0$ → $y^\star=(2,1)$, value $10$ (**strong duality**). Shadow prices: budget $y_1=2$, cap $y_2=1$ (verified: +1 budget → +2 response; +1 cap → +1 response). Both constraints bind ⇒ both dual prices positive (complementary slackness).
- **WE1 simplex:** slacks $s_1,s_2$: $x_1+x_2+s_1=4,\ x_1+s_2=2$. Start BFS $(0,0,4,2)$. Pivot 1: $x_1$ enters (reduced cost), min-ratio $\min(4/1,2/1)=2$ ⇒ $s_2$ leaves, BFS $(2,0,2,0)$, value $6$. Pivot 2: $x_2$ enters, min-ratio on $s_1$ row ⇒ $s_1$ leaves, BFS $(2,2,0,0)$, value $10$ — optimal (all reduced costs $\ge0$).
- **WE3 PWL allocation** (B=5): channel 1 segment slopes $[5,3,1]$, channel 2 $[4,2,1]$, each width $2$. Fill steepest first: ch1-seg1 ($5$, 2 units), ch2-seg1 ($4$, 2 units), ch1-seg2 ($3$, 1 unit) — budget exhausted. $b^\star=(3,2)$, total response $5\cdot2+3\cdot1+4\cdot2=21$, **budget shadow price $\lambda=3$** (slope of the partially filled segment). All filled segments have slope $\ge3$, all unfilled $\le3$ — the PWL equal-marginal/water-filling condition.
- **PWL → smooth convergence** (ties to Ch11): $r_i(b)=a_i\sqrt b$, $a=(2,1)$, $B=9$, smooth optimum $b^\star=(7.2,1.8)$, $\lambda=1/\sqrt{7.2}=0.3727$, response $6.7082$. PWL LP with $K$ secant segments per channel over $[0,9]$: $K{=}3\to b=(6,3),\lambda=0.367$; $K{=}10\to b=(7.2,1.8),\lambda=0.362$; $K{=}200\to b=(7.2,1.8),\lambda=0.372,\text{resp }6.7082$. Allocation converges fast; the dual price (a derivative) trails toward $0.3727$.

---

## Task 1: Front matter + Motivation

**Files:** Modify `parts/04-optimization/02-linear-programming.qmd` (keep H1 + anchors line; write `## Motivation`).

- [ ] **Step 1: Write Motivation (~400–520 words).**
  - Pick up Ch11's hand-off: when the objective and constraints are all **affine**, Ch11's convex theory becomes the most complete corner of optimization — linear programming. The feasible set is a **polyhedron**, an optimum sits at a **vertex**, and Ch11's Lagrangian duality sharpens into **LP duality**.
  - Driving question: *what does optimization look like when everything is linear?* Answer: you can find the global optimum **with certainty** (a vertex), an algorithm (**simplex**) walks there, and every constraint earns a **shadow price** — the marginal value of relaxing it.
  - The MMM hook, stated but deferred: a budget problem with a **piecewise-linear** response is an LP, and its budget shadow price is exactly the equal-marginal-return $\lambda$ Ch11 derived by hand — so this chapter rediscovers Ch11's capstone through duality. (Promise the PWL bridge; don't open with it.)
  - Honesty up front: LP owns the *linear/concave-separable* world exactly; smooth and non-concave response is Ch13's job.
  - Roadmap: standard form → polyhedral geometry → Fundamental Theorem → simplex → duality → shadow prices → the PWL budget LP → the boundary where LP stops.
- [ ] **Step 2: Lint + commit.** `git -C <WT> add … && git -C <WT> commit -m "feat(ch12): Motivation — optimization when everything is linear"`.

---

## Task 2: Theory rungs 1–2 — Standard form & polyhedral geometry (vertex⇔BFS proof)

**Files:** Modify chapter (begin `## Theory & Proofs`; rungs 1–2).

- [ ] **Step 1: Open `## Theory & Proofs`** with a short ladder paragraph (standard form → geometry → Fundamental Theorem → simplex → duality → shadow prices → PWL payoff → boundary).

- [ ] **Step 2: Rung 1 — LP standard form.**
  - Standard form $\min c^\top x$ s.t. $Ax=b,\ x\ge0$; canonical inequality form $\min c^\top x$ s.t. $Ax\ge b,\ x\ge0$. Maximization is $-\min(-c^\top x)$.
  - **Slack/surplus variables** convert inequalities to equalities (e.g. $a^\top x\le b \Rightarrow a^\top x + s = b,\ s\ge0$); every LP has a standard form.
  - The feasible set $\{x:Ax=b,x\ge0\}$ is a polyhedron (intersection of hyperplanes and the nonnegative orthant — cite Ch11 Rung 1).

- [ ] **Step 3: Rung 2 — Polyhedral geometry: vertices and basic feasible solutions (proof).**
  - **Vertex / extreme point:** a feasible $x$ that is not a convex combination of two other feasible points.
  - **Basic feasible solution (BFS):** for $Ax=b$ with $A\in\mathbb{R}^{m\times n}$ rank $m$, pick $m$ linearly independent columns (a **basis** $B$), set the $n-m$ nonbasic variables to $0$, solve $x_B=A_B^{-1}b$; if $x_B\ge0$ it is a BFS.
  - **Theorem.** $x$ is a vertex of $\{x:Ax=b,x\ge0\}$ $\iff$ $x$ is a basic feasible solution.
  - **Proof (∎, both directions).** (BFS ⇒ vertex) if a BFS $x$ were $\tfrac12(u+v)$ for feasible $u\ne v$, then $u,v$ agree with $x$ ($=0$) on nonbasic coordinates (nonnegativity forces the averaged zeros to be zero in both), so $u,v$ have the same support as $x$ on the basic columns; $A_B$ invertible ⇒ $u_B=v_B=A_B^{-1}b$ ⇒ $u=v$, contradiction. (vertex ⇒ BFS) if $x$ is feasible with support $S=\{j:x_j>0\}$ and the columns $\{A_j:j\in S\}$ were linearly dependent, a nonzero $d$ supported on $S$ with $Ad=0$ lets $x\pm\varepsilon d$ stay feasible for small $\varepsilon$, exhibiting $x$ as a midpoint — so a vertex has independent active columns, which extend to a basis, making $x$ a BFS. $\blacksquare$
  - Significance: **finitely many** bases ⇒ finitely many vertices; the next rung says one of them is optimal, and simplex searches among them.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch12): Theory rungs 1-2 — standard form and vertices (BFS proof)"`.

---

## Task 3: Theory rungs 3–4 — Fundamental Theorem (proof) & the simplex method

**Files:** Modify chapter (rungs 3–4).

- [ ] **Step 1: Rung 3 — Fundamental Theorem of LP (full proof).**
  - **Theorem.** For an LP in standard form: (i) if it has a feasible solution it has a basic feasible one; (ii) if it has an optimal solution it has one at a vertex (BFS). Trichotomy: every LP is **infeasible**, **unbounded**, or has a **finite optimum attained at a vertex**.
  - **Proof (∎).** Take an optimal $x^\star$ with minimal support. If its positive columns are dependent, build a feasible direction $d$ ($Ad=0$, supported on the support) along which the objective $c^\top d$ is (WLOG) $\le0$; move $x^\star+\theta d$ until a coordinate hits $0$ — the new point is feasible, no worse (if $c^\top d<0$ strictly better, contradicting optimality, so $c^\top d=0$ and it is also optimal) and has smaller support. Repeat until the support columns are independent — a BFS, hence a vertex. Finitely many vertices ⇒ comparing objective values gives an optimal vertex. (If the objective decreases without bound along a feasible ray, the LP is unbounded.) $\blacksquare$
  - This is *why* LP is tractable: the search reduces from a continuum to finitely many vertices.

- [ ] **Step 2: Rung 4 — The simplex method.**
  - The constructive engine of the Fundamental Theorem: start at a BFS; express the objective in terms of nonbasic variables to get **reduced costs** $\bar c_N = c_N - A_N^\top (A_B^{-1})^\top c_B$; if all $\bar c_N\ge0$, the vertex is optimal (no improving direction); otherwise pick an **entering** variable with $\bar c_j<0$, and a **leaving** variable by the **min-ratio test** $\min_i \{(A_B^{-1}b)_i/(A_B^{-1}A_j)_i : (\cdot)_i>0\}$ to stay feasible; **pivot** to the adjacent BFS. Repeat.
  - Geometric reading: walking along **edges** of the polyhedron, downhill, corner to corner, until no edge descends.
  - Termination: finitely many vertices ⇒ terminates if cost strictly decreases; **degeneracy** can stall (cycling), prevented by **Bland's rule** or lexicographic pivoting (stated, cite `@bertsimas1997`).
  - Complexity reality (one sentence): exponential worst case (Klee–Minty), excellent in practice; **interior-point/barrier** methods give polynomial guarantees and are the modern default for large LPs — and their smooth-KKT view foreshadows Ch13.
  - Include **one worked pivot step** inline (full sequence is Worked Example 1): show the first pivot of WE1.
- [ ] **Step 3: Lint + commit.** One new `\blacksquare` (Rung 3). `-m "feat(ch12): Theory rungs 3-4 — Fundamental Theorem (proof) and simplex"`.

---

## Task 4: Theory rungs 5–6 — LP duality (weak & strong duality + complementary slackness)

**Files:** Modify chapter (rungs 5–6).

- [ ] **Step 1: Rung 5 — Duality construction & weak duality (proof).**
  - The **primal–dual pair**: primal $\min c^\top x$ s.t. $Ax\ge b,\ x\ge0$ ↔ dual $\max b^\top y$ s.t. $A^\top y\le c,\ y\ge0$. Derive the dual as Ch11's **Lagrangian dual**: $L(x,y)=c^\top x - y^\top(Ax-b)$, $\min_{x\ge0}L$ finite only when $A^\top y\le c$, giving the dual objective $b^\top y$. The dual variables ARE the constraint multipliers of Ch11.
  - **Theorem (weak duality) + proof (∎).** For primal-feasible $x$ and dual-feasible $y$: $b^\top y \le (Ax)^\top y = x^\top(A^\top y) \le c^\top x$, using $y\ge0,\ Ax\ge b$ for the first and $x\ge0,\ A^\top y\le c$ for the last. So every dual value lower-bounds every primal value. Corollaries: if $b^\top y=c^\top x$ both are optimal; if the primal is unbounded below the dual is infeasible.

- [ ] **Step 2: Rung 6 — Strong duality & complementary slackness (proof).**
  - **Theorem (strong duality).** If the primal has a finite optimum, so does the dual, and their optimal values are **equal** (zero gap).
  - **Proof (∎).** Present as the LP specialization of Ch11's KKT/Slater strong-duality result: an LP is a convex program with affine constraints, every feasible LP trivially satisfies a constraint qualification (affine constraints need no Slater), so Ch11's strong-duality theorem applies and the duality gap is zero. (Alternatively/in addition, note the constructive route: at an optimal basis the simplex reduced-cost condition $\bar c_N\ge0$ yields a dual-feasible $y=(A_B^{-1})^\top c_B$ with $b^\top y=c^\top x^\star$, exhibiting equal values.) Cite `@bertsimas1997` for the self-contained Farkas-lemma proof.
  - **Theorem (complementary slackness).** Primal-feasible $x$ and dual-feasible $y$ are **both optimal** $\iff$ $y_i(a_i^\top x - b_i)=0$ for all $i$ **and** $x_j(c_j-(A^\top y)_j)=0$ for all $j$. **Proof (∎):** weak-duality chain is tight $\iff$ both inequalities are equalities $\iff$ each product term vanishes. This is the LP face of Ch11's $\lambda_i g_i=0$: a constraint with slack has dual price $0$; a positive dual price forces its constraint to bind.
- [ ] **Step 3: Lint + commit.** Three new `\blacksquare` (weak, strong, complementary slackness). `-m "feat(ch12): Theory rungs 5-6 — LP duality (weak/strong + complementary slackness)"`.

---

## Task 5: Theory rungs 7–9 — Shadow prices, the PWL budget LP, and the boundary

**Files:** Modify chapter (rungs 7–9, completing Theory).

- [ ] **Step 1: Rung 7 — Shadow prices (the economic reading).**
  - The optimal dual variable $y_i^\star$ is the **shadow price** of constraint $i$: within the range where the optimal basis is unchanged, $\partial p^\star/\partial b_i = y_i^\star$ — the marginal change in the optimal objective per unit relaxation of the right-hand side (state; short derivation from $p^\star=b^\top y^\star$ and basis stability). By complementary slackness a **slack** constraint has shadow price $0$, a **binding** one a positive price. For a budget constraint, the shadow price is the **marginal response per extra dollar of budget** — precisely Ch11's $\lambda$.

- [ ] **Step 2: Rung 8 — The MMM payoff: budget allocation as an LP (the PWL bridge).**
  - Approximate each concave response $r_i$ by a **piecewise-linear concave** function: split spend $b_i=\sum_k s_{ik}$, $s_{ik}\in[0,w_{ik}]$, with **decreasing** slopes $m_{i1}>m_{i2}>\cdots$ (concavity ⇒ decreasing marginal returns). Then
$$
\max \sum_{i,k} m_{ik}\,s_{ik}\quad \text{s.t.}\quad \sum_{i,k}s_{ik}=B,\ \ 0\le s_{ik}\le w_{ik},\ \ (\text{+ per-channel min/max caps})
$$
is a genuine **LP**.
  - **Why concavity makes it an LP that "just works":** because slopes decrease, the LP optimum **fills the steepest segments first** (a lower-slope segment is never filled while a steeper one has room — else swapping budget improves the objective). The **budget dual price equals the slope of the marginal (last-funded) segment** — the common marginal return $\lambda$, **Ch11's equal-marginal-returns capstone rediscovered by LP duality**. Worked-in numbers: WE3 gives $b^\star=(3,2)$, $\lambda=3$; refining a $\sqrt{}$-response PWL converges to Ch11's $(7.2,1.8)$, $\lambda\to0.3727$ (forward-reference Code/Worked).
  - Layered **caps/minimums** get their own shadow prices (marginal value of relaxing a cap) — multi-constraint LP duality in action.

- [ ] **Step 3: Rung 9 — The honest rung: where LP stops.**
  - PWL-as-LP is an **approximation**; refining segments converges to the smooth concave optimum (error shrinks with segment width — state).
  - But it rests entirely on **concavity**: "fill steepest first" is optimal *only* because slopes decrease. For an **S-shaped** response (Ch11's $n>1$ Hill, convex toe), a PWL fit has **increasing-then-decreasing** slopes; filling steepest-first is wrong, the LP relaxation no longer matches the true problem, and choosing whether to fund the low-slope toe at all is a **combinatorial (integer) / nonlinear** decision — not an LP. (One sentence on integer programming / branch-and-bound as the discrete tool.)
  - The boundary: **LP owns the linear and concave-separable world exactly; the smooth, possibly non-concave world is Ch13's SLSQP.** Bridge explicitly — Ch13 attacks the nonlinear KKT system iteratively (as Ch11 foreshadowed) and runs from multiple starts because the response can be non-concave. End the Theory section here.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch12): Theory rungs 7-9 — shadow prices, PWL budget LP, and the boundary"`.

---

## Task 6: Worked Examples

**Files:** Modify chapter (`## Worked Examples`). Reproduce verified anchors exactly.

- [ ] **Step 1: Example 1 — a 2-variable LP, solved twice.** $\max 3x_1+2x_2$ s.t. $x_1+x_2\le4,\ x_1\le2,\ x\ge0$. (a) **Graphically:** list the four vertices $(0,0),(2,0),(2,2),(0,4)$, evaluate the objective $0,6,10,8$, optimum at vertex $(2,2)$, value $10$ — the Fundamental Theorem concrete. (b) **By simplex:** standard form with slacks $s_1,s_2$; pivot 1 ($x_1$ enters, $s_2$ leaves, min-ratio $2$, value $6$); pivot 2 ($x_2$ enters, $s_1$ leaves, value $10$); all reduced costs $\ge0$ ⇒ optimal $(2,2)$.
- [ ] **Step 2: Example 2 — the dual and complementary slackness.** Dual $\min 4y_1+2y_2$ s.t. $y_1+y_2\ge3,\ y_1\ge2,\ y\ge0$. Solve: both primal vars positive ⇒ dual constraints bind ($y_1+y_2=3,\ y_1=2$) ⇒ $y^\star=(2,1)$, value $10$ = primal (**strong duality**). Complementary slackness: both primal constraints bind ⇒ both dual prices positive; interpret $y_1=2$ as "+1 budget → +2 response," $y_2=1$ as "+1 cap → +1 response."
- [ ] **Step 3: Example 3 — PWL-concave budget allocation.** Two channels, slopes ch1 $[5,3,1]$, ch2 $[4,2,1]$, widths $2$, budget $B=5$. Fill steepest first: ch1-seg1 ($5$, full $2$), ch2-seg1 ($4$, full $2$), ch1-seg2 ($3$, $1$ of $2$) — budget gone. $b^\star=(3,2)$, total $10+3+8=21$, **budget shadow price $\lambda=3$** (the marginal segment's slope). Check the water-filling condition: all funded segment slopes $\ge3$, all unfunded $\le3$. Note this is Ch11's equal-marginal rule in PWL clothing, and that refining the segmentation (Code Tie-in) converges to Ch11's smooth $(7.2,1.8)$, $\lambda\to0.3727$.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch12): Worked Examples — LP+simplex, dual, PWL allocation"`.

---

## Task 7: Code Tie-in

**Files:** Modify chapter (`## Code Tie-in`, single `{python}` cell, NumPy/Matplotlib/SciPy).

- [ ] **Step 1: Build the cell.** In order:
  1. Imports (`numpy`, `matplotlib.pyplot`, `from scipy.optimize import linprog`).
  2. **Solve WE1** with `linprog` (minimize $-(3x_1+2x_2)$, `A_ub=[[1,1],[1,0]]`, `b_ub=[4,2]`, `method="highs"`); print optimal $x=(2,2)$, value $10$; read **dual values** from `res.ineqlin.marginals` (negate sign convention) → shadow prices $(2,1)$; assert they match. Verify **strong duality** (primal value = $b^\top y$) and **complementary slackness** numerically.
  3. **PWL budget LP (WE3):** slopes `[5,3,1,4,2,1]`, widths `[2]*6`, `A_eq=[ones]`, `b_eq=[5]`, segment bounds; solve; recover $b=(3,2)$; read `res.eqlin.marginals` → $\lambda=3$; assert.
  4. **PWL → smooth convergence (the payoff):** response $r_i=a_i\sqrt b$, $a=(2,1)$, $B=9$; build PWL with $K$ secant segments per channel; solve for $K\in\{3,10,50,200\}$; print allocation and budget dual price each; show allocation → $(7.2,1.8)$ and $\lambda\to0.3727$ (Ch11's value); `assert` the $K{=}200$ allocation is within tol of $(7.2,1.8)$ and the response within tol of $6.7082$.
  5. **Figures (each `plt.show()`):** (a) the WE1 feasible polygon with vertices and the optimal corner marked; (b) the PWL approximations of $a_1\sqrt b$ for increasing $K$ overlaid on the smooth curve (showing convergence); (c) the budget dual price $\lambda$ vs $K$ approaching the dashed line $0.3727$.
- [ ] **Step 2: Verify headless.** Extract to `/tmp/ch12_cell.py`; `cd <WT> && MPLBACKEND=Agg python3 /tmp/ch12_cell.py`; iterate until asserts pass and prints match. **Pin printed numbers into surrounding prose.**
- [ ] **Step 3: Lint + commit.** Even `$`; one ```{python} cell. `-m "feat(ch12): Code Tie-in — linprog, duality, PWL convergence to Ch11's lambda"`.

---

## Task 8: Exercises (written directly by the controller)

**Files:** Modify chapter (`## Exercises`). Four tiers, self-contained. **No `psmallmatrix`.**

- [ ] **Step 1: Write C/B/P/A.**
  - **C (3):** why an LP optimum can always be taken at a vertex (and what that buys algorithmically); what a zero vs positive shadow price means (complementary slackness); why the PWL-as-LP trick needs concavity and fails for an S-curve.
  - **B (3):** solve a small 2-var LP graphically (fresh numbers, vertices + optimum); write and solve its dual, verify strong duality + read complementary slackness; a small PWL fill-steepest-first allocation with a budget shadow price. **NumPy/SciPy-verify each fresh number.**
  - **P (2–3):** prove weak duality; prove vertex ⇒ BFS (or BFS ⇒ vertex); prove that for a concave PWL the "fill steepest first" allocation is optimal (exchange argument) / equivalently that the budget dual price equals the marginal segment slope.
  - **A (2–3):** solve an LP with `linprog` and extract shadow prices; build a PWL approximation of a concave response and show the LP optimum converges to the smooth optimum as segments refine; (stretch) show a non-concave PWL where fill-steepest-first gives a suboptimal/incorrect allocation.
- [ ] **Step 2: Lint + commit.** `-m "feat(ch12): Exercises — C/B/P/A tiers"`.

---

## Task 9: Summary (written directly by the controller)

**Files:** Modify chapter (`## Summary`, auto-included; inline math only).

- [ ] **Step 1: Key concepts + Key identities (inline math only — no `$$`).**
  - **Concepts:** LP standard form; polyhedron & vertices = basic feasible solutions; Fundamental Theorem (optimum at a vertex); simplex (vertex-hopping, reduced costs, min-ratio, Bland's rule); weak & strong duality; complementary slackness; shadow prices; the PWL-concave budget LP and the boundary where LP stops.
  - **Identities (inline):** standard form $\min c^\top x$ s.t. $Ax=b,x\ge0$; primal–dual $\min c^\top x,\ Ax\ge b,x\ge0$ ↔ $\max b^\top y,\ A^\top y\le c,y\ge0$; weak duality $b^\top y\le c^\top x$; strong duality $b^\top y^\star=c^\top x^\star$; complementary slackness $y_i(a_i^\top x-b_i)=0,\ x_j(c_j-(A^\top y)_j)=0$; shadow price $\partial p^\star/\partial b_i=y_i^\star$; PWL budget dual price $=\lambda$ (Ch11's equal-marginal return).
- [ ] **Step 2: Lint + commit.** `-m "feat(ch12): Summary — key concepts and identities"`.

---

## Task 10: Appendix solutions block

**Files:** Modify `appendix/solutions.qmd` (insert before the final `:::`).

- [ ] **Step 1: Append `## Chapter 12 — Linear Programming`** with full worked solutions to every C/B/P/A exercise (NumPy/SciPy-verified numerics; code cells use `linprog`), inside the `show-solutions` block, before the file's last `:::`. **No `psmallmatrix`.** No inline links back.
- [ ] **Step 2: Lint + commit.** Even `$`; closes with single `:::`; Ch12 heading present (em-dash); appendix code cells run headless. `-m "feat(ch12): appendix solutions for Linear Programming"`.

---

## Task 11: Final review + PR + CI

- [ ] **Step 1: Chapter lint.** H1 `# Linear Programming`; six template H2 in order; no bare `\begin{align}`; even `$`; **no `psmallmatrix` anywhere**; `\blacksquare` count ≥5 (vertex⇔BFS, Fundamental Theorem, weak, strong, complementary slackness); citations resolve (`@bertsimas1997`/`@boyd2004`/`@nocedal2006`).
- [ ] **Step 2: Appendix lint.** Even `$`; closes with single `:::`; `## Chapter 12 — Linear Programming` present; no `psmallmatrix`.
- [ ] **Step 3: Re-run chapter + appendix code cells headless** — asserts pass, prints match pinned numbers.
- [ ] **Step 4: Through-line check.** Ch11 KKT→duality continuity; PWL budget LP rediscovers Ch11's λ; shadow-price reading; honest boundary bridges to Ch13.
- [ ] **Step 5: Push + PR.** `git -C <WT> push -u origin ch12-linear-programming`; `gh pr create` (title "Chapter 12 — Linear Programming").
- [ ] **Step 6: Watch CI `quarto render` (HTML **and PDF**) to green** — the `psmallmatrix`-class bug only surfaces in the PDF stage, so confirm BOTH. Report PR link + status. Then the user merges.

---

## Self-Review (controller, before dispatching Task 1)

- **Spec coverage:** standard form ✓(T2) · vertices=BFS ✓(T2) · Fundamental Theorem ✓(T3) · simplex ✓(T3) · weak duality ✓(T4) · strong duality + complementary slackness ✓(T4) · shadow prices ✓(T5) · PWL budget LP ✓(T5) · LP boundary→Ch13 ✓(T5) · 4–5 proofs ✓.
- **Placeholder scan:** none (PR body at T11).
- **Anchor consistency:** WE1 $(2,2)/10$, dual $(2,1)/10$, shadow prices $(2,1)$; WE3 $(3,2)/\lambda{=}3$; convergence $(7.2,1.8)/\lambda\to0.3727/6.7082$ — all SciPy-verified, used identically in Worked Examples (T6) and Code (T7); Exercises (T8) use fresh, separately-verified numbers.
- **PDF safety:** psmallmatrix banned in conventions and rechecked at T1/T8/T10/T11; CI watch covers HTML **and** PDF.
