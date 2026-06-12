# Chapter 11 — Convexity Theory Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 11 — Convexity Theory, the opening chapter of Part IV: convex sets and functions, the first/second-order conditions, Jensen, the local=global keystone, Lagrange/KKT, the equal-marginal-returns capstone, and an honest rung on S-curve non-convexity — anchored on the Hill saturation curve and the MMM budget-allocation problem, shipped as a PR whose CI `quarto render` passes.

**Architecture:** Single `{python}` cell (NumPy/Matplotlib/SciPy). Six full proofs (first-order, second-order, Jensen, local=global, KKT sufficiency, equal-marginal-returns). Fixed chapter template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary).

**Tech Stack:** Quarto `.qmd`, KaTeX math, Python 3 (NumPy, Matplotlib, SciPy `MPLBACKEND=Agg`), `references.bib`.

---

## Conventions (enforce in every task)

- **Worktree:** all work in `~/.config/superpowers/worktrees/mmm_foundational_text/ch11-convexity`. **Use explicit `git -C <worktree>` for every git command** (a bare `cd` may not apply and would send commits to main, where a hook blocks them). Identity already set (`jlh530i`).
- **KaTeX:** use `aligned` inside `$$…$$`; **never** bare `\begin{align}`. `$$` delimiters on their own lines. Even `$` count per file (check with `python3 -c "print(open(F).read().count(chr(36))%2)"` → 0). Inline math only in Summary.
- **Chapter file:** `parts/04-optimization/01-convexity.qmd`. Keep H1 `# Convexity Theory` and the canonical-anchors italic line. Replace the stub body.
- **Citations:** `@boyd2004` (anchor), `@nocedal2006`, `@gelman2013` (convex-posterior link). No new bib keys.
- **Verification:** the single `{python}` cell runs top-to-bottom under `MPLBACKEND=Agg python3`; figures end `plt.show()`. SciPy available. Real gate: CI `quarto render` (HTML+PDF) green on PR.
- **Voice/structure exemplars (read-only):** `parts/03-sampling/04-model-checking.qmd` (Ch10, immediate predecessor) for rung/proof style; `parts/02-regression-bayes/02-bayesian-inference.qmd` (Ch5) for the Hill curve notation $f(x;K,n)=x^n/(K^n+x^n)$.

## Verified numeric anchors (NumPy/SciPy/SymPy-checked)

- **Hill inflection** (Rung 4 / Worked c): for $f(x;K,n)=x^n/(K^n+x^n)$, $f''=0$ at $x^\star = K\left(\frac{n-1}{n+1}\right)^{1/n}$. For $n=2,K=3$: $x^\star = \sqrt{3} \approx 1.7321$. For $n=1$: $f''(x) = -2K/(K+x)^3 < 0$ on $x>0$ (concave).
- **Equal-marginal-returns allocation** (Capstone / Worked b): responses $r_i(b)=a_i\sqrt{b}$, so $r_i'(b)=a_i/(2\sqrt{b})$; equalizing gives $b_i^\star = B\,a_i^2/\sum_j a_j^2$. For $a=(2,1),\,B=9$: $b^\star=(7.2,\,1.8)$, common marginal $\lambda = 1/\sqrt{7.2} \cdot \tfrac12 \cdot 2 = 0.3727$, optimal total $2\sqrt{7.2}+\sqrt{1.8}=6.7082$ vs equal-split $3\sqrt{4.5}=6.3640$ (gain $0.3442$). NumPy optimizer confirms $b=(7.200,1.800)$.
- **Non-convex multi-start** (Rung 10 / Code): two channels, $\max\,[\,\text{hill}(b_1)+\text{hill}(b_2)\,]$ s.t. $b_1+b_2=6$, Hill $n=2,K=3$. Local optimizer lands on the **symmetric** interior optimum $b=(3,3)$ with total $1.000$ from central starts, but on a **corner** $b=(0,6)$ or $(6,0)$ with total $0.800$ from extreme starts — multiple local optima, the non-convexity signature. (hill(6)$=36/45=0.8$; hill(3)+hill(3)$=0.5+0.5=1.0$.)
- **Convexity certificate** (Worked a / Code): quadratic with $Q=\begin{psmallmatrix}2&0.5\\0.5&1\end{psmallmatrix}$ has eigenvalues $\{0.7929,2.2071\}>0$ ⇒ Hessian PD ⇒ strictly convex.

---

## Task 1: Front matter + Motivation

**Files:** Modify `parts/04-optimization/01-convexity.qmd` (keep H1 + anchors line; write `## Motivation`).

- [ ] **Step 1: Write Motivation (~400–550 words).**
  - Open Part IV: Parts I–III built and *checked* a posterior; now we *act* on it — allocate a fixed budget across channels to maximize response. That is an optimization problem, and whether it is solvable hinges on one property: **convexity**.
  - The driving question: *why is one optimization problem solvable and another a minefield?* A convex problem has no false summits — any local optimum is the global one, so a downhill algorithm cannot be fooled. A non-convex one is a landscape of traps.
  - Where it shows up in an MMM: the budget-allocation problem $\max \sum_i r_i(b_i)$ s.t. $\mathbf{1}^\top b = B$. If each response $r_i$ is **concave** (diminishing returns — the saturation curve of Ch5), the problem is convex and has a clean, unique answer with a beautiful structure (equal marginal returns). If the curves are **S-shaped** (a convex "toe" at low spend), the problem is non-convex and the optimizer can land on a bad allocation — which is exactly why Ch13 reaches for multi-start SLSQP.
  - Roadmap: convex sets → convex functions → the conditions that certify convexity → the keystone (local=global) → constrained optimality (Lagrange/KKT) → the MMM capstone (equal marginal returns) → the honest reckoning with non-convexity.
  - Position Ch11 as the theory hub: Ch12 (LP) and Ch13 (SLSQP) draw on the optimality conditions proved here.
- [ ] **Step 2: Lint + commit.** Even `$`, no bare align. `git -C <WT> add … && git -C <WT> commit -m "feat(ch11): Motivation — why convexity makes optimization solvable"`.

---

## Task 2: Theory rungs 1–2 — Convex sets and convex functions

**Files:** Modify chapter (begin `## Theory & Proofs`; rungs 1–2).

- [ ] **Step 1: Open `## Theory & Proofs` with a short ladder paragraph** (sets → functions → conditions → optimality, the path to the two keystones).

- [ ] **Step 2: Rung 1 — Convex sets.**
  - Definition: $C$ convex if $x,y\in C,\ \lambda\in[0,1] \Rightarrow \lambda x+(1-\lambda)y\in C$ (the segment stays inside).
  - Examples: hyperplane $\{x:a^\top x=b\}$, halfspace $\{x:a^\top x\le b\}$, norm ball, **polyhedron** $\{x:Ax\le b\}$, and the **budget feasible set** $\mathcal{B}=\{b\in\mathbb{R}^n: b\ge 0,\ \mathbf{1}^\top b\le B\}$ (a simplex/polyhedron).
  - **Proposition + proof:** the intersection of convex sets is convex (pick two points in the intersection, the segment lies in each set, hence in the intersection — `\blacksquare`). Corollary: $\mathcal{B}$ is convex (intersection of halfspaces).

- [ ] **Step 3: Rung 2 — Convex functions.**
  - Definition (requires convex domain): $f(\lambda x+(1-\lambda)y)\le \lambda f(x)+(1-\lambda)f(y)$ — the chord lies above the graph. Strict convexity with strict inequality.
  - **Epigraph** characterization: $f$ convex $\iff \text{epi}\,f=\{(x,t):f(x)\le t\}$ is a convex set (state; one-line proof or cite `@boyd2004`).
  - **Concave** $=-$convex; a diminishing-returns response is concave. Introduce the Hill curve $f(x;K,n)=x^n/(K^n+x^n)$ as the motivating concave response (for $n=1$) — maximizing a concave function is equivalent to minimizing its convex negative, so "concave maximization" and "convex minimization" are the same theory.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch11): Theory rungs 1-2 — convex sets and functions"`.

---

## Task 3: Theory rungs 3–4 — First- and second-order conditions (two proofs)

**Files:** Modify chapter (rungs 3–4).

- [ ] **Step 1: Rung 3 — First-order condition (full proof).**
  - **Theorem.** For differentiable $f$ on a convex domain, $f$ is convex $\iff f(y)\ge f(x)+\nabla f(x)^\top(y-x)$ for all $x,y$ — the graph lies above each tangent.
  - **Proof (∎):** ($\Rightarrow$) from convexity, $f(x+\lambda(y-x))\le f(x)+\lambda(f(y)-f(x))$; rearrange to $\frac{f(x+\lambda(y-x))-f(x)}{\lambda}\le f(y)-f(x)$, let $\lambda\to0^+$, the left side $\to \nabla f(x)^\top(y-x)$. ($\Leftarrow$) apply the tangent inequality at $z=\lambda x+(1-\lambda)y$ to both $x$ and $y$, take the $\lambda$-weighted combination, the gradient terms cancel, yielding the chord inequality. Use `aligned` inside `$$`.
  - Geometric reading: supporting hyperplane; at a stationary point this becomes $f(y)\ge f(x)$ globally (foreshadow keystone).

- [ ] **Step 2: Rung 4 — Second-order condition (full proof) + Hill application.**
  - **Theorem.** For twice-differentiable $f$ on an open convex domain, $f$ convex $\iff \nabla^2 f(x)\succeq 0$ (Hessian PSD) everywhere.
  - **Proof (∎):** restrict to the line $g(t)=f(x+tv)$; then $g''(t)=v^\top\nabla^2 f(x+tv)\,v$. $f$ convex $\iff$ every such 1-D restriction $g$ is convex $\iff g''\ge0 \iff v^\top\nabla^2 f\,v\ge0$ for all $v$ $\iff \nabla^2 f\succeq0$ (reuse Ch1 PSD characterization). Include the 1-D fact $g$ convex $\iff g''\ge0$.
  - **Hill application:** compute $f''(x;K,n)$; show for $n=1$, $f''(x)=-2K/(K+x)^3<0$ on $x>0$ (strictly concave); for $n>1$ the curve has an **inflection** at $x^\star=K\left(\frac{n-1}{n+1}\right)^{1/n}$ — convex below, concave above (sets up Rung 10). State the $n=2,K=3$ value $x^\star=\sqrt3\approx1.7321$.
- [ ] **Step 3: Lint + commit.** Two new `\blacksquare`. `-m "feat(ch11): Theory rungs 3-4 — first/second-order conditions (proofs)"`.

---

## Task 4: Theory rungs 5–6 — Operations preserving convexity; Jensen (proof)

**Files:** Modify chapter (rungs 5–6).

- [ ] **Step 1: Rung 5 — Operations that preserve convexity.**
  - State the toolkit: (i) nonnegative weighted sums $\sum w_i f_i$ ($w_i\ge0$); (ii) composition with an affine map $f(Ax+b)$; (iii) pointwise maximum/supremum $\max_i f_i$; (iv) scalar composition rules ($h\circ g$ convex when, e.g., $h$ convex nondecreasing and $g$ convex).
  - **Proofs** for (i) and (ii) (direct from the chord definition — `\blacksquare` on the combined proposition). State (iii)–(iv) with `@boyd2004`.
  - Payoff: certify built-up objectives without re-deriving — e.g. the **negative log-posterior of Ch5's ridge model**, $\tfrac{1}{2\sigma^2}\|y-X\beta\|^2 + \tfrac{1}{2\tau^2}\|\beta\|^2$, is convex in $\beta$ as a nonnegative sum of convex quadratics (so MAP estimation is a convex problem — connect back to Part II).

- [ ] **Step 2: Rung 6 — Jensen's inequality (full proof).**
  - **Theorem.** For convex $f$ and random $X$ (with finite mean), $f(\mathbb{E}[X])\le\mathbb{E}[f(X)]$ (reverse for concave).
  - **Proof (∎):** supporting-hyperplane at $\mu=\mathbb{E}[X]$: $f(X)\ge f(\mu)+\nabla f(\mu)^\top(X-\mu)$ (Rung 3); take expectations, the linear term vanishes since $\mathbb{E}[X-\mu]=0$. (Optionally the finite-sample induction.)
  - MMM reading: response is concave, so by Jensen $r(\mathbb{E}[\text{spend}])\ge\mathbb{E}[r(\text{spend})]$ — **steady spend beats erratic spend** of the same average; the gap is the curvature, a quantitative diminishing-returns statement.
- [ ] **Step 3: Lint + commit.** One new `\blacksquare` (operations) + one (Jensen). `-m "feat(ch11): Theory rungs 5-6 — convexity-preserving operations and Jensen (proof)"`.

---

## Task 5: Theory rung 7 — KEYSTONE: local = global (proof)

**Files:** Modify chapter (rung 7, keystone).

- [ ] **Step 1: Rung 7 — local minimum = global minimum.**
  - **Theorem.** Let $f$ be convex on a convex set $C$. Then every local minimizer of $f$ on $C$ is a global minimizer. If $f$ is differentiable, $\nabla f(x^\star)=0$ (more generally $\nabla f(x^\star)^\top(y-x^\star)\ge0$ for all $y\in C$) is **sufficient** for global optimality. If $f$ is strictly convex, the minimizer is unique.
  - **Proof (∎):** suppose $x^\star$ is a local but not global min, so $\exists\,z\in C$ with $f(z)<f(x^\star)$. For $\lambda\in(0,1]$, convexity gives $f(\lambda z+(1-\lambda)x^\star)\le\lambda f(z)+(1-\lambda)f(x^\star)<f(x^\star)$. Taking $\lambda$ small puts a point arbitrarily close to $x^\star$ with strictly smaller value — contradicting local minimality. Sufficiency of $\nabla f(x^\star)=0$ via the first-order condition (Rung 3): $f(y)\ge f(x^\star)+\nabla f(x^\star)^\top(y-x^\star)=f(x^\star)$. Strict convexity ⇒ uniqueness (two distinct minimizers would give a strictly lower midpoint).
  - Significance: this is why convex optimization *works* — and why detecting convexity (Rungs 3–4) is the first move on any optimization problem.
- [ ] **Step 2: Lint + commit.** New `\blacksquare`. `-m "feat(ch11): Theory rung 7 — local=global (keystone proof)"`.

---

## Task 6: Theory rung 8 — Lagrange & KKT (sufficiency proof)

**Files:** Modify chapter (rung 8).

- [ ] **Step 1: Rung 8 — Constrained optimality.**
  - The convex program: $\min f(x)$ s.t. $g_i(x)\le0\ (i=1..m)$, $h_j(x)=0\ (j=1..p)$, with $f,g_i$ convex and $h_j$ affine.
  - **Lagrangian** $L(x,\lambda,\nu)=f(x)+\sum_i\lambda_i g_i(x)+\sum_j\nu_j h_j(x)$, multipliers $\lambda\ge0$, $\nu$ free.
  - **KKT conditions** at $(x^\star,\lambda^\star,\nu^\star)$: stationarity $\nabla f(x^\star)+\sum_i\lambda_i^\star\nabla g_i(x^\star)+\sum_j\nu_j^\star\nabla h_j(x^\star)=0$; primal feasibility $g_i(x^\star)\le0,\,h_j(x^\star)=0$; dual feasibility $\lambda_i^\star\ge0$; **complementary slackness** $\lambda_i^\star g_i(x^\star)=0$.
  - **Theorem (sufficiency).** For a convex program, any point satisfying KKT is a **global** optimum.
  - **Proof (∎):** define $\tilde L(x)=f(x)+\sum_i\lambda_i^\star g_i(x)+\sum_j\nu_j^\star h_j(x)$. It is convex in $x$ (nonneg $\lambda^\star$ times convex $g_i$, plus affine), and KKT stationarity says $\nabla\tilde L(x^\star)=0$, so by the keystone (Rung 7) $x^\star$ minimizes $\tilde L$. Then for any feasible $y$: $f(y)\ge f(y)+\sum_i\lambda_i^\star g_i(y)+\sum_j\nu_j^\star h_j(y)=\tilde L(y)\ge\tilde L(x^\star)=f(x^\star)$, where the first $\ge$ uses $\lambda_i^\star\ge0,\,g_i(y)\le0,\,h_j(y)=0$, and the last equality uses complementary slackness. Hence $x^\star$ is globally optimal.
  - **Necessity & Slater (stated):** if a strictly feasible point exists (Slater's condition), KKT is also *necessary*; strong duality holds. Cite `@boyd2004`, `@nocedal2006`. Note this is the hand-off: Ch12 uses KKT/duality on polyhedra, Ch13's SLSQP solves these conditions iteratively.
- [ ] **Step 2: Lint + commit.** New `\blacksquare`. `-m "feat(ch11): Theory rung 8 — Lagrange & KKT (sufficiency proof)"`.

---

## Task 7: Theory rungs 9–10 — Capstone (equal marginal returns) + non-convexity

**Files:** Modify chapter (rungs 9–10).

- [ ] **Step 1: Rung 9 — CAPSTONE: equal marginal returns (full proof via KKT).**
  - The budget problem: $\max_{b\ge0}\sum_i r_i(b_i)$ s.t. $\mathbf{1}^\top b=B$, each $r_i$ concave, increasing, differentiable. Equivalently $\min -\sum_i r_i(b_i)$ (convex).
  - **Theorem.** At the optimum there is a common $\lambda$ (the shadow price of budget) such that every funded channel ($b_i^\star>0$) satisfies $r_i'(b_i^\star)=\lambda$, and every unfunded channel ($b_i^\star=0$) satisfies $r_i'(0)\le\lambda$. (Equal marginal returns / water-filling.)
  - **Proof (∎):** form the Lagrangian with equality multiplier $\nu$ for the budget and $\mu_i\ge0$ for $b_i\ge0$; KKT stationarity gives $-r_i'(b_i^\star)+\nu-\mu_i=0$; complementary slackness $\mu_i b_i^\star=0$. For $b_i^\star>0$: $\mu_i=0\Rightarrow r_i'(b_i^\star)=\nu\equiv\lambda$. For $b_i^\star=0$: $\mu_i\ge0\Rightarrow r_i'(0)=\nu-\mu_i\le\lambda$. Concavity makes KKT sufficient (Rung 8) so this is the global optimum.
  - Worked-in number: $a=(2,1),B=9$, $r_i=a_i\sqrt b$ ⇒ $b^\star=(7.2,1.8)$, $\lambda=0.3727$, beating the equal split (forward-reference Worked Examples).
  - Interpretation: budget flows until marginal returns equalize; a channel with steeper early returns gets more — *this is the principle Ch13's SLSQP converges to numerically.*

- [ ] **Step 2: Rung 10 — The honest rung: when the MMM problem is not convex.**
  - The Hill curve with $n>1$ is **S-shaped**: convex on the toe $x<x^\star$, concave above (Rung 4's inflection). So a real saturation response is **not globally concave**.
  - Consequence: $\sum_i r_i(b_i)$ is **not concave**, the budget problem is **non-convex**, the keystone's guarantee is void — **multiple local optima** exist, and complementary slackness/KKT identifies *stationary points* that may be local, not global.
  - Concrete: two channels, $n=2,K=3$, $B=6$ — a local optimizer reaches the symmetric optimum $(3,3)$ (total $1.0$) from central starts but a **corner** $(0,6)$ (total $0.8$) from extreme starts (forward-reference Code Tie-in).
  - What practitioners do: **multi-start** local optimization, careful initialization, sometimes convex relaxations or fitting concave surrogates; accept a *certified-local* rather than *certified-global* solution.
  - Bridge: Ch12 (LP) is the fully-convex special case where geometry is completely understood; Ch13 (SLSQP) is the constrained local solver run from multiple starts precisely because real response curves break convexity. End the Theory section here.
- [ ] **Step 3: Lint + commit.** New `\blacksquare` (capstone). `-m "feat(ch11): Theory rungs 9-10 — equal marginal returns (capstone) and non-convexity"`.

---

## Task 8: Worked Examples

**Files:** Modify chapter (`## Worked Examples`).

- [ ] **Step 1: Example 1 — convexity certificate by Hessian.** Take the quadratic spend-cost $c(b)=b^\top Q b$ with $Q=\begin{psmallmatrix}2&0.5\\0.5&1\end{psmallmatrix}$ (symmetric). Compute eigenvalues $\{0.7929,2.2071\}$ (both $>0$) ⇒ $Q\succ0$ ⇒ strictly convex by Rung 4. Contrast: locate the Hill inflection $x^\star=K(\tfrac{n-1}{n+1})^{1/n}$; for $n=2,K=3$, $x^\star=\sqrt3\approx1.7321$, so $f$ is convex on $[0,1.7321)$ and concave above — *not* globally convex.
- [ ] **Step 2: Example 2 — equal-marginal-returns allocation by hand.** Two channels, $r_i(b)=a_i\sqrt b$, $a=(2,1)$, budget $B=9$. Set $r_1'(b_1)=r_2'(b_2)$: $\frac{a_1}{2\sqrt{b_1}}=\frac{a_2}{2\sqrt{b_2}}\Rightarrow \frac{b_1}{b_2}=\frac{a_1^2}{a_2^2}=4$, with $b_1+b_2=9$ ⇒ $b^\star=(7.2,1.8)$. Common marginal $\lambda=\frac{a_1}{2\sqrt{7.2}}=0.3727$. Total response $2\sqrt{7.2}+\sqrt{1.8}=6.7082$, vs equal split $3\sqrt{4.5}=6.3640$ — the optimal allocation gains $0.3442$ by sending more to the steeper channel.
- [ ] **Step 3: Example 3 — an S-curve breaks convexity.** Hill $n=2,K=3$: show $f''$ changes sign at $x^\star=\sqrt3$ (convex toe, concave top). Argue the two-channel $B=6$ problem then has the symmetric stationary point $(3,3)$ AND the corner solutions $(0,6)/(6,0)$ as KKT points with different objective values ($1.0$ vs $0.8$), so KKT no longer certifies global optimality.
- [ ] **Step 4: Lint + commit.** `-m "feat(ch11): Worked Examples — convexity certificate, allocation, S-curve"`.

---

## Task 9: Code Tie-in

**Files:** Modify chapter (`## Code Tie-in`, single `{python}` cell).

- [ ] **Step 1: Build the cell (NumPy/Matplotlib/SciPy).** In order:
  1. Imports; `np.random` seed if needed.
  2. **Hill curves:** define `hill(x,K,n)`; plot concave ($n=1$) vs S-shaped ($n=2$) for $K=3$; mark the inflection $x^\star=\sqrt3$ for $n=2$.
  3. **Convexity certificate:** $Q=[[2,0.5],[0.5,1]]$; print `np.linalg.eigvalsh(Q)` ($\approx[0.793,2.207]$); assert all $>0$. Also a chord-test helper confirming $-$Hill (n=1) is convex on a grid.
  4. **Concave budget allocation (capstone, numeric):** responses $r_i(b)=a_i\sqrt b$, $a=(2,1)$, $B=9$; solve with `scipy.optimize.minimize` (SLSQP, equality constraint $\sum b=B$, bounds $b\ge0$); print $b^\star$ (expect $\approx(7.2,1.8)$), the per-channel marginals $a_i/(2\sqrt{b_i^\star})$ (expect equal $\approx0.373$), and confirm the response $6.7082>6.3640$ equal split. `assert np.allclose(marginals[0],marginals[1])` and `assert b_star ≈ (7.2,1.8)`.
  5. **Non-convex multi-start (Rung 10):** two channels, Hill $n=2,K=3$, $B=6$; minimize $-(\text{hill}(b_1)+\text{hill}(B-b_1))$ from several starts $\{0.1,1,3,5,5.9\}$; print the reached $b_1$ and total for each; show central starts → $(3,3)$ total $1.0$, extreme starts → corner total $0.8$. `assert len(set(round(total,3))) > 1` (multiple optima found).
  6. **Figures (each `plt.show()`):** Hill concave-vs-S plot with inflection; a bar/line of the allocation marginals equalizing; and (optional) the non-convex objective $b_1\mapsto$ total over $[0,6]$ showing two humps + corners.
- [ ] **Step 2: Verify headless.** Extract to `/tmp/ch11_cell.py`; `cd <WT> && MPLBACKEND=Agg python3 /tmp/ch11_cell.py`. Iterate until asserts pass and prints match. **Pin the printed numbers into surrounding prose.**
- [ ] **Step 3: Lint + commit.** Even `$`; one ```{python} cell. `-m "feat(ch11): Code Tie-in — convexity, equal-marginal allocation, non-convex multi-start"`.

---

## Task 10: Exercises (written directly by the controller)

**Files:** Modify chapter (`## Exercises`). Four tiers, self-contained (no inline solution links).

- [ ] **Step 1: Write C/B/P/A.**
  - **C (3):** why local=global makes a downhill method trustworthy; why a concave response ⇒ convex budget problem; what is lost when curves are S-shaped (no global certificate).
  - **B (3):** verify a given Hessian is PSD (convexity); solve a 2-channel equal-marginal allocation with a fresh concave form (e.g. $r_i(b)=a_i\log(1+b)$ or fresh $a$, $B$); locate a Hill inflection for given $n,K$. **NumPy-verify each fresh number before writing.**
  - **P (2–3):** prove the first-order condition (one direction) OR the keystone local=global; prove Jensen from the supporting hyperplane; prove that a nonnegative weighted sum of convex functions is convex.
  - **A (2–3):** implement a chord/Hessian convexity check; solve an $n$-channel concave allocation and confirm equal marginals; run a multi-start on an S-curve problem and report the spread of optima.
- [ ] **Step 2: Lint + commit.** `-m "feat(ch11): Exercises — C/B/P/A tiers"`.

---

## Task 11: Summary (written directly by the controller)

**Files:** Modify chapter (`## Summary`, auto-included; inline math only).

- [ ] **Step 1: Key concepts + Key identities.**
  - **Key concepts:** convex set; convex/concave function; first-order (tangent below) and second-order (Hessian PSD) conditions; operations preserving convexity; Jensen; local=global keystone; Lagrangian/KKT; equal-marginal-returns capstone; S-curve non-convexity and multi-start.
  - **Key identities (inline only):** chord inequality $f(\lambda x+(1-\lambda)y)\le\lambda f(x)+(1-\lambda)f(y)$; first-order $f(y)\ge f(x)+\nabla f(x)^\top(y-x)$; second-order $\nabla^2 f\succeq0$; Jensen $f(\mathbb{E}X)\le\mathbb{E}f(X)$; KKT stationarity + complementary slackness $\lambda_i g_i=0$; equal marginal returns $r_i'(b_i^\star)=\lambda$; Hill inflection $x^\star=K(\tfrac{n-1}{n+1})^{1/n}$.
- [ ] **Step 2: Lint + commit.** Inline math only (no `$$`). `-m "feat(ch11): Summary — key concepts and identities"`.

---

## Task 12: Appendix solutions block

**Files:** Modify `appendix/solutions.qmd` (insert before the final `:::`).

- [ ] **Step 1: Append `## Chapter 11 — Convexity Theory`** with full worked solutions to every C/B/P/A exercise (NumPy-verified numerics), inside the `show-solutions` block, before the file's last `:::`. No inline links back.
- [ ] **Step 2: Lint + commit.** Even `$`; file still closes with single `:::`; Ch11 heading present (em-dash). `-m "feat(ch11): appendix solutions for Convexity Theory"`.

---

## Task 13: Final review + PR + CI

- [ ] **Step 1: Chapter lint.** H1 `# Convexity Theory`; six template H2 in order; no bare `\begin{align}`; even `$`; `\blacksquare` count ≥6 (first-order, second-order, operations, Jensen, keystone, KKT, capstone — i.e. 6–7); citations resolve (`@boyd2004`/`@nocedal2006`/`@gelman2013`).
- [ ] **Step 2: Appendix lint.** Even `$`; closes with single `:::`; `## Chapter 11 — Convexity Theory` present.
- [ ] **Step 3: Re-run the code cell headless** — asserts pass, prints match pinned numbers.
- [ ] **Step 4: Through-line check.** Part IV opener framing; Hill/Ch5 continuity; Lagrange/KKT positioned as the hub for Ch12–13; non-convexity bridge to Ch13.
- [ ] **Step 5: Push + PR.** `git -C <WT> push -u origin ch11-convexity`; `gh pr create` (title "Chapter 11 — Convexity Theory").
- [ ] **Step 6: Watch CI `quarto render` to green** (HTML+PDF; deploy skipped on PR). Report PR link + status. Then the user merges.

---

## Self-Review (controller, before dispatching Task 1)

- **Spec coverage:** convex sets ✓(T2) · convex functions/epigraph ✓(T2) · first-order ✓(T3) · second-order + Hill ✓(T3) · operations ✓(T4) · Jensen ✓(T4) · keystone local=global ✓(T5) · Lagrange/KKT ✓(T6) · capstone equal-marginal ✓(T7) · non-convexity S-curve ✓(T7) · six proofs ✓.
- **Placeholder scan:** none (PR body filled at T13).
- **Anchor consistency:** Hill inflection $\sqrt3$, allocation $b^\star=(7.2,1.8)/\lambda=0.3727$/response $6.7082$ vs $6.3640$, multi-start $1.0$ vs $0.8$, $Q$ eigenvalues $\{0.793,2.207\}$ — all NumPy/SymPy-verified, used identically in Worked Examples (T8) and Code (T9); Exercises (T10) use fresh, separately-verified numbers.
