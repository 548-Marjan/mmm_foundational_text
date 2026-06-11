# Chapter 7 — Markov Chains Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/03-sampling/01-markov-chains.qmd` with the opening chapter of Part III — the theory of Markov chains that makes MCMC work, anchored on customer-journey attribution — and add matching appendix solutions.

**Architecture:** One Quarto `.qmd` chapter authored section-by-section against the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus a `## Chapter 7 — Markov Chains` block appended to the shared `appendix/solutions.qmd` (gated by `show-solutions`). Finite state spaces throughout, with general-state analogues stated and cited at the end. The chapter reuses Chapter 1's eigenvalue/spectral decomposition directly and closes the loop on Chapter 6's deferred promise (the Gibbs cycle is a Markov chain). The single runnable Python cell is verified headless before commit; the real render gate is CI (`quarto render`), Quarto not being installed locally.

**Tech Stack:** Quarto book, KaTeX (HTML) / LaTeX (PDF) math, Python code cell (numpy, scipy, matplotlib pinned in `requirements.txt`). Canonical references via `references.bib` keys: `@robert2004`, `@norris1997`, `@gelman2013`, `@strang2016`, `@axler2015`.

---

## Conventions (apply to every task)

- **Math/KaTeX:** Use `aligned` inside `$$ ... $$`, never bare `\begin{align}`. Keep `$$` delimiters on their own lines so display-math counts stay even. Inline math with single `$`. Use the **row-vector convention** $\pi_{t+1}=\pi_t P$ consistently in prose and code; $P$ is **row-stochastic**; $\pi$ is a **left** eigenvector of $P$ for eigenvalue $1$.
- **Voice:** Match Chapters 1–6 (read `parts/02-regression-bayes/03-hierarchical-regression.qmd` for the immediate predecessor). Rigorous, prose-led, each theory subsection ("rung") ending by tying back to the through-line / an MMM payoff. End every full proof with `$\blacksquare$`; use `> **Theorem (name).** ...` blockquotes for named results.
- **Reuse Chapter 1 explicitly:** eigenvalues, eigenvectors, spectral/eigen-decomposition, the spectral theorem. Cite/refer rather than re-deriving.
- **Citations:** `[@robert2004]`, `[@norris1997]`, `[@gelman2013]`, `[@strang2016]`, `[@axler2015]` only (these keys exist in `references.bib`). Do not invent keys.
- **Verification:** Verify every numeric claim against NumPy before committing. The code cell must run top-to-bottom under `MPLBACKEND=Agg python3` and end figures with `plt.show()`. No pytest exists; the build is the test.
- **Commits:** One commit per task, message prefix `feat(ch07): ...`. Run `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"` once if commits fail with "Author identity unknown".
- **No inline solution links** in exercises; solutions live only in the appendix, gated by `::: {.content-visible when-meta="show-solutions"}`.

## Verified anchor numbers (use these exact values)

- **Worked (a) customer-journey attribution.** Transient states Start, Display, Search; absorbing Convert, Null. Transitions: Start→Display $0.5$, Start→Search $0.5$; Display→Search $0.3$, Display→Convert $0.2$, Display→Null $0.5$; Search→Convert $0.4$, Search→Null $0.6$. Absorbing block form $P=\big[\begin{smallmatrix}Q&R\\0&I\end{smallmatrix}\big]$ with transient order (Start, Display, Search):
  $$Q=\begin{bmatrix}0&0.5&0.5\\0&0&0.3\\0&0&0\end{bmatrix},\quad R=\begin{bmatrix}0&0\\0.2&0.5\\0.4&0.6\end{bmatrix}\ (\text{columns Convert, Null}).$$
  Fundamental matrix $N=(I-Q)^{-1}$; absorption probabilities $B=NR$. **Conversion probability from Start $= 0.36$** (Search alone $0.40$, Display alone $0.32$). **Removal effect** (delete a channel; journeys into it are lost to Null): removing **Search** drops conversion to $0.10$ → effect $(0.36-0.10)/0.36\approx 0.722$; removing **Display** drops it to $0.20$ → effect $(0.36-0.20)/0.36\approx 0.444$. Search carries the larger attribution credit.
- **Worked (b) 2-state convergence.** $P=\begin{bmatrix}0.7&0.3\\0.4&0.6\end{bmatrix}$. Eigenvalues $\lambda_1=1$, $\lambda_2=\operatorname{tr}P-1=0.3$. Stationary $\pi=(4/7,3/7)\approx(0.5714,0.4286)$ (solve $\pi P=\pi$, $0.3\pi_1=0.4\pi_2$). **Detailed balance holds:** $\pi_1 P_{12}=\tfrac47\cdot0.3=\tfrac{12}{70}$ and $\pi_2 P_{21}=\tfrac37\cdot0.4=\tfrac{12}{70}$ — reversible. Convergence $\lVert\pi_t-\pi\rVert\sim 0.3^t$; spectral gap $1-|\lambda_2|=0.7$ (fast mixing).

---

## File Structure

- **Modify (replace body):** `parts/03-sampling/01-markov-chains.qmd` — the chapter. Keep the H1 `# Markov Chains` and the `*Canonical anchors: Robert & Casella; Norris.*` line; replace the stub callout and all empty section bodies with full content.
- **Modify (append):** `appendix/solutions.qmd` — insert a `## Chapter 7 — Markov Chains` block immediately before the closing `:::` of the `show-solutions` visible div (currently the last line, line 1008). Match the format of the existing `## Chapter 6 — Hierarchical Models` block.

---

### Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd:1-11`

- [ ] **Step 1: Clean the stub and write the Motivation**

Keep line 1 `# Markov Chains` and line 3 `*Canonical anchors: Robert & Casella; Norris.*`. Delete the `::: {.callout-note ...}` stub block. Replace the `## Motivation` placeholder body with a full section; leave the other section headings/placeholders for later tasks.

The Motivation (~3 paragraphs, Chapter-6 density) must:
- Open on the **two faces of Markov chains in MMM**: (i) marketers already model the **customer journey** across channel touchpoints as a Markov chain and attribute conversions via the **removal effect** (drop a channel, recompute the conversion probability, the decrease is that channel's credit); (ii) Chapter 6 left an unpaid promise — cycling the Gibbs full conditionals was *claimed* to converge to the posterior, with the justifying theorem deferred to here.
- State that both rest on the same questions: when does a chain settle into a unique long-run distribution, how fast, and when do **time-averages along one trajectory** equal **expectations** under that distribution?
- Pose the driving question verbatim (emphasized): *"Run the Gibbs sampler for a thousand steps and average — why should that average be the posterior mean, and how long is long enough?"*
- Preview the **flip**: Part II asked "given a chain, find its stationary distribution"; MCMC inverts it to "given a target distribution (the posterior), engineer a chain that converges to it." Name the three theorems the chapter will prove (detailed balance ⟹ stationarity; convergence via the spectral gap; the ergodic theorem) and that they respectively become MCMC's construction tool, convergence guarantee, and estimation license. Anchor to the Part II posterior $p(\beta\mid y)$ as the eventual target $\pi$.

- [ ] **Step 2: Commit**

```bash
cd /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch7-markov-chains
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): front matter and Motivation"
```

---

### Task 2: Theory & Proofs — rungs 1–3 (chains, classification, detailed balance)

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (the `## Theory & Proofs` section, first three rungs)

Write `## Theory & Proofs` with a 1–2 sentence ladder lead, then three `###` subsections. **Stop before rung 4** (next task). Prose-led, each rung closes on the through-line.

- [ ] **Step 1: Rung 1 — Markov chains and transition matrices**

Finite state space $S=\{1,\dots,m\}$. The **Markov property**:
$$
P(X_{t+1}=j\mid X_t=i, X_{t-1},\dots)=P(X_{t+1}=j\mid X_t=i)=P_{ij}.
$$
Define the **stochastic (row-stochastic) matrix** $P$ ($P_{ij}\ge0$, $\sum_j P_{ij}=1$). State **Chapman–Kolmogorov**: the $n$-step transition matrix is $P^n$, i.e. $P(X_{t+n}=j\mid X_t=i)=(P^n)_{ij}$ (prove in one line by induction/matrix multiplication). Distribution evolution as left-multiplication: if $\pi_t$ is the row vector of state probabilities then $\pi_{t+1}=\pi_t P$ and $\pi_t=\pi_0 P^t$. Introduce the **customer-journey chain**: states Start, Display, Search, and absorbing Convert, Null, with the concrete transitions from the anchor numbers; write its $P$.

- [ ] **Step 2: Rung 2 — Classification**

Define: state $j$ **accessible** from $i$ ($i\to j$) if $(P^n)_{ij}>0$ some $n$; **communication** ($i\leftrightarrow j$); **communicating classes**; **irreducible** = single class. **Period** $d(i)=\gcd\{n:(P^n)_{ii}>0\}$; **aperiodic** = all periods 1. **Recurrent vs transient**; **absorbing state** ($P_{ii}=1$). Put the customer-journey chain in **absorbing block form**
$$
P=\begin{bmatrix}Q&R\\0&I\end{bmatrix},
$$
$Q$ transient→transient, $R$ transient→absorbing. State the **fundamental matrix** $N=(I-Q)^{-1}=\sum_{k\ge0}Q^k$ (expected visit counts) and that **absorption probabilities** are $B=NR$ — the engine of the attribution example. Note these classification conditions (irreducible + aperiodic) are exactly the hypotheses of the convergence theorem in rung 4.

- [ ] **Step 3: Rung 3 — Stationary distributions & detailed balance (PROOF)**

Define the **stationary distribution**: a probability row vector $\pi$ with $\pi P=\pi$ — a **left eigenvector** of $P$ for eigenvalue $1$. State (citing **Perron–Frobenius**, `[@strang2016]`) that an irreducible finite chain has a **unique** strictly positive stationary $\pi$. Define **reversibility / detailed balance**: $\pi_i P_{ij}=\pi_j P_{ji}$ for all $i,j$.
> **Theorem (detailed balance ⟹ stationarity).** If $\pi$ satisfies $\pi_i P_{ij}=\pi_j P_{ji}$ for all $i,j$, then $\pi P=\pi$.

**Proof.** Sum the balance equation over $i$:
$$
(\pi P)_j=\sum_i \pi_i P_{ij}=\sum_i \pi_j P_{ji}=\pi_j\sum_i P_{ji}=\pi_j,
$$
using row-stochasticity $\sum_i P_{ji}=1$. Hence $\pi P=\pi$. End `$\blacksquare$`. Emphasize the MCMC payoff: detailed balance is a **local, checkable** condition (one equation per pair of states), so engineering a chain reversible w.r.t. a target $\pi$ is far easier than solving the global system $\pi P=\pi$ — this is precisely how Metropolis–Hastings (Chapter 8) will be built.

- [ ] **Step 4: Commit**

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Theory rungs 1-3 — chains, classification, detailed balance"
```

---

### Task 3: Theory & Proofs — rungs 4–6 (convergence keystone, ergodic theorem, the MCMC flip)

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (append three subsections to `## Theory & Proofs`, before `## Worked Examples`)

- [ ] **Step 1: Rung 4 — Convergence to stationarity (KEYSTONE proof, spectral gap)**

> **Theorem (convergence).** For an irreducible, aperiodic finite chain, $P^n\to\mathbf 1\pi$ entrywise (every row tends to $\pi$), so $\pi_t=\pi_0 P^t\to\pi$ from any initial $\pi_0$.

**Proof (diagonalizable case, via the Chapter 1 spectral decomposition).** By Perron–Frobenius, irreducibility + aperiodicity make $\lambda_1=1$ a **simple** eigenvalue and the **strictly dominant** one: all other eigenvalues satisfy $|\lambda_k|<1$. The right eigenvector for $\lambda_1=1$ is $\mathbf 1$ (since $P\mathbf 1=\mathbf 1$) and the left eigenvector is the stationary $\pi$ (since $\pi P=\pi$), normalized so $\pi\mathbf 1=1$. Writing the spectral decomposition (Chapter 1) $P=\sum_{k=1}^m \lambda_k\, r_k\ell_k^\top$ with $r_1=\mathbf 1$, $\ell_1=\pi^\top$,
$$
P^n=\sum_{k=1}^m \lambda_k^n\, r_k\ell_k^\top=\mathbf 1\pi+\sum_{k=2}^m \lambda_k^n\, r_k\ell_k^\top.
$$
Because $|\lambda_k|<1$ for $k\ge2$, the sum vanishes as $n\to\infty$, leaving $P^n\to\mathbf 1\pi$; and $\pi_0 P^n\to\pi_0\mathbf 1\pi=\pi$ since $\pi_0\mathbf 1=1$. The **rate** is set by the **subdominant eigenvalue** $\lambda_2$ (largest $|\lambda_k|$, $k\ge2$):
$$
\lVert\pi_t-\pi\rVert\le C\,|\lambda_2|^{\,t}.
$$
End `$\blacksquare$`. Define the **spectral gap** $\gamma=1-|\lambda_2|$ (larger gap ⇒ faster mixing) and the **mixing time** as the number of steps to bring $\lVert\pi_t-\pi\rVert$ below a tolerance. Note the general (non-diagonalizable / Jordan-block, or coupling) proof gives the same geometric rate; cite `[@norris1997]`. Reuses Chapter 1 eigen-decomposition directly. MMM payoff: the spectral gap is exactly what determines "how long is long enough" for an MCMC run.

- [ ] **Step 2: Rung 5 — The ergodic theorem (PROOF, finite case)**

> **Theorem (ergodic theorem for finite chains).** For an irreducible finite chain with stationary distribution $\pi$ and any bounded $f:S\to\mathbb R$, from any start,
> $$\frac1n\sum_{t=1}^n f(X_t)\xrightarrow{n\to\infty}\mathbb E_\pi[f]=\sum_{x}f(x)\pi_x\quad\text{almost surely.}$$

**Proof (renewal / return-time argument).** Fix a state $s$ (recurrent, since the chain is finite and irreducible). Let $0\le T_0<T_1<\dots$ be successive visit times to $s$. By the **strong Markov property**, the trajectory segments ("cycles") between consecutive visits are i.i.d. Let $\tau_k=T_k-T_{k-1}$ be the cycle lengths and $S_k=\sum_{t=T_{k-1}}^{T_k-1}f(X_t)$ the cycle sums; the pairs $(\tau_k,S_k)$ are i.i.d. with finite means (positive recurrence gives $\mathbb E[\tau]=1/\pi_s<\infty$). After $N$ completed cycles the time average is the ratio
$$
\frac{\sum_{t<T_N}f(X_t)}{T_N}=\frac{\frac1N\sum_{k=1}^N S_k}{\frac1N\sum_{k=1}^N \tau_k}\xrightarrow{\text{LLN}}\frac{\mathbb E[S]}{\mathbb E[\tau]}.
$$
The **renewal-reward identity** evaluates the limit: the expected number of visits to state $x$ within one $s$-cycle is $\pi_x/\pi_s=\pi_x\,\mathbb E[\tau]$, so $\mathbb E[S]=\sum_x f(x)\pi_x\,\mathbb E[\tau]$ and the ratio is $\sum_x f(x)\pi_x=\mathbb E_\pi[f]$. (Number of steps not aligned to a cycle boundary is negligible as $n\to\infty$.) End `$\blacksquare$`. This is the **law of large numbers for Markov chains** — the result that licenses estimating a posterior expectation $\mathbb E_{p(\beta\mid y)}[f]$ from a **single long MCMC run**, with no independence between draws required.

- [ ] **Step 3: Rung 6 — The flip: from analysis to MCMC**

Reverse the program of rungs 1–5. So far: given $P$, deduce $\pi$ and long-run behavior. **MCMC inverts it**: given a **target** $\pi$ — the Chapter 5 / Chapter 6 posterior $p(\beta\mid y)$, known only up to its normalizing constant — **construct** a transition rule $P$ such that (i) $\pi$ is stationary, achieved by enforcing **detailed balance** against $\pi$ (rung 3, a local condition needing only ratios $\pi_j/\pi_i$, so the unknown normalizer cancels); (ii) $P$ is irreducible and aperiodic, so **convergence** (rung 4) drives $\pi_t\to\pi$; then (iii) the **ergodic theorem** (rung 5) makes the trajectory average $\frac1n\sum f(X_t)$ an estimator of $\mathbb E_\pi[f]$. Show that Chapter 6's **Gibbs cycle** — resampling each block from its full conditional — leaves the posterior invariant, so it is exactly such a chain (the deferred Chapter 6 promise, now paid). Name **Metropolis–Hastings** as Chapter 8's general construction. **State** (citing `[@robert2004]`, `[@norris1997]`) the general-state-space analogues that continuous posteriors require: transition **kernels** $P(x,\mathrm dy)$, detailed balance against a density, **Harris recurrence**, and the general ergodic theorem — same logic, measure-theoretic statement. End the Theory section on this baton-pass to Chapter 8.

- [ ] **Step 4: Commit**

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Theory rungs 4-6 — convergence, ergodic theorem, the MCMC flip"
```

---

### Task 4: Worked Examples

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (the `## Worked Examples` section)

Match Chapter 6's worked-example style: bold step mini-headings, every number shown, closing interpretation.

- [ ] **Step 1: Example (a) — Customer-journey attribution by hand**

Use the anchor chain (transient Start, Display, Search; absorbing Convert, Null) with $Q,R$ as in the anchor numbers. **Conversion probability from Start:** solve by backward substitution — from Search, $P(\text{Convert})=0.4$; from Display, $0.2+0.3\cdot0.4=0.32$; from Start, $0.5\cdot0.32+0.5\cdot0.4=0.36$ (equivalently $B=NR$, the $(I-Q)^{-1}R$ Convert column). **Removal effect:** removing **Search** (its inbound traffic is lost to Null) gives Display→Convert $0.2$, Start $=0.5\cdot0.2+0.5\cdot0=0.10$; effect $(0.36-0.10)/0.36\approx0.722$. Removing **Display**: Start $=0.5\cdot0+0.5\cdot0.4=0.20$; effect $(0.36-0.20)/0.36\approx0.444$. Interpret: Search receives the larger conversion credit; this is a genuine Markov-chain MMM attribution, computed with rung 2's fundamental matrix.

- [ ] **Step 2: Example (b) — Convergence by hand (2-state chain)**

$P=\begin{bmatrix}0.7&0.3\\0.4&0.6\end{bmatrix}$. **Stationary distribution:** solve $\pi P=\pi$ → $0.3\pi_1=0.4\pi_2$ with $\pi_1+\pi_2=1$ → $\pi=(4/7,3/7)\approx(0.571,0.429)$. **Eigenvalues:** $\lambda_1=1$, $\lambda_2=\operatorname{tr}P-1=1.3-1=0.3$. **Detailed balance:** $\pi_1 P_{12}=\tfrac47(0.3)=\tfrac{12}{70}$ equals $\pi_2 P_{21}=\tfrac37(0.4)=\tfrac{12}{70}$ — the chain is reversible. **Convergence:** by rung 4, $\lVert\pi_t-\pi\rVert\sim|\lambda_2|^t=0.3^t$; the spectral gap $1-0.3=0.7$ is large, so the chain mixes in a handful of steps (e.g. $0.3^5\approx0.0024$).

- [ ] **Step 3: Verify both, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
Q=np.array([[0,.5,.5],[0,0,.3],[0,0,0.]]); R=np.array([[0,0],[.2,.5],[.4,.6]])
N=np.linalg.inv(np.eye(3)-Q); B=N@R
print('(a) conv from Start',round(B[0,0],4))
P=np.array([[.7,.3],[.4,.6]]); w,v=np.linalg.eig(P.T)
pi=np.real(v[:,np.argmin(abs(w-1))]); pi/=pi.sum()
print('(b) lam2',round(sorted(np.real(np.linalg.eigvals(P)))[0],4),'pi',np.round(pi,4),
      'DB',round(pi[0]*.3,4),round(pi[1]*.4,4))"
```
Expected: `(a) conv from Start 0.36`; `(b) lam2 0.3 pi [0.5714 0.4286] DB 0.1714 0.1714`.

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Worked Examples — attribution and convergence"
```

---

### Task 5: Code Tie-in

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (the `## Code Tie-in` section)

- [ ] **Step 1: Write a single self-contained ```` ```{python} ```` cell**

Model the structure on Chapter 6's Code Tie-in. Intro paragraph, one ```` ```{python} ```` block (numpy + matplotlib, seeded `rng = np.random.default_rng(...)`), closing paragraph quoting actual printed numbers. The cell, with small helpers, must:
1. **Attribution chain:** build $Q,R$; compute the fundamental matrix $N=(I-Q)^{-1}$, absorption probs $B=NR$, print the conversion probability from Start ($0.36$); implement a `removal_effect(channel)` that deletes a transient state (its inbound mass goes to Null) and recompute, printing the removal effects for Display and Search ($\approx0.444$, $\approx0.722$).
2. **Stationary two ways:** for an irreducible aperiodic $P$ (e.g. the 2-state $[[0.7,0.3],[0.4,0.6]]$ or a 3-state irreducible chain), compute $\pi$ as the normalized left eigenvector for eigenvalue 1 **and** by power iteration $\pi_{t+1}=\pi_t P$ until convergence; **assert** they match (`np.allclose`).
3. **Convergence experiment:** iterate $\pi_t$ from a few starting distributions; plot $\lVert\pi_t-\pi\rVert_1$ on a log-$y$ axis; overlay the line $|\lambda_2|^t$ and confirm the empirical decay slope matches $\log|\lambda_2|$ (print the measured vs predicted rate).
4. **Ergodic-theorem demo:** simulate one long trajectory (sample states by the row distributions); plot the running time-average of an indicator/`f` converging to $\mathbb E_\pi[f]=\sum_x f(x)\pi_x$; print the final time-average vs the analytic value.
5. **Metropolis preview:** given a *target* discrete $\pi^\star$ on a few states, build a Metropolis chain (proposal + accept ratio $\min(1,\pi^\star_j/\pi^\star_i)$ for a symmetric proposal — note the normalizer cancels), run it, and show the empirical visit frequencies converge to $\pi^\star$ (`np.allclose` within tolerance). A concrete bridge to Chapter 8.

End figures with `plt.show()`. Closing paragraph quotes the conversion prob, the two removal effects, the matched convergence rate, and the Metropolis-recovers-$\pi^\star$ confirmation.

- [ ] **Step 2: Extract the cell and run headless**

Save the cell body to `/tmp/ch7_code.py` and run:
```bash
MPLBACKEND=Agg python3 /tmp/ch7_code.py
```
Expected: runs top-to-bottom; prints conversion $0.36$, removal effects $\approx0.444/0.722$, matched convergence rate, and the Metropolis target recovery. Make the prose match the actual output.

- [ ] **Step 3: Commit**

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Code Tie-in — attribution, convergence, ergodic, Metropolis preview"
```

---

### Task 6: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (the `## Exercises` section with the four `###` tier headings)

Match Chapter 6's exercise style and heading text exactly: `### C -- Conceptual / Reading Comprehension`, `### B -- By Hand`, `### P -- Prove It`, `### A -- Applied / Code`. Each problem fully self-contained (state all matrices/data inline). **No links to solutions.**

- [ ] **Step 1: Tier C (3 problems)**
- **C1.** State the Markov property precisely and explain what it does and does not assume (memorylessness given the present); give an MMM example where it is reasonable and one where it is questionable.
- **C2.** Explain irreducibility and aperiodicity in plain terms and why **both** are needed for $P^n\to\mathbf 1\pi$; explain what the spectral gap $1-|\lambda_2|$ and the mixing time control for an MCMC run ("how long is long enough").
- **C3.** Explain why **detailed balance ⟹ stationarity** is the central MCMC design tool — why a local, ratio-only condition is easier to engineer than $\pi P=\pi$, and how it lets MCMC target an unnormalized posterior. Then state why the ergodic theorem licenses single-trajectory estimation.

- [ ] **Step 2: Tier B (3 problems)**
- **B1.** Given a 2-state $P=\begin{bmatrix}0.9&0.1\\0.5&0.5\end{bmatrix}$: compute $\pi$ by solving $\pi P=\pi$, the eigenvalue $\lambda_2$, and verify detailed balance. (Clean: $\pi=(5/6,1/6)$, $\lambda_2=0.4$.)
- **B2.** A 3-state absorbing chain (one transient state $T$, absorbing $A,B$) with $T\to A\,0.3$, $T\to B\,0.2$, $T\to T\,0.5$: compute the fundamental scalar $N=(1-0.5)^{-1}=2$ and absorption probabilities to $A$ and $B$ ($0.6$, $0.4$).
- **B3.** Evolve a distribution: given $P=\begin{bmatrix}0.7&0.3\\0.4&0.6\end{bmatrix}$ and $\pi_0=(1,0)$, compute $\pi_1,\pi_2$ by hand and the distance $\lVert\pi_2-\pi\rVert$ to the stationary $\pi=(4/7,3/7)$; compare to the bound $|\lambda_2|^2=0.09$.

- [ ] **Step 3: Tier P (3 problems)**
- **P1.** Prove detailed balance ⟹ stationarity (sum the balance equation; use row-stochasticity).
- **P2.** Prove that for a diagonalizable irreducible aperiodic finite chain $P^n\to\mathbf 1\pi$ with $\lVert\pi_t-\pi\rVert\le C|\lambda_2|^t$, using the spectral decomposition and Perron–Frobenius ($\lambda_1=1$ simple and strictly dominant).
- **P3.** Prove the ergodic theorem for a finite irreducible chain via the renewal/return-time argument (i.i.d. cycles between returns to a fixed state; LLN on the ratio of cycle sums to cycle lengths; the renewal-reward identity giving $\mathbb E_\pi[f]$). *Alternatively/additionally:* prove the stationary distribution of an irreducible finite chain is unique.

- [ ] **Step 4: Tier A (2 problems)**
- **A1.** Attribution + removal effect: build a customer-journey absorbing chain (≥3 touchpoints), compute conversion probabilities via $N=(I-Q)^{-1}$, and report removal-effect attribution per channel; sanity-check that removal effects rank channels sensibly.
- **A2.** Convergence + ergodic + Metropolis: (i) for an irreducible aperiodic $P$, plot $\lVert\pi_t-\pi\rVert$ vs $|\lambda_2|^t$ and confirm the rate; (ii) simulate one trajectory and show the time-average of an $f$ converges to $\mathbb E_\pi[f]$; (iii) build a Metropolis chain for a chosen discrete target $\pi^\star$ and show empirical frequencies converge to $\pi^\star$.

- [ ] **Step 5: Commit**

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Exercises — C/B/P/A tiers"
```

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/03-sampling/01-markov-chains.qmd` (append `## Summary` as the final section)

- [ ] **Step 1: Write `## Summary`**

Match Chapter 6's Summary format exactly: one-line lead, then **Key concepts** (bold lead-ins) and **Key identities** (**inline math only**). Cover:
- **Key concepts:** the Markov property and row-stochastic $P$; irreducibility + aperiodicity as the convergence hypotheses; stationary distribution as the eigenvalue-1 left eigenvector; detailed balance as the MCMC design tool; convergence governed by the spectral gap / subdominant eigenvalue (mixing time); the ergodic theorem licensing single-run estimation; the fundamental matrix for absorption/attribution; the MCMC flip (construct a chain whose stationary distribution is the posterior).
- **Key identities (inline):**
  - Markov / Chapman–Kolmogorov: $(P^n)_{ij}=P(X_{t+n}=j\mid X_t=i)$, evolution $\pi_t=\pi_0 P^t$.
  - Stationarity: $\pi P=\pi$; detailed balance: $\pi_i P_{ij}=\pi_j P_{ji}$.
  - Convergence: $P^n\to\mathbf 1\pi$, $\lVert\pi_t-\pi\rVert\le C|\lambda_2|^t$, spectral gap $1-|\lambda_2|$.
  - Ergodic theorem: $\frac1n\sum_{t=1}^n f(X_t)\to\mathbb E_\pi[f]$.
  - Absorption: $N=(I-Q)^{-1}$, $B=NR$.

- [ ] **Step 2: Commit**

```bash
git add parts/03-sampling/01-markov-chains.qmd
git commit -m "feat(ch07): Summary — key concepts and identities"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (insert before the closing `:::` at line 1008)

- [ ] **Step 1: Append the Chapter 7 solutions block**

Immediately before the final `:::` that closes the `show-solutions` visible div, add a `## Chapter 7 — Markov Chains` heading followed by full solutions to every exercise from Task 6 (C1–C3, B1–B3, P1–P3, A1–A2). Match the format of the existing `## Chapter 6 — Hierarchical Models` block: em-dash sub-headings `### C — …` etc., bold labels (`**C1.**`), full worked math for B and P, plain ```` ```python ```` reference sketches (NOT `{python}`) plus a qualitative takeaway for A. Each proof ends `$\blacksquare$`. Verify every numeric answer against NumPy (B1 $\pi=(5/6,1/6)$, $\lambda_2=0.4$; B2 $N=2$, absorption $0.6/0.4$; B3 $\pi_1=(0.7,0.3)$, $\pi_2=(0.61,0.39)$, $\lVert\pi_2-\pi\rVert$ vs $0.09$).

- [ ] **Step 2: Verify by-hand answers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
# B1
P=np.array([[.9,.1],[.5,.5]]); w,v=np.linalg.eig(P.T); pi=np.real(v[:,np.argmin(abs(w-1))]); pi/=pi.sum()
print('B1 pi',np.round(pi,4),'lam2',round(sorted(np.real(np.linalg.eigvals(P)))[0],4))
# B3
P2=np.array([[.7,.3],[.4,.6]]); p0=np.array([1.,0]); p1=p0@P2; p2=p1@P2
pi2=np.array([4/7,3/7]); print('B3 p1',np.round(p1,4),'p2',np.round(p2,4),'dist',round(np.abs(p2-pi2).sum(),4))"
```
Expected: `B1 pi [0.8333 0.1667] lam2 0.4`; `B3 p1 [0.7 0.3] p2 [0.61 0.39] dist ...` (a small number $\le 0.18$).

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch07): appendix solutions for Markov Chains exercises"
```

---

### Task 9: Final whole-chapter review

**Files:**
- Review: `parts/03-sampling/01-markov-chains.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: KaTeX / structure lint**

```bash
cd /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch7-markov-chains
grep -n 'begin{align}' parts/03-sampling/01-markov-chains.qmd || echo "no bare align — good"
grep -c '^\$\$' parts/03-sampling/01-markov-chains.qmd   # expect even
grep -nE '^\## ' parts/03-sampling/01-markov-chains.qmd  # expect the 6 template headings in order
grep -nE '^# ' parts/03-sampling/01-markov-chains.qmd    # expect '# Markov Chains'
grep -c '^:::$' appendix/solutions.qmd   # expect 2 (file still closes cleanly)
```
Confirm no bare `\begin{align}`, even `$$` count, the six template headings in order, and citations only the valid keys.

- [ ] **Step 2: Re-run the code cell headless**

```bash
MPLBACKEND=Agg python3 /tmp/ch7_code.py && echo "CODE CELL OK"
```
Expected: `CODE CELL OK`, conversion $0.36$, removal effects $\approx0.444/0.722$, matched convergence rate, Metropolis target recovered.

- [ ] **Step 3: Through-line check (manual)**

Confirm every theory subsection ties back to the attribution / MCMC-flip through-line, the three named proofs are present and correct, Chapter 1's spectral results are reused by reference, and the Chapter 6 Gibbs-as-Markov-chain promise is explicitly paid in rung 6. Fix gaps inline and commit if changed.

---

## Self-Review (completed by plan author)

- **Spec coverage:** Motivation (T1); rungs 1–3 incl. detailed-balance proof (T2); rungs 4–6 incl. convergence keystone, ergodic theorem, MCMC flip (T3); both worked examples (T4); code tie-in with the 5 required behaviors (T5); C/B/P/A exercises (T6); auto-included Summary (T7); appendix solutions (T8); render/lint gate (T9). All spec success criteria mapped.
- **Placeholder scan:** none — every task carries the actual formulas and NumPy-verified anchor numbers.
- **Consistency:** row-vector convention $\pi_{t+1}=\pi_t P$ and row-stochastic $P$ throughout; the attribution chain $Q,R$ and the 2-state $P$ are identical in T2/T4/T5; the three proofs (detailed balance, spectral convergence, ergodic) are stated identically across T3/T6/T7; worked/appendix numbers are NumPy-verified. The fundamental matrix $N=(I-Q)^{-1}$ defined in T2 is reused unchanged in T4/T5/T7/T8.
- **Note:** the real render gate is CI `quarto render` (HTML + PDF) on the PR — Quarto is not installed locally, so local verification is the Python cell + lint checks above.
