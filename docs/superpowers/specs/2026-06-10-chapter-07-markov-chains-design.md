# Part III, Ch 1 — Markov Chains: Design Spec

**Date:** 2026-06-10
**Branch:** `ch7-markov-chains`
**File to author:** `parts/03-sampling/01-markov-chains.qmd`
**Canonical anchors:** Robert & Casella (2004) *Monte Carlo Statistical Methods*; Norris (1997) *Markov Chains*; Gelman et al. *BDA3* (2013) for the MCMC pointer.

## Goal

Replace the scaffolded stub with the opening chapter of Part III. It builds the theory of
Markov chains that makes MCMC work, and closes the loop on Chapter 6, which ended by claiming
that cycling the Gibbs full conditionals converges to the posterior — the theorem proved here.
Scope is **MMM-targeted core**: every topic earns its place by a downstream use in MCMC
(Chapters 8–10) or by the genuine MMM technique of **customer-journey attribution**. Deliverable
this cycle is a **full chapter draft**: prose, proofs, worked examples, a runnable Python tie-in,
full C/B/P/A exercises with appendix solutions, and the closing **Summary** (auto-included).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

## Through-line (anchor)

Markov chains appear in MMM **twice**: as an existing attribution tool and as the engine of MCMC.
The recurring object is the **customer-journey chain** — a customer moves across channel
touchpoints (Display → Search → Email → …) until absorbed into **Convert** or **Null** — a finite
Markov chain marketers already use to attribute conversion credit. The chapter develops the theory
on this concrete chain, then performs the central **flip**:

- Part II's question was *analysis*: given a chain $P$, find its stationary distribution $\pi$ and
  long-run behavior.
- **MCMC inverts it to engineering**: given a target $\pi$ — the Chapter 5 / Chapter 6 posterior —
  **construct** a chain $P$ whose stationary distribution is $\pi$, make it irreducible and
  aperiodic, run it, and read off samples.

Detailed balance is the construction tool, convergence guarantees you reach $\pi$, and the ergodic
theorem guarantees time-averages estimate posterior expectations. Chapter 6's Gibbs cycle **is**
exactly such a chain; this chapter proves why running it works, and hands the explicit construction
(Metropolis–Hastings) to Chapter 8.

## Section-by-section content

### Motivation
Open on the two faces of Markov chains in MMM. First, marketers already model the **customer
journey** as a chain over channel touchpoints and attribute conversions via the **removal effect**
(drop a channel, recompute the conversion probability, the decrease is that channel's credit).
Second, Chapter 6 left an unpaid promise: cycling the Gibbs full conditionals was claimed to
converge to the posterior, with the justifying theorem deferred. Both rest on the same theory —
when does a chain settle into a unique long-run distribution, how fast, and when do time-averages
along one trajectory equal expectations under that distribution? Driving question (verbatim
emphasis): *"Run the Gibbs sampler for a thousand steps and average — why should that average be
the posterior mean, and how long is long enough?"* This chapter answers it and then flips the
machinery into the engine of MCMC.

### Theory & Proofs
A ladder; finite state spaces throughout, with general-state analogues stated at the end.

1. **Markov chains and transition matrices.** Finite state space $S$, the **Markov property**
   $P(X_{t+1}=j\mid X_t=i,\text{past})=P(X_{t+1}=j\mid X_t=i)=P_{ij}$, the **stochastic matrix** $P$
   (rows nonnegative, sum to 1), $n$-step transitions via **Chapman–Kolmogorov** $P^{(n)}=P^n$, and
   distribution evolution as a left-multiplication $\pi_{t+1}=\pi_t P$ (row-vector convention).
   Introduce the running **customer-journey chain** with transient touchpoints and absorbing
   Convert/Null states, and a concrete $P$.
2. **Classification.** Communicating classes; **irreducibility** (single class); **period** and
   **aperiodicity**; **recurrence/transience**; **absorbing states** and the transient/absorbing
   block form $P=\big[\begin{smallmatrix}Q&R\\0&I\end{smallmatrix}\big]$. State that absorption
   probabilities solve $(I-Q)^{-1}R$ (the fundamental matrix $N=(I-Q)^{-1}$), the engine of the
   attribution example. These conditions are the hypotheses of the convergence theorem.
3. **Stationary distributions & detailed balance** (**proof**). Define the stationary distribution
   $\pi P=\pi$ ($\pi$ a left eigenvector of $P$ for eigenvalue $1$). State existence and uniqueness
   for irreducible finite chains via **Perron–Frobenius** (a unique strictly positive $\pi$).
   Define **reversibility / detailed balance** $\pi_i P_{ij}=\pi_j P_{ji}$. **Prove**: detailed
   balance $\Rightarrow$ $\pi$ stationary (sum the balance equation over $i$ and use $\sum_i P_{ji}=1$).
   Emphasize this is the MCMC design lever — engineering detailed balance against a target is far
   easier than solving $\pi P=\pi$ directly.
4. **Convergence to stationarity** (**KEYSTONE proof, Perron–Frobenius / spectral gap**). For an
   irreducible, aperiodic finite chain, $P^n\to\mathbf{1}\pi$ (every row converges to $\pi$), so
   $\pi_t\to\pi$ from any start. Prove via the **spectral decomposition** (Chapter 1): $P$ has a
   simple eigenvalue $1$ (Perron–Frobenius) with all other $|\lambda|<1$; expand the initial
   distribution in the eigenbasis; the components along $|\lambda|<1$ decay geometrically, leaving
   the eigenvalue-1 (stationary) component. The **subdominant eigenvalue** $\lambda_2$ sets the rate:
   $\lVert\pi_t-\pi\rVert\le C\,|\lambda_2|^t$. Define the **spectral gap** $1-|\lambda_2|$ and the
   **mixing time**. Keep the proof at the diagonalizable case; cite the general (Jordan / coupling)
   argument. Reuses Chapter 1 eigen-decomposition directly.
5. **The ergodic theorem** (**proof, finite case**). For an irreducible finite chain with stationary
   $\pi$, the **time average** converges to the **space average**:
   $\frac1n\sum_{t=1}^n f(X_t)\to\mathbb E_\pi[f]$ almost surely, for any bounded $f$ and any start.
   Prove at the finite level via the **renewal / return-time** argument (returns to a fixed state
   split the trajectory into i.i.d. cycles; apply the ordinary LLN to cycle sums; positive
   recurrence gives finite mean return time). This is the law of large numbers for Markov chains —
   the result that licenses estimating a posterior expectation from a **single long run**.
6. **The flip: from analysis to MCMC.** Reverse the program. Given a target $\pi$ (the posterior),
   **construct** a $P$ with $\pi$ stationary by enforcing **detailed balance** (rung 3); ensure
   irreducibility + aperiodicity so **convergence** (rung 4) drives $\pi_t\to\pi$; then the **ergodic
   theorem** (rung 5) makes the trajectory average an estimator of $\mathbb E_\pi[f]$. Show Chapter
   6's Gibbs cycle leaves the posterior invariant, so it is precisely such a chain. Name
   **Metropolis–Hastings** as Chapter 8's general recipe. **State** the general-state-space analogues
   — transition **kernels** $P(x,\mathrm dy)$, detailed balance against a density, Harris recurrence
   and the general ergodic theorem — that continuous posteriors require, citing `[@robert2004]`,
   `[@norris1997]`. End the Theory section on this baton-pass to Chapter 8.

### Worked Examples
- (a) **Customer-journey attribution by hand.** Transient states Start, Display, Search; absorbing
  Convert, Null. Transition probabilities: Start → Display $0.5$, Start → Search $0.5$;
  Display → Search $0.3$, Display → Convert $0.2$, Display → Null $0.5$; Search → Convert $0.4$,
  Search → Null $0.6$. Compute the **conversion probability from Start** via the fundamental matrix
  (or backward substitution): $0.36$. Then compute **removal-effect attribution** — remove Search
  (its inbound journeys are lost to Null) and recompute: conversion drops to $0.10$, a removal
  effect $(0.36-0.10)/0.36\approx 0.722$ for Search; removing Display gives $0.20$, removal effect
  $\approx 0.444$. (All NumPy-verified.) A genuine MMM attribution computation.
- (b) **Convergence by hand.** Two-state chain $P=\big[\begin{smallmatrix}0.7&0.3\\0.4&0.6\end{smallmatrix}\big]$.
  Eigenvalues $1$ and $\lambda_2=0.3$ (trace $-1$); stationary distribution $\pi=(4/7,3/7)\approx(0.571,0.429)$;
  verify **detailed balance** $\pi_1 P_{12}=\pi_2 P_{21}=12/70$ (the chain is reversible); show
  $\lVert\pi_t-\pi\rVert$ decays at rate $|\lambda_2|^t=0.3^t$ — a spectral gap of $0.7$, fast mixing.

### Code Tie-in
NumPy, runnable, headless (`MPLBACKEND=Agg python3`), figures end `plt.show()`:
- Build the attribution chain; compute the stationary distribution of an irreducible variant **two
  ways** — left eigenvector of $P$ for eigenvalue 1, and power iteration $\pi_{t+1}=\pi_t P$ until
  convergence — and confirm they agree. For the absorbing attribution chain, compute conversion
  probabilities via $N=(I-Q)^{-1}$ and the **removal-effect** attribution for each channel.
- **Convergence experiment:** pick an irreducible aperiodic $P$; iterate from several starting
  distributions; plot $\lVert\pi_t-\pi\rVert_1$ on a log axis decaying at slope $\log|\lambda_2|$;
  confirm the empirical rate matches the subdominant eigenvalue (spectral gap).
- **Ergodic-theorem demo:** simulate one long trajectory; show the running time-average of an
  $f(X_t)$ converges to $\mathbb E_\pi[f]$.
- **Metropolis preview:** construct (by detailed balance) a small chain targeting a *given* discrete
  $\pi$ and show $\pi_t\to\pi$ — a concrete bridge to Chapter 8.

### Exercises
- **C — Conceptual:** the Markov property; irreducible/aperiodic in plain terms; why detailed
  balance $\Rightarrow$ stationarity is the MCMC design tool; what the spectral gap / mixing time
  controls; why the ergodic theorem licenses single-trajectory Monte-Carlo estimation.
- **B — By hand:** evolve a distribution $\pi_t P^n$; solve $\pi P=\pi$ for a 2–3 state chain; verify
  detailed balance; compute absorption/conversion probabilities for a small absorbing chain (and a
  removal effect); compute $\lambda_2$ and the implied convergence rate.
- **P — Prove it:** detailed balance $\Rightarrow$ stationarity; finite irreducible aperiodic
  $\Rightarrow P^n\to\mathbf 1\pi$ with geometric rate $|\lambda_2|$ (spectral / Perron–Frobenius,
  diagonalizable case); the ergodic theorem for finite irreducible chains (renewal-cycle argument),
  or uniqueness of the stationary distribution for an irreducible $P$.
- **A — Applied / code:** attribution chain + removal-effect attribution; convergence-rate
  experiment ($\lVert\pi_t-\pi\rVert$ vs $|\lambda_2|^t$); ergodic-theorem time-average demo; a
  Metropolis-on-discrete-$\pi$ preview converging to a chosen target.

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions` (no inline links),
heading `## Chapter 7 — Markov Chains` to match siblings (em-dash; `### C — …` etc.).

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only) covering the Markov
property, stochastic matrix and Chapman–Kolmogorov, irreducibility/aperiodicity, stationary
$\pi P=\pi$, detailed balance $\pi_i P_{ij}=\pi_j P_{ji}$, Perron–Frobenius convergence with the
spectral gap / mixing rate, the ergodic theorem (time average $\to$ expectation), the fundamental
matrix $N=(I-Q)^{-1}$ for absorption, and the MCMC flip.

## Rigor level

Opening chapter of Part III — **Norris / Robert–Casella style on finite state spaces**: clean linear
algebra, every result cashed out on the customer-journey chain, no measure theory. **Prove in full**:
(1) detailed balance $\Rightarrow$ stationarity; (2) convergence to stationarity for irreducible
aperiodic finite chains via Perron–Frobenius and the spectral decomposition (diagonalizable case),
with the geometric rate $|\lambda_2|$; (3) the ergodic theorem for finite irreducible chains via the
renewal/return-time argument. **State** (with citations) Perron–Frobenius in full generality and the
general-state-space analogues (transition kernels, Harris recurrence, the general ergodic theorem).
Reuse Chapter 1 (eigenvalues, spectral decomposition, the spectral theorem) explicitly. Anchor every
abstraction to the attribution chain and the Gibbs-sampler-as-Markov-chain story.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned` inside `$$`, never
  bare `\begin{align}`; `$$` delimiters on their own lines so display-math counts stay even.
- Citations via existing `references.bib` keys: `@robert2004`, `@norris1997`, `@gelman2013` (and
  `@strang2016`/`@axler2015` if a Chapter-1 spectral fact is invoked). **Do not invent keys.**
- Code original, minimal, runnable against pinned `requirements.txt` (numpy, scipy, matplotlib);
  verify every numeric claim with NumPy before commit; the cell must run headless under
  `MPLBACKEND=Agg python3` and end figures with `plt.show()`. Use a row-vector convention
  $\pi_{t+1}=\pi_t P$ consistently in prose and code.
- Match the established voice of Ch 1–6. Reuse Chapter 1's eigenvalue/spectral notation and the
  Part II posterior notation when describing the MCMC flip (the target $\pi$ is the posterior
  $p(\beta\mid y)$). Row-stochastic $P$; left eigenvectors for $\pi$.
- Git identity in the fresh worktree: `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"`.

## Out of scope (YAGNI)

Measure-theoretic general-state-space machinery (state and cite, do not develop); continuous-time
Markov chains and generators; Markov decision processes and reinforcement learning; coupling-from-
the-past / perfect sampling; quantitative mixing-time theory beyond the spectral gap (conductance,
cutoff); spectral theory for non-reversible/non-diagonalizable chains beyond a cited statement; and
the actual MCMC algorithms — Metropolis–Hastings and Gibbs are Chapter 8, HMC/NUTS Chapter 9,
convergence diagnostics Chapter 10.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR (Quarto not
  installed locally).
- Every theory subsection ties back to the attribution chain / MCMC-flip through-line.
- The three named proofs (detailed balance $\Rightarrow$ stationarity; finite convergence via
  spectral gap; the ergodic theorem) are present and correct, and Chapter 1's spectral results are
  reused by reference.
- The code tie-in runs top-to-bottom headless: stationary distribution two ways agree, removal-effect
  attribution computed, the convergence rate matches $|\lambda_2|$, the time-average converges, and
  the Metropolis preview reaches its target.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–6.
