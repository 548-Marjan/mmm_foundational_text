# Part III, Ch 2 — Markov Chain Monte Carlo: Design Spec

**Date:** 2026-06-11
**Branch:** `ch8-mcmc`
**File to author:** `parts/03-sampling/02-mcmc.qmd`
**Canonical anchors:** Robert & Casella (2004) *Monte Carlo Statistical Methods*; Gelman et al. *BDA3* (2013).

## Goal

Replace the scaffolded stub with the second chapter of Part III. Chapter 7 proved that an
irreducible, aperiodic chain satisfying **detailed balance** with respect to a target $\pi$
converges to $\pi$ and that its time-averages estimate $\mathbb E_\pi[f]$ — but left the
**construction** of such a chain open. This chapter builds the two workhorse algorithms,
**Metropolis–Hastings** and **Gibbs**, proves they target $\pi$, and pays off with a working
sampler for a **non-conjugate MMM posterior** that Chapter 5's "wall" rendered intractable in
closed form. Scope is **MMM-targeted core**. Deliverable this cycle is a **full chapter draft**:
prose, proofs, worked examples, a runnable Python tie-in, full C/B/P/A exercises with appendix
solutions, and the closing **Summary** (auto-included).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

## Through-line (anchor)

Chapter 7 ended with a recipe stated but not built: *to sample a posterior $\pi=p(\beta\mid y)$,
construct a chain that is irreducible, aperiodic, and reversible with respect to $\pi$; then run
it.* This chapter supplies the construction. Two anchors carry it:

1. The **correlated bivariate Gaussian** $\mathcal N(0,\Sigma)$, $\Sigma=\big[\begin{smallmatrix}1&\rho\\\rho&1\end{smallmatrix}\big]$ — the clean,
   fully analytic vehicle. Its Gibbs full conditionals are one-dimensional Gaussians, and the
   $\rho\to1$ slow-mixing pathology is the concrete motivation for Chapter 9's Hamiltonian Monte
   Carlo.
2. The **non-conjugate MMM posterior** that Chapter 5 declared intractable (a sign-constrained or
   saturated channel coefficient — a posterior known only up to its normalizing constant). MCMC
   samples it directly, delivering the payoff: *the posterior we could not write down, we can now
   sample.*

The recurring slogan: $\pi$ enters every accept/reject decision **only through ratios**
$\pi(y)/\pi(x)$, so the unknown normalizer cancels — which is exactly why MCMC works on posteriors
we cannot normalize. Gibbs is revealed as the Chapter 6 hierarchical sampler, now formally
justified; random-walk Metropolis is shown to mix slowly in high dimensions, handing the baton to
Chapter 9.

## Section-by-section content

### Motivation
Open on the gap Chapter 7 left: the ergodic theorem guarantees that *if* we can build a reversible
chain for $\pi$ then averaging its trajectory estimates posterior expectations — but Chapter 7
never built one. Meanwhile Chapter 5 ended at a wall: introduce a non-negativity prior or an
adstock/saturation nonlinearity and the posterior loses its closed form; even its normalizing
constant $\int p(y\mid\beta)p(\beta)\,d\beta$ is an intractable integral. Driving question (verbatim
emphasis): *"The posterior is a formula you can evaluate at any $\beta$ but cannot integrate, normalize,
or sample from directly — so how do you draw from it?"* The answer is to manufacture a Markov chain
whose stationary distribution is that posterior, using only pointwise evaluations of the unnormalized
density. Name the two algorithms and the efficiency questions (acceptance rate, autocorrelation,
effective sample size) the chapter develops, and that random-walk Metropolis's limitations motivate
Chapter 9.

### Theory & Proofs
A ladder, each rung an MMM payoff:

1. **The Monte Carlo principle and why MCMC.** Monte Carlo estimation $\mathbb E_\pi[f]\approx\frac1n\sum_i f(X_i)$
   with the CLT error $\sigma_f/\sqrt n$ — **dimension-free** in the number of samples, the reason
   sampling beats quadrature in high dimensions. The obstacle: we cannot draw i.i.d. from a complex
   posterior. The resolution: generate **dependent** draws from a Markov chain with stationary
   distribution $\pi$; Chapter 7's ergodic theorem certifies the time-average still converges to
   $\mathbb E_\pi[f]$. The task reduces to **building the chain**.
2. **Metropolis–Hastings** (**KEYSTONE proof**). The algorithm: from state $x$, draw a proposal
   $y\sim q(\cdot\mid x)$ and accept ($X_{t+1}=y$) with probability
   $$\alpha(x,y)=\min\!\left(1,\ \frac{\pi(y)\,q(x\mid y)}{\pi(x)\,q(y\mid x)}\right),$$
   else stay ($X_{t+1}=x$). **Prove** the resulting transition kernel satisfies **detailed balance**
   $\pi(x)P(x,y)=\pi(y)P(y,x)$ with respect to $\pi$ (hence $\pi$ is stationary, by Chapter 7 rung 3),
   via the standard $\min$ identity $\pi(x)q(y\mid x)\alpha(x,y)=\min(\pi(x)q(y\mid x),\pi(y)q(x\mid y))$,
   which is symmetric in $x\leftrightarrow y$. Emphasize $\pi$ enters only through the **ratio**
   $\pi(y)/\pi(x)$, so the normalizer cancels. Special cases: **symmetric proposal** ($q(y\mid x)=q(x\mid y)$)
   gives the **Metropolis** rule $\alpha=\min(1,\pi(y)/\pi(x))$; **random-walk** Metropolis
   ($y=x+\varepsilon$); the **independence sampler**. Note irreducibility/aperiodicity hold for any
   proposal with positive density on the support (cite `[@robert2004]`), so Chapter 7's convergence
   and ergodic theorems apply.
3. **The Gibbs sampler** (**proof of invariance**). Cycle through coordinate blocks, replacing
   $x_i$ by a draw from its **full conditional** $p(x_i\mid x_{-i})$. **Prove** each full-conditional
   update leaves $\pi$ invariant (if $x\sim\pi$ then after the update $x\sim\pi$, since the conditional
   times the marginal of the rest reconstitutes the joint), and show **Gibbs is Metropolis–Hastings
   with proposal $=$ the full conditional and acceptance probability identically 1** (the MH ratio
   collapses to 1). This is exactly **Chapter 6's hierarchical sampler**, now formally justified.
   Systematic-scan vs random-scan; requirement that the full conditionals be available in closed form.
4. **Efficiency: acceptance, autocorrelation, effective sample size.** The random-walk step size sets
   a tradeoff: too small ⇒ high acceptance but tiny moves and high autocorrelation; too large ⇒ most
   proposals rejected, the chain stalls. Define the **lag-$k$ autocorrelation** $\rho_k$ of the chain,
   the integrated autocorrelation, and the **effective sample size**
   $$\text{ESS}=\frac{n}{1+2\sum_{k\ge1}\rho_k},$$
   with the MCMC standard error $\sigma_f/\sqrt{\text{ESS}}$. State (and lightly motivate, citing
   `[@robert2004]`, `[@gelman2013]`) the heuristic optimal acceptance rates (~0.234 for high-dimensional
   random-walk Metropolis, ~0.44 in one dimension) without a full proof. (Formal multi-chain convergence
   diagnostics — $\hat R$, rank-normalized ESS — are deferred to Chapter 10.)
5. **Burn-in, mixing, and the limits — the bridge to HMC.** **Burn-in / warm-up**: discard the
   pre-convergence transient (Chapter 7's geometric approach to $\pi$). Why random-walk Metropolis
   **mixes slowly**: in high dimensions a blind random walk must take $O(d)$ tiny steps to traverse the
   distribution, and on **correlated** or **funnel-shaped** posteriors (the MMM posterior has strongly
   correlated coefficients) the lag-1 autocorrelation of a Gibbs/RW coordinate approaches 1 as the
   correlation $\rho\to1$ — concretely $\rho^2$ for the bivariate-Gaussian Gibbs chain. This is the
   structural inefficiency that **Hamiltonian Monte Carlo** (Chapter 9) removes by using the gradient of
   $\log\pi$ to propose distant, high-acceptance moves. End the Theory section on this baton-pass.

### Worked Examples
- (a) **Metropolis acceptance by hand.** Unnormalized target $\pi(x)\propto e^{-x^2/2}$ (a standard
  normal with the $1/\sqrt{2\pi}$ deliberately dropped to show it is unneeded), symmetric random-walk
  proposal. Compute the acceptance probability $\alpha=\min(1,e^{-(y^2-x^2)/2})$ for three concrete
  moves: $x=0\to y=1$ gives $\alpha=e^{-1/2}\approx0.607$; $x=1\to y=0.5$ gives $\alpha=1$ (a move to
  higher density is always accepted); $x=0.5\to y=2$ gives $\alpha\approx0.153$. Interpret: uphill moves
  are automatic, downhill moves are accepted probabilistically in proportion to the density ratio, and
  the normalizer never appears.
- (b) **Gibbs for the correlated bivariate normal.** Target $\mathcal N(0,\Sigma)$ with $\Sigma=\big[\begin{smallmatrix}1&\rho\\\rho&1\end{smallmatrix}\big]$,
  $\rho=0.8$. Derive the full conditionals $x\mid y\sim\mathcal N(\rho y,\,1-\rho^2)=\mathcal N(0.8y,0.36)$
  and symmetrically for $y\mid x$. Carry two Gibbs sweeps by hand from a stated start. State the classic
  mixing fact: the lag-1 autocorrelation of either coordinate equals $\rho^2=0.64$, so as $\rho\to1$ the
  chain crawls — the motivation for Chapter 9. (NumPy-verified: conditional variance $0.36$; simulated
  recovery of $\mathrm{corr}=0.80$; empirical lag-1 autocorrelation $\approx0.64$.)

### Code Tie-in
NumPy/SciPy, runnable, headless (`MPLBACKEND=Agg python3`), figures end `plt.show()`:
- **Random-walk Metropolis** on the correlated 2-D Gaussian: implement the sampler; show trace plots; sweep
  the proposal step size (too small / good / too large) and tabulate **acceptance rate** and **ESS** for each,
  exhibiting the tuning tradeoff; confirm the sample mean/covariance match the target.
- **Gibbs** on the same target: implement the closed-form conditional sweeps; confirm recovery of $\Sigma$
  (correlation $\approx\rho$); measure the lag-1 autocorrelation and confirm it tracks $\rho^2$; show mixing
  degrade as $\rho\to0.99$.
- **The payoff — sample a non-conjugate posterior.** A 1-D MMM posterior with **no closed form**: e.g. a
  channel-effect coefficient with a Gaussian likelihood times a **half-normal (non-negativity) prior**, or a
  logistic/saturation link — evaluable only up to its normalizing constant. Run Metropolis on the unnormalized
  log-posterior and show the sampler recovers the (truncated/non-Gaussian) posterior, overlaying the MCMC
  histogram on a fine-grid normalized reference. This is the posterior Chapter 5 could not handle.

### Exercises
- **C — Conceptual:** the Monte Carlo principle and why *dependent* MCMC draws still estimate
  $\mathbb E_\pi[f]$; the MH acceptance ratio and why the normalizer cancels (so unnormalized posteriors are
  fine); Gibbs as the special case of MH with acceptance 1, and its relation to Chapter 6; what acceptance
  rate, autocorrelation, and ESS measure; why random-walk MH struggles in high/correlated dimensions (→ HMC).
- **B — By hand:** compute MH acceptance probabilities for a given unnormalized target and proposal (incl. an
  asymmetric proposal where the $q$-ratio does not cancel); derive a Gibbs full conditional for the bivariate
  normal (and/or a 2-block Gaussian); compute ESS from given autocorrelations.
- **P — Prove it:** Metropolis–Hastings satisfies detailed balance with respect to $\pi$ (the $\min$ identity);
  the Gibbs full-conditional update leaves $\pi$ invariant; the symmetric-proposal Metropolis acceptance is the
  special case $\alpha=\min(1,\pi(y)/\pi(x))$.
- **A — Applied / code:** implement random-walk Metropolis for a posterior, tune the step size and report
  acceptance and ESS; implement Gibbs for the bivariate normal and show mixing degrade with $\rho$; sample a
  non-conjugate posterior with MH and validate against a grid reference.

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions` (no inline links),
heading `## Chapter 8 — Markov Chain Monte Carlo` to match siblings (em-dash; `### C — …` etc.).

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only) covering the Monte Carlo estimator and
its CLT error, the Metropolis–Hastings acceptance rule and its detailed-balance property, the
normalizer-cancellation, Gibbs as full-conditional sampling (= MH with acceptance 1, = the Chapter 6 sampler),
the acceptance/autocorrelation/ESS efficiency triple, burn-in, and the random-walk slow-mixing limit pointing to
HMC.

## Rigor level

Second chapter of Part III — **Robert–Casella / BDA3 style**: algorithmic but rigorous, every result cashed out
on the bivariate-Gaussian and non-conjugate-posterior anchors, no measure theory beyond what Chapter 7 already
stated. **Prove in full**: (1) the Metropolis–Hastings kernel satisfies detailed balance with respect to $\pi$
(keystone, via the $\min$ identity), and derive the symmetric Metropolis special case; (2) the Gibbs
full-conditional update leaves $\pi$ invariant, and Gibbs $=$ MH with acceptance 1. **State** (with citations) the
ergodic guarantees inherited from Chapter 7 (irreducibility/aperiodicity of MH for positive-density proposals) and
the heuristic optimal acceptance rates; do not prove the ~0.234 result. Reuse Chapter 7 (detailed balance ⇒
stationarity, the ergodic theorem) and Chapter 6 (the hierarchical Gibbs sampler) explicitly.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned` inside `$$`, never bare
  `\begin{align}`; `$$` delimiters on their own lines so display-math counts stay even.
- Citations via existing `references.bib` keys: `@robert2004`, `@gelman2013`, `@hoff2009` (and `@gelmanhill2006`
  where the Chapter 6 Gibbs link is invoked). **Do not invent keys** (there is no dedicated Metropolis/Hastings/Geman
  entry — attribute the algorithms in prose and cite `@robert2004`/`@gelman2013`).
- Code original, minimal, runnable against pinned `requirements.txt` (numpy, scipy, matplotlib); verify every numeric
  claim with NumPy before commit; the cell must run headless under `MPLBACKEND=Agg python3` and end figures with
  `plt.show()`. Seed all RNGs.
- Match the established voice of Ch 1–7. Reuse Chapter 7's notation ($\pi$, transition kernel, detailed balance
  $\pi(x)P(x,y)=\pi(y)P(y,x)$, the ergodic theorem) and Chapter 5/6 posterior notation ($p(\beta\mid y)$, the
  unnormalized $p(y\mid\beta)p(\beta)$). The target the sampler is built for *is* the posterior.
- Git identity in the fresh worktree: `git config user.name "jlh530i" && git config user.email "jlh530i@gmail.com"`.
  Commit from the worktree with explicit `git -C <worktree>` to avoid a cwd pitfall.

## Out of scope (YAGNI)

Hamiltonian Monte Carlo and NUTS (Chapter 9); formal multi-chain convergence diagnostics — $\hat R$, rank-normalized
/ bulk-tail ESS, trace-rank plots (Chapter 10, though ESS is introduced here); reversible-jump and trans-dimensional
MCMC; adaptive-MCMC theory and its ergodicity caveats; slice sampling, particle MCMC/SMC, parallel tempering; the full
proof of optimal scaling (~0.234); and quantitative mixing-time bounds beyond the autocorrelation/ESS heuristics.

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR (Quarto not installed locally).
- Every theory subsection ties back to the construct-a-chain-for-the-posterior / efficiency / HMC-bridge through-line.
- The two named proofs (MH detailed balance; Gibbs invariance $=$ MH-with-acceptance-1) are present and correct, with
  the symmetric-Metropolis special case derived.
- The code tie-in runs top-to-bottom headless: RW-Metropolis tunes acceptance/ESS on the bivariate Gaussian, Gibbs
  recovers $\Sigma$ and shows $\rho$-dependent mixing, and Metropolis samples the non-conjugate posterior matching a
  grid reference.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–7.
