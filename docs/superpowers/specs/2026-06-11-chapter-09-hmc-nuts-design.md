# Part III, Ch 3 — Hamiltonian Monte Carlo & NUTS: Design Spec

**Date:** 2026-06-11
**Branch:** `ch9-hmc-nuts`
**File to author:** `parts/03-sampling/03-hmc-nuts.qmd`
**Canonical anchors:** Neal (2011) *MCMC using Hamiltonian dynamics*; Betancourt (2017) *A Conceptual Introduction to HMC*; Hoffman & Gelman (2014) *The No-U-Turn Sampler*.

## Goal

Replace the scaffolded stub with the third chapter of Part III. Chapter 8 ended with random-walk
Metropolis crawling along the correlated MMM posterior ridge. This chapter builds **Hamiltonian
Monte Carlo** — which augments the parameters with momentum and follows the gradient of $\log\pi$
to make distant, high-acceptance proposals — and **NUTS**, which removes the trajectory-length
knob. The payoff is efficient sampling of the high-dimensional, correlated posterior (the engine
inside Stan and PyMC), and HMC's canonical hard case is exactly Chapter 6's hierarchical funnel.
Scope is **MMM-targeted core**. Deliverable this cycle is a **full chapter draft**: prose, proofs,
worked examples, a runnable Python tie-in, full C/B/P/A exercises with appendix solutions, and the
closing **Summary** (auto-included).

Template (fixed): Motivation → Theory & Proofs → Worked Examples → Code Tie-in →
Exercises (C/B/P/A) → Summary.

## Through-line (anchor)

Chapter 8's random-walk Metropolis explores a correlated posterior by blind, isotropic steps that
must be tiny to be accepted, so it crawls along the ridge produced by collinear channel spends.
HMC fixes the *direction* problem at the source: augment the position $q=\beta$ with a **momentum**
$p$, define a potential energy $U(q)=-\log\pi(q)$ and kinetic energy $K(p)=\tfrac12 p^\top M^{-1}p$,
and simulate **frictionless physics** on the landscape. A particle launched with a random momentum
glides along the contours of the posterior, traveling far in a single proposal while staying in
high-probability regions — because energy is conserved. The recurring demonstration: on a strongly
correlated Gaussian (the Chapter 8 hard case), HMC delivers an effective sample size **orders of
magnitude** larger than random-walk Metropolis for the same number of iterations. NUTS removes the
last manual tuning knob (how long to simulate), and the chapter closes on HMC's own failure mode —
**divergences** on the hierarchical **funnel** (Chapter 6) — cured by reparameterization, handing
the diagnostics to Chapter 10.

## Section-by-section content

### Motivation
Open on Chapter 8's diagnosis: random-walk Metropolis is *uninformed* — its proposal is symmetric
in every direction and ignores the shape of $\pi$, so on a correlated posterior it makes mostly
cross-ridge proposals that are rejected and tiny along-ridge steps that are accepted, collapsing the
effective sample size. The missing information is the **gradient** $\nabla\log\pi$. Driving question
(verbatim emphasis): *"You can compute the gradient of the log-posterior — so why propose blindly?
Use it to glide along the distribution instead of stumbling across it."* Introduce the physical idea
(roll a frictionless puck on the surface $-\log\pi$, given a random kick) and name the deliverables:
Hamiltonian dynamics and why its conservation laws make a long trajectory a valid proposal; the
leapfrog integrator; the HMC accept step; NUTS as the automation that made HMC the default in modern
probabilistic programming; and the funnel as the geometry that still defeats a naive HMC. Anchor to
the correlated MMM posterior $p(\beta\mid y)$.

### Theory & Proofs
A ladder; full proofs of the conservation/symplectic/invariance results, NUTS conceptual-but-precise.

1. **Augmenting with momentum.** Recap the random-walk failure. Introduce the **augmented target**:
   position $q\in\mathbb R^d$ with **potential energy** $U(q)=-\log\pi(q)$ (up to a constant), an
   auxiliary **momentum** $p\in\mathbb R^d$ with **kinetic energy** $K(p)=\tfrac12 p^\top M^{-1}p$
   (mass matrix $M$, default $I$), and the **Hamiltonian** $H(q,p)=U(q)+K(p)$. The joint density is
   $\pi(q,p)\propto e^{-H(q,p)}=e^{-U(q)}e^{-K(p)}$, so $q$ and $p$ are independent and the **marginal
   over $p$ is exactly $\pi(q)$** — sampling $(q,p)$ and discarding $p$ samples the posterior. The plan:
   move through $(q,p)$-space along trajectories that keep $H$ (nearly) constant, so proposals stay in
   high-probability regions.
2. **Hamiltonian dynamics** (**proof**). **Hamilton's equations** $\dot q=\partial H/\partial p=M^{-1}p$,
   $\dot p=-\partial H/\partial q=-\nabla U(q)$. Prove the three properties that make a trajectory a valid
   proposal: (i) **conservation of energy**, $\frac{d}{dt}H(q,p)=0$ along solutions (full proof via the
   chain rule — the cross terms cancel); (ii) **volume preservation** (Liouville's theorem — the flow has
   divergence-free velocity field, so phase-space volume is preserved and the proposal needs **no Jacobian
   correction**); (iii) **time-reversibility** (negating $p$ and running backward retraces the path). State
   that exact Hamiltonian flow with momentum negation is therefore a measure-preserving, reversible
   involution — a perfect Metropolis proposal with acceptance 1. The catch: we cannot solve the equations
   exactly.
3. **The leapfrog integrator** (**proof**). Discretize with step size $\varepsilon$ via the **leapfrog
   (Störmer–Verlet)** scheme: a half-step momentum update, a full-step position update, a half-step
   momentum update,
   $$p_{1/2}=p-\tfrac{\varepsilon}{2}\nabla U(q),\quad q'=q+\varepsilon M^{-1}p_{1/2},\quad p'=p_{1/2}-\tfrac{\varepsilon}{2}\nabla U(q').$$
   **Prove** leapfrog is (i) **volume-preserving** — each of the three sub-steps is a **shear** that changes
   one block using only the other, so its Jacobian is triangular with unit diagonal, determinant 1, and the
   composition has Jacobian 1; and (ii) **reversible** — negating $p$, applying leapfrog, and negating $p$
   again returns the start. Energy is **not** exactly conserved: the local error is $O(\varepsilon^3)$ per
   step and the trajectory error $O(\varepsilon^2)$, but it does not drift (symplecticity), so $\Delta H$
   stays small and bounded. Contrast Euler, which is neither volume-preserving nor stable. This $O(\varepsilon^2)$
   energy error is exactly what the Metropolis step in Rung 4 corrects.
4. **The HMC algorithm** (**KEYSTONE proof**). One HMC step: (a) **resample momentum** $p\sim\mathcal N(0,M)$;
   (b) run $L$ **leapfrog** steps from $(q,p)$ to a proposal $(q^\*,p^\*)$; (c) **accept** $q^\*$ with
   probability $\min\!\big(1,\,e^{-(H(q^\*,p^\*)-H(q,p))}\big)=\min(1,e^{-\Delta H})$, else keep $q$.
   **Prove** this leaves $\pi$ invariant: the momentum resample draws $p$ from its exact conditional
   $\mathcal N(0,M)$ (a Gibbs step preserving the joint); the leapfrog-plus-momentum-negation proposal is a
   **volume-preserving, reversible involution** (Rung 3), so the Metropolis–Hastings acceptance reduces to
   $\min(1,e^{-\Delta H})$ with **no Jacobian and no proposal-density ratio**, and satisfies detailed balance
   with respect to the joint $e^{-H}$ (Chapter 8's MH theorem applied to the deterministic proposal). Marginalizing
   $p$ leaves $\pi(q)$ invariant. Because leapfrog nearly conserves $H$, $\Delta H\approx 0$ and the acceptance
   probability is **near 1 even for distant proposals** — the decisive contrast with random walk. State the two
   tuning knobs: step size $\varepsilon$ (set by a target acceptance, ~0.65–0.8) and trajectory length $L$
   (number of leapfrog steps). Cite `[@neal2011]`.
5. **NUTS: removing the trajectory-length knob** (conceptual but precise). Choosing $L$ is delicate: too small
   and HMC degenerates toward a random walk; too large and the trajectory makes a **U-turn**, wasting
   computation looping back. The **No-U-Turn Sampler** builds the trajectory adaptively, **doubling** its length
   (forward or backward in time) until the endpoints begin to approach each other — the **U-turn criterion**
   $(q^+-q^-)\cdot p^-<0$ or $(q^+-q^-)\cdot p^+<0$ (the trajectory has started to retrace) — then samples a
   state from the trajectory in a way that preserves reversibility and detailed balance. State the recursive
   tree-doubling and the slice/multinomial selection at the conceptual level; **cite** `[@hoffman2014]` for the
   full algorithm and its correctness proof. Describe **step-size adaptation by dual averaging** during warm-up,
   targeting a chosen acceptance statistic. Note: NUTS + dual-averaged $\varepsilon$ is what **Stan and PyMC run
   by default**, which is why modern MMM is fit this way.
6. **Divergences, the funnel, and the bridge to Chapter 10.** HMC's distinctive failure mode is the **divergent
   transition**: in a region of high curvature the fixed step size $\varepsilon$ makes leapfrog unstable, $\Delta H$
   blows up, the proposal is rejected, and — crucially — the sampler **cannot reach** part of the distribution,
   biasing the result. The canonical hard geometry is **Neal's funnel**, $v\sim\mathcal N(0,3^2)$,
   $x_i\mid v\sim\mathcal N(0,e^{v})$ — exactly the shape of a **hierarchical variance parameter** and its
   coefficients (Chapter 6): the neck of the funnel has curvature that no single $\varepsilon$ handles. The fix is
   **non-centered reparameterization** — sample $\tilde x_i\sim\mathcal N(0,1)$ and set $x_i=\tilde x_i e^{v/2}$, so
   the geometry the sampler sees is decoupled and flat. Divergences are not a nuisance but a **trustworthy alarm**
   that the sampler is missing mass. This hands off to Chapter 10, which develops the formal convergence diagnostics
   ($\hat R$, effective sample size, divergence counts) that detect these pathologies in practice. Cite `[@betancourt2017]`.

### Worked Examples
- (a) **Hamiltonian dynamics for a 1-D Gaussian, by hand.** Target $\pi(q)\propto e^{-q^2/2}$, so $U(q)=\tfrac12 q^2$,
  $\nabla U=q$, and with $M=1$ Hamilton's equations are $\dot q=p$, $\dot p=-q$ — the **harmonic oscillator** with
  solution $q(t)=q_0\cos t+p_0\sin t$, $p(t)=p_0\cos t-q_0\sin t$. Verify $H=\tfrac12(q^2+p^2)$ is conserved (it is
  $\tfrac12(q_0^2+p_0^2)$ for all $t$). From $(q_0,p_0)=(1,0)$: at $t=\pi/2$ the state is $(0,-1)$, at $t=\pi$ it is
  $(-1,0)$ — a single trajectory carries the particle clear across the distribution, the geometric reason HMC
  decorrelates so fast, versus a random-walk step of size $O(\varepsilon)$.
- (b) **One leapfrog step, by hand.** Same target, $\varepsilon=0.3$, start $(q,p)=(1,0)$. Half-step momentum
  $p_{1/2}=0-\tfrac{0.3}{2}\cdot 1=-0.15$; full-step position $q'=1+0.3\cdot(-0.15)=0.955$; half-step momentum
  $p'=-0.15-\tfrac{0.3}{2}\cdot 0.955=-0.29325$. Energies $H=\tfrac12(1)=0.5$ and $H'=\tfrac12(0.955^2+0.29325^2)=0.49901$,
  so $\Delta H=-0.00099$ — third-order small, the controlled error leapfrog leaves for the Metropolis step. The
  acceptance probability $\min(1,e^{-\Delta H})=1$ (the proposal slightly *lowered* the energy). NumPy-verified.

### Code Tie-in
NumPy, runnable, headless (`MPLBACKEND=Agg python3`), figures end `plt.show()`:
- **HMC on a strongly correlated 2-D Gaussian** (the Chapter 8 hard case, $\rho=0.95$). Implement the potential
  $U(q)=\tfrac12 q^\top\Sigma^{-1}q$ and its gradient, the **leapfrog** integrator, and the full HMC step. Run HMC and
  random-walk Metropolis for the same number of iterations on the **same target**, and print the **effective sample
  size** of each (HMC should dominate by ~2 orders of magnitude — the verified anchor: HMC ESS ≈ 4000 vs RWM ESS ≈ 24
  on $\rho=0.95$). The money plot: side-by-side trace/scatter showing HMC traversing the ridge while RWM crawls.
- **Energy conservation / acceptance:** plot $\Delta H$ per trajectory (small, mean ≈ 0) and report the HMC acceptance
  rate (~0.8–0.9); show acceptance degrading as $\varepsilon$ grows too large.
- **Neal's funnel and divergences:** sample the funnel $v\sim\mathcal N(0,3^2)$, $x\mid v\sim\mathcal N(0,e^{v})$ with a
  fixed-$\varepsilon$ HMC; flag **divergences** (trajectories where $\Delta H$ exceeds a large threshold) and show they
  cluster in the neck; then apply the **non-centered reparameterization** and show the divergences vanish and the
  $v$-marginal is recovered. This is HMC's failure mode and its cure, on Chapter 6's geometry.

### Exercises
- **C — Conceptual:** why augment with momentum and how marginalizing $p$ recovers $\pi$; the roles of energy
  conservation, volume preservation, and reversibility in making a long trajectory a valid proposal; why leapfrog and not
  Euler; what NUTS automates and what the U-turn criterion detects; what a divergence signals and why the funnel is hard.
- **B — By hand:** solve Hamilton's equations for the 1-D Gaussian (harmonic oscillator) and verify energy conservation;
  carry out one leapfrog step and compute $\Delta H$ and the acceptance; compute the HMC acceptance $\min(1,e^{-\Delta H})$
  for a given $\Delta H$.
- **P — Prove it:** energy conservation $\frac{d}{dt}H=0$ along Hamilton's equations; leapfrog is volume-preserving (the
  shear/unit-Jacobian argument) and reversible; the HMC accept step $\min(1,e^{-\Delta H})$ leaves $e^{-H}$ invariant
  (detailed balance for the deterministic, volume-preserving, reversible proposal).
- **A — Applied / code:** implement HMC for a correlated Gaussian and compare ESS to random-walk Metropolis; study
  acceptance and $\Delta H$ as functions of $\varepsilon$; reproduce divergences on Neal's funnel and cure them by
  non-centered reparameterization.

Solutions go in the shared `appendix/solutions.qmd`, gated by `show-solutions` (no inline links), heading
`## Chapter 9 — Hamiltonian Monte Carlo & NUTS` to match siblings (em-dash; `### C — …` etc.).

### Summary
Auto-included: **Key concepts** + **Key identities** (inline math only) covering the augmented target and the
$e^{-H}$ joint, Hamilton's equations, energy conservation / volume preservation / reversibility, the leapfrog
scheme, the HMC accept rule $\min(1,e^{-\Delta H})$ and why acceptance is near 1, NUTS and the U-turn criterion,
and divergences / the funnel / non-centered reparameterization.

## Rigor level

Third chapter of Part III — **Neal / Betancourt style**: physical intuition married to full proofs of the core
results. **Prove in full**: (1) energy conservation $\frac{d}{dt}H=0$ along Hamilton's equations; (2) leapfrog
volume preservation (unit-Jacobian shear composition) and reversibility; (3) HMC leaves $\pi$ invariant (the
acceptance $\min(1,e^{-\Delta H})$ is the Metropolis–Hastings rule for a volume-preserving reversible involution,
satisfying detailed balance with respect to $e^{-H}$). **State and cite** Liouville's theorem in its general form,
the symplectic-integrator error theory ($O(\varepsilon^2)$ global, no drift), NUTS and its correctness
(`[@hoffman2014]`), and the optimal HMC/NUTS acceptance targets. Reuse Chapter 8 (the Metropolis–Hastings
detailed-balance theorem, ESS) and Chapter 6 (the hierarchical funnel) explicitly. Keep the differential geometry
to the intuition level (cite `[@betancourt2017]`); no measure theory beyond Chapter 7.

## Conventions & constraints

- Quarto `.qmd`, KaTeX (HTML) / LaTeX (PDF); number-sections on. Use `aligned` inside `$$`, never bare
  `\begin{align}`; `$$` delimiters on their own lines so display-math counts stay even.
- Citations via existing `references.bib` keys: `@neal2011`, `@betancourt2017`, `@hoffman2014`, `@gelman2013`,
  `@robert2004`. **Do not invent keys.**
- Code original, minimal, runnable against pinned `requirements.txt` (numpy, scipy, matplotlib); verify every
  numeric claim with NumPy before commit; the cell must run headless under `MPLBACKEND=Agg python3` and end figures
  with `plt.show()`. Seed all RNGs.
- Match the established voice of Ch 1–8. Reuse Chapter 8's notation (target $\pi$, the Metropolis acceptance and
  detailed balance, ESS) and the posterior notation $p(\beta\mid y)$; here position $q=\beta$, momentum $p$,
  Hamiltonian $H=U+K$ with $U=-\log\pi$.
- Git identity in the fresh worktree, and **commit from the worktree with explicit `git -C <worktree>`** to avoid a
  cwd pitfall (a bare `cd` in a compound command may not take effect, sending the commit to `main` where a hook
  blocks it): `git -C <wt> config user.name "jlh530i" && git -C <wt> config user.email "jlh530i@gmail.com"`.

## Out of scope (YAGNI)

Riemannian-manifold HMC / RHMC and position-dependent mass matrices; symplectic-integrator theory beyond the leapfrog
scheme; the full NUTS recursive tree-building, slice-sampling auxiliary, and correctness proof (state conceptually,
cite); the differential-geometric foundations of HMC beyond intuition; formal optimal-acceptance / optimal-scaling
derivations (state the targets); dual-averaging internals beyond a one-paragraph description; and the formal
convergence diagnostics $\hat R$ and effective-sample-size estimation theory (Chapter 10 — though divergences and a
practical ESS computation appear here).

## Success criteria

- Chapter renders cleanly (`quarto render`, HTML + PDF) — verified in CI on the PR (Quarto not installed locally).
- Every theory subsection ties back to the glide-don't-stumble / efficiency / funnel through-line.
- The three named proofs (energy conservation; leapfrog volume preservation + reversibility; HMC invariance via the
  $\min(1,e^{-\Delta H})$ rule) are present and correct, with NUTS developed conceptually-but-precisely and cited.
- The code tie-in runs top-to-bottom headless: HMC's ESS dominates random-walk Metropolis on the correlated Gaussian,
  energy conservation / acceptance are shown, and the funnel's divergences are reproduced and then cured by
  non-centered reparameterization.
- All four exercise tiers populated with matching appendix solutions.
- A Summary section closes the chapter, consistent with Ch 1–8.
