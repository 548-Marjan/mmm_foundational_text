# Chapter 9 — Hamiltonian Monte Carlo & NUTS Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/03-sampling/03-hmc-nuts.qmd` with the third chapter of Part III — Hamiltonian Monte Carlo and NUTS, the gradient-based samplers that fix Chapter 8's slow random-walk mixing — and add matching appendix solutions.

**Architecture:** One Quarto `.qmd` chapter authored section-by-section against the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus a `## Chapter 9 — Hamiltonian Monte Carlo & NUTS` block appended to the shared `appendix/solutions.qmd` (gated by `show-solutions`). The chapter reuses Chapter 8 (the Metropolis–Hastings detailed-balance theorem, ESS) and Chapter 6 (the hierarchical funnel) directly. The single runnable Python cell is verified headless before commit; the real render gate is CI (`quarto render`), Quarto not being installed locally.

**Tech Stack:** Quarto book, KaTeX math, Python code cell (numpy, scipy, matplotlib pinned in `requirements.txt`). Canonical references via `references.bib` keys: `@neal2011`, `@betancourt2017`, `@hoffman2014`, `@gelman2013`, `@robert2004`.

---

## Conventions (apply to every task)

- **Math/KaTeX:** Use `aligned` inside `$$ ... $$`, never bare `\begin{align}`. Keep `$$` delimiters on their own lines so display-math counts stay even. Inline math with single `$`. Notation: position $q=\beta\in\mathbb R^d$; momentum $p$; potential $U(q)=-\log\pi(q)$; kinetic $K(p)=\tfrac12 p^\top M^{-1}p$ (default $M=I$); Hamiltonian $H(q,p)=U(q)+K(p)$; joint $\pi(q,p)\propto e^{-H}$. Reuse Chapter 8's $\pi$, Metropolis acceptance, detailed balance, ESS.
- **Voice:** Match Chapters 1–8 (read `parts/03-sampling/02-mcmc.qmd`, the immediate predecessor). Rigorous, prose-led, each theory subsection ("rung") ending by tying back to the through-line / an MMM payoff. End every full proof with `$\blacksquare$`; use `> **Theorem (name).** ...` blockquotes for named results.
- **Reuse explicitly:** Chapter 8 (MH detailed-balance theorem, ESS), Chapter 6 (the hierarchical funnel) — cite/refer rather than re-deriving.
- **Citations:** `[@neal2011]`, `[@betancourt2017]`, `[@hoffman2014]`, `[@gelman2013]`, `[@robert2004]` only (these keys exist). Do not invent keys.
- **Verification:** Verify every numeric claim against NumPy before committing. The code cell must run top-to-bottom under `MPLBACKEND=Agg python3` and end figures with `plt.show()`. Seed all RNGs. No pytest; the build is the test.
- **Commits:** One commit per task, message prefix `feat(ch09): ...`. **This is a git worktree — commit with explicit path** `git -C /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts ...` (a bare `cd` in a compound command may not take effect, sending the commit to `main` where a hook blocks it). Set identity once with `git -C <wt> config ...`.
- **No inline solution links** in exercises; solutions live only in the appendix, gated by `::: {.content-visible when-meta="show-solutions"}`.

## Verified anchor numbers (use these exact values)

- **Worked (a) — 1-D Gaussian Hamiltonian dynamics.** $\pi(q)\propto e^{-q^2/2}$, $U=\tfrac12 q^2$, $\nabla U=q$, $M=1$. Hamilton's equations $\dot q=p,\ \dot p=-q$ are the **harmonic oscillator** with solution $q(t)=q_0\cos t+p_0\sin t$, $p(t)=p_0\cos t-q_0\sin t$; $H=\tfrac12(q^2+p^2)=\tfrac12(q_0^2+p_0^2)$ is conserved for all $t$. From $(q_0,p_0)=(1,0)$: $t=\tfrac\pi2\Rightarrow(0,-1)$; $t=\pi\Rightarrow(-1,0)$ — one trajectory crosses the whole distribution.
- **Worked (b) — one leapfrog step.** Same target, $\varepsilon=0.3$, start $(q,p)=(1,0)$: $p_{1/2}=0-\tfrac{0.3}{2}(1)=-0.15$; $q'=1+0.3(-0.15)=0.955$; $p'=-0.15-\tfrac{0.3}{2}(0.955)=-0.29325$. $H=0.5$, $H'=\tfrac12(0.955^2+0.29325^2)=0.49901$, $\Delta H=-0.00099$ (third-order small), acceptance $\min(1,e^{-\Delta H})=1$.
- **Code money plot — HMC vs RWM on a correlated Gaussian.** Target $\mathcal N(0,\Sigma)$, $\Sigma=\big[\begin{smallmatrix}1&0.95\\0.95&1\end{smallmatrix}\big]$. In a verified run (seed 0, 4000 iters, HMC $\varepsilon=0.25,L=20$; RWM step 0.4): **HMC ESS ≈ 4000 (essentially independent), RWM ESS ≈ 24** — roughly a 150–170× gain; HMC acceptance ≈ 0.88, RWM ≈ 0.52. Exact ESS is seed-dependent; the code must **compute both and assert HMC ESS is at least, say, 10× RWM ESS**, not hardcode 4000.
- **Neal's funnel** (Rung 6 / code): $v\sim\mathcal N(0,3^2)$, $x_i\mid v\sim\mathcal N(0,e^{v})$. Fixed-$\varepsilon$ HMC produces **divergences** (trajectories with $|\Delta H|$ exceeding a large threshold, e.g. $>1000$, or NaN) clustered in the neck ($v$ very negative); **non-centered reparameterization** $\tilde x_i\sim\mathcal N(0,1),\ x_i=\tilde x_i e^{v/2}$ removes them and recovers the $v$-marginal $\mathcal N(0,9)$.

---

## File Structure

- **Modify (replace body):** `parts/03-sampling/03-hmc-nuts.qmd` — the chapter. Keep the H1 `# Hamiltonian Monte Carlo & NUTS` and the `*Canonical anchors: Neal 2011; Betancourt 2017; Hoffman & Gelman 2014.*` line; replace the stub callout and all empty section bodies with full content.
- **Modify (append):** `appendix/solutions.qmd` — insert a `## Chapter 9 — Hamiltonian Monte Carlo & NUTS` block immediately before the closing `:::` of the `show-solutions` visible div (currently the last line, line 1564). Match the format of the existing `## Chapter 8 — Markov Chain Monte Carlo` block.

---

### Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd:1-11`

- [ ] **Step 1: Clean the stub and write the Motivation**

Keep line 1 `# Hamiltonian Monte Carlo & NUTS` and line 3 `*Canonical anchors: Neal 2011; Betancourt 2017; Hoffman & Gelman 2014.*`. Delete the `::: {.callout-note ...}` stub block. Replace the `## Motivation` placeholder body with a full section; leave the other section headings/placeholders for later tasks.

The Motivation (~3 paragraphs, Chapter-8 density) must:
- Open on Chapter 8's diagnosis: random-walk Metropolis is **uninformed** — its proposal is symmetric in all directions and ignores the shape of $\pi$, so on a correlated posterior it makes mostly-rejected cross-ridge proposals and tiny accepted along-ridge steps, collapsing the effective sample size. The information it ignores is the **gradient** $\nabla\log\pi$.
- Pose the driving question verbatim (emphasized): *"You can compute the gradient of the log-posterior — so why propose blindly? Use it to glide along the distribution instead of stumbling across it."*
- Introduce the physical idea (a frictionless puck on the surface $U=-\log\pi$, given a random momentum kick, conserving energy as it glides) and name the deliverables: Hamiltonian dynamics and why its conservation laws make a long trajectory a valid proposal; the leapfrog integrator; the HMC accept step; **NUTS** as the automation that made HMC the default in Stan and PyMC; and the **funnel** (Chapter 6's hierarchical geometry) as the case that still defeats a naive HMC. Anchor to the correlated MMM posterior $p(\beta\mid y)$ and foreshadow the money plot (HMC's effective sample size dwarfs random walk's on the same target).

- [ ] **Step 2: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" config user.name "jlh530i"; git -C "$WT" config user.email "jlh530i@gmail.com"
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): front matter and Motivation"
```

---

### Task 2: Theory & Proofs — rungs 1–3 (augmentation, Hamiltonian dynamics, leapfrog)

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (the `## Theory & Proofs` section, first three rungs)

Write `## Theory & Proofs` with a 1–2 sentence ladder lead, then three `###` subsections. **Stop before rung 4** (next task). Prose-led, each rung closes on the through-line.

- [ ] **Step 1: Rung 1 — Augmenting with momentum**

Recap the random-walk failure (Chapter 8). Introduce the augmented target: position $q\in\mathbb R^d$, **potential energy** $U(q)=-\log\pi(q)$ (up to a constant); auxiliary **momentum** $p\in\mathbb R^d$, **kinetic energy** $K(p)=\tfrac12 p^\top M^{-1}p$ (mass matrix $M$, default $I$); the **Hamiltonian** $H(q,p)=U(q)+K(p)$. The **joint density**
$$
\pi(q,p)\propto e^{-H(q,p)}=e^{-U(q)}\,e^{-K(p)},
$$
factorizes, so $q\sim\pi$ and $p\sim\mathcal N(0,M)$ are **independent** and the **marginal over $p$ is exactly $\pi(q)$**: sampling $(q,p)$ and discarding $p$ samples the posterior. State the plan: move through $(q,p)$-space along trajectories that keep $H$ (nearly) constant, so a long proposal stays in high-probability regions — the cure for the random walk's blindness.

- [ ] **Step 2: Rung 2 — Hamiltonian dynamics (PROOF of conservation)**

State **Hamilton's equations**:
$$
\dot q=\frac{\partial H}{\partial p}=M^{-1}p,\qquad \dot p=-\frac{\partial H}{\partial q}=-\nabla U(q).
$$
> **Theorem (energy conservation).** Along any solution of Hamilton's equations, $H(q(t),p(t))$ is constant: $\frac{d}{dt}H=0$.

**Proof.** By the chain rule and Hamilton's equations,
$$
\frac{dH}{dt}=\sum_i\left(\frac{\partial H}{\partial q_i}\dot q_i+\frac{\partial H}{\partial p_i}\dot p_i\right)=\sum_i\left(\frac{\partial H}{\partial q_i}\frac{\partial H}{\partial p_i}-\frac{\partial H}{\partial p_i}\frac{\partial H}{\partial q_i}\right)=0.
$$
$\blacksquare$
Then **state** (citing `[@neal2011]`) the other two properties: **volume preservation** (Liouville — the flow's velocity field $(\dot q,\dot p)$ has zero divergence, so phase-space volume is preserved and a dynamics proposal needs **no Jacobian correction**), and **time-reversibility** (negating $p$ and integrating backward retraces the path). Conclude: exact Hamiltonian flow followed by momentum negation is a volume-preserving, reversible involution that conserves $H$ — a Metropolis proposal accepted with probability 1. The obstacle: the equations have no closed-form solution for a general $U$, forcing numerical integration (Rung 3).

- [ ] **Step 3: Rung 3 — The leapfrog integrator (PROOF of volume preservation + reversibility)**

Define the **leapfrog (Störmer–Verlet)** step with step size $\varepsilon$:
$$
p_{1/2}=p-\tfrac{\varepsilon}{2}\nabla U(q),\qquad q'=q+\varepsilon M^{-1}p_{1/2},\qquad p'=p_{1/2}-\tfrac{\varepsilon}{2}\nabla U(q').
$$
> **Theorem (leapfrog is volume-preserving and reversible).** Each leapfrog step has Jacobian determinant 1, and the map becomes an involution when composed with momentum negation.

**Proof (volume preservation).** Each of the three sub-steps is a **shear**: the momentum half-steps change $p$ using only $q$ (so $\partial p'/\partial p=I$, and the Jacobian is block-triangular $\big[\begin{smallmatrix}I&0\\ *&I\end{smallmatrix}\big]$ with unit determinant), and the position step changes $q$ using only $p_{1/2}$ (Jacobian $\big[\begin{smallmatrix}I& *\\0&I\end{smallmatrix}\big]$, unit determinant). The full step is their composition, so its Jacobian determinant is $1\cdot1\cdot1=1$. **Reversibility:** negating $p$, applying the leapfrog step, and negating $p$ again returns $(q,p)$ — the half-step structure is symmetric in time. $\blacksquare$
Then discuss **energy error**: leapfrog does **not** conserve $H$ exactly, but because it is **symplectic** the error does not drift — the global error in $H$ is $O(\varepsilon^2)$ and stays bounded over long trajectories (state, cite `[@neal2011]`). Contrast **Euler**, which is neither volume-preserving nor stable (energy spirals out). This small, bounded $\Delta H$ is precisely what the Metropolis correction of Rung 4 cleans up. Reuse Worked Example (b)'s numbers informally if helpful (a single step from $(1,0)$ at $\varepsilon=0.3$ gives $\Delta H\approx-0.001$).

- [ ] **Step 4: Verify the leapfrog numbers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
eps=0.3; q,p=1.0,0.0
ph=p-0.5*eps*q; qn=q+eps*ph; pn=ph-0.5*eps*qn
print('p_half',ph,'q*',round(qn,4),'p*',round(pn,5),'dH',round(0.5*(qn*qn+pn*pn)-0.5,5))"
```
Expected: `p_half -0.15 q* 0.955 p* -0.29325 dH -0.00099`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Theory rungs 1-3 — momentum augmentation, Hamiltonian dynamics, leapfrog"
```

---

### Task 3: Theory & Proofs — rungs 4–6 (HMC keystone, NUTS, divergences/funnel)

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (append three subsections to `## Theory & Proofs`, before `## Worked Examples`)

- [ ] **Step 1: Rung 4 — The HMC algorithm (KEYSTONE proof)**

State the algorithm. One HMC step from $q$:
1. **Resample momentum** $p\sim\mathcal N(0,M)$.
2. Run $L$ **leapfrog** steps from $(q,p)$ to a proposal $(q^\*,p^\*)$.
3. **Accept** $q^\*$ with probability $\min\!\big(1,\,e^{-(H(q^\*,p^\*)-H(q,p))}\big)=\min(1,e^{-\Delta H})$; otherwise keep $q$.
> **Theorem (HMC leaves $\pi$ invariant).** The HMC transition has $\pi(q)$ as a stationary distribution.

**Proof.** Two sub-steps. (i) **Momentum resampling** draws $p$ from its exact conditional $\mathcal N(0,M)$ under the joint $\pi(q,p)\propto e^{-H}$ (since $p\perp q$); this is a Gibbs step and leaves the joint invariant. (ii) **Leapfrog proposal + accept.** Let $T$ be "run $L$ leapfrog steps, then negate momentum." By Rung 3, $T$ is **volume-preserving** ($|\det \partial T|=1$) and an **involution** ($T\circ T=\mathrm{id}$, by reversibility). For a deterministic, volume-preserving, involutive proposal the Metropolis–Hastings acceptance (Chapter 8) reduces to
$$
\alpha=\min\!\left(1,\ \frac{\pi(T(q,p))}{\pi(q,p)}\,\big|\det \partial T\big|\right)=\min\!\left(1,\ \frac{e^{-H(q^\*,p^\*)}}{e^{-H(q,p)}}\right)=\min(1,e^{-\Delta H}),
$$
with **no Jacobian factor and no proposal-density ratio** — and this kernel satisfies detailed balance with respect to the joint $e^{-H}$ (Chapter 8's MH theorem). Negating $p^\*$ at the end is irrelevant to $q$ and is undone by the next resample. Both sub-steps preserve $e^{-H}$, so its $q$-marginal $\pi(q)$ is invariant. $\blacksquare$
Land the payoff: because leapfrog **nearly conserves $H$**, $\Delta H\approx 0$ and the acceptance probability is **near 1 even for a proposal $L\varepsilon$ away** — the decisive contrast with random walk, whose acceptance forces $O(1/\sqrt d)$ steps. State the two tuning knobs: $\varepsilon$ (set by a target acceptance, optimal ~0.65–0.8, cite `[@neal2011]`) and $L$ (trajectory length). Note that the gradient $\nabla U=-\nabla\log\pi$ is what the sampler needs — available for any differentiable log-posterior, e.g. via automatic differentiation in Stan/PyMC.

- [ ] **Step 2: Rung 5 — NUTS: removing the trajectory-length knob (conceptual but precise)**

Motivate: choosing $L$ is delicate — too small and HMC degenerates toward a random walk; too large and the trajectory makes a **U-turn**, looping back and wasting computation (and can land *back near the start*). The **No-U-Turn Sampler** builds the trajectory **adaptively**, **doubling** its length forward or backward in time (chosen randomly) and stopping when the trajectory starts to retrace — the **U-turn criterion**: with $q^-,p^-$ and $q^+,p^+$ the two ends, stop when
$$
(q^+-q^-)\cdot p^-<0\quad\text{or}\quad(q^+-q^-)\cdot p^+<0,
$$
i.e. continuing would reduce the distance between the ends. State that NUTS then **samples a state from the trajectory** (via a slice or multinomial scheme over the visited points) in a way that **preserves reversibility and detailed balance**, and that the recursive **tree-doubling** keeps the cost proportional to the trajectory length. **Cite `[@hoffman2014]`** for the full recursive algorithm and its correctness proof — state it is conceptual here. Describe **step-size adaptation by dual averaging** during warm-up (tune $\varepsilon$ to hit a target acceptance statistic, then freeze it). Close: **NUTS with dual-averaged $\varepsilon$ is what Stan and PyMC run by default**, which is why modern MMM is fit with it — the analyst supplies the model, the sampler tunes itself.

- [ ] **Step 3: Rung 6 — Divergences, the funnel, and the bridge to Chapter 10**

HMC's distinctive failure mode is the **divergent transition**: where the posterior has **high curvature**, a fixed step size $\varepsilon$ makes leapfrog **unstable**, $\Delta H$ blows up (or goes to infinity / NaN), the proposal is rejected, and the sampler **cannot enter** that region — biasing the result by silently missing mass. The canonical hard geometry is **Neal's funnel**:
$$
v\sim\mathcal N(0,3^2),\qquad x_i\mid v\sim\mathcal N(0,e^{v}),
$$
exactly the shape of a **hierarchical variance parameter** $v$ (a log-scale) and its coefficients $x_i$ — Chapter 6's model: when $v$ is very negative the $x_i$ are squeezed into a sharp neck whose curvature no single $\varepsilon$ can handle. The fix is **non-centered reparameterization**: sample $\tilde x_i\sim\mathcal N(0,1)$ and set $x_i=\tilde x_i\,e^{v/2}$, so the sampler explores the **decoupled, flat** $(\tilde x,v)$ geometry and the funnel disappears. Stress that **divergences are a trustworthy alarm**, not a nuisance — they flag exactly the regions the sampler is failing to explore. Hand off to **Chapter 10**, which develops the formal convergence diagnostics ($\hat R$, effective sample size, divergence counts) that detect these pathologies in practice. Cite `[@betancourt2017]`.

- [ ] **Step 4: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Theory rungs 4-6 — HMC algorithm, NUTS, divergences and the funnel"
```

---

### Task 4: Worked Examples

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (the `## Worked Examples` section)

Match Chapter 8's worked-example style: bold step mini-headings, every number shown, closing interpretation.

- [ ] **Step 1: Example (a) — Hamiltonian dynamics for a 1-D Gaussian, by hand**

Target $\pi(q)\propto e^{-q^2/2}$, so $U(q)=\tfrac12 q^2$, $\nabla U=q$, $M=1$. Hamilton's equations $\dot q=p$, $\dot p=-q$ are the **harmonic oscillator**; verify the solution $q(t)=q_0\cos t+p_0\sin t$, $p(t)=p_0\cos t-q_0\sin t$ satisfies them (differentiate). Show $H=\tfrac12(q^2+p^2)$ is conserved: substitute to get $\tfrac12(q_0^2+p_0^2)$ independent of $t$. From $(q_0,p_0)=(1,0)$: at $t=\tfrac\pi2$ the state is $(0,-1)$, at $t=\pi$ it is $(-1,0)$ — a single trajectory carries the particle clear across the distribution, the geometric reason HMC decorrelates so fast, versus a random-walk step of size $O(\varepsilon)$. Interpret: the momentum converts potential to kinetic energy and back, sweeping out a circle in phase space at constant $H$.

- [ ] **Step 2: Example (b) — One leapfrog step, by hand**

Same target, $\varepsilon=0.3$, start $(q,p)=(1,0)$. Half-step momentum $p_{1/2}=0-\tfrac{0.3}{2}(1)=-0.15$; full-step position $q'=1+0.3(-0.15)=0.955$; half-step momentum $p'=-0.15-\tfrac{0.3}{2}(0.955)=-0.29325$. Energies $H=\tfrac12(1^2+0^2)=0.5$ and $H'=\tfrac12(0.955^2+0.29325^2)=0.49901$, so $\Delta H=-0.00099$ — third-order small in $\varepsilon$, the controlled error leapfrog leaves for the Metropolis step. The acceptance probability $\min(1,e^{-\Delta H})=\min(1,e^{0.00099})=1$ (the step slightly *lowered* the energy, so it is accepted with certainty). Interpret: even a coarse $\varepsilon$ tracks the true orbit closely, which is why HMC can take long, cheap, high-acceptance trajectories.

- [ ] **Step 3: Verify, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
# (a) conservation on exact orbit
for t in [np.pi/2, np.pi]:
    q=np.cos(t); p=-np.sin(t); print('(a) t',round(t,3),'state',round(q,4),round(p,4),'H',round(0.5*(q*q+p*p),4))
# (b) leapfrog
eps=0.3; q,p=1.,0.; ph=p-0.5*eps*q; qn=q+eps*ph; pn=ph-0.5*eps*qn
print('(b)',round(ph,4),round(qn,4),round(pn,5),'dH',round(0.5*(qn*qn+pn*pn)-0.5,5))"
```
Expected: `(a) t 1.571 state 0.0 -1.0 H 0.5`, `(a) t 3.142 state -1.0 -0.0 H 0.5`; `(b) -0.15 0.955 -0.29325 dH -0.00099`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Worked Examples — Hamiltonian dynamics and a leapfrog step"
```

---

### Task 5: Code Tie-in

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (the `## Code Tie-in` section)

- [ ] **Step 1: Write a single self-contained ```` ```{python} ```` cell**

Model the structure on Chapter 8's Code Tie-in. Intro paragraph, one ```` ```{python} ```` block (numpy + matplotlib, seeded `rng = np.random.default_rng(0)`), closing paragraph quoting actual printed numbers. The cell, with small helpers, must:
1. **HMC vs random-walk Metropolis on a strongly correlated 2-D Gaussian** ($\Sigma=[[1,0.95],[0.95,1]]$). Implement `U(q)=0.5*q@Siginv@q`, `gradU(q)=Siginv@q`, a `leapfrog(q,p,eps,L)` step, and a full `hmc(n,eps,L)` step (resample $p\sim\mathcal N(0,I)$, $L$ leapfrog steps, accept on `log(u) < -(H_new-H_old)`). Also implement `rwm(n,step)`. Run both for the same `n` (e.g. 4000) on the **same target**; compute the **ESS** of the first coordinate for each (helper `ess(chain)` = `n/(1+2*sum(rho_k))`, summing autocorrelations until the first negative). **Print HMC ESS, RWM ESS, and both acceptance rates, and `assert ess_hmc > 10*ess_rwm`** (the verified anchor is ~4000 vs ~24; HMC $\varepsilon=0.25,L=20$, RWM step ~0.4 — tune so HMC acc ~0.8–0.9 and the assert holds). The money plot: scatter/trace of both chains showing HMC covering the ridge while RWM crawls.
2. **Energy conservation / acceptance vs step size.** For the same target, plot the per-iteration $\Delta H$ for HMC (small, centered near 0) and sweep $\varepsilon\in\{0.1,0.25,0.6,1.2\}$ printing acceptance — show it degrading as $\varepsilon$ grows.
3. **Neal's funnel and divergences.** Sample the funnel $v\sim\mathcal N(0,3^2)$, $x\mid v\sim\mathcal N(0,e^{v})$ (1 or a few $x$). Run a fixed-$\varepsilon$ HMC on the **centered** parameterization (potential from the joint $-\log p(v,x)$); count **divergences** (iterations where $|\Delta H|$ exceeds a large threshold, e.g. 1000, or is non-finite) and record where they occur in $v$ (the neck, $v$ negative). Then run the **non-centered** reparameterization $\tilde x\sim\mathcal N(0,1)$, $x=\tilde x e^{v/2}$ and show the divergence count drops to ~0 and the recovered $v$-marginal mean/variance match $\mathcal N(0,9)$ (print). Plot centered vs non-centered $(v,x)$ scatters.

End each figure with `plt.show()`. The closing paragraph quotes: HMC vs RWM ESS and acceptance, the $\Delta H$/acceptance-vs-$\varepsilon$ behavior, and the funnel divergence counts before/after reparameterization.

- [ ] **Step 2: Extract the cell and run headless**

Save the cell body to `/tmp/ch9_code.py` and run:
```bash
MPLBACKEND=Agg python3 /tmp/ch9_code.py
```
Expected: runs top-to-bottom, the `ess_hmc > 10*ess_rwm` assert passes, divergences drop after reparameterization. Make the prose match the actual printed numbers.

- [ ] **Step 3: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Code Tie-in — HMC vs RWM ESS, energy conservation, the funnel"
```

---

### Task 6: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (the `## Exercises` section with the four `###` tier headings)

Match Chapter 8's exercise style and heading text exactly: `### C -- Conceptual / Reading Comprehension`, `### B -- By Hand`, `### P -- Prove It`, `### A -- Applied / Code`. Each problem fully self-contained. **No links to solutions.**

- [ ] **Step 1: Tier C (3 problems)**
- **C1.** Why augment the parameter with a momentum, and how does marginalizing $p$ from the joint $e^{-H}$ recover the posterior $\pi(q)$? Explain the roles of the potential $U=-\log\pi$ and kinetic $K=\tfrac12 p^\top M^{-1}p$ energies, and why a trajectory that conserves $H$ stays in high-probability regions.
- **C2.** State the three properties of Hamiltonian flow — energy conservation, volume preservation, reversibility — and explain what each one buys for the validity of the proposal (no Jacobian, detailed balance, acceptance near 1). Then explain why **leapfrog**, not Euler, is the right discretization, and what role the Metropolis $\min(1,e^{-\Delta H})$ step plays.
- **C3.** Explain what **NUTS** automates and what the **U-turn criterion** detects, and why this matters for a practitioner. Then explain HMC's **divergence** failure mode: what a divergent transition is, why **Neal's funnel** (a hierarchical variance and its coefficients) provokes it, how **non-centered reparameterization** fixes it, and why divergences are a trustworthy diagnostic.

- [ ] **Step 2: Tier B (3 problems)**
- **B1.** For $\pi(q)\propto e^{-q^2/2}$ ($U=\tfrac12 q^2$, $M=1$): write Hamilton's equations, verify that $q(t)=q_0\cos t+p_0\sin t$, $p(t)=p_0\cos t-q_0\sin t$ solves them, and show $H=\tfrac12(q^2+p^2)$ is conserved. Starting from $(q_0,p_0)=(2,0)$, give the state at $t=\pi/2$ and $t=\pi$. (Clean: $(0,-2)$ and $(-2,0)$; $H=2$ throughout.)
- **B2.** One leapfrog step for the same target with $\varepsilon=0.5$ from $(q,p)=(0,1)$: compute $p_{1/2}$, $q'$, $p'$, and $\Delta H$. (Clean: $\nabla U(0)=0$ so $p_{1/2}=1$; $q'=0.5$; $p'=1-0.25\cdot0.5=0.875$; $H=0.5$, $H'=\tfrac12(0.25+0.765625)=0.5078$, $\Delta H\approx0.0078$.)
- **B3.** Given a proposed HMC move with $\Delta H=0.4$, compute the acceptance probability $\min(1,e^{-\Delta H})$. Repeat for $\Delta H=-0.2$ (accepted with certainty) and $\Delta H=5$ (a near-divergence). (Clean: $e^{-0.4}\approx0.670$; $1$; $e^{-5}\approx0.0067$.)

- [ ] **Step 3: Tier P (3 problems)**
- **P1.** Prove **energy conservation**: along solutions of Hamilton's equations $\dot q=\partial H/\partial p$, $\dot p=-\partial H/\partial q$, show $\frac{d}{dt}H=0$ (chain rule; the cross terms cancel).
- **P2.** Prove the **leapfrog step is volume-preserving** (Jacobian determinant 1) by writing each of the three sub-steps as a shear with unit-determinant triangular Jacobian and composing, and prove it is **reversible** (momentum negation conjugates the step to its inverse).
- **P3.** Prove that the **HMC accept step leaves $e^{-H}$ invariant**: argue that "leapfrog then negate momentum" is a volume-preserving involution, so the Metropolis–Hastings acceptance reduces to $\min(1,e^{-\Delta H})$ with no Jacobian, and invoke Chapter 8's detailed-balance theorem; conclude the $q$-marginal $\pi(q)$ is invariant.

- [ ] **Step 4: Tier A (2 problems)**
- **A1.** Implement HMC (with a leapfrog integrator) and random-walk Metropolis for a strongly correlated Gaussian; compare effective sample size at matched iteration counts and show HMC's advantage grows with the correlation; study how acceptance and mean $|\Delta H|$ depend on the step size $\varepsilon$.
- **A2.** Reproduce **Neal's funnel** $v\sim\mathcal N(0,3^2)$, $x\mid v\sim\mathcal N(0,e^v)$; run fixed-$\varepsilon$ HMC on the centered parameterization and count divergences; then apply the **non-centered** reparameterization and show the divergences vanish and the $v$-marginal $\mathcal N(0,9)$ is recovered.

- [ ] **Step 5: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Exercises — C/B/P/A tiers"
```

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/03-sampling/03-hmc-nuts.qmd` (append `## Summary` as the final section)

- [ ] **Step 1: Write `## Summary`**

Match Chapter 8's Summary format exactly: one-line lead, then **Key concepts** (bold lead-ins) and **Key identities** (**inline math only**). Cover:
- **Key concepts:** augment with momentum so the joint is $e^{-H}$ and its $q$-marginal is $\pi$; Hamiltonian dynamics conserves energy, preserves volume, and is reversible — making a long trajectory a valid proposal; leapfrog is the symplectic discretization that keeps volume preservation and reversibility with $O(\varepsilon^2)$ bounded energy error; HMC accepts with $\min(1,e^{-\Delta H})$, near 1 because energy is nearly conserved, so it glides across correlated posteriors where random walk crawls; NUTS automates the trajectory length via the U-turn criterion and is the Stan/PyMC default; divergences flag high-curvature regions (the hierarchical funnel), cured by non-centered reparameterization.
- **Key identities (inline):**
  - Augmented target: $H(q,p)=U(q)+K(p)$, $U=-\log\pi$, $K=\tfrac12 p^\top M^{-1}p$, joint $\propto e^{-H}$.
  - Hamilton's equations: $\dot q=M^{-1}p$, $\dot p=-\nabla U(q)$; conservation $\frac{d}{dt}H=0$.
  - Leapfrog: $p_{1/2}=p-\tfrac{\varepsilon}{2}\nabla U(q)$, $q'=q+\varepsilon M^{-1}p_{1/2}$, $p'=p_{1/2}-\tfrac{\varepsilon}{2}\nabla U(q')$; Jacobian determinant 1.
  - HMC acceptance: $\min(1,e^{-\Delta H})$, $\Delta H=H(q^\*,p^\*)-H(q,p)$.
  - U-turn criterion: stop when $(q^+-q^-)\cdot p^-<0$ or $(q^+-q^-)\cdot p^+<0$.
  - Non-centered funnel: $\tilde x\sim\mathcal N(0,1)$, $x=\tilde x\,e^{v/2}$.

- [ ] **Step 2: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add parts/03-sampling/03-hmc-nuts.qmd
git -C "$WT" commit -m "feat(ch09): Summary — key concepts and identities"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (insert before the closing `:::` at line 1564)

- [ ] **Step 1: Append the Chapter 9 solutions block**

Immediately before the final `:::` that closes the `show-solutions` visible div, add a `## Chapter 9 — Hamiltonian Monte Carlo & NUTS` heading followed by full solutions to every exercise from Task 6 (C1–C3, B1–B3, P1–P3, A1–A2). Match the format of the existing `## Chapter 8 — Markov Chain Monte Carlo` block: em-dash sub-headings `### C — …` etc., bold labels (`**C1.**`), full worked math for B and P, plain ```` ```python ```` reference sketches (NOT `{python}`) plus a qualitative takeaway for A. Each proof ends `$\blacksquare$`. Verify every numeric answer against NumPy (B1: states $(0,-2),(-2,0)$, $H=2$; B2: $p_{1/2}=1,q'=0.5,p'=0.875,\Delta H\approx0.0078$; B3: $0.670,\,1,\,0.0067$).

- [ ] **Step 2: Verify by-hand answers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
# B2
eps=0.5; q,p=0.,1.; ph=p-0.5*eps*q; qn=q+eps*ph; pn=ph-0.5*eps*qn
print('B2', ph, qn, pn, 'dH', round(0.5*(qn*qn+pn*pn)-0.5*(0+1),4))
# B3
print('B3', round(np.exp(-0.4),4), min(1,np.exp(0.2)), round(np.exp(-5),4))"
```
Expected: `B2 1.0 0.5 0.875 dH 0.0078`; `B3 0.6703 1 0.0067`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
git -C "$WT" add appendix/solutions.qmd
git -C "$WT" commit -m "feat(ch09): appendix solutions for HMC & NUTS exercises"
```

---

### Task 9: Final whole-chapter review

**Files:**
- Review: `parts/03-sampling/03-hmc-nuts.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: KaTeX / structure lint**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch9-hmc-nuts
grep -n 'begin{align}' "$WT"/parts/03-sampling/03-hmc-nuts.qmd || echo "no bare align — good"
grep -c '^\$\$' "$WT"/parts/03-sampling/03-hmc-nuts.qmd   # expect even
grep -nE '^\## ' "$WT"/parts/03-sampling/03-hmc-nuts.qmd  # expect the 6 template headings in order
grep -nE '^# ' "$WT"/parts/03-sampling/03-hmc-nuts.qmd    # expect '# Hamiltonian Monte Carlo & NUTS'
grep -c '^:::$' "$WT"/appendix/solutions.qmd   # expect 2 (file still closes cleanly)
```
Confirm no bare `\begin{align}`, even `$$` count, the six template headings in order, and citations only the valid keys.

- [ ] **Step 2: Re-run the code cell headless**

```bash
MPLBACKEND=Agg python3 /tmp/ch9_code.py && echo "CODE CELL OK"
```
Expected: `CODE CELL OK`, `ess_hmc > 10*ess_rwm` passes, funnel divergences drop after reparameterization.

- [ ] **Step 3: Through-line check (manual)**

Confirm every theory subsection ties back to the glide-don't-stumble / efficiency / funnel through-line, the three named proofs (energy conservation; leapfrog volume preservation + reversibility; HMC invariance) are present and correct, NUTS is conceptual-but-precise with the U-turn criterion, and Chapter 8 (MH detailed balance, ESS) and Chapter 6 (the funnel) are reused by reference. Fix gaps inline and commit if changed.

---

## Self-Review (completed by plan author)

- **Spec coverage:** Motivation (T1); rungs 1–3 incl. conservation + leapfrog proofs (T2); rungs 4–6 incl. the HMC keystone proof, NUTS, divergences/funnel (T3); both worked examples (T4); code tie-in with the 3 required behaviors incl. the funnel cure (T5); C/B/P/A exercises (T6); auto-included Summary (T7); appendix solutions (T8); render/lint gate (T9). All spec success criteria mapped.
- **Placeholder scan:** none — every task carries the actual formulas and NumPy-verified anchor numbers.
- **Consistency:** $H=U+K$, $U=-\log\pi$, Hamilton's equations, the leapfrog triple, the acceptance $\min(1,e^{-\Delta H})$, and the U-turn criterion appear identically across T2/T3/T4/T6/T7. Worked numbers are NumPy-verified; the HMC-vs-RWM ESS gap is computed-and-asserted in code (not hardcoded), since exact ESS is seed-dependent.
- **Note:** the real render gate is CI `quarto render` (HTML + PDF) on the PR — Quarto is not installed locally, so local verification is the Python cell + lint checks above.
