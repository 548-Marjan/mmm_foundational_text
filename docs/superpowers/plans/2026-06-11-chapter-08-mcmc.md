# Chapter 8 — Markov Chain Monte Carlo Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the scaffolded stub `parts/03-sampling/02-mcmc.qmd` with the second chapter of Part III — the Metropolis–Hastings and Gibbs algorithms that build the chains Chapter 7 proved sufficient — and add matching appendix solutions.

**Architecture:** One Quarto `.qmd` chapter authored section-by-section against the fixed template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), plus a `## Chapter 8 — Markov Chain Monte Carlo` block appended to the shared `appendix/solutions.qmd` (gated by `show-solutions`). The chapter reuses Chapter 7 (detailed balance ⇒ stationarity, the ergodic theorem) and Chapter 6 (the hierarchical Gibbs sampler) directly. Two anchors: the correlated bivariate Gaussian (clean, analytic, motivates HMC) and a non-conjugate 1-D posterior with no closed form (the payoff). The single runnable Python cell is verified headless before commit; the real render gate is CI (`quarto render`), Quarto not being installed locally.

**Tech Stack:** Quarto book, KaTeX (HTML) / LaTeX (PDF) math, Python code cell (numpy, scipy, matplotlib pinned in `requirements.txt`). Canonical references via `references.bib` keys: `@robert2004`, `@gelman2013`, `@hoff2009`, `@gelmanhill2006`.

---

## Conventions (apply to every task)

- **Math/KaTeX:** Use `aligned` inside `$$ ... $$`, never bare `\begin{align}`. Keep `$$` delimiters on their own lines so display-math counts stay even. Inline math with single `$`. Reuse Chapter 7 notation: target $\pi$, transition kernel $P(x,y)$, detailed balance $\pi(x)P(x,y)=\pi(y)P(y,x)$, the ergodic theorem; and Chapter 5/6 posterior notation $p(\beta\mid y)\propto p(y\mid\beta)p(\beta)$.
- **Voice:** Match Chapters 1–7 (read `parts/03-sampling/01-markov-chains.qmd`, the immediate predecessor). Rigorous, prose-led, each theory subsection ("rung") ending by tying back to the through-line / an MMM payoff. End every full proof with `$\blacksquare$`; use `> **Theorem (name).** ...` blockquotes for named results.
- **Citations:** `[@robert2004]`, `[@gelman2013]`, `[@hoff2009]`, `[@gelmanhill2006]` only (these keys exist in `references.bib`). There is NO dedicated Metropolis/Hastings/Geman entry — attribute the algorithms by name in prose and cite `@robert2004`/`@gelman2013`. Do not invent keys.
- **Verification:** Verify every numeric claim against NumPy before committing. The code cell must run top-to-bottom under `MPLBACKEND=Agg python3` and end figures with `plt.show()`. Seed all RNGs. No pytest exists; the build is the test.
- **Commits:** One commit per task, message prefix `feat(ch08): ...`. **This is a git worktree — commit with an explicit path** `git -C /Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc ...` to avoid a cwd pitfall (a bare `cd` in a compound command may not take effect, sending the commit to `main` where a hook blocks it). Set identity once: `git -C <wt> config user.name "jlh530i" && git -C <wt> config user.email "jlh530i@gmail.com"`.
- **No inline solution links** in exercises; solutions live only in the appendix, gated by `::: {.content-visible when-meta="show-solutions"}`.

## Verified anchor numbers (use these exact values)

- **Worked (a) Metropolis acceptance by hand.** Unnormalized target $\pi(x)\propto e^{-x^2/2}$ (standard normal with the $1/\sqrt{2\pi}$ dropped), symmetric random-walk proposal so $\alpha=\min(1,e^{-(y^2-x^2)/2})$:
  - $x=0\to y=1$: $\alpha=e^{-1/2}\approx0.6065$.
  - $x=1\to y=0.5$: ratio $e^{-(0.25-1)/2}=e^{0.375}>1$, so $\alpha=1$ (downhill in $x^2$ ⇒ uphill in density ⇒ always accept).
  - $x=0.5\to y=2$: $\alpha=e^{-(4-0.25)/2}=e^{-1.875}\approx0.1534$.
- **Worked (b) Gibbs for the correlated bivariate normal.** $\mathcal N(0,\Sigma)$, $\Sigma=\big[\begin{smallmatrix}1&\rho\\\rho&1\end{smallmatrix}\big]$, $\rho=0.8$. Full conditionals $x\mid y\sim\mathcal N(\rho y,\,1-\rho^2)=\mathcal N(0.8y,\,0.36)$ and symmetrically $y\mid x\sim\mathcal N(0.8x,\,0.36)$. **Classic mixing fact:** the lag-1 autocorrelation of either coordinate equals $\rho^2=0.64$; as $\rho\to1$ it $\to1$ (slow mixing). NumPy-verified: simulated recovery $\mathrm{corr}\approx0.80$, empirical lag-1 autocorrelation $\approx0.64$.
- **Code non-conjugate payoff.** A 1-D logistic (binary-conversion) posterior: $y_i\sim\text{Bernoulli}(\sigma(\beta x_i))$ with $\sigma(z)=1/(1+e^{-z})$, prior $\beta\sim\mathcal N(0,\tau^2)$, $\tau^2=4$. The log-posterior $\ell(\beta)=\sum_i[y_i\beta x_i-\log(1+e^{\beta x_i})]-\beta^2/(2\tau^2)$ has **no closed form** and an intractable normalizer. RW-Metropolis on $\ell$ matches a fine-grid normalized reference (in one verified run: grid mean $1.574$, MH mean $1.577$, acceptance $\approx0.68$). The exact numbers are seed/data-dependent — the code must **compute the grid reference and assert** the MH mean/sd match it within tolerance, not hardcode $1.574$.

---

## File Structure

- **Modify (replace body):** `parts/03-sampling/02-mcmc.qmd` — the chapter. Keep the H1 `# Markov Chain Monte Carlo` and the `*Canonical anchors: Robert & Casella; BDA3.*` line; replace the stub callout and all empty section bodies with full content.
- **Modify (append):** `appendix/solutions.qmd` — insert a `## Chapter 8 — Markov Chain Monte Carlo` block immediately before the closing `:::` of the `show-solutions` visible div (currently the last line, line 1272). Match the format of the existing `## Chapter 7 — Markov Chains` block.

---

### Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd:1-11`

- [ ] **Step 1: Clean the stub and write the Motivation**

Keep line 1 `# Markov Chain Monte Carlo` and line 3 `*Canonical anchors: Robert & Casella; BDA3.*`. Delete the `::: {.callout-note ...}` stub block. Replace the `## Motivation` placeholder body with a full section; leave the other section headings/placeholders for later tasks.

The Motivation (~3 paragraphs, Chapter-7 density) must:
- Open on the gap Chapter 7 left: the ergodic theorem guarantees that *if* we can build an irreducible, aperiodic chain reversible with respect to $\pi$ then its trajectory average estimates $\mathbb E_\pi[f]$ — but Chapter 7 never constructed such a chain. Meanwhile Chapter 5 hit a wall: a non-negativity prior or an adstock/saturation nonlinearity destroys the closed-form posterior, and even its normalizing constant $\int p(y\mid\beta)p(\beta)\,d\beta$ is an intractable integral.
- Pose the driving question verbatim (emphasized): *"The posterior is a formula you can evaluate at any $\beta$ but cannot integrate, normalize, or sample from directly — so how do you draw from it?"*
- State the resolution: manufacture a Markov chain whose stationary distribution is that posterior, using only **pointwise evaluations of the unnormalized density** — because (foreshadowing Rung 2) the accept/reject decision depends on $\pi$ only through ratios $\pi(y)/\pi(x)$, in which the unknown normalizer cancels. Name the two algorithms (Metropolis–Hastings, the general recipe; Gibbs, the coordinate-wise sampler that *is* Chapter 6's hierarchical sampler), the efficiency questions (acceptance rate, autocorrelation, effective sample size), and that random-walk Metropolis's slow mixing in high/correlated dimensions motivates Chapter 9's Hamiltonian Monte Carlo. Anchor to $p(\beta\mid y)$.

- [ ] **Step 2: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" config user.name "jlh530i"; git -C "$WT" config user.email "jlh530i@gmail.com"
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): front matter and Motivation"
```

---

### Task 2: Theory & Proofs — rungs 1–2 (Monte Carlo principle + Metropolis–Hastings keystone)

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (the `## Theory & Proofs` section, first two rungs)

Write `## Theory & Proofs` with a 1–2 sentence ladder lead, then two `###` subsections. **Stop before rung 3** (next task). Prose-led, each rung closes on the through-line.

- [ ] **Step 1: Rung 1 — The Monte Carlo principle and why MCMC**

Monte Carlo estimation: for i.i.d. $X_i\sim\pi$, $\hat\mu_n=\frac1n\sum_{i=1}^n f(X_i)$ estimates $\mathbb E_\pi[f]$, with (by the CLT) standard error $\sigma_f/\sqrt n$ where $\sigma_f^2=\operatorname{Var}_\pi(f)$ — a rate **independent of the dimension** of the state space, the structural reason sampling beats deterministic quadrature for high-dimensional posteriors. The obstacle: we cannot draw i.i.d. from a complex posterior (no inverse-CDF, no rejection envelope in high dimension). The resolution: generate **dependent** draws $X_1,X_2,\dots$ from a Markov chain whose stationary distribution is $\pi$; **Chapter 7's ergodic theorem** certifies $\frac1n\sum_t f(X_t)\to\mathbb E_\pi[f]$ regardless of the dependence. So the entire problem reduces to **constructing a chain with stationary distribution $\pi$** — the task of rungs 2–3.

- [ ] **Step 2: Rung 2 — Metropolis–Hastings (KEYSTONE proof)**

Present the algorithm: from the current state $x$, draw a **proposal** $y\sim q(\cdot\mid x)$ and **accept** it ($X_{t+1}=y$) with probability
$$
\alpha(x,y)=\min\!\left(1,\ \frac{\pi(y)\,q(x\mid y)}{\pi(x)\,q(y\mid x)}\right),
$$
otherwise **reject** and stay ($X_{t+1}=x$). For $x\ne y$ the transition density is $P(x,y)=q(y\mid x)\,\alpha(x,y)$ (plus a rejection mass at $x$).
> **Theorem (Metropolis–Hastings is reversible).** The Metropolis–Hastings kernel satisfies detailed balance with respect to $\pi$: $\pi(x)P(x,y)=\pi(y)P(y,x)$ for all $x\ne y$. Hence (Chapter 7, rung 3) $\pi$ is stationary for the chain.

**Proof.** For $x\ne y$, multiply through and use the **min identity**:
$$
\pi(x)P(x,y)=\pi(x)q(y\mid x)\,\alpha(x,y)=\min\!\big(\pi(x)q(y\mid x),\ \pi(y)q(x\mid y)\big),
$$
because $a\cdot\min(1,b/a)=\min(a,b)$ for $a,b>0$. The right-hand side is **symmetric** under swapping $x\leftrightarrow y$; carrying out the same computation for $\pi(y)P(y,x)=\pi(y)q(x\mid y)\alpha(y,x)=\min(\pi(y)q(x\mid y),\pi(x)q(y\mid x))$ gives the identical expression. Therefore $\pi(x)P(x,y)=\pi(y)P(y,x)$. The rejection (diagonal) part is trivially in balance. By Chapter 7's detailed-balance theorem, $\pi$ is stationary. $\blacksquare$

Then emphasize the decisive point: $\pi$ enters $\alpha$ **only through the ratio $\pi(y)/\pi(x)$**, so any unknown normalizing constant cancels — Metropolis–Hastings runs on an **unnormalized** posterior $p(y\mid\beta)p(\beta)$, exactly the situation Chapter 5's wall created. Derive the special cases: a **symmetric proposal** $q(y\mid x)=q(x\mid y)$ collapses the ratio to the **Metropolis** rule
$$
\alpha(x,y)=\min\!\left(1,\ \frac{\pi(y)}{\pi(x)}\right);
$$
the **random-walk Metropolis** ($y=x+\varepsilon$, $\varepsilon$ symmetric) and the **independence sampler** ($q(y\mid x)=q(y)$) are named instances. State (cite `[@robert2004]`) that for any proposal with positive density on the support the chain is irreducible and aperiodic, so Chapter 7's convergence and ergodic theorems apply. Close on the through-line: we have, at last, the construction Chapter 7 promised — a chain that targets any posterior we can evaluate pointwise.

- [ ] **Step 3: Verify the worked-(a) acceptance numbers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
f=lambda x: np.exp(-x*x/2)
for x,y in [(0,1),(1,0.5),(0.5,2)]:
    print(x,'->',y,'alpha=',round(min(1,f(y)/f(x)),4))"
```
Expected: `0 -> 1 alpha= 0.6065`, `1 -> 0.5 alpha= 1.0`, `0.5 -> 2 alpha= 0.1534`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Theory rungs 1-2 — Monte Carlo principle and Metropolis-Hastings"
```

---

### Task 3: Theory & Proofs — rungs 3–5 (Gibbs invariance, efficiency/ESS, mixing → HMC)

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (append three subsections to `## Theory & Proofs`, before `## Worked Examples`)

- [ ] **Step 1: Rung 3 — The Gibbs sampler (proof of invariance)**

Present the algorithm: partition the parameter into blocks $x=(x_1,\dots,x_d)$ and, in turn, replace each $x_i$ by a fresh draw from its **full conditional** $p(x_i\mid x_{-i})$ (all other blocks held fixed). One full pass is a **sweep**.
> **Theorem (Gibbs leaves $\pi$ invariant).** The full-conditional update of block $i$ — propose $x_i'\sim p(\cdot\mid x_{-i})$, keep $x_{-i}'=x_{-i}$ — preserves $\pi$: if $x\sim\pi$ then the updated state is distributed as $\pi$.

**Proof.** Marginalize the joint over the old $x_i$ and apply the conditional. The density of the updated state $(x_i',x_{-i})$ is
$$
\int \pi(x_i,x_{-i})\,p(x_i'\mid x_{-i})\,dx_i = p(x_i'\mid x_{-i})\int \pi(x_i,x_{-i})\,dx_i = p(x_i'\mid x_{-i})\,\pi(x_{-i}) = \pi(x_i',x_{-i}),
$$
using that $\int\pi(x_i,x_{-i})\,dx_i=\pi(x_{-i})$ is the marginal and $p(x_i'\mid x_{-i})\pi(x_{-i})=\pi(x_i',x_{-i})$ reconstitutes the joint. $\blacksquare$

Then show **Gibbs is Metropolis–Hastings with acceptance 1**: take the proposal $q(x'\mid x)=p(x_i'\mid x_{-i})$ (with $x_{-i}'=x_{-i}$). The MH ratio is
$$
\frac{\pi(x')\,q(x\mid x')}{\pi(x)\,q(x'\mid x)}=\frac{\pi(x_i',x_{-i})\,p(x_i\mid x_{-i})}{\pi(x_i,x_{-i})\,p(x_i'\mid x_{-i})}=1,
$$
since $\pi(x_i',x_{-i})=p(x_i'\mid x_{-i})\pi(x_{-i})$ and $\pi(x_i,x_{-i})=p(x_i\mid x_{-i})\pi(x_{-i})$ cancel exactly — so every Gibbs proposal is accepted. Identify this as **Chapter 6's hierarchical sampler**, now formally justified: each block update there ($\beta_g$, $\mu_\beta$, $\Sigma_\beta$, $\sigma^2$ from their conjugate conditionals) is precisely a Gibbs step, and the sweep leaves the joint posterior invariant. Note systematic-scan vs random-scan, and that Gibbs requires the full conditionals in closed form (when they are not, fall back to a Metropolis step within the sweep — "Metropolis-within-Gibbs").

- [ ] **Step 2: Rung 4 — Efficiency: acceptance, autocorrelation, effective sample size**

The catch with dependent draws: they carry **less information** than i.i.d. ones. For random-walk Metropolis the **step size** sets a tradeoff — too small ⇒ nearly every proposal accepted but the chain inches along, producing highly autocorrelated draws; too large ⇒ most proposals land in the low-density tail and are rejected, so the chain stalls at one point. Define the **lag-$k$ autocorrelation** $\rho_k=\operatorname{Corr}(f(X_t),f(X_{t+k}))$ and the **effective sample size**
$$
\text{ESS}=\frac{n}{1+2\sum_{k\ge1}\rho_k},
$$
the number of i.i.d. draws carrying equivalent information; the **MCMC standard error** of $\hat\mu_n$ is then $\sigma_f/\sqrt{\text{ESS}}$ (contrast the i.i.d. $\sigma_f/\sqrt n$ of rung 1 — positive autocorrelation inflates the denominator's deficit). State (cite `[@robert2004]`, `[@gelman2013]`), without proof, the heuristic **optimal acceptance rates**: about $0.234$ for high-dimensional random-walk Metropolis and about $0.44$ in one dimension. Note that formal multi-chain convergence diagnostics ($\hat R$, rank-normalized ESS) are the subject of Chapter 10; here ESS is the single-chain efficiency measure that the tuning of the code tie-in optimizes.

- [ ] **Step 3: Rung 5 — Burn-in, mixing, and the limits (bridge to HMC)**

**Burn-in / warm-up:** the chain starts out of equilibrium and approaches $\pi$ geometrically (Chapter 7's spectral-gap rate), so the initial transient is discarded before averaging. **Why random-walk Metropolis mixes slowly:** in $d$ dimensions a blind isotropic random walk, to keep acceptance from collapsing, must scale its step like $1/\sqrt d$, so it needs $O(d)$ steps to traverse the distribution — and on a **strongly correlated** posterior it is worse still. Make this concrete with the bivariate-Gaussian Gibbs chain of Worked Example (b): the lag-1 autocorrelation of a coordinate is $\rho^2$, which $\to1$ as the correlation $\rho\to1$, so near-degenerate posteriors crawl. The MMM posterior — channel coefficients that move together because their spends are collinear (Chapter 4) — is exactly this regime. The structural fix is to stop proposing blindly and instead **follow the gradient of $\log\pi$** to make distant, high-acceptance moves: **Hamiltonian Monte Carlo**, Chapter 9. End the Theory section on this baton-pass (and note Chapter 10 supplies the diagnostics that detect slow mixing in practice).

- [ ] **Step 4: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Theory rungs 3-5 — Gibbs, efficiency/ESS, mixing and the HMC bridge"
```

---

### Task 4: Worked Examples

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (the `## Worked Examples` section)

Match Chapter 7's worked-example style: bold step mini-headings, every number shown, closing interpretation.

- [ ] **Step 1: Example (a) — Metropolis acceptance by hand**

Unnormalized target $\pi(x)\propto e^{-x^2/2}$ (note explicitly: the standard normal with $1/\sqrt{2\pi}$ dropped, to demonstrate the normalizer is unneeded), symmetric random-walk proposal, so by the Metropolis rule $\alpha=\min(1,\,e^{-(y^2-x^2)/2})$. Compute three moves: $x=0\to y=1$ gives $\alpha=e^{-1/2}\approx0.607$; $x=1\to y=0.5$ gives ratio $e^{0.375}>1$ so $\alpha=1$ (a move to higher density is always accepted); $x=0.5\to y=2$ gives $\alpha=e^{-1.875}\approx0.153$. Interpret: uphill moves are automatic, downhill moves accepted in proportion to the density ratio, the normalizer never appears, and the long-run visiting frequencies reproduce $\pi$ by the detailed-balance proof of Rung 2.

- [ ] **Step 2: Example (b) — Gibbs for the correlated bivariate normal**

Target $\mathcal N(0,\Sigma)$, $\Sigma=\big[\begin{smallmatrix}1&\rho\\\rho&1\end{smallmatrix}\big]$, $\rho=0.8$. **Derive** the full conditionals from the bivariate-normal conditioning formula: $x\mid y\sim\mathcal N(\rho y,\,1-\rho^2)=\mathcal N(0.8\,y,\,0.36)$ and symmetrically $y\mid x\sim\mathcal N(0.8\,x,\,0.36)$. Carry **two Gibbs sweeps** by hand from a stated start (e.g. $(x_0,y_0)=(2,2)$: draw $x_1$ centered at $0.8\cdot2=1.6$, then $y_1$ centered at $0.8\,x_1$, etc. — show the conditional means and variances; the random draws can be stated symbolically or with one fixed normal draw). State the **mixing fact**: the lag-1 autocorrelation of either coordinate is $\rho^2=0.64$, so as $\rho\to1$ the chain mixes ever more slowly — the concrete motivation for Chapter 9. Tie to MMM: collinear channel spends produce exactly such correlated posteriors.

- [ ] **Step 3: Verify, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
rho=0.8; print('cond var',1-rho**2,'mean coef',rho)
rng=np.random.default_rng(0); n=200000; x=y=0.0; s=(1-rho**2)**0.5; xs=np.empty(n)
for t in range(n):
    x=rho*y+s*rng.standard_normal(); y=rho*x+s*rng.standard_normal(); xs[t]=x
x0=xs-xs.mean(); print('lag1 autocorr',round(np.dot(x0[:-1],x0[1:])/np.dot(x0,x0),4),'rho^2',rho**2)"
```
Expected: `cond var 0.36 mean coef 0.8`, `lag1 autocorr ~0.64 rho^2 0.64`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Worked Examples — Metropolis acceptance and Gibbs for the bivariate normal"
```

---

### Task 5: Code Tie-in

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (the `## Code Tie-in` section)

- [ ] **Step 1: Write a single self-contained ```` ```{python} ```` cell**

Model the structure on Chapter 7's Code Tie-in. Intro paragraph, one ```` ```{python} ```` block (numpy/scipy/matplotlib, seeded `rng = np.random.default_rng(0)`), closing paragraph quoting actual printed numbers. The cell, with small helpers, must:
1. **Random-walk Metropolis on the correlated 2-D Gaussian** ($\Sigma=[[1,0.8],[0.8,1]]$, log-target $-\tfrac12 z^\top\Sigma^{-1}z$). Implement `rw_metropolis(logpi, x0, step, n)`. Sweep step sizes (e.g. `0.1` too small, `1.2` good, `8.0` too large); for each print the **acceptance rate** and a simple **ESS estimate** (via the summed autocorrelations of the first coordinate, truncating at the first negative $\rho_k$). Confirm the good-step sample mean $\approx0$ and sample covariance $\approx\Sigma$. Trace plot of one coordinate for each step size.
2. **Gibbs on the same target.** Implement the closed-form conditional sweeps $x\mid y\sim\mathcal N(\rho y,1-\rho^2)$. Confirm recovered correlation $\approx0.8$; measure the lag-1 autocorrelation ($\approx\rho^2=0.64$); rerun at $\rho=0.99$ and show the autocorrelation jump toward 1 (slow mixing) — print both. (Optional: ESS comparison Gibbs vs RW-Metropolis.)
3. **The payoff — sample a non-conjugate posterior.** A 1-D logistic (binary-conversion) posterior: generate $n=40$ points $x_i\sim\mathcal N(0,1)$, $y_i\sim\text{Bernoulli}(\sigma(\beta_{\text{true}}x_i))$ with $\beta_{\text{true}}=1.5$; prior $\beta\sim\mathcal N(0,4)$; log-posterior $\ell(\beta)=\sum_i[y_i\beta x_i-\log(1+e^{\beta x_i})]-\beta^2/8$ — **no closed form**. Run random-walk Metropolis on $\ell$; build a **fine-grid normalized reference** of the posterior; **assert** the MH posterior mean and sd match the grid values within tolerance (e.g. `atol=0.05` for the mean) and print both. Plot the MH histogram with the grid posterior overlaid. This is the posterior Chapter 5 could not handle.

End figures with `plt.show()`. Closing paragraph quotes: the acceptance/ESS tradeoff across step sizes, the Gibbs autocorrelation at $\rho=0.8$ vs $0.99$, and the MH-vs-grid agreement on the logistic posterior.

- [ ] **Step 2: Extract the cell and run headless**

Save the cell body to `/tmp/ch8_code.py` and run:
```bash
MPLBACKEND=Agg python3 /tmp/ch8_code.py
```
Expected: runs top-to-bottom, all asserts pass (good-step covariance ≈ Σ; Gibbs corr ≈ 0.8; MH logistic mean ≈ grid mean). Make the prose match the actual printed numbers.

- [ ] **Step 3: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Code Tie-in — RW-Metropolis tuning, Gibbs mixing, non-conjugate posterior"
```

---

### Task 6: Exercises (C / B / P / A)

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (the `## Exercises` section with the four `###` tier headings)

Match Chapter 7's exercise style and heading text exactly: `### C -- Conceptual / Reading Comprehension`, `### B -- By Hand`, `### P -- Prove It`, `### A -- Applied / Code`. Each problem fully self-contained (state all data inline). **No links to solutions.**

- [ ] **Step 1: Tier C (3 problems)**
- **C1.** The Monte Carlo principle: why the estimator error $\sigma_f/\sqrt n$ is dimension-free, why we resort to *dependent* MCMC draws, and why (Chapter 7's ergodic theorem) the time-average still converges despite the dependence.
- **C2.** The Metropolis–Hastings acceptance ratio: write it out, explain why $\pi$ enters only through $\pi(y)/\pi(x)$ so the normalizer cancels, and why that is exactly what lets MCMC target an unnormalized posterior $p(y\mid\beta)p(\beta)$. State when the proposal ratio $q(x\mid y)/q(y\mid x)$ may be dropped.
- **C3.** Gibbs as a special case: explain why Gibbs is Metropolis–Hastings with acceptance probability 1, how it relates to Chapter 6's hierarchical sampler, and what "Metropolis-within-Gibbs" does when a full conditional has no closed form. Then explain what acceptance rate, autocorrelation, and ESS measure, and why random-walk Metropolis mixes slowly on correlated posteriors (pointing to HMC).

- [ ] **Step 2: Tier B (3 problems)**
- **B1.** Metropolis acceptance by hand: for the unnormalized target $\pi(x)\propto e^{-x^2/2}$ with symmetric proposal, compute $\alpha$ for $x=0\to y=1.5$, $x=2\to y=1$, and $x=-1\to y=-2$. (Clean: $e^{-1.125}\approx0.3247$; ratio $>1$ so $1$; $e^{-1.5}\approx0.2231$.)
- **B2.** An **asymmetric** proposal where the $q$-ratio does **not** cancel: target $\pi(x)\propto e^{-x}$ on $x>0$ (unnormalized Exp(1)); proposal $q(y\mid x)$ a fixed log-normal or a stated asymmetric density. Given specific $q(y\mid x)$ and $q(x\mid y)$ values and a move, compute the full Hastings ratio $\min(1,[\pi(y)q(x\mid y)]/[\pi(x)q(y\mid x)])$, emphasizing the proposal correction.
- **B3.** Effective sample size: given a chain of $n=10000$ draws with lag autocorrelations $\rho_1=0.6,\rho_2=0.36,\rho_3=0.216,\dots$ (geometric $\rho_k=0.6^k$), compute the integrated autocorrelation $1+2\sum_{k\ge1}0.6^k=1+2\cdot\frac{0.6}{0.4}=4$ and hence $\text{ESS}=10000/4=2500$; interpret.

- [ ] **Step 3: Tier P (3 problems)**
- **P1.** Prove the Metropolis–Hastings kernel satisfies detailed balance with respect to $\pi$ (the $\min$ identity $a\min(1,b/a)=\min(a,b)$ and its $x\leftrightarrow y$ symmetry).
- **P2.** Prove the Gibbs full-conditional update leaves $\pi$ invariant (marginalize then re-condition), and prove that Gibbs equals Metropolis–Hastings with acceptance probability identically 1.
- **P3.** Derive the symmetric-proposal **Metropolis** special case $\alpha=\min(1,\pi(y)/\pi(x))$ from the general Hastings ratio, and derive the full conditional $x\mid y\sim\mathcal N(\rho y,1-\rho^2)$ of the bivariate normal $\mathcal N(0,\big[\begin{smallmatrix}1&\rho\\\rho&1\end{smallmatrix}\big])$ from the conditioning formula.

- [ ] **Step 4: Tier A (2 problems)**
- **A1.** Implement random-walk Metropolis for a target of your choice (e.g. the correlated 2-D Gaussian or a Bayesian-regression posterior); tune the step size across a grid; plot acceptance rate and ESS versus step size and identify the efficient regime; validate the sample moments against the known target.
- **A2.** (i) Implement Gibbs for the bivariate normal and show the lag-1 autocorrelation rising toward 1 as $\rho\to1$. (ii) Sample a **non-conjugate** 1-D posterior (e.g. logistic/binary-conversion with a Gaussian prior, or a half-normal-prior channel effect) with Metropolis and validate the MCMC histogram against a fine-grid normalized reference.

- [ ] **Step 5: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Exercises — C/B/P/A tiers"
```

---

### Task 7: Summary (auto-included)

**Files:**
- Modify: `parts/03-sampling/02-mcmc.qmd` (append `## Summary` as the final section)

- [ ] **Step 1: Write `## Summary`**

Match Chapter 7's Summary format exactly: one-line lead, then **Key concepts** (bold lead-ins) and **Key identities** (**inline math only**). Cover:
- **Key concepts:** Monte Carlo estimation and its dimension-free $1/\sqrt n$ error; MCMC turns the i.i.d. requirement into a chain with stationary distribution $\pi$ (valid by Chapter 7's ergodic theorem); Metropolis–Hastings constructs a reversible chain from any pointwise-evaluable unnormalized $\pi$; the normalizer cancels in the acceptance ratio; Gibbs samples the full conditionals and is MH with acceptance 1 (= the Chapter 6 sampler); efficiency is governed by acceptance rate, autocorrelation, and ESS; random-walk mixing degrades on correlated/high-dimensional posteriors, motivating HMC.
- **Key identities (inline):**
  - Monte Carlo: $\hat\mu_n=\frac1n\sum_i f(X_i)$, standard error $\sigma_f/\sqrt n$ (i.i.d.) or $\sigma_f/\sqrt{\text{ESS}}$ (MCMC).
  - MH acceptance: $\alpha(x,y)=\min\big(1,\frac{\pi(y)q(x\mid y)}{\pi(x)q(y\mid x)}\big)$; symmetric case $\min(1,\pi(y)/\pi(x))$.
  - Detailed balance of MH: $\pi(x)q(y\mid x)\alpha(x,y)=\min(\pi(x)q(y\mid x),\pi(y)q(x\mid y))$, symmetric in $x\leftrightarrow y$.
  - Gibbs: draw $x_i\sim p(x_i\mid x_{-i})$; acceptance $\equiv1$.
  - Bivariate-normal conditional: $x\mid y\sim\mathcal N(\rho y,1-\rho^2)$; coordinate lag-1 autocorrelation $\rho^2$.
  - ESS: $\text{ESS}=n/(1+2\sum_{k\ge1}\rho_k)$.

- [ ] **Step 2: Commit**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add parts/03-sampling/02-mcmc.qmd
git -C "$WT" commit -m "feat(ch08): Summary — key concepts and identities"
```

---

### Task 8: Appendix solutions

**Files:**
- Modify: `appendix/solutions.qmd` (insert before the closing `:::` at line 1272)

- [ ] **Step 1: Append the Chapter 8 solutions block**

Immediately before the final `:::` that closes the `show-solutions` visible div, add a `## Chapter 8 — Markov Chain Monte Carlo` heading followed by full solutions to every exercise from Task 6 (C1–C3, B1–B3, P1–P3, A1–A2). Match the format of the existing `## Chapter 7 — Markov Chains` block: em-dash sub-headings `### C — …` etc., bold labels (`**C1.**`), full worked math for B and P, plain ```` ```python ```` reference sketches (NOT `{python}`) plus a qualitative takeaway for A. Each proof ends `$\blacksquare$`. Verify every numeric answer against NumPy (B1: $0.3247$, $1$, $0.2231$; B3: integrated autocorrelation $4$, ESS $2500$).

- [ ] **Step 2: Verify by-hand answers, then commit**

```bash
MPLBACKEND=Agg python3 -c "
import numpy as np
f=lambda x: np.exp(-x*x/2)
print('B1', round(min(1,f(1.5)/f(0)),4), min(1,f(1)/f(2))>=1, round(min(1,f(-2)/f(-1)),4))
print('B3 iact', 1+2*sum(0.6**k for k in range(1,200)), 'ESS', 10000/(1+2*(0.6/0.4)))"
```
Expected: `B1 0.3247 True 0.2231`; `B3 iact 4.0 ESS 2500.0`.

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
git -C "$WT" add appendix/solutions.qmd
git -C "$WT" commit -m "feat(ch08): appendix solutions for Markov Chain Monte Carlo exercises"
```

---

### Task 9: Final whole-chapter review

**Files:**
- Review: `parts/03-sampling/02-mcmc.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: KaTeX / structure lint**

```bash
WT=/Users/jameshenson/.config/superpowers/worktrees/mmm_foundational_text/ch8-mcmc
grep -n 'begin{align}' "$WT"/parts/03-sampling/02-mcmc.qmd || echo "no bare align — good"
grep -c '^\$\$' "$WT"/parts/03-sampling/02-mcmc.qmd   # expect even
grep -nE '^\## ' "$WT"/parts/03-sampling/02-mcmc.qmd  # expect the 6 template headings in order
grep -nE '^# ' "$WT"/parts/03-sampling/02-mcmc.qmd    # expect '# Markov Chain Monte Carlo'
grep -c '^:::$' "$WT"/appendix/solutions.qmd   # expect 2 (file still closes cleanly)
```
Confirm no bare `\begin{align}`, even `$$` count, the six template headings in order, and citations only the valid keys.

- [ ] **Step 2: Re-run the code cell headless**

```bash
MPLBACKEND=Agg python3 /tmp/ch8_code.py && echo "CODE CELL OK"
```
Expected: `CODE CELL OK`, all asserts pass.

- [ ] **Step 3: Through-line check (manual)**

Confirm every theory subsection ties back to the construct-a-chain / efficiency / HMC-bridge through-line, the two named proofs (MH detailed balance; Gibbs invariance = MH-with-acceptance-1) are present and correct, and Chapter 7 (detailed balance, ergodic theorem) and Chapter 6 (the hierarchical Gibbs sampler) are reused by reference. Fix gaps inline and commit if changed.

---

## Self-Review (completed by plan author)

- **Spec coverage:** Motivation (T1); rungs 1–2 incl. the MH keystone proof (T2); rungs 3–5 incl. Gibbs invariance proof, efficiency/ESS, mixing→HMC (T3); both worked examples (T4); code tie-in with the 3 required behaviors incl. the non-conjugate payoff (T5); C/B/P/A exercises (T6); auto-included Summary (T7); appendix solutions (T8); render/lint gate (T9). All spec success criteria mapped.
- **Placeholder scan:** none — every task carries the actual formulas and NumPy-verified anchor numbers.
- **Consistency:** the MH acceptance $\alpha=\min(1,\pi(y)q(x\mid y)/(\pi(x)q(y\mid x)))$, the symmetric special case, the Gibbs invariance/acceptance-1 identity, and the bivariate-normal conditional $\mathcal N(\rho y,1-\rho^2)$ with autocorrelation $\rho^2$ appear identically across T2/T3/T4/T6/T7. Worked/appendix numbers are NumPy-verified. The non-conjugate logistic posterior numbers are computed-and-asserted in code (not hardcoded), since they are seed/data-dependent.
- **Note:** the real render gate is CI `quarto render` (HTML + PDF) on the PR — Quarto is not installed locally, so local verification is the Python cell + lint checks above.
