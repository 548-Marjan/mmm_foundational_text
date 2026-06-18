# Chapter 22 — The Prior Store & the Calibration–Optimization Loop Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author Chapter 22 — Part VI's closing capstone — turning Chapter 21's single-shot calibration into an evidence-accumulating prior store, and proving the central decision-theoretic payoff: the Chapter 18 EVPI gap decays monotonically to zero as ridge-aimed studies accumulate (the loop closes).

**Architecture:** Replace the stub body of `parts/06-causal-calibration/04-prior-store.qmd` with the six-heading template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary). Append a Chapter 22 solutions block to `appendix/solutions.qmd`. Add one BibTeX entry. Verify via headless code run + CI `quarto render`, ship as a PR.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), `numpy`/`scipy`/`matplotlib` only, `references.bib`.

---

## Conventions (enforced every task)

- **HOUSE RULES (critical):** NEVER name PyMC-Marketing or any MMM/PPL/sampler/causal library (no PyMC, Stan, Orbit, numpyro, pyro, lightweight_mmm, causalimpact). Everything is a generic **method**, not a product. `numpy`/`scipy`/`matplotlib` only.
- **KaTeX:** `aligned` inside `$$ … $$` only — never bare `\begin{align}`. `$$` delimiters on their own lines. **Even** `$$` count per file. NEVER `\begin{psmallmatrix}`.
- Proofs end with `$\blacksquare$`.
- "Key identities" in the Summary must be a **bulleted list**.
- Keep H1 `# The Prior Store & the Calibration–Optimization Loop` and the anchors line `*Canonical anchors: power priors (Ibrahim & Chen); cumulative meta-analysis.*`. Remove the stub callout.
- Citation keys: `@ibrahim2000` (exists), `@gelman2013` (exists), `@dersimonian1986` (added in T9).
- Commit identity `jlh530i` / `jlh530i@gmail.com`. Branch → PR → user merges; never self-merge; never commit to main.

## NumPy-verified anchors (all confirmed before planning)

- **WE1 sequential = batch:** prior precision $4$; study 1 ($s_1=0.2$) precision $25$ → $29$; study 2 ($s_2=0.25$) precision $16$ → $45$. Batch $4+25+16=45$. Posterior std $1/\sqrt{45}=0.149$, mean $2.0$.
- **WE2 conflict pooling:** $\hat g_A=2.0$, $\hat g_B=2.8$, each $s=0.2$ (precision $25$). Fixed-effect mean $2.4$, std $1/\sqrt{50}=0.141$. Cochran's $Q=8$ on $1$ df. DerSimonian–Laird $\tau^2=(8-1)/25=0.28$. Random-effects weight $1/(0.04+0.28)=3.125$ each, pooled std $0.4$ ($\approx2.8\times$ wider). Conflict $z=0.8/\sqrt{0.08}=2.83$.
- **WE3 EVPI staircase:** ridge std $\sigma_u=1/\sqrt{\Lambda_{0,u}+k\bar\kappa}$ with $\Lambda_{0,u}=1/0.45^2=4.938$, $\bar\kappa=3\Lambda_{0,u}=14.815$; EVPI over $k=0..5$ = $0.119, 0.030, 0.018, 0.013, 0.010, 0.008$ (monotone non-increasing, $O(1/k)$). Uses the Chapter 18 Monte-Carlo EVPI (abar $(2,1)$, $B=9$, $v$-std $0.08$, seed `default_rng(22)`).

## Critical files

- `parts/06-causal-calibration/04-prior-store.qmd` — the chapter (currently a stub; replace body).
- `appendix/solutions.qmd` — append Chapter 22 block before the final `:::` (currently last line; file is 2902 lines; one `content-visible` gated div).
- `references.bib` — add `@dersimonian1986`.
- Voice/structure exemplars (read-only): `parts/06-causal-calibration/03-advanced-calibration.qmd` (Ch21, immediate predecessor; owns the calibration factor + EVPI heal), `parts/05-mmm-modeling/04-budget-optimization.qmd` (Ch18; the EVPI machinery to reuse), `parts/02-regression-bayes/03-hierarchical-regression.qmd` (Ch6; partial pooling reused in Rung 3).
- Authoritative design: `docs/superpowers/specs/2026-06-18-chapter-22-prior-store-design.md`.

---

### Task 1: Front matter + Motivation

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`

- [ ] **Step 1: Strip the stub.** Keep lines 1–3 (H1 + blank + anchors line). Delete the stub callout and the placeholder italic lines. Keep the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises` with its four `### C/B/P/A` subheads, `## Summary`).

- [ ] **Step 2: Write Motivation (3–4 paragraphs).** Driving question: *"You will run not one experiment but many, over years, of varying quality and freshness. How do you accumulate them correctly — and does the decision actually converge?"* Cover:
  - Chapter 21 folded *one* calibration estimate into the posterior and healed the ridge once. Real practice runs many experiments over time, of varying credibility and freshness; the question is how to accumulate them **correctly** and whether the decision **converges**.
  - The prior store is the answer: an append-only, versioned ledger of calibration likelihood factors keyed by the Chapter 20/21 taxonomy. Preview the four mechanisms (sequential = batch; power-prior discounting/decay; pooling conflicts; monotone EVPI decay) and the climax — the Chapter 18 EVPI gap driven to zero.
  - The full software build is Part VII; this chapter is the mathematics of evidence accumulation.
  - Name the anchors: $S(x;\theta)=\theta\sqrt x$ (truth $\theta=2$) and the two-channel $(a_1,a_2)$ ridge from Chapters 18/21.
  - Do NOT name any library.

- [ ] **Step 3: Verify** even `$$` count, headings present/ordered. Commit: `feat(ch22): motivation`.

---

### Task 2: Theory & Proofs — Rungs 1–3 (mechanism)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`. Match Ch21 rung prose density and proof rigor.

- [ ] **Step 1: Rung 1 — The store as a product of likelihood factors (sequential = batch).** Recursion $p_j\propto p_{j-1}L_j$. **Theorem (cumulative meta-analysis):**
$$
p_k(\theta)\;\propto\;p_0(\theta)\prod_{j=1}^{k}L_j(\hat g_j\mid\theta).
$$
  **Proof:** induction on $k$ using associativity of multiplication and conditional independence of experiments given $\theta$. End `$\blacksquare$`. State the correctness guarantee (order-independent, append without re-fit) and the MMM reading (append-only ledger; live posterior is its product).

- [ ] **Step 2: Rung 2 — Power-prior discounting and temporal decay.** Accumulated posterior $p_0\prod_j L_j^{w_j}$, $w_j\in[0,1]$, $w_j=(\text{design credibility})\times\rho^{\,t_{\text{now}}-t_j}$, $\rho<1$. $w=0$ = validation-tier (never enters likelihood). Note briefly the normalized-vs-unnormalized power-prior subtlety (Ibrahim & Chen): raising a likelihood to a power changes its normalizing constant, which matters only when $w$ is unknown; with $w$ fixed it is a clean tempered likelihood. Temporal decay = a forgetting store tracking a drifting landscape. Cite `[@ibrahim2000]`.

- [ ] **Step 3: Rung 3 — Hierarchical pooling of conflicting studies.** Fixed-effect pooling $\hat g_{\text{FE}}=\sum_j w_j\hat g_j/\sum_j w_j$, variance $1/\sum_j w_j$ — overconfident under disagreement. Random-effects: $w_j^\star=1/(s_j^2+\tau^2)$ with DerSimonian–Laird $\tau^2=\max(0,(Q-(k-1))/(\sum w_j-\sum w_j^2/\sum w_j))$, Cochran's $Q=\sum_j w_j(\hat g_j-\hat g_{\text{FE}})^2$. Connect to Chapter 6 partial pooling / shrinkage ($\tau^2$ = group-level variance, studies exchangeable around a population functional). Conflict detection (large $Q$) = the Chapter 21 validation-tier alarm: pool, don't multiply. Cite `[@dersimonian1986]`.

- [ ] **Step 4: Verify** even `$$` count, no bare `\begin{align}`/`psmallmatrix`, Rung 1 proof ends `$\blacksquare$`. Commit: `feat(ch22): theory rungs 1-3 (sequential=batch, power prior, pooling)`.

---

### Task 3: Theory & Proofs — Rungs 4–5 (keystone + handoff)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`

- [ ] **Step 1: Rung 4 — KEYSTONE: Monotone EVPI decay (the loop closes).** Precision accumulation
$$
\Lambda_k=\Lambda_0+\sum_{j=1}^{k}w_j\,s_j^{-2}\,c_j c_j^\top\;\succeq\;\Lambda_{k-1}.
$$
  **Lemma (Loewner monotonicity):** $\Lambda_k\succeq\Lambda_{k-1}\Rightarrow\Sigma_k\preceq\Sigma_{k-1}$ (variance non-increasing along every direction, deterministically). **Theorem:** the Chapter 18 EVPI is increasing in the ridge-direction variance $\sigma_u^2=u^\top\Sigma_k u$; hence $\text{EVPI}_k$ is non-increasing, and with studies aimed at $u$, $\sigma_u^2=1/(\Lambda_{0,u}+k\bar\kappa)=O(1/k)\to0$ so $\text{EVPI}_k\to0$ at $O(1/k)$. **Proof:** Loewner gives the ordering (cite that $A\succeq B\succ0\Rightarrow B^{-1}\succeq A^{-1}$, the operator-antitone inverse); EVPI monotone in $\sigma_u^2$ from the Chapter 18 optimizer's-curse gap (vanishing as the curve is pinned); the $1/k$ rate from the recursion. End `$\blacksquare$`. State plainly: the wound of Chapters 18–19 driven to zero by the accumulating store.

- [ ] **Step 2: Rung 5 — The store as a data product (brief) and handoff to Part VII.** Schema = the Chapter 20/21 taxonomy: each record = (estimand class secant/tangent/validation, target spend / constraint direction $c$, measurement variance $s^2$, power weight $w$, timestamp). Append-only + versioned (a study is superseded, never edited; the posterior at any past version is Rung 1's recursion replayed). Conflicts pooled not multiplied (Rung 3). Keep brief — the system (storage, migrations, governance, the re-optimize-on-update pipeline) is Part VII (Chapters 23–26). Close on the full loop and Part VI.

- [ ] **Step 3: Verify** even `$$` count, `@dersimonian1986`/`@ibrahim2000` used, Rung 4 proof ends `$\blacksquare$`. Commit: `feat(ch22): theory rungs 4-5 (keystone EVPI monotonicity + data-product handoff)`.

---

### Task 4: Worked Examples (WE1–WE3)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`. Use the EXACT anchor numbers above.

- [ ] **Step 1: WE1 — Two studies fold in; sequential = batch.** $\theta\sim\mathcal N(2,0.5^2)$ (precision $4$); study 1 over $[4,9]$, $s_1=0.2$ (precision $25$) → $29$; study 2 over $[9,16]$ where $S(16)-S(9)=\theta$, $s_2=0.25$ (precision $16$) → $45$; batch $4+25+16=45$ identical; std $0.149$, mean $2.0$. Store appended study 2 without retouching study 1.

- [ ] **Step 2: WE2 — Conflicting studies: fixed-effect vs random-effects.** $\hat g_A=2.0$, $\hat g_B=2.8$, each $s=0.2$; differ by $0.8$ vs combined std $0.283$ ($z\approx2.83$). Fixed-effect: mean $2.4$, std $0.141$ (overconfident). Random-effects: $Q=8$ on $1$ df, $\tau^2=0.28$, re-weighted variance $3.125$/study, pooled std $0.4$ ($\approx2.8\times$ wider). Large $Q$ = conflict alarm; pool not multiply.

- [ ] **Step 3: WE3 — The EVPI staircase (climax).** From the Chapter 18 ridge (contrast std $0.45$, EVPI $\approx0.12$), accumulate equal-precision ridge-aimed studies; $\sigma_u^2=1/(\Lambda_{0,u}+k\bar\kappa)$; EVPI staircase $0.119\to0.030\to0.018\to0.013\to0.010\to0.008$, monotone, $O(1/k)\to0$. The wound closing as a limit; the decision converging on the full-information optimum.

- [ ] **Step 4: Verify** even `$$` count; numbers match anchors. Commit: `feat(ch22): worked examples WE1-WE3`.

---

### Task 5: Code Tie-in (single runnable cell)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`

- [ ] **Step 1: Write a single ```{python}``` cell** under `## Code Tie-in`, with a short prose preface (Ch21-style). `numpy` + `matplotlib`; seed `np.random.default_rng(22)`; figures end `plt.show()`. Three blocks, every numeric claim `assert`ed:
  1. **Sequential = batch (WE1):** fold two studies into $\theta$ sequentially ($4\to29\to45$); compute batch precision ($4+25+16$); `assert` equal and std $\approx0.149$, mean $\approx2.0$.
  2. **Random-effects pooling (WE2):** fixed-effect mean/std ($2.4$, $0.141$); Cochran's $Q$ ($8$); DerSimonian–Laird $\tau^2$ ($0.28$); random-effects std ($0.4$); `assert` RE std $> 2.5\times$ FE std. Figure: the two study estimates with FE and RE pooled intervals.
  3. **Monotone EVPI decay (WE3):** the Chapter 18 Monte-Carlo EVPI as ridge-aimed studies accumulate ($k=0..5$); build the staircase; `assert` the sequence is non-increasing (each $\le$ previous + tiny tol) and the tail $<0.02$. Figure: EVPI vs number of studies (staircase) with an $O(1/k)$ reference envelope.

- [ ] **Step 2: Extract and run headless.** Copy the cell body to `/tmp/ch22_code.py`, run `MPLBACKEND=Agg python3`. Expected: all asserts pass; prints show $45=45$ / std $0.149$, $\tau^2=0.28$ / RE std $0.4$, the monotone staircase $0.12\to\cdots\to\approx0.008$. Fix until clean.

- [ ] **Step 3: Verify** even `$$` count; no banned libraries (`grep -inwE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm'` returns nothing). Commit: `feat(ch22): code tie-in (sequential=batch, pooling, EVPI staircase)`.

---

### Task 6: Exercises (C / B / P / A — self-contained)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`

- [ ] **Step 1: Write the four exercise blocks** under the existing `### C/B/P/A` subheads. No inline solution links.
  - **C:** (C1) Why does sequential updating equal batch, and why does that make an append-only store correct? (C2) Why must conflicting studies be *pooled* not *multiplied*? (C3) Why does temporal decay make the store track a drifting landscape, and why is EVPI decay monotone in the Gaussian setting but only in-expectation in general?
  - **B:** (B1) WE1 sequential-vs-batch precision ($4,29,45$) and posterior std. (B2) WE2 fixed-effect mean/std, Cochran's $Q$, DerSimonian–Laird $\tau^2$, random-effects std. (B3) Read two steps of the EVPI staircase off the ridge-variance recursion $\sigma_u^2=1/(\Lambda_{0,u}+k\bar\kappa)$.
  - **P:** (P1) Prove sequential = batch (cumulative meta-analysis) by induction. (P2) Prove $\Lambda_k\succeq\Lambda_{k-1}\Rightarrow\Sigma_k\preceq\Sigma_{k-1}$ and hence that the ridge-direction variance and the EVPI are non-increasing in $k$, with the $O(1/k)$ rate under ridge-aimed studies.
  - **A:** (A1) Add a *stale* study with decay weight $w=\rho^{\Delta t}$ and show the live posterior down-weights it. (A2) Drive a conflicting study into the store and show the random-effects pooled variance widens while the EVPI staircase stalls until the conflict resolves.

- [ ] **Step 2: Verify** even `$$` count. Commit: `feat(ch22): exercises C/B/P/A`.

---

### Task 7: Summary (auto-included)

**Files:** Modify `parts/06-causal-calibration/04-prior-store.qmd`

- [ ] **Step 1: Write `## Summary`** with a one-paragraph wrap, then **bulleted** "Key concepts" and **bulleted** "Key identities" (inline math). Identities (each a bullet):
  - cumulative posterior $p_k\propto p_0\prod_j L_j^{w_j}$
  - power weight $w_j=\text{cred}\times\rho^{\,t_{\text{now}}-t_j}$
  - random-effects weight $w_j^\star=1/(s_j^2+\tau^2)$ with DerSimonian–Laird $\tau^2$
  - precision accumulation $\Lambda_k=\Lambda_0+\sum_j w_j s_j^{-2}c_jc_j^\top$
  - monotone decay $\text{EVPI}_k\downarrow0$ at rate $O(1/k)$
  Close tying back to Chapter 18 (the EVPI gap driven to zero) and forward to Part VII (the store built as a data product).

- [ ] **Step 2: Verify** Key identities is a bulleted list; even `$$` count. Commit: `feat(ch22): summary`.

---

### Task 8: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Append `## Chapter 22 — The Prior Store & the Calibration–Optimization Loop`** immediately before the final `:::` (currently the last line). Full C/B/P/A:
  - **C1–C3:** sequential=batch ⟹ order-independent append-only correctness; conflicting studies multiplied would falsely compound into overconfidence, so pool with $\tau^2$; temporal decay continuously down-weights old evidence so the posterior tracks drift; Gaussian precision adds deterministically (Loewner) so EVPI is pathwise monotone, whereas in general a surprising study can transiently raise posterior variance so only the *expected* posterior EVPI is non-increasing.
  - **B1:** precision $4\to29\to45$, std $1/\sqrt{45}=0.149$, mean $2.0$. **B2:** FE mean $2.4$ std $0.141$; $Q=8$, df $1$; $\tau^2=(8-1)/25=0.28$; RE std $0.4$. **B3:** $\sigma_u^2(0)=0.2025$ (EVPI $\approx0.119$), $\sigma_u^2(1)=1/(4.938+14.815)=0.0506$ (std $0.225$, EVPI $\approx0.030$).
  - **P1:** induction — base $k=1$ is one update; step multiplies the $k$-th conditionally-independent factor onto the unrolled product. End `$\blacksquare$`. **P2:** $\Lambda_k=\Lambda_{k-1}+w_k s_k^{-2}c_kc_k^\top\succeq\Lambda_{k-1}$ (PSD rank-one increment); inverse is operator-antitone so $\Sigma_k\preceq\Sigma_{k-1}$; thus $u^\top\Sigma_k u$ non-increasing; EVPI increasing in that variance ⟹ non-increasing; with $u^\top c_j\ne0$ and equal increments $\sigma_u^2=1/(\Lambda_{0,u}+k\bar\kappa)=O(1/k)$. End `$\blacksquare$`.
  - **A1, A2:** describe expected behavior — a stale study's weight $\rho^{\Delta t}\to0$ so its precision contribution vanishes and the posterior reverts toward the undeweighted product; a conflicting study raises $Q$ and $\tau^2$, widening the pooled variance, so the ridge variance (and EVPI) stops falling — the staircase stalls — until the conflict is resolved (more studies, or a credibility down-weight).

- [ ] **Step 2: Verify** even `$$` count in `appendix/solutions.qmd`; block before final `:::`; proofs end `$\blacksquare$`. Commit: `feat(ch22): appendix solutions (C/B/P/A)`.

---

### Task 9: Bibliography + final review + PR

**Files:** Modify `references.bib`; verify `parts/06-causal-calibration/04-prior-store.qmd`, `appendix/solutions.qmd`

- [ ] **Step 1: Add the BibTeX entry** to `references.bib`:

```bibtex
@article{dersimonian1986,
  author  = {DerSimonian, Rebecca and Laird, Nan},
  title   = {Meta-Analysis in Clinical Trials},
  journal = {Controlled Clinical Trials},
  year    = {1986},
  volume  = {7},
  number  = {3},
  pages   = {177--188}
}
```

- [ ] **Step 2: Lint pass on both files.** Confirm: even `$$` count in the chapter and in the appendix; no bare `\begin{align}`; no `\begin{psmallmatrix}`; six template headings present and in order; H1 and anchors line intact; stub callout gone; all proofs end `$\blacksquare$`; Key identities bulleted; only `@ibrahim2000`/`@gelman2013`/`@dersimonian1986` cited and all resolve; no banned library names (`grep -inwE 'pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm|marketing'`).

- [ ] **Step 3: Re-run the extracted code cell headless** one final time (`MPLBACKEND=Agg python3`); all asserts pass.

- [ ] **Step 4: Verify local CI hooks** (`git config --get core.hooksPath`; this repo ships no tracked hook suite, so the default is expected — CI render is the gate). Commit: `feat(ch22): random-effects meta-analysis reference (DerSimonian & Laird 1986)`.

- [ ] **Step 5: Push the branch and open a PR** against `main` (title `feat: Chapter 22 — The Prior Store & the Calibration–Optimization Loop`; body summarizing the keystone + the EVPI staircase + Part VI completion). Start the background merge-watcher: `until gh pr view <N> --json state -q .state | grep -qE 'MERGED|CLOSED'; do sleep 120; done`.

- [ ] **Step 6: Watch CI `quarto render` (HTML + PDF) to a green conclusion** — the real gate. If it fails, read the log, fix, push, re-watch. Report the green PR; the user merges.

---

## Self-Review (done at plan time)

- **Spec coverage:** Rungs 1–5 → Tasks 2–3; WE1–3 → Task 4; code 3 blocks → Task 5; exercises C/B/P/A → Task 6; summary → Task 7; appendix → Task 8; bib `@dersimonian1986` → Task 9. All spec sections mapped.
- **Anchor consistency:** WE1 ($4/29/45$, std $0.149$), WE2 (FE $2.4/0.141$, $Q=8$, $\tau^2=0.28$, RE std $0.40$, $z=2.83$), WE3 staircase ($0.119,0.030,0.018,0.013,0.010,0.008$) are identical across Theory, Worked Examples, Code, and Appendix tasks — all NumPy-verified above.
- **No placeholders:** every task carries the actual math/numbers/commands.
- **House rules:** library ban, KaTeX rules, even `$$`, bulleted Key identities, `$\blacksquare$`, branch→PR→user-merges all enforced per task.
