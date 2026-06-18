# Chapter 22 — The Prior Store & the Calibration–Optimization Loop: Design Spec

**Date:** 2026-06-18
**Part:** VI — Causal Grounding & Calibration (fourth chapter; the closing capstone)
**File:** `parts/06-causal-calibration/04-prior-store.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter closes Part VI. Chapter 19 named the structural wound (observational MMM cannot
identify the response curve under confounding); Chapter 20 produced interventional, taxonomy-tagged
estimands; Chapter 21 folded *one* estimate into the posterior as a calibration likelihood factor and
showed it heals the ridge. Chapter 22 turns that one-shot calibration into a **system that accumulates
evidence over time** — the **prior store** — and proves the decision-theoretic payoff: as ridge-aimed
calibration studies accumulate, the Chapter 18 EVPI gap decays **monotonically to zero**. The loop the
book has been building — **intervene → identify → calibrate → re-optimize** — provably converges.

The full software-engineering build of the store (real schema, versioning, governance, pipelines) is
**deferred to Part VII** (Chapters 23–26). This chapter is Part VI's *mathematical* capstone: the
evidence-accumulation mathematics, with enough engineering framing to motivate the data product.

**Driving question (locked):** *"You will run not one experiment but many, over years, of varying
quality and freshness. How do you accumulate them correctly — and does the decision actually
converge?"*

**Center of gravity (locked, via AskUserQuestion):** the **EVPI monotonicity** keystone — the loop
closes. The mechanism rungs (sequential = batch; power-prior discounting; hierarchical pooling) support
it; EVPI monotone decay + convergence is the single climax.

**Anchors:** power priors (Ibrahim & Chen 2000); cumulative meta-analysis; random-effects pooling
(DerSimonian & Laird 1986). Same surrogate $S(x;\theta)=\theta\sqrt{x}$ (truth $\theta=2$) and
two-channel $(a_1,a_2)$ ridge as Chapters 18 and 21.

---

## 2. Scope discipline

- Does **not** re-derive the calibration likelihood factor or the secant/tangent/validation taxonomy
  (Chapter 21 owns them); it *accumulates* such factors.
- Does **not** build the store as a real software artifact (schema migrations, storage, access control,
  CI) — Part VII. Rung 5 gives the **conceptual** schema/versioning framing only and hands off.
- **House rule (critical):** explain everything as a **generic** technique. **Never** name
  PyMC-Marketing or any MMM/PPL/sampler/causal library. It is a method, not a product.
  `numpy`/`scipy`/`matplotlib` only.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their own lines;
  even `$$` count.

---

## 3. Theory & Proofs — the rungs (5)

**Rung 1 — The store as a product of likelihood factors (sequential = batch).**
The prior store is, mathematically, the running product of calibration factors. With observational
posterior $p_0(\theta)\propto p(\theta)L_{\text{obs}}$ and a stream of experiments $\hat g_1,\dots,\hat
g_k$, the Bayesian recursion is $p_j(\theta)\propto p_{j-1}(\theta)\,L_j(\hat g_j\mid\theta)$.
**Theorem (cumulative meta-analysis / sequential = batch):** unrolling the recursion,
$$
p_k(\theta)\;\propto\;p_0(\theta)\prod_{j=1}^{k}L_j(\hat g_j\mid\theta),
$$
so folding studies in **one at a time** (sequentially, as the store receives them) yields exactly the
**batch** posterior that conditions on all $k$ at once. **Proof:** induction on $k$ using associativity
of multiplication and conditional independence of the experiments given $\theta$ (each randomization is
a separate data source). This is the **correctness guarantee** of an incremental store: order does not
matter, and a study can be appended without re-fitting from scratch. **MMM reading:** the store is an
append-only ledger of likelihood factors; the live posterior is its product.

**Rung 2 — Power-prior discounting and temporal decay.**
Not every study deserves full weight. The **power prior** raises each factor to a fractional power,
$L_j^{\,w_j}$ with $w_j\in[0,1]$, so the accumulated posterior is $p_0(\theta)\prod_j L_j^{\,w_j}$.
Compose $w_j=(\text{design credibility})\times\rho^{\,t_{\text{now}}-t_j}$ with decay rate $\rho<1$:
recent high-quality studies approach $w=1$, stale or weak ones approach $0$ (a validation-tier study
sets $w=0$ and never enters the likelihood). Note briefly the **normalized-vs-unnormalized** power-prior
subtlety (Ibrahim & Chen): raising a likelihood to a power changes its normalizing constant, which
matters when $w$ itself is treated as unknown; with $w$ fixed (the store's design choice) the factor is
a clean tempered likelihood. Temporal decay makes the store **forgetting**: old evidence is continuously
down-weighted, so the live posterior tracks a drifting media landscape. [@ibrahim2000]

**Rung 3 — Hierarchical pooling of conflicting studies.**
When several studies measure the *same* functional but disagree by more than their stated noise, naive
**fixed-effect** pooling (precision-weighted average, $\hat g_{\text{FE}}=\sum_j w_j\hat g_j/\sum_j
w_j$, variance $1/\sum_j w_j$) is **overconfident** — it treats disagreement as if it did not exist.
**Random-effects** pooling adds a between-study variance $\tau^2$, re-weighting by $w_j^\star=1/(s_j^2+
\tau^2)$ and inflating the pooled variance accordingly. State the DerSimonian–Laird moment estimator
$\tau^2=\max\!\big(0,(Q-(k-1))/(\sum w_j-\sum w_j^2/\sum w_j)\big)$ with Cochran's
$Q=\sum_j w_j(\hat g_j-\hat g_{\text{FE}})^2$. **Connect to Chapter 6:** this is exactly partial
pooling / shrinkage — the studies are exchangeable draws around a population functional, and $\tau^2$ is
the group-level variance. **Conflict detection** is the Chapter 21 **validation tier** operationalized:
a large $Q$ (heterogeneity) is the posterior-predictive alarm that the store's factors are mutually
inconsistent and must be pooled, not multiplied. [@dersimonian1986]

**Rung 4 — KEYSTONE: Monotone EVPI decay (the loop closes).**
In the linear-Gaussian setting, each power-weighted calibration factor on a linear functional
$c_j^\top\theta$ contributes precision $w_j s_j^{-2}c_j c_j^\top$, so after $k$ studies
$$
\Lambda_k=\Lambda_0+\sum_{j=1}^{k}w_j\,s_j^{-2}\,c_j c_j^\top\;\succeq\;\Lambda_{k-1}.
$$
**Lemma (Loewner monotonicity):** $\Lambda_k\succeq\Lambda_{k-1}\Rightarrow\Sigma_k\preceq\Sigma_{k-1}$,
so the posterior variance along **every** direction is non-increasing — *deterministically*, not merely
in expectation. **Theorem (monotone EVPI decay + convergence):** the Chapter 18 EVPI is an increasing
function of the ridge-direction variance $\sigma_{u}^2=u^\top\Sigma_k u$; therefore $\text{EVPI}_k$ is
**non-increasing** in $k$, and if studies are repeatedly aimed at the ridge direction $u$ (so
$\sum_j w_j s_j^{-2}(u^\top c_j)^2\to\infty$), then $\sigma_u^2=O(1/k)\to 0$ and $\text{EVPI}_k\to 0$ at
rate $O(1/k)$. **Proof:** Loewner monotonicity gives the variance ordering; EVPI monotone-increasing in
$\sigma_u^2$ (from the Chapter 18 optimizer's-curse expression, where the gap vanishes as the curve is
pinned) gives the decay; the $1/k$ rate follows from $\sigma_u^2=1/(\Lambda_{0,u}+k\bar\kappa)$ for
equal per-study ridge precision $\bar\kappa$. **This is the loop closing:** the wound named in Chapter 18
and reopened structurally in Chapter 19 is driven to zero by the accumulating store.

**Rung 5 — The store as a data product (engineering framing, brief) and the handoff to Part VII.**
The mathematics above dictates the store's design. Its **schema is the Chapter 20/21 taxonomy**: each
record carries the estimand class (secant / tangent / validation), the target spend interval or
operating point (the constraint direction $c$), the measurement variance $s^2$, the power weight $w$
(credibility × decay), and a timestamp. The store is **append-only and versioned** — a study is never
edited, only superseded — so the live posterior is reproducible from the ledger at any past version
(Rung 1's recursion replayed to that point). Conflicting records are **pooled, not multiplied**
(Rung 3). Keep this brief: the *mechanism* is here; the *system* — storage, migrations, governance,
the pipeline that re-optimizes on update — is built in Part VII (Chapters 23–26). Close on the full loop
and on Part VI: the book has gone from "the optimizer needs a curve the data cannot identify" to "a
versioned ledger of experiments drives the decision-relevant uncertainty monotonically to zero."

---

## 4. Worked Examples

**WE1 — Two studies fold in; sequential equals batch (numeric).**
Scale model $S(x;\theta)=\theta\sqrt x$, observational posterior $\theta\sim\mathcal N(2,0.5^2)$
(precision $4$). Study 1: lift over $[4,9]$ measures $\theta$ with $s_1=0.2$ (precision $25$) →
posterior precision $4+25=29$. Study 2: lift over $[9,16]$, where $S(16;\theta)-S(9;\theta)=\theta(4-3)=
\theta$, measured with $s_2=0.25$ (precision $16$) → posterior precision $29+16=45$. The **batch**
posterior conditioning on both at once has precision $4+25+16=45$ — identical. Posterior std
$1/\sqrt{45}=0.149$, mean held at $2.0$ (studies agree). The store appended Study 2 without re-touching
Study 1: sequential = batch, in arithmetic.

**WE2 — Conflicting studies: fixed-effect overconfidence vs random-effects honesty.**
Two studies measure the scale $\theta$ but disagree: $\hat g_A=2.0$ and $\hat g_B=2.8$, each with
$s=0.2$ (precision $25$). They differ by $0.8$ against a combined std $\sqrt{0.04+0.04}=0.283$ — a
$z\approx2.83$ conflict. **Fixed-effect** pooling gives mean $2.4$ with std $1/\sqrt{50}=0.141$:
tight, but straddling a real $0.8$ gap — overconfident. **Random-effects:** Cochran's $Q=8$ on
$k-1=1$ df; DerSimonian–Laird $\tau^2=(8-1)/25=0.28$; re-weighted variance $1/(0.04+0.28)=3.125$ per
study, pooled std $0.4$ — nearly $3\times$ wider, honestly reflecting the disagreement. The large $Q$
is the conflict alarm: the store **pools** these two records rather than multiplying their factors.

**WE3 — The EVPI staircase (the loop closes, the climax).**
Start at the Chapter 18 ridge: contrast std $0.45$, EVPI $\approx0.12$. Accumulate ridge-aimed
calibration studies of equal precision; the ridge variance follows $\sigma_u^2=1/(\Lambda_{0,u}+
k\bar\kappa)$ and the EVPI traces a monotone staircase down:
$0.119 \to 0.030 \to 0.018 \to 0.013 \to 0.010 \to 0.008$ over studies $k=0,\dots,5$ — non-increasing
and converging toward $0$ at the $O(1/k)$ rate of Rung 4. State plainly: this is the wound of Chapters
18–19 closing as a *limit*, the decision converging on the full-information optimum as the store fills.

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + `matplotlib`; verified headless; seed
`np.random.default_rng(22)`). It:
1. **Sequential = batch (Rung 1 / WE1):** fold two studies into $\theta$ sequentially ($4\to29\to45$)
   and compute the batch precision ($4+25+16=45$); `assert` they are equal and the std is $0.149$.
2. **Random-effects pooling (Rung 3 / WE2):** compute fixed-effect mean/std ($2.4$, $0.141$), Cochran's
   $Q$ ($8$), DerSimonian–Laird $\tau^2$ ($0.28$), and the random-effects std ($0.4$); `assert` the
   random-effects interval is $\sim3\times$ wider than fixed-effect. Figure: the two pooled intervals
   against the two study estimates.
3. **Monotone EVPI decay (Rung 4 / WE3):** the Chapter 18 Monte-Carlo EVPI as ridge-aimed studies
   accumulate; build the staircase $0.12\to0.03\to\cdots$; `assert` the sequence is non-increasing and
   the tail is below a small tolerance (converging). Figure: EVPI vs number of studies (staircase), with
   the $O(1/k)$ envelope.
Every numeric claim asserted.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why does sequential updating equal batch, and why does that make an append-only store correct?
  Why must conflicting studies be *pooled* rather than *multiplied*? Why does temporal decay make the
  store track a drifting landscape? Why is EVPI decay monotone in the Gaussian setting but only
  in-expectation in general?
- **B:** WE1 sequential-vs-batch precision ($4,29,45$); WE2 fixed-effect mean/std, Cochran's $Q$,
  DerSimonian–Laird $\tau^2$, random-effects std; read two steps of the EVPI staircase off the ridge
  variance recursion.
- **P:** (i) Prove sequential = batch (cumulative meta-analysis) by induction. (ii) Prove the Loewner
  monotonicity $\Lambda_k\succeq\Lambda_{k-1}\Rightarrow\Sigma_k\preceq\Sigma_{k-1}$ and hence that the
  ridge-direction variance — and the EVPI — is non-increasing in $k$, with the $O(1/k)$ rate under
  ridge-aimed studies.
- **A:** Extend the Code Tie-in — (i) add a *stale* study with temporal decay weight $w=\rho^{\Delta t}$
  and show the live posterior down-weights it; (ii) drive a conflicting study into the store and show
  the random-effects pooled variance widens while the EVPI staircase stalls until the conflict is
  resolved.

---

## 7. Appendix solutions

Append `## Chapter 22 — The Prior Store & the Calibration–Optimization Loop` to
`appendix/solutions.qmd`, **in chapter order** (after the Chapter 21 block), inside the existing
`content-visible` gated div. Full C/B/P/A; the P-block carries both proofs (sequential = batch by
induction; Loewner monotonicity ⟹ monotone EVPI decay + $O(1/k)$ rate).

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities: cumulative
posterior $p_k\propto p_0\prod_j L_j^{w_j}$; power weight $w_j=\text{cred}\times\rho^{t_{\text{now}}-
t_j}$; random-effects weight $w_j^\star=1/(s_j^2+\tau^2)$ with DerSimonian–Laird $\tau^2$; precision
accumulation $\Lambda_k=\Lambda_0+\sum_j w_j s_j^{-2}c_jc_j^\top$; monotone decay
$\text{EVPI}_k\downarrow0$ at $O(1/k)$. Close tying back to Chapter 18 (the EVPI gap, now driven to
zero) and forward to Part VII (the store as a built data product).

---

## 9. Cross-references

- **Chapter 6 (Hierarchical Models):** partial pooling / shrinkage reused for random-effects study
  pooling.
- **Chapter 16 (Building & Fitting):** the Fisher-information ridge the accumulating store fills in.
- **Chapter 18 (Budget Optimization):** the EVPI gap, now driven monotonically to zero.
- **Chapter 19–21 (Part VI):** the wound, the interventional estimands, and the single-shot calibration
  the store accumulates.
- **Part VII (Chapters 23–26):** the store built as a real software artifact.

---

## 10. Bibliography

Reuse `@ibrahim2000` (power priors — added in Chapter 21) and `@gelman2013` (hierarchical / decision —
already present). Add:
- `@dersimonian1986` — DerSimonian & Laird, "Meta-Analysis in Clinical Trials" (Controlled Clinical
  Trials 1986) — random-effects pooling and the $\tau^2$ moment estimator.
The plan confirms exact BibTeX and that all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor NumPy-verified (done: WE1 $4\to29\to45$, std $0.149$; WE2 FE $2.4/0.141$, $Q=8$,
  $\tau^2=0.28$, RE std $0.40$, conflict $z=2.83$; WE3 staircase
  $0.119,0.030,0.018,0.013,0.010,0.008$ monotone).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`).
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no `\begin{align}`,
  no `psmallmatrix`, valid citation keys, no banned library names.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
