# Funnel-Aware Bayesian MMM — Background for Teaching the Fundamental Tools

*Internal learning notes. Source material: three PyMC-Labs blog posts —*
- *["Funnel-Aware MMM"](https://www.pymc-labs.com/blog-posts/funnel-aware-mmm) (Part I — business case)*
- *["Full-Funnel MMM Optimization"](https://www.pymc-labs.com/blog-posts/full-funnel-mmm-optimization) (optimization layer)*
- *["Extending the Funnel-Aware Bayesian MMM"](https://www.pymc-labs.com/blog-posts/extending-funnel-aware-bayesian-marketing-mix-model) (production extensions)*

*The blog's promised **Part II** ("technical architecture: censored demand modeling, the two-likelihood design, convergence engineering") does not appear to be published yet — it is referenced only as narrative text, and candidate URLs 404 to the site's SPA shell. The concrete equations and code it would contain are, however, available in PyMC-Marketing's **official documentation notebooks**, which are the canonical implementation of this exact method:*
- *["Introduction: Measuring Upper-Funnel Impact"](https://www.pymc-marketing.io/en/latest/notebooks/mmm/mmm_intro_upper_funnel.html) — the minimal two-model version.*
- *["Measuring Upper-Funnel Impact with PyMC-Marketing"](https://www.pymc-marketing.io/en/stable/notebooks/mmm/mmm_upper_funnel_causal_approach.html) — the advanced version with the full data-generating process, backdoor criterion, and naive-model comparison.*
- *Everything grounded in those notebooks is captured in **Appendix A** below (exact equations, priors, and API).*

The three **blog** posts are deliberately conceptual — they **name** the tools but rarely write the equations. These notes pair each named tool with its standard formulation so the methods can actually be taught and implemented. Where a formula is *my supplement* rather than something stated in the posts, it is flagged **[supplied]**; where it is confirmed verbatim by the official docs, it is flagged **[docs]**.

---

## 0. The one-paragraph summary

A standard MMM regresses a single outcome (leads/revenue) on transformed channel spends and treats every channel as an independent input. A **funnel-aware** MMM instead encodes the *causal ordering* of the funnel: upper-funnel spend (video, demand-gen) creates latent demand that only later shows up as lower-funnel activity (branded / non-branded search), which finally converts to leads. Two facts break the standard model and motivate every tool below: (1) upper-funnel effects are **mediated** through lower-funnel channels, so ignoring the path biases upper-funnel contribution downward; and (2) lower-funnel channels hit **daily budget caps**, so observed spend is a *censored* view of true demand. The engagement described (Nürnberger Versicherung, German insurance) reported **27%+ CPL reduction** and continued use into 2026.

---

## 1. The funnel as a causal DAG (mediation)

**What the posts apply.** They model the funnel as a chain, not a flat regression:

```
upper-funnel spend ──direct──────────────────────────────► leads
        │                                                    ▲
        └──indirect──► lower-funnel demand ──► (caps) ──► lower-funnel spend
```

Upper-funnel channels have **two effects**: a *direct* effect on leads (standard MMM machinery, and broader than last-touch — it also lifts organic search), and an *indirect* effect **mediated** through lower-funnel spend.

**The fundamental tool: mediation analysis on a DAG.** This is the classic mediator decomposition. With treatment `X` (upper spend), mediator `M` (lower demand/spend), outcome `Y` (leads):

- **[supplied]** Total effect = Direct effect + Indirect (mediated) effect.
- In a linear sketch: `M = a·X + …`, `Y = c'·X + b·M + …` ⇒ indirect effect `= a·b`, direct effect `= c'`, total `= c' + a·b`.
- In the MMM the relationships are nonlinear (adstock+saturation on each edge), so the decomposition is done by *simulation on the graph* rather than a product of coefficients — but the conceptual split (direct vs. mediated) is exactly mediation analysis.

**Why it matters for teaching.** Everything else (censoring, the `do` operator, joint optimization) exists to estimate this DAG faithfully. Teach the DAG first; the tools are answers to "how do I estimate each edge honestly?"

**Key assumption:** temporal/causal ordering — upper-funnel precedes and drives lower-funnel; no major unmeasured confounder of the mediator→outcome edge.

---

## 2. Foundational media transforms: adstock & saturation

**What the posts apply.** "Adstock (decay over time) and saturation (diminishing returns) are applied to **both** the upper-to-lower relationship and the lower-to-target relationship — each with independently estimated parameters." Concrete intuitions given:
- Upper-funnel video → leads: **long adstock** (awareness lingers), **moderate saturation**.
- Video → search: **short adstock** (queries happen fast), **steep saturation**.

**The fundamental tools (standard forms) [supplied]:**

*Adstock / carryover* — geometric adstock with retention `θ ∈ [0,1)` (optionally with lag `L` and delay `d`):
```
x*_t = Σ_{l=0}^{L} θ^l · x_{t-l}   (normalized so weights sum to 1)
```
A longer effective memory ⇔ larger `θ`. "Long decay" = `θ` near 1; "short decay" = `θ` small.

*Saturation / diminishing returns* — Hill or logistic curve. Hill form:
```
S(x) = x^s / (κ^s + x^s)          (s = shape/slope, κ = half-saturation point)
```
"Steep saturation" ⇔ small `κ` (saturates at low spend); "moderate saturation" ⇔ larger `κ`.

**Teaching point:** the novelty is not the transforms themselves — it's that **each edge of the funnel gets its own adstock+saturation parameters**. The upper→lower edge and lower→outcome edge decay and saturate differently, and both are estimated jointly.

---

## 3. Censored likelihood for budget caps (the core innovation)

**What the posts apply.** Lower-funnel channels have **daily budget caps**. When a campaign "hits its daily cap by noon," extra demand created upstream never converts — observed spend understates true demand. Modeling this with an ordinary likelihood makes the model interpret capped days as *saturation* ("returns flattened"), systematically **underestimating** channel effectiveness. The fix: a **censored likelihood** that treats capped days as *right-censored* observations of a latent demand.

Stated mechanics: "when a paid search campaign hits its daily cap, the model estimates what demand *would have been* without that constraint"; posterior says "demand was at least €500, probably higher, with a distribution shaped by what the model knows about the upper-funnel drivers." Also: periods before a channel became active are **masked**.

**The fundamental tool: censored / Tobit regression.** Observed spend is
```
y_obs,t = min(D_t, cap_t)
```
where `D_t` is latent demand (driven by upper-funnel adstock/saturation + baseline). The likelihood splits by whether the cap binds:
```
uncensored (y_obs < cap):   contribute  f(y_obs | μ_t, σ)          -- the density
right-censored (y_obs = cap): contribute  P(D_t ≥ cap_t) = 1 − F(cap_t | μ_t, σ)  -- the tail mass
```
with `f`, `F` the pdf/cdf of the assumed demand distribution (e.g. Normal/LogNormal). **[supplied]** This is exactly the Tobit type-I likelihood; the censoring threshold `cap_t` **varies by day** and comes from domain knowledge (paid-search team input). In PyMC this is `pm.Censored(...)`.

**Consequence stated in the posts:** using censored likelihoods "revealed **higher upper-funnel estimates**" than earlier iterations — because demand suppressed by caps is no longer misread as weak upstream effect or as saturation.

**Teaching sequence:** Tobit/censored regression → why `min(demand, cap)` is right-censoring → why ignoring it confounds *constraint* with *saturation*.

---

## 4. The two-likelihood (joint) architecture

**What the posts apply.** A single joint model with **two likelihoods**:
1. **Spend-to-spend** (upper-funnel drives lower-funnel demand) — observed at **daily** granularity, fine-grained, information-rich.
2. **Spend-to-leads** (lower-funnel → conversions) — coarser, noisier.

The daily lower-funnel likelihood carries **strong information about the upper-funnel adstock/saturation parameters** — but that information enters *through the censoring constraint*, not through direct observation. This creates **posterior geometry challenges**: "two likelihoods can pull parameters in different directions, inducing multimodality or slow mixing," requiring careful priors and "convergence engineering."

**The fundamental tools:**
- **Joint likelihood factorization [supplied]:** `p(data | θ) = p(leads | lower, θ) · p(lower_spend | upper, caps, θ)`. Sharing parameters across factors is what lets daily spend data identify upstream effects.
- **Posterior geometry / identifiability:** ridges and multimodality when two likelihoods constrain shared parameters. This is why the priors and reparameterizations (§6) matter and why MCMC diagnostics (§9) are non-negotiable.

---

## 5. Time-varying baselines via HSGP

**What the posts apply.** "**HSGP (Hilbert Space Gaussian Process)** allows the baseline spend of lower-funnel channels to fluctuate over time, capturing seasonality and organic growth." The intercept is not constant — it's a smooth latent function of time.

**The fundamental tool: Gaussian-process trend, HSGP approximation [supplied].**
- A GP prior on the baseline: `baseline(t) ~ GP(0, k(t,t'))` with e.g. a Matérn/exponential-quadratic kernel controlling smoothness/lengthscale.
- Exact GPs cost `O(n³)`. **HSGP** approximates the kernel with a truncated basis (eigenfunctions of the Laplacian on a bounded domain):
```
f(t) ≈ Σ_{j=1}^{m} β_j · φ_j(t),   β_j ~ N(0, S(ω_j))
```
where `φ_j` are fixed sinusoidal basis functions and `S` is the kernel's spectral density evaluated at eigenvalues `ω_j`. Cost drops to roughly `O(n·m)` with `m ≪ n` basis terms.
- **Teaching point:** HSGP = "a GP you can afford," used here as a *flexible time-varying intercept* so that seasonality/organic drift isn't misattributed to paid channels.

---

## 6. Discrete campaign-change multipliers + Dirichlet identification

**What the posts apply.** Channel effectiveness shifts **abruptly** at known dates (e.g. creative switches to performance-focused on **July 1**). Modeled as **discrete per-segment multipliers**, not a smooth trend: segment `i` effectiveness `= β · m_i`. Naively, `N` free multipliers with a baseline `β` are **unidentified** — you can scale all `m_i` up and `β` down with no change in fit (a **posterior ridge**). Fix: a **Dirichlet prior** forcing the multipliers onto a simplex (they sum to a constant), removing one degree of freedom and turning them into **relative shares**.

**The fundamental tool: relative-share reparameterization for identification [supplied].**
- Ridge/aliasing: `(m_i, β) → (λ m_i, β/λ)` leaves the likelihood invariant.
- Constraint `Σ_i m_i = C` (Dirichlet: `m ~ Dirichlet(α)`, then scaled) removes the aliasing. Multipliers now mean "campaign B is 80–120% as effective as campaign A" — a prior you can actually elicit from a channel manager.
- Same *family* of trick as softmax/sum-to-one constraints and centered/non-centered reparameterization: **fix the scale to kill the ridge.**

**Decision-not-metric principle (stated):** this extension does **not** improve held-out MAPE, yet is kept because the model exists to answer a **causal** question ("did switching to performance creative improve ROI?"), not to minimize forecast error. *"Build for the decision, not the metric."*

---

## 7. Latent-quality × observable-reach decomposition (influencers)

**What the posts apply.** Influencer spend impact splits into **observable reach** (from follower count) and **latent quality** (learned from outcomes). Two 500K-follower creators can have identical reach but different conversion power (Influencer A converts 30% better than B). Reach **saturates** in follower count (2M followers ≠ 2× the impact of 1M).

**The fundamental tool: latent-variable decomposition [supplied].**
```
effect_k = r(followers_k) · q_k
r(F)  = saturating function of follower count  (sublinear, jointly estimated)
q_k   = latent per-influencer quality, learned from spend×outcome history
```
- `r(·)` is a saturation curve (same §2 idea) on an *audience-size* input.
- `q_k` is a latent parameter with a prior (e.g. `q_k ~ LogNormal`) identified only through the outcome likelihood — a mini measurement model. This is how the model "learns heterogeneity without explicit quality ratings."

---

## 8. Causal interventions: the `do` operator

**What the posts apply.** Optimization is done by **intervening** on the fitted causal graph: "PyMC-Marketing's `BudgetOptimizer` was extended to jointly optimize across both funnel layers using **PyMC's `do` operator** to simultaneously replace channel spend and lower-funnel budget caps in the computational graph." The optimizer then "traces how upper-funnel spend flows through the latent demand for lower-funnel channels (accounting for censoring) and ultimately into leads."

**The fundamental tool: do-calculus / graph surgery [supplied].**
- `do(X = x)` = *replace* the node `X` with a fixed value, severing its incoming edges, then propagate through the rest of the structural model to read off the interventional distribution `P(Y | do(X=x))`.
- Practically (PyMC's `do`): swap the spend/cap nodes in the compute graph for candidate values and forward-simulate the posterior of leads. This is what makes the optimizer *causal* (evaluating counterfactual allocations) rather than merely correlational.
- **Teaching bridge:** connects Pearl's intervention calculus to a concrete software operation on a probabilistic program.

---

## 9. Full-funnel budget optimization

**What the posts apply.** Two decision problems:
1. **Budget split:** given total monthly budget €X, split across {demand-gen, video, performance-max, paid-search} **and** set daily **caps** for paid-search brand/non-brand. Decision variables = channel spends **and** lower-funnel caps.
2. **Marginal cap analysis:** "If I increase video spend by 20%, by how much should I raise the paid-search caps?" — because extra upstream demand is **wasted** if the downstream cap binds by noon.

The optimizer "jointly searches the space of all combinations to **maximize expected leads (or minimize cost-per-lead)** subject to the total budget constraint," balancing **direct and indirect** effects across both layers.

**The fundamental tool: constrained nonlinear optimization [supplied].**
```
maximize   E[ Leads( s, c ) ]           over channel spends s and caps c
subject to Σ_j s_j ≤ B                  (total budget)
           s_j ≥ 0,  c_k ≥ 0            (box constraints / per-channel bounds)
```
where `Leads(s,c)` is the *causal* response surface obtained via `do(·)` on the fitted graph (adstock→saturation→censoring→leads).
- **Marginal ROAS / KKT intuition:** at the optimum, the budget flows so that the **marginal return per euro is equalized across channels** (equal to the constraint's shadow price / Lagrange multiplier `λ`):
```
∂E[Leads]/∂s_j = λ   for every channel used interior to its bounds
```
This is the standard water-filling result of constrained allocation — the same logic taught in classical MMM budget optimization, but here the response surface *includes the mediated funnel path and the censoring correction.*
- **Why caps are decision variables:** raising a cap relaxes a downstream censoring constraint, so upstream marginal value depends on the cap — the two must be optimized **jointly**, not sequentially.

---

## 10. Time-slice cross-validation with dynamic feature sets

**What the posts apply.** Channels **enter mid-series** (a new platform appears in later months). Standard time-series CV breaks: a late channel `D` has **no signal** in early folds (zero training spend) yet can't simply be dropped from the final model. Solution: **time-slice CV** where, per fold, channels/controls with **zero training-set spend are excluded from both train and test**, while the **full model trained on all data retains** late-entry channels. (Contributed back to open source.)

**The fundamental tool: rolling-origin CV with a dynamic feature mask [supplied].**
- Expanding/rolling train window; each fold's design matrix is masked to only-active features, so evaluation isn't contaminated by all-zero columns.
- Separates the two jobs cleanly: **CV** measures generalization on features that existed; the **production fit** uses every feature including late entrants.

---

## 11. Bayesian inference & workflow scaffolding

**What the posts apply / imply.**
- **PyMC-Marketing / PyMC** as the stack; full-posterior Bayesian inference (specific MCMC not named — practically NUTS/HMC).
- **`mu_effects`** — a *pluggable-component* framework: Funnel class, campaign multipliers, influencer effects, HSGP intercept each built and tested **independently**, then composed ("modularity over monolith").
- **YAML-configured models** — every prior/assumption inspectable; supports reproducibility and CI/CD ownership by the client team.
- **Convergence engineering** — required because the two-likelihood geometry is hard; diagnostics implied (R̂, ESS, divergences), plus **posterior-predictive checks** and **component-wise validation** before composing.
- **Hierarchical priors** — a planned extension to share strength across comparable channels with sparse data.

**Fundamental tools to teach [supplied]:** MCMC/NUTS basics; convergence diagnostics `R̂ ≈ 1`, effective sample size, divergent transitions; posterior predictive checks; prior predictive checks and prior elicitation; hierarchical/partial-pooling priors; non-centered reparameterization for hard geometries.

---

## 12. Dependency graph for teaching (suggested order)

```
1. Bayesian regression + MCMC/NUTS + diagnostics        (prereq)
2. MMM basics: adstock (§2) + saturation (§2)           (prereq)
3. Causal DAGs + mediation (direct vs indirect) (§1)    ← the frame
4. Censored / Tobit likelihood (§3)                     ← core innovation
5. Joint / two-likelihood models + identifiability (§4) 
6. Time-varying effects: GP → HSGP (§5)
7. Reparameterization for identification: Dirichlet shares (§6)
8. Latent-variable measurement models (§7)
9. do-calculus / interventions on a fitted graph (§8)
10. Constrained optimization + KKT/marginal-ROAS (§9)
11. Rolling-origin CV with dynamic features (§10)
12. Bayesian workflow: modular components, PPCs, YAML   (wraps it all)
```

**The through-line to emphasize:** items 4–10 are each an answer to one honesty requirement of the DAG in item 3 —
- *"observed spend ≠ true demand when caps bind"* → censoring (§3),
- *"daily data should inform upstream effects"* → joint likelihood (§4),
- *"the baseline drifts"* → HSGP (§5),
- *"effectiveness jumps at campaign changes"* → Dirichlet multipliers (§6),
- *"reach ≠ quality"* → latent decomposition (§7),
- *"I want counterfactual allocations"* → `do` operator (§8) + constrained optimization (§9),
- *"channels come and go"* → time-slice CV (§10).

---

## 13. Honest gaps in the source material

The three **blog** posts are conceptual; they **omit** the following. Most are now filled by the official docs notebooks (**Appendix A**); the rest remain genuine gaps:
- ~~No explicit adstock/saturation functional forms~~ → **resolved [docs]**: Geometric Adstock + Michaelis–Menten, forms in §A.3.
- ~~No written censoring distribution~~ → **resolved [docs]**: censored **Normal**, `lower=0`, via `pymc_extras.prior.Censored` (§A.4). *Note:* the docs censor at zero (no negatives); the **blog's** censoring is at the **budget cap** (right-censoring) — same mechanism, different threshold. Both matter; teach zero-censoring first, cap-censoring second.
- **HSGP** kernel/lengthscale/basis-count: still unspecified (blog-only feature, not in these two notebooks).
- **Optimization solver / `BudgetOptimizer` / `do` operator internals**: still not written out numerically (blog-only; the docs show *manual* g-computation counterfactuals instead — §A.6).
- **Dirichlet campaign multipliers, influencer reach×quality, time-slice CV, two-likelihood spend-to-spend**: blog-only extensions, not in these base notebooks.
- No posterior summaries, R̂/ESS values, or sensitivity analyses in any source (the advanced notebook does note the all-in model **diverges** during sampling — a diagnostic signal in itself).
- The blog's **Part II** appears unpublished; **Appendix A** (official docs) is the best available substitute for its technical content.

---

## Appendix A — Concrete reference implementation (from PyMC-Marketing official docs)

*Everything here is **[docs]** — quoted or closely paraphrased from the two upper-funnel notebooks named at the top. This is the runnable, equation-level version of §§1–3 and §8. It uses a clean synthetic 4-channel funnel, which is the ideal teaching example.*

### A.1 The canonical DAG (concrete, 4 channels)

```
        U1        U2      U3       U4
        │          │        │        │
        ▼          ▼        ▼        ▼
holidays──► X1 ──► X2 ──┐
              └──► X3 ──┴──► X4 ──► Y (new users)
                                    ▲
              holidays, exogenous ──┘
```
- **X1** = upper-funnel awareness (video/display). **X2, X3** = mid-funnel (social, remarketing) — the **mediators**. **X4** = lower-funnel brand search — the **sole proximal driver** of Y.
- **Critical structural fact:** `X1 → (X2,X3) → X4 → Y`. **X1 has *no* direct edge to Y.** Its entire effect is mediated. This is what standard "throw-everything-in-one-regression" MMM gets wrong.
- The DAG comes from **domain knowledge**, not the data (causal discovery is out of scope).

### A.2 Structural causal equations (the data-generating process) **[docs]**

Additive at every node; nonlinearity *only* at the final X4→Y map:
```
X₁,t = μ₁,t + β₀·Eₜ                      (E = event signal)
X₂,t = μ₂,t + β₁₂·X₁,t
X₃,t = μ₃,t + β₁₃·X₁,t
X₄,t = μ₄,t + β₂₄·X₂,t + β₃₄·X₃,t
Yₜ   = f(X₄,t; θ) + Trendₜ + Eventsₜ + εₜ
```
`μ_{j,t}` are exogenous (random-walk) components; `f` = adstock∘saturation. In the sim: `β₁₂=.02, β₁₃=.03, β₂₄=.04, β₃₄=.05`, so the **true total X1→X4 coefficient** is `β₁₂·β₂₄ + β₁₃·β₃₄` — the sum of path products (this is the mediation-formula check).

### A.3 The two media transforms, with exact forms **[docs]**

- **Geometric adstock**, `l_max = 24` (up to 24 days carryover), decay `α`:
  `x*_t = Σ_{l=0}^{24} α^l x_{t-l}` (optionally normalized). Prior: **`α ~ Beta(1,1)`** (uniform on [0,1]).
- **Michaelis–Menten saturation** (this is the Hill curve with shape `s=1`):
  `S(x) = α · x / (λ + x)` — `α` = max effect, `λ` = half-saturation point (`S(λ)=α/2`).
  Priors: **`λ ~ Gamma(μ=2, σ=2)`**, **`α ~ Gamma(μ=1, σ=1)`**.
- API: `GeometricAdstock(l_max=24, priors=...)`, `MichaelisMentenSaturation(priors=...)`, composed inside `MMM(...)`.

### A.4 Censored likelihood (concrete) **[docs]**

```python
model_config = {"likelihood": Censored(
    Prior("Normal", sigma=Prior("HalfNormal", sigma=1), dims="date"),
    lower=0)}
```
- Here censoring is at **zero** (outcome can't be negative) via `pymc_extras.prior.Censored`. Same machinery the blog points at the **budget cap** for right-censored demand — teach this as the general tool, then move the threshold from `0` to `cap_t`.

### A.5 The two-block estimation strategy **[docs]** — the heart of the method

Because X4 is a **mediator**, the backdoor criterion says: to get X1's total effect, do **not** put X1 and X4 in one regression. Instead fit two models:

**Block 1 — Outcome model (X4 → Y).** Full MMM with adstock+saturation+censored Normal, `channel_columns=["impressions_x4"]`, controls = events + trend. X4 alone "blocks all backdoor paths for the rest of the variables," so this cleanly recovers `f(X4)`.

**Block 2 — Mediator model (X1 → X4).** A *linear, instantaneous* model on **first differences**:
```
ΔX₄,t = α + β·ΔX₁,t + ηₜ
```
`NoAdstock(l_max=1)`, `NoSaturation(priors={"beta": Gamma(μ=0.7, σ=0.4)})`, target = `impressions_x4_diff`, feature = `impressions_x1_diff`.
- **Why first differences:** X1, X4 are near-`I(1)` (random walks). Levels-on-levels `X₄=βX₁+ε` is **spurious regression** (drifting `β̂`, invalid inference) unless cointegrated. Differencing → stationary `I(0)` relation → OLS consistently estimates the **total instantaneous effect** `β`. (`α≠0` because level drifts leave a nonzero mean in differences.)
- **Why no adstock/saturation here:** the modeling choice is "X1's effect on X4 doesn't itself carry over or saturate — all the dynamics live in X4→Y." Keeps the mediator link interpretable.

### A.6 Why the naive "all-in" model fails (the teachable derivation) **[docs]**

Fit `Y = α + β̂₁X₁ + β̂₄X₄ + ε`. The normal equations force orthogonal residuals:
```
E[X₁(Y − α − β̂₁X₁ − β̂₄X₄)] = 0
E[X₄(Y − α − β̂₁X₁ − β̂₄X₄)] = 0
```
These depend on **Cov(X₁,X₄)**. Since X1→X4 makes them collinear, and Y depends **nonlinearly** on X4 (so `f(X4)` isn't fully captured by a linear `β̂₄X₄`), the unmodeled nonlinear residual gets **redistributed onto the correlated X1**:
```
β̂₁ ∝ Cov(X₁, Y | X₄)  >  0     even though the true direct X₁→Y effect is ZERO.
```
Result (confirmed in the notebook): the all-in model **over-credits** upper funnel *and* **diverges** during sampling. Conditioning on a mediator **blocks** the very path you want to measure — the classic "don't control for a mediator" lesson.

### A.7 Backdoor criterion (formal statement) **[docs]**

A set `Z` satisfies the backdoor criterion for `X→Y` if: (1) no node in `Z` is a **descendant** of `X`; and (2) `Z` blocks every path into `X` (backdoor/confounding paths). Then
```
E[Y | do(X=x)] = ∫ E[Y | X=x, Z=z] dP(z).
```
Here: **adjust for** trend/events (common causes of X1 and Y); **never adjust for** X2/X3/X4 (they're descendants/mediators of X1).

### A.8 The counterfactual: g-computation / do-operator **[docs]**

To answer "what if X1 = 0?", propagate the intervention through the fitted chain (this *is* the `do` operator, done by hand):
1. Set `impressions_x1 = 0`, recompute `impressions_x1_diff`.
2. **Block 2** predicts counterfactual X4 under `do(X1=0)`: `X4_cf = X4 − (contribution of X1 to X4)`.
3. Feed `X4_cf` into **Block 1** to predict `Y | do(X1=0)`.
4. **Causal effect** `= Y|do(X1=1) − Y|do(X1=0)`, with full posterior uncertainty (HDI).
This is the **mediation formula / g-computation**: intervene, re-simulate each edge, difference the outcomes. The notebook shows the recovered effect matches the simulated ground truth in- and out-of-sample.

### A.9 Inference & API cheat-sheet **[docs]**

- Sampler: NUTS via **`nuts_sampler="numpyro"`** (JAX); `tune≈1000–1500, draws≈500–1000, chains=4, target_accept 0.84–0.9`.
- Core classes: `MMM(date_column, target_column, channel_columns, control_columns, adstock, saturation, model_config, sampler_config)`; transforms `GeometricAdstock / MichaelisMentenSaturation / NoAdstock / NoSaturation`; priors via `pymc_extras.prior.Prior` and `Censored`.
- `.build_model()` → `.add_original_scale_contribution_variable()` → `.fit()` → `.sample_posterior_predictive()`; diagnostics/plots via **ArviZ** (`az.plot_hdi`, `az.plot_posterior`), DGP built symbolically in **PyTensor** (`pytensor.xtensor`).

### A.10 Minimal teaching arc using this example

1. Draw the DAG (§A.1) → assert X1 has no direct Y edge.
2. Fit the **all-in** model → show it over-credits X1 *and* diverges (§A.6).
3. Derive *why* via normal equations + Cov(X1,X4) (§A.6) and the backdoor criterion (§A.7).
4. Fix it with **two blocks** (§A.5): outcome X4→Y, mediator ΔX1→ΔX4.
5. Recover X1's causal effect by **g-computation** `do(X1=0)` (§A.8), show it matches truth.
6. *Then* layer the blog's production concerns on top: cap-censoring (§3), two-likelihood spend-to-spend (§4), HSGP baselines (§5), and joint budget optimization (§9).

---

## Appendix B — The machinery behind the blog-only extensions (from PyMC-Marketing docs)

*There is **no influencer notebook** in PyMC-Marketing — the reach × latent-quality influencer model (§7) is a **bespoke, unreleased** component built for the Nürnberger engagement. But it is assembled from three documented pieces below: the **custom-component framework** (how you plug in any new effect), **HSGP time-varying parameters** (§5's tool, and the natural home for latent quality drifting over time), and **time-slice CV** (§10 / blog Part III). All **[docs]**. Notebooks pulled: `mmm_components`, `mmm_tvp_example`, `mmm_time_slice_cross_validation`.*

### B.1 Custom components — how a bespoke effect (e.g. influencer) is built

The `MMM` class is assembled from **swappable components**, and you can compose your own inside a raw `pm.Model`. This is the mechanism the blog's `mu_effects` framework uses to add influencer/campaign effects.

- Every transform is an object with **configurable priors** and an **`.apply(x)`** method that (a) creates the parameter RVs and (b) applies the transform. Example: `MichaelisMentenSaturation(priors={...})`, then `saturation.apply(data)` inside a model context. Introspect with `.function_priors`, `.default_priors`; preview with `.sample_prior()` → `.sample_curve()` → `.plot_curve()`.
- **Priors are where the modeling lives.** "Change the priors, change the model." A parameter can be **pooled** (common across channels), **independent** (`dims="channel"`), or **hierarchical** (a prior whose hyperparameter is itself a `Prior` with its own dims):
  ```python
  hierarchical_lam = Prior("HalfNormal", sigma=Prior("HalfNormal", sigma=1), dims="channel")
  # multi-level: alpha hierarchical across channels AND geos
  hierarchical_alpha = Prior("Gamma",
      mu=Prior("HalfNormal", sigma=1, dims="geo"),
      sigma=Prior("HalfNormal", sigma=1, dims="geo"),
      dims=("channel", "geo"))
  ```
  This is the exact tool for the blog's *"hierarchical priors for comparable channels with sparse data"* — partial pooling shares strength.
- **How the influencer effect maps onto this [supplied bridge]:** define a custom component whose contribution is `r(followers)·q_k` (§7) — `r` a saturation object on follower count, `q_k` a `Prior("LogNormal", dims="influencer")` latent quality — and add its `.apply(...)` output into the model's `mu` alongside `media_contributions`. Broadcasting rule (stated in the notebook): the parameter `dims` just need to be **broadcastable** with the data's dims.
- The generic assembly pattern (verbatim shape): `mu = intercept + media_contributions + fourier_trend`, then `Normal("target", mu=mu, sigma=sigma, observed=...)`. Events can be added as **Gaussian-bump basis** functions (`GaussianBasis`, `HalfGaussianBasis`, `AsymmetricGaussianBasis` via `EventEffect`) — the same Gaussian-pulse events used in the DGP of Appendix A.
- Seasonality: `YearlyFourier(n_order=...)` / `MonthlyFourier` (Fourier basis), also prior-configurable and hierarchical.

### B.2 Time-varying parameters via HSGP (concrete — fills §5)

Turn on with a single flag: **`MMM(..., time_varying_intercept=True)`**.

- **Structural form [docs]:** the time-varying intercept is a **Hilbert-Space GP constrained to mean μ=1**, then **multiplied by a scalar baseline intercept**. So the GP models the **percentage deviation** from baseline over time: if baseline = 1000 units/week, the GP is the multiplicative wiggle around it. (Same idea extends to `time_varying_media` for time-varying channel effectiveness.)
- **Configuration object:** `HSGPKwargs(m=500, L=188, eta_lam=5.0, ls_mu=5.0, ls_sigma=10.0, cov_func=None)`:
  - **`m`** = number of basis functions (default **200**) — more = higher-frequency fidelity, slower sampling.
  - **`L`** = half-width of the GP domain box `[−L, L]`; rule of thumb in the notebook `L ≈ (1+c)·(n_dates/2)` with a small padding factor `c≈0.2–0.3`.
  - **`ls_mu`, `ls_sigma`** = **length-scale** prior mean/sd (in time units) — the single most important knob. Default mean ≈ 2 years; **shorten it to catch short-timescale events** (the notebook refits with `ls_mu=52.18` weeks = 1 year to capture event spikes).
  - **`eta_lam`** = amplitude prior; **`cov_func`** = optional custom kernel, and can be the **sum of two kernels** to capture two distinct timescales simultaneously.
- **When to use it (the notebook's verdict — teach this honestly):**
  - Seasonality → *works, but a **Fourier basis** is better/cheaper.*
  - Steady upward trend → *works poorly; use a **linear control** instead.*
  - **Unexpected events / unexplained variance no other term captures → excellent** — this is the real use case. It's a **catch-all** for temporal variance.
  - GPs **extrapolate poorly**: out-of-sample the intercept reverts to prior mean and uncertainty balloons (fine, since MMM's job is *in-sample* contribution decomposition).
- Prior: baseline intercept forced positive, e.g. `Prior("Normal", mu=0, sigma=0.1, transform="sigmoid")`. Sampler here: **`nutpie`** with JAX backend.

### B.3 Time-slice cross-validation (concrete — fills §10, blog Part III)

API: **`TimeSliceCrossValidator(n_init, forecast_horizon, date_column, step_size)`**, then `.run(X, y, sampler_config=..., yaml_path=...)`.

- **Expanding-window** scheme (simulates production: training set grows over time). Number of folds:
  ```
  n_splits = y.size − n_init − forecast_horizon + 1
  ```
  (`n_init` = initial train size, e.g. 163; `forecast_horizon` = e.g. 12 weeks = 3 months; `step_size` grows the window each fold.)
- **Data-leakage discipline (emphasized):** no future data may touch training — *including preprocessing*. Cost-share features must be recomputed **per fold**; global/trend features must be leakage-safe (they use a simple incrementing `t`, +1 per slice, rather than a globally-fit trend).
- **What it evaluates — two things, not just accuracy:**
  1. **Parameter stability** over folds via `cv.plot.param_stability(var_names=["adstock_alpha", "saturation_beta", "saturation_lam"])` — good models don't have parameters (hence ROAS) that lurch fold-to-fold.
  2. **Probabilistic out-of-sample accuracy** via **CRPS** (Continuous Ranked Probability Score) — `cv.plot.crps()` — the proper scoring rule for a model that predicts a *distribution*, not a point. Diverging train-vs-test CRPS trends flag **overfitting**.
- **Diagnostics:** `results` is an ArviZ `InferenceData`; check `results.sample_stats["diverging"].sum()` for divergences.
- **Connection to the blog's late-entry-channel problem (§10):** this is the harness into which the blog's dynamic-feature masking plugs — per fold, channels with zero training-set spend are excluded from that fold's model, while the full production fit keeps them.

### B.4 What is still genuinely unavailable in any public source

Even after the docs, these remain **blog-only, no open code/equations**: the **spend-to-spend two-likelihood** joint fit (§4), **cap-level (right-)censoring** of demand (the docs only censor at zero, §A.4), **Dirichlet campaign-change multipliers** (§6), the **reach×quality influencer component** itself (§7 — only its ingredients are documented, §B.1), and the **`BudgetOptimizer` + `do`-operator joint budget/cap optimization** numerics (§9 — the docs show manual g-computation instead, §A.8). These are the genuine "read Part II / engage PyMC-Labs" gaps.
