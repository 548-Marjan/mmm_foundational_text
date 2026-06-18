# Book TODO

Editorial review of *The Mathematics and Engineering of Marketing Mix Modeling*,
captured as an actionable backlog. Reviewed from three perspectives: a PhD-level
editor (content, pedagogy, proofs), a would-be reader, and a quantitative-marketing
professional.

**Headline finding (updated):** chapters 1–22 are now fully written — Parts I–IV (the
prerequisite mathematics), Part V (the MMM itself: DGP, building & fitting, DLM /
state-space, budget optimization), Part VI (causal grounding & calibration: causal
foundations, quasi-experimental designs, advanced calibration, the prior store), and
the first chapter of Part VII (CS foundations). The foundations' proofs held up under
line-by-line scrutiny; the Part V–VI chapters still want the verification pass in P0b.
What remains unwritten is **chapters 23–26**: Part VII's software-architecture,
data-engineering, and testing/reliability chapters, and the Part VIII capstone — all
still ~31–33-line stubs.

---

## P0 — Remaining unwritten chapters

- [x] **Write Part V (Marketing Mix Modeling: Modeling), ch. 14–17.** Done — DGP /
      adstock / saturation, building & fitting, DLM / state-space, budget optimization.
      The payload the preface promises ("build a full Bayesian MMM end-to-end").
- [x] **Write Part VI (Causal Grounding & Calibration), ch. 18–21.** Done — the
      *differentiator*: potential-outcomes / do-calculus foundations, quasi-experimental
      designs, advanced calibration, and the prior store / calibration–optimization
      loop. The **identification** mandate is met (the chapters confront the
      observational-spend confound head-on). The potential-outcomes foundation and the
      Kalman-filter derivation still want the foundations'-standard verification pass —
      tracked under P0b, the credibility requirement your notes flag.
- [ ] **Write the rest of Part VII (SWE & CS), ch. 23–25, and Part VIII (Capstone),
      ch. 26.** Ch. 22 (CS foundations) is written; **ch. 23 (software architecture),
      ch. 24 (data engineering), ch. 25 (testing & reliability)** remain ~33-line stubs,
      as does the **ch. 26 capstone**. This is the payload that remains.
- [ ] The preface's end-to-end promise is now substantially met through Part VI;
      before claiming the *complete* book (and resolving the `[link TBD]` CTA), finish
      the Part VII SWE chapters and the capstone.

## P0a — Authoring note: Chapter 16 (DLM / State-Space MMM) — DLM vs. BTVC arc

**Status:** ch. 16 (`parts/05-mmm-modeling/03-dlm.qmd`) is now written (~528 lines).
The checklist below was the authoring spec; reconcile each item against the shipped
chapter rather than treating it as open work — most or all has likely been addressed.

Captured from a review session comparing the foundational MMM literature
(Jin et al. 2017), **BTVC** (Ng, Wang & Dai, Uber 2021,
[arXiv 2106.03322](https://arxiv.org/pdf/2106.03322)), and Bayes R² (Gelman et al.
2018). Use this as the spine for ch. 16's time-varying-coefficient section —
present DLM/Kalman as the principled baseline, then BTVC as the production-driven
reframing. Follows the chapter template (intuition → theory/proofs → tradeoffs →
code tie-in).

- [ ] **Frame the shared problem.** After a log-log transform the MMM is additive
      with time-varying coefficients:
      `ln(ŷ_t) = l_t + s_t + Σ_p ln(x_{t,p}) β_{t,p}` (trend, seasonality,
      channel coefficients `β_{t,p}`). The "natural" fit is a state-space model —
      DLM or Kalman filter. Derive the Kalman filter at the foundations' rigor
      standard (this is also a credibility requirement flagged under P0/Part VI).

- [ ] **Why BTVC argues DLMs break down (teach all three points).**
  - *Caveat 1 — DLM + MCMC is too costly at MMM scale.* A DLM models the
    coefficient as a latent state evolving every time step (random walk), so
    inference samples a state for **every regressor × every time step** — a
    high-dimensional problem under MMM's "small *n*, large *p*" reality.
  - *Caveat 2 — the Kalman filter is fast but too rigid.* Closed-form updates buy
    speed but forbid the customizations MMM needs: **non-negative** spend
    coefficients (Gaussian conjugacy can't impose sign constraints) and a
    **t-distributed / robust** noise process for outlier resistance.
  - *Deeper structural reason (the unifying point).* A DLM **couples the parameter
    count to the number of time steps** — the recursive `β_t | β_{t-1}` state
    forces a sequential, `T`-dimensional inference problem. That creates the
    dilemma: pay for MCMC over a huge state space (caveat 1) or buy speed with
    Kalman's restrictive assumptions (caveat 2). You can't get scalable *and*
    flexible in the state-space framing.

- [ ] **BTVC's solution (kernel-smoothed latent knots).** Stop modeling the
      coefficient as a sequential stochastic process; represent the coefficient
      *curve* as a smooth function of a small set of latent **knots** (GAM /
      kernel-regression view). For each regressor place `J ≪ T` knots `b_{j,p}` at
      times `t_j`:
      `β_{t,p} = Σ_j w_j(t) b_{j,p}`, `w_j(t) = k(t,t_j) / Σ_i k(t,t_i)`; matrix
      form `β = K b`. Gaussian kernel for regression coefficients
      (`k(t,t_j;ρ) = exp(−(t−t_j)²/2ρ²)`, bandwidth `ρ` sets smoothness);
      triangular-like kernel for trend/seasonality.
  - *Bayesian hierarchy.* Laplace prior on adjacent trend/seasonality knots
    (change points, à la Prophet). Two-layer **folded-normal** hierarchy on
    regression knots — `μ_reg ~ N⁺(μ_pool, σ²_pool)`, `b_reg ~ N⁺(μ_reg, σ²_reg)`
    — giving **positivity** for free and **shrinkage toward the channel's grand
    mean over sparse/zero-spend periods**. Experiment lift-test results ingestible
    as channel priors. Posteriors via **Stochastic Variational Inference (SVI)**.
  - *Show the fix mapping explicitly:* params decoupled from `T` (`J` knots, kernel
    interpolates) → solves caveat 1; static hierarchical regression over knots +
    SVI → solves the MCMC cost; folded-normal priors → sign constraints Kalman
    couldn't do; free Bayesian likelihood → robust/custom noise + ingested priors.

- [ ] **Land the one-liner contrast.** DLM asks *"how does the coefficient drift
      step to step?"* (sequential, recursive, `T`-dimensional); BTVC asks *"what
      smooth curve through a few latent knots best explains the data?"* (smoothing,
      basis-function, `J`-dimensional). The smoothness a DLM gets from the
      random-walk transition, BTVC gets from kernel bandwidth `ρ`; the
      regularization a DLM gets from state-noise variance, BTVC gets from
      hierarchical pooling. Flag this as the **Advanced/Extensions track** (ties to
      the P1 "mark the Advanced track" item) — core arc shouldn't depend on it.

## P0b — Pre-publication verification (AI-assisted authorship) & the free-publish decision

This book was authored with heavy AI assistance (spec → plan → subagent cycle). That
is a legitimate way to write — synthesis, architecture, and exposition are real
authorial work — but it creates specific obligations that must be cleared before the
book goes public under an author's name. The method is fine *if the author has
actually verified it and can defend any page when an expert pushes back*; it is not
fine as unread output. The items below are the gate.

- [ ] **Spot-check the keystone proofs by hand.** Prioritize the load-bearing
      results that carry the book's claim to rigor: the omitted-variable-bias theorem
      and back-door adjustment (ch. 18/19), the DiD-identifies-ATT-under-parallel-trends
      and IV/Wald-is-a-secant proofs (ch. 19/20), the calibration ridge-compression
      result (ch. 20), and the monotone-EVPI-decay convergence result (ch. 21). AI math
      is fluently wrong in hard-to-spot ways; a textbook's implicit promise is that
      these hold.
- [ ] **Verify every citation resolves to a real source making the claimed point.**
      Audit `references.bib` and every "Canonical anchors" line against the in-text use.
      AI invents or misattributes references; one fabricated citation in a
      rigor-branded book torches credibility with the exact expert audience. (Overlaps
      the P1 Perron–Frobenius and Hill cross-reference fixes — fold those in.)
- [ ] **Render the full build and run every code tie-in; confirm the numbers match.**
      The prose asserts specific figures (e.g. the 3.5-vs-2.0 confounding result, the
      EVPI dropping ≈0.12 → ≈0.03). Execute the cells and check the asserts actually
      pass and produce the quoted values.
- [ ] **Decide on an AI-assistance disclosure line and resolve author identity.** A
      one-line note ("drafted with AI assistance, reviewed and verified by the author")
      raises credibility rather than lowering it and protects against a reader inferring
      it and feeling misled. Resolve the `"548 / Marjan"` author placeholder (see P1).
- [ ] **Publish-format decision: lean free + open on the author's site.** The book's
      real payoff is reputation/discovery in a niche where public material is shallow,
      not direct sales — free HTML maximizes that, with an optional paid PDF/print or
      tip jar as a hedge. Free also lowers the stakes on residual imperfection. Gate
      this decision on the verification items above being cleared, not on format.

## P1 — Consistency / correctness fixes (mechanical, do now)

- [ ] **README is a structural version behind.** It describes a **7-part** layout
      (SWE = Part VI, Capstone = Part VII, no separate causal-calibration part) and a
      Status line claiming only ch. 1–2 are written. The live book is **8 parts**
      (causal calibration = Part VI; SWE = VII; Capstone = VIII) with **12** chapters
      written. Rewrite README structure table + Status to match `_quarto.yml`/`index.qmd`.
- [ ] **Hill-curve cross-reference error.** The convexity chapter (ch. 11) repeatedly
      cites "the Hill saturation curve introduced in Chapter 5." Chapter 5 is *Bayesian
      Inference* (ridge/MAP) — correctly cited elsewhere in the same chapter. The Hill /
      saturation curve actually belongs to the MMM DGP (ch. 14, a stub). Fix the
      reference; note the deeper issue that Part IV's central running example depends on
      a concept only *defined* in an unwritten chapter.
- [ ] **Stub contradicts the locked scope decision.** `parts/05-mmm-modeling/01-mmm-dgp.qmd`
      lists "*Canonical anchors: Jin et al. 2017; pymc-marketing.*" The positioning note
      and project memory record a deliberate decision to stay framework-agnostic and
      **never name PyMC-Marketing**. Scrub `pymc-marketing` from every stub anchor now.
- [ ] **Perron–Frobenius mis-cited.** Ch. 7 attributes stationary-distribution
      existence/uniqueness to `[@strang2016]`. Re-cite to **Norris (1997)**, already in
      the bibliography.
- [ ] **Placeholders.** `_quarto.yml` author is `"548 / Marjan"`; preface CTA is
      `[link TBD]`. Resolve before any paid publication.
- [ ] **Mark the Advanced/Extensions track.** The strategy note calls for DLM + QE
      calibration to be explicitly flagged as an advanced track the core arc does not
      depend on. No such marking exists in `_quarto.yml` yet.

## P2 — Proof / rigor polish (written chapters)

- [ ] **Markov convergence proof (ch. 7)** is given only for the diagonalizable case,
      with Jordan blocks waved to Norris. For a "without hand-waving" book, add the
      one-line reason the non-diagonalizable case still vanishes
      (`t^k |λ|^t → 0` for `|λ| < 1`) rather than only citing it.
- [ ] **Independent expert verification pass on Parts I–III derivations** before any
      paid publication — flagged in the positioning note as a pre-ship gate. Proofs
      checked in this review (normal equations, spectral theorem, PD, MLE=LS, LLN,
      detailed balance / convergence / ergodic, convex first/second-order, KKT
      sufficiency, equal-marginal-returns, Beta–Binomial & Gaussian conjugacy) all read
      as correct, but a rigor-branded book wants a second human signature.

## P3 — Pedagogy / reader experience

- [ ] **Prerequisite copy undersells the real bar.** "Who this book is for" leads with
      "Calculus 2, linear algebra, calculus-based statistics," then admits the rigor
      reaches "advanced master's / first-year PhD." The chapters confirm the latter
      (conjugate-transpose eigenvalue arguments, biorthogonal spectral decompositions,
      renewal–reward proofs appear early). Set expectations honestly up front.
- [ ] **Lighten the motivation sections.** Several (esp. Markov chains) are walls of
      text — one sentence runs ~120 words. Add a 2–3 sentence hook + a bulleted "by the
      end you can…" before the dense exposition. McElreath (your stated pedagogy
      benchmark) is far more whitespace-generous.

## What is already working (keep)

- Graduate-level, correct proofs across the written chapters — the "without
  hand-waving" bet is being met where content exists.
- The "rungs" ladder, per-chapter rhythm, genuinely hand-worked examples, four-tier
  (C/B/P/A) exercises, and the 2,263-line solutions appendix.
- Practitioner bridges done right: Markov-chain **removal-effect attribution** (and the
  correct note that removal effects need not sum to 1), the **equal-marginal-returns /
  water-filling** result with shadow-price reading, and the honest **S-curve
  non-convexity** reckoning that motivates multi-start SLSQP.
- Hierarchical/geo pooling already exists (ch. 6) — table stakes for credible MMM.
