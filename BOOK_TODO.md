# Book TODO

Editorial review of *The Mathematics and Engineering of Marketing Mix Modeling*,
captured as an actionable backlog. Reviewed from three perspectives: a PhD-level
editor (content, pedagogy, proofs), a would-be reader, and a quantitative-marketing
professional.

**Headline finding:** chapters 1–12 (Parts I–IV — the prerequisite mathematics) are
fully written and genuinely strong; the proofs hold up under line-by-line scrutiny.
Chapters 13–26 (Parts V–VIII — the MMM itself, causal calibration, software
engineering, and the capstone) are ~33-line stubs. Everything that makes this *an MMM
book* rather than an applied-math book is still scaffolding.

---

## P0 — The book's actual subject is unwritten

- [ ] **Write Part V (Marketing Mix Modeling: Modeling), ch. 14–17.** DGP / adstock /
      saturation, building & fitting, DLM / state-space, budget optimization. This is
      the payload the preface promises ("build a full Bayesian MMM end-to-end").
- [ ] **Write Part VI (Causal Grounding & Calibration), ch. 18–21.** This is the
      *differentiator* — the "calibration step most treatments skip" and the "where it
      quietly misleads you" promise. Must confront **identification** (you cannot pin a
      saturation curve / true incrementality from observational, collinear,
      seasonally-confounded spend). Ensure the potential-outcomes foundation (QE) and
      the Kalman-filter derivation (DLM) land at the foundations' rigor standard — your
      own notes call these a credibility requirement.
- [ ] **Write Part VII (SWE & CS), ch. 22–25 and Part VIII (Capstone), ch. 26.**
- [ ] Until the above exist, the preface's "*this preface is a free sample… the
      complete book is available at [link TBD]*" framing overpromises. Either gate the
      claim or finish the payload.

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
