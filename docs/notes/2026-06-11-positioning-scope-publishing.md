# Strategy Notes — Positioning, Scope & Publishing

**Date:** 2026-06-11
**Context:** Conversation reviewing the book (read only from the published site
<https://jlhf80.github.io/mmm_foundational_text/>) and working through identity, competitive
positioning, scope, and publishing/revenue strategy. Captured here for reference.

---

## What the book is

*The Mathematics and Engineering of Marketing Mix Modeling* — a comprehensive textbook bridging
quantitative mathematics and software engineering to teach modern Marketing Mix Modeling (MMM),
"from standard quantitative prerequisites all the way to advanced practice, without hand-waving."

- **Audience:** readers with Calc 2, linear algebra, calculus-based statistics. Rigor escalates
  undergraduate → advanced graduate.
- **Structure (7 parts):** I Mathematical Foundations · II Regression & Bayesian Inference ·
  III Computation/Sampling (MCMC, HMC) · IV Optimization · V MMM Synthesis (adstock, saturation,
  fitting, budget optimization) · VI Software Engineering & CS · VII Capstone (end-to-end on
  synthetic data).
- **Per-chapter format:** intuition → theory → worked example → runnable code (PyMC/SciPy) →
  tiered exercises → summary. Pragmatic on proofs: "Where a result is genuinely illuminating we
  prove it; where a proof is mechanical we cite it."

## Identity decision: keep it theory-forward

Confirmed direction — **theory-first is the moat.** The math-first spine differentiates it from
the many applied "run PyMC-Marketing" tutorials. Implications:
- The "without hand-waving" promise is the whole bet — proofs/derivations must be genuinely
  rigorous and correct (a hand-wave in a theory book is the flaw readers catch).
- Code tie-ins stay subordinate: demonstrations of the theory, not standalone how-tos.
- Keep the prerequisite bar where it is; don't soften it (dilutes the rigor).

## Competitive landscape

**MMM-specific:** the market is thin and almost empty where this book sits. Existing material is
either *theory-free & practical* (vendor/CMO e-books — Think with Google, Forvio, Baxter's
*Handbook*; blog tutorials — 1749.io, Ryan O'Sullivan, ELIYA; library docs — PyMC-Marketing,
Robyn, Google Meridian) or *practice-free & narrow* (papers — Jin et al. 2017; 2024–25 arXiv on
time-varying effects, NNN, DeepCausalMMM). **No rigorous, ground-up MMM textbook exists.**

**Spiritual analogs (the real benchmark):** rigorous Bayesian textbooks with code —
- Gelman et al., *Bayesian Data Analysis (BDA3)* — rigor gold standard, general, code an afterthought.
- McElreath, *Statistical Rethinking* — pedagogy benchmark (intuition-before-formalism, single arc);
  this book's chapter structure mirrors it.
- Martin/Kumar/Lao, *Bayesian Modeling and Computation in Python* — closest technically (PyMC-native).

**Positioning statement:** *"Statistical Rethinking's pedagogy + BDA's rigor, narrowed to MMM, and
extended through optimization and software engineering."* No direct competitor. Differentiators:
(1) domain-specific end-to-end vertical slice; (2) optimization + SWE as first-class material (no
analog in the Bayesian canon); (3) math taught for the sake of building the system.
**Watch-item:** readers who know the canon will hold Parts I–III to the BDA/Rethinking standard —
those foundations must be unimpeachable.

## Scope decision: DLM + QE calibration → advanced track inside Book 1

**Add both topics** (Bayesian DLM / time-varying coefficients, and quasi-experimental design for
calibration):
- **DLM** — frontier topic (time-varying effects; cf. "Your MMM is Broken" 2024), natural extension
  of the regression → state-space ladder, aligns with `model/dlm_stub.py` in the backend. Scope as
  *"time-varying MMM"* (DLM as bridge), not a full state-space text. Owes a real Kalman-filter
  derivation.
- **QE calibration** — close to a *credibility requirement*. Addresses MMM's core weakness
  (identification/causality); calibrating with experiments (geo lift, RDD, ITS, diff-in-diff) is
  what separates credible MMM from curve-fitting. Forces causal-inference foundations (potential
  outcomes). Aligns with backend `experiments/` + `calibration.py`.

**Placement:** an explicitly-marked **"Advanced / Extensions" track** the core arc does NOT depend
on. DLM chapter in Part V (state-space seeded in Part III); QE-calibration chapter near the capstone.

**Do NOT make a separate Volume II now.** DLM + QE don't form a coherent standalone book (different
themes; would need to re-lay or hard-depend on Book 1's foundations); forking before Book 1 ships
risks shipping neither; the DLM/QE reader is the same reader as Book 1's. Revisit Volume II only
**post-ship**, and only if (a) Book 1 has shipped with a proven audience, (b) advanced material has
grown past ~120–150 pp and unbalances the book, and (c) it has broadened into a coherent theme
(e.g. "Causal & Dynamic MMM": + hierarchical/geo, full causal inference, nonstationary effects).

## Publishing (made with Claude)

AI assistance is not a blocker. Three things to handle:
1. **Own the correctness** — highest-value step for a rigor-branded book. Foundational/derivation
   chapters need genuine human/expert verification before publishing. (General info, not legal advice.)
2. **Copyright** — US Copyright Office guidance (2023–25): pure AI output isn't copyrightable, but
   human selection/arrangement/editing/original writing is. Register claiming human authorship;
   disclaim AI-assisted portions.
3. **Disclosure** — KDP asks at upload (answer "AI-assisted" given heavy human authorship);
   traditional publishers increasingly add AI clauses/warranties — read carefully.

**Paths:** (1) stay free + open (already a Quarto site on GitHub Pages — reach/reputation engine);
(2) self-publish PDF/EPUB via Quarto → Leanpub/KDP (revenue, full rights); (3) traditional/academic
publisher (CRC, Springer, Manning, No Starch — credibility, but AI clauses + lower royalty).

**Recommended:** finish + verify the math → keep free web edition as funnel → self-publish polished
PDF on Leanpub + KDP paperback. Decide EPUB/Kindle only after checking worst equation page on-device.

### Self-publish mechanics (Quarto → Leanpub/KDP)

- **Math + EPUB don't mix well.** PDF (LaTeX) renders math beautifully → lead with PDF. EPUB3 MathML
  support is inconsistent on e-readers; don't ship an EPUB whose worst equation page you haven't
  inspected on a real device.
- **`_quarto.yml`:** add `pdf` (documentclass `scrbook`, 6×9 geometry with larger inner/gutter
  margin, embed fonts, `keep-tex`) and optionally `epub` (cover image). Keep existing `html`.
- **KDP print specs:** 6×9 trim; gutter scales with page count (KDP margin table); embed all fonts
  (`pdffonts` → emb yes); no bleed unless images to edge; black text only.
- **Render:** `quarto render --to pdf`; verify with `pdffonts`/`pdfinfo` (page size 432×648 pt).
- **Cover:** KDP needs full wraparound PDF (spine width from KDP Cover Calculator); Leanpub/ebook
  needs front image only (~1600×2560).
- **ISBN:** KDP gives a free one (KDP = publisher of record); buy own only to be named publisher
  (separate ISBNs for print vs ebook).
- **Leanpub:** use "Bring Your Own Book" (BYOB) to upload the Quarto PDF; variable pricing; free
  updates to buyers; ~80% royalty.
- **KDP:** Create → Paperback; upload interior + wraparound cover PDF; use online previewer; answer
  AI-disclosure as "AI-assisted."

### Leanpub revenue expectations

- **Royalty ~80%** (vs KDP ebook 35–70%, traditional 10–15%). At $39 list ≈ $31/copy. Sensible price
  band for a rigorous technical book: **$29–$49**; variable pricing often lifts the average paid.
- **Scenarios @ $39:** no audience 30–80 copies (~$1–2.5k); modest niche 200–400 (~$6–12k); strong
  launch 800–1,500 (~$25–47k); outlier 3,000+ (~$90k+).
- **The binding constraint is audience, not the royalty rate.** Leanpub has near-zero organic
  discovery — revenue ≈ reach × conversion. The free web edition is both funnel and competition
  (need a paid reason: polished PDF, print, support-the-author, gated extras like solutions/datasets).
- **Highest-leverage move:** build the email list *before* launch, using the free edition as the
  magnet. Realistic without pre-built audience: low four figures year one; five figures achievable
  with real distribution (MMM is a narrow but monied professional niche).

---

## Open next steps (not yet done)

- [ ] Set up `_quarto.yml` PDF/EPUB formats; get a clean 6×9 render + `pdffonts`/`pdfinfo` checks.
- [ ] Brainstorm chapter design for the DLM and QE-calibration advanced chapters (scope, prereq
      dependencies, exact placement in the 7-part structure).
- [ ] Independent/expert verification pass on Parts I–III derivations before any paid publication.
- [ ] Pre-launch funnel/email-list plan using the free site.
