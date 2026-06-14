# DLM vs. BTVC — Structural Decision

**Date:** 2026-06-14
**Context:** Conversation comparing time-varying-coefficient methodology across the
foundational MMM literature — Jin et al. 2017 (static + adstock/saturation), BTVC
(Ng, Wang & Dai, Uber 2021, [arXiv 2106.03322](https://arxiv.org/pdf/2106.03322)),
and Bayes R² (Gelman et al. 2018). Records the editorial call on whether to keep
**both** DLM and BTVC in the book, and how to position them. Companion to the
Chapter 16 authoring note (`BOOK_TODO.md`, section **P0a**) and the
[Part V/VI reorg plan](2026-06-12-part-v-vi-reorg-plan.md).

---

## Decision

**Keep both DLM and BTVC — but rank them.** DLM/state-space is the **core,
required baseline**; BTVC is **one Advanced/Extensions endpoint**, clearly flagged
and skippable. They are *not* two co-equal mandatory methods.

The single escalating spine for the dynamic-MMM material:

> **static (Jin et al.) → dynamic via DLM → dynamic via BTVC**

each step motivated by a concrete limitation of the step before it.

## Why keep both (the contrast *is* the lesson)

1. **Different transferable skills, not two recipes.** DLM teaches sequential
   state-space reasoning (recursive filtering, the random walk as structural prior,
   Gaussian conjugacy). BTVC teaches basis-function / kernel-smoothing reasoning
   (a function as a few latent knots; parameter count decoupled from series length;
   flexibility from priors, not from a transition equation). A reader who sees only
   one mistakes a *method* for *the* way to make coefficients vary.
2. **DLM is the pedagogically necessary baseline.** BTVC's whole motivation is a
   critique of DLM. You can't teach *why* BTVC exists without first making the
   reader feel the DLM dilemma. Cut DLM and BTVC becomes an unmotivated recipe —
   the "call this library function" failure mode the preface positions against.
3. **It is the book's thesis made concrete.** The preface promises "where it
   quietly misleads you." The textbook-canonical answer (state-space) has real
   failure modes under production constraints; a different framing resolves them.
   That escalation is the pitch, demonstrated.
4. **DLM is reusable infrastructure.** The Kalman derivation pays off again in the
   Part VI causal-calibration material and anywhere principled uncertainty
   propagation is needed. Not a dead-end detour.

## The real risk (why "rank," not just "include")

The danger is **scope creep and a muddied core path**, not redundancy. If both are
presented as co-equal full methods, we risk bloating Part V, implying the reader
must master both before shipping anything, and burying the single end-to-end MMM
the preface actually promises.

## Resolution / structure

- **Core path (everyone):** static-coefficient Bayesian MMM → DLM/state-space as
  the natural way to make coefficients move, *with its failure modes made
  explicit* (MCMC cost; Kalman rigidity on sign constraints and robust noise; the
  structural parameter-vs-`T` coupling). This alone is enough to build and ship a
  working dynamic MMM.
- **Advanced track (flagged, skippable):** BTVC as the production-driven reframing
  that resolves the DLM dilemma (kernel-smoothed latent knots `β = K b`, `J ≪ T`;
  folded-normal two-layer hierarchy for positivity + shrinkage; SVI). Mark it so
  the core arc does not depend on it.
- **Keep it a clean escalation, not a three-way pileup.** Jin et al. is the Part V
  backbone; DLM and BTVC each enter motivated by the prior step's limitation. Push
  BTVC-specific digressions (SVI internals, Orbit specifics) to a reference or
  appendix rather than the main text.

## Cross-references / follow-ups

- Implements the existing **P1** "mark the Advanced/Extensions track" item in
  `BOOK_TODO.md` for the DLM + QE-calibration material.
- Chapter 16 authoring detail lives in `BOOK_TODO.md` **P0a**.
- Affects `_quarto.yml`: the Advanced-track flagging still does not exist there yet
  (see P1) — apply when Ch. 16 is written.
