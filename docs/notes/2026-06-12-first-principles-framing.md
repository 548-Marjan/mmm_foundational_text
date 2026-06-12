# First-Principles Framing

**Date:** 2026-06-12
**Context:** Conversation clarifying what "from first principles" actually means and assessing
whether *The Mathematics and Engineering of Marketing Mix Modeling* earns that label. Captured
here as positioning reference. Companion to
[`2026-06-11-positioning-scope-publishing.md`](2026-06-11-positioning-scope-publishing.md).

---

## What "from first principles" means

Borrowed from philosophy and physics (Aristotle's *archē*; the way Feynman or Musk talk about
reasoning). A **first principle** is a foundational truth that can't be deduced from anything more
basic — an axiom, a definition, a law taken as bedrock.

To reason **from first principles** means you:

1. **Decompose** a problem down to those bedrock truths, and
2. **Rebuild** understanding upward by deduction — *not* by analogy, authority, or convention.

The contrast is **reasoning by analogy** ("we do it this way because that's how it's done / because
the library does it"). First-principles reasoning refuses to take an intermediate result on faith:
you either derive it, or you explicitly flag that you're borrowing it.

### Nuances (the phrase gets abused as a slogan)

- **It's relative to a chosen floor.** Nobody re-derives Peano arithmetic to do regression. "First
  principles" always means *first principles relative to a declared starting point*. The honesty is
  in **declaring the floor** and then not smuggling in unearned results above it.
- **It's about the chain of justification, not difficulty.** A book can be hard and still *not* be
  from first principles (asserts advanced results, drills mechanics). A book can be elementary and
  *deeply* first-principles (earns every step).
- **The tell is "why is this true?"** at every link. If the answer is ever "because the
  textbook/library says so" — that link isn't from first principles.

## Does the book earn the label?

**Largely yes — and it's explicit about its floor, which is the honest version of the claim.**

| First-principles trait | The book |
|---|---|
| **Declares its floor** | Yes — Calc 2, linear algebra, calculus-based stats are stated prerequisites. Doesn't pretend to derive calculus. |
| **Builds upward by deduction** | Yes — "foundations → applications," "each part builds strictly on the ones before it." MMM at the end is a *synthesis* of Parts I–IV, not an opening assertion. |
| **Refuses hand-waving** | Explicitly: "without hand-waving and without assuming you already know Bayesian modeling, MCMC, or constrained optimization." |
| **Proves vs. cites, deliberately** | "Where a result is genuinely illuminating we prove it; where a proof is mechanical we cite it." The correct discipline — derive what illuminates, declare what you borrow. |
| **Doesn't reason by analogy to tooling** | Strongly — house rule is *explain calibration generically, never name PyMC-Marketing.* Teaches the mechanism, not "do what the library does." This is the essence of first-principles vs. analogy. |

### "First principles" with an asterisk (fine — every such book has one)

- The floor is **Calc 2 + linalg + stats**, not ZFC or real analysis. So it's "first principles
  *for a quantitatively-prepared reader*," not "from nothing." The subtitle already says this.
- "Mechanical proofs are cited" means a few links are *borrowed by reference* rather than rebuilt —
  a deliberate, declared exception, which is the intellectually honest way to do it.

## Bottom line / positioning

Fair to call it a from-first-principles MMM text — arguably its central differentiator. The
accurate framing isn't "from zero" but:

> **From a declared mathematical floor, every step to advanced MMM is earned by derivation rather
> than asserted or borrowed from a library.**

That's a stronger, more defensible claim than the bare slogan.

**Naming thought:** the current subtitle ("self-contained path from Calc 2… to advanced MMM
practice") *describes* first-principles without using the phrase. If the phrase is wanted as
positioning, lean on the **"every step earned, library-agnostic"** angle — that's what genuinely
backs it up, and it differentiates from tool-tutorial MMM books that reason purely by analogy to
one vendor's API.
