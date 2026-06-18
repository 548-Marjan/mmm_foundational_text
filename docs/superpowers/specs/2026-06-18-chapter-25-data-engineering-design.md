# Chapter 25 — Data Engineering: Design Spec

**Date:** 2026-06-18
**Part:** VII — Software Engineering & Computer Science (third chapter)
**File:** `parts/07-swe-cs/03-data-engineering.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter builds the Chapter 22 prior store as a real, versioned data product, in the idiom of
Kleppmann's *Designing Data-Intensive Applications* [@kleppmann2017]. It is where the back half of the
book converges: Chapter 22's **sequential = batch** theorem becomes the data-engineering correctness
guarantee, and Chapter 24's **replay-the-ledger** becomes time-travel reproducibility. Chapter 24 defined
the store's *interface* and its *position* in the architecture; this chapter defines its *implementation*
as an append-only event log whose fold is the live posterior.

**Center of gravity (locked):** the prior store as an **append-only event log**, with the live posterior
as a **deterministic fold** over it. The keystone is that a **deduplicated** log with a
**commutative/associative** fold yields a posterior that is reproducible, reorder-safe, and retry-safe —
and the sharp insight is that *commutativity comes for free* (Chapter 22) while *idempotency must be
engineered* (dedup by event id), because Bayesian updates are not idempotent.

**Driving question (locked):** *"Experiments arrive over years, out of order, sometimes twice, under a
schema that keeps changing. How do you store them so the live posterior is always correct, reproducible,
and safe to recompute?"*

**Styling (locked, user directive):** match Part IV (`parts/04-optimization/01-convexity.qmd`) — rungs
as `### Rung N — Title` H3 headings, airy step-by-step proofs (each step its own short paragraph leading
into a `$$` display on its own lines, `\begin{aligned}` for multi-line).

---

## 2. Scope discipline

- Implements the Chapter 22 prior store / Chapter 24 store interface; introduces no new modeling.
- Does **not** re-derive the calibration factor (Chapter 21) or sequential = batch (Chapter 22) — it
  *uses* them. Does **not** cover testing mechanics (Chapter 26) beyond noting the invariants worth
  property-testing.
- **House rule (critical):** `numpy` + Python standard library + `matplotlib` only. No database/ORM/
  streaming library; the event log, fold, and dedup are implemented directly. Never name PyMC-Marketing
  or any MMM/PPL/sampler library.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their own lines;
  even `$$` count.

---

## 3. Theory & Proofs — the rungs (5; keystone at Rung 3)

**Rung 1 — The store as an append-only event log.**
Event sourcing: rather than store the posterior and mutate it, store the **immutable, append-only log**
of the events that produced it — one record per calibration experiment, carrying an `event_id`, the
estimand class (secant / tangent / validation, the Chapter 20/21 taxonomy), the constraint direction
$c$, the measurement variance $s^2$, the power weight $w$, and a timestamp. The live posterior is
**derived state**: a fold (reduction) over the log, never the source of truth. State why source-of-truth
is the log, not the posterior: the log is auditable, replayable, and survives schema and code changes,
while the posterior is a cache that can always be rebuilt.

**Rung 2 — The fold and its order-independence.**
Define the fold: $\Lambda_k = \text{fold}(\oplus,\ \Lambda_0,\ \text{events}[1{:}k])$. In the
linear-Gaussian store the combine operation adds a precision contribution:
$$
\Lambda_k = \Lambda_0 + \sum_{j=1}^{k} w_j\, s_j^{-2}\, c_j c_j^\top.
$$
**Proposition (order-independence).** Because each event contributes an additive (positive-semidefinite)
term and addition is commutative and associative, the fold is invariant to the order in which events are
applied. **Proof:** the sum $\sum_j w_j s_j^{-2} c_j c_j^\top$ is unchanged under any permutation of its
terms (commutativity and associativity of matrix addition); this is the matrix form of Chapter 22's
sequential = batch result. $\blacksquare$
**Consequences:** the store may ingest events **out of order**, and **replay-to-version-$k$** (fold the
first $k$ events by timestamp) is a well-defined, reproducible query — the data-engineering form of
Chapter 24's "replay the ledger." Anchor: with prior precision $4$ and contributions $25$ then $16$, the
fold in order and reversed both give $45$.

**Rung 3 — Proof P (KEYSTONE): idempotency must be engineered; dedup makes at-least-once safe.**
Real delivery is **at-least-once**: a producer that does not hear an acknowledgement retries, so the log
can receive the **same** event twice. A Bayesian update is **not idempotent** — re-applying a factor
adds its precision again:
$$
\Lambda_0 + s^{-2}cc^\top + s^{-2}cc^\top \neq \Lambda_0 + s^{-2}cc^\top,
$$
so a duplicated event double-counts the evidence (anchor: $4 + 25 + 25 + 16 = 70$, not $45$).
**Theorem.** If the log is **deduplicated by `event_id`** (each id folded at most once), then processing
the log under at-least-once delivery yields exactly the posterior of exactly-once delivery. **Proof:**
dedup makes the fold a function of the **set** of distinct event ids rather than the delivered sequence
(a multiset); two deliveries of the same id collapse to one application, so the realized fold equals the
fold over the de-duplicated set, which is the exactly-once result. The de-duplicated update is therefore
idempotent at the log level: re-delivering a processed id leaves the posterior unchanged. $\blacksquare$
State the distinction plainly and as the chapter's thesis: **commutativity/associativity is free** (Rung
2, inherited from the mathematics of the model), but **idempotency is an engineering obligation** —
without an explicit dedup key the model's own non-idempotence corrupts the posterior under ordinary
retry. Anchor: duplicate without dedup $\to 70$; with dedup $\to 45$.

**Rung 4 — Schema evolution and compatibility.**
The store outlives any single schema: new estimand types, new fields, deprecated columns. Define
**backward compatibility** (new readers can read old events) and **forward compatibility** (old readers
tolerate new fields), the pair that lets producers and consumers deploy independently. The secant /
tangent / validation taxonomy is the **stable schema core** (carried since Chapters 20–21 and made an
interface in Chapter 24); evolution adds optional fields around it rather than changing it. Treat a
**migration** as a pure function applied to the log (or to the fold), versioned alongside the data, so
that any historical version is still reconstructible. Keep this rung brief and DDIA-grounded.

**Rung 5 — Synthesis: the store as a data product.**
Assemble the pieces: an append-only, deduplicated, schema-versioned event log whose fold yields the live
posterior, served behind the Chapter 24 interface and consumed by the Chapter 18 optimizer. The three
guarantees line up with the three rungs — reorder-safe (Rung 2), retry-safe (Rung 3), reproducible
(replay-to-version, Rungs 2–4) — and together they make the Chapter 22 prior store a genuine
data-intensive product rather than a mutable cache. Close forward: Chapter 26 turns these invariants
(order-independence, idempotency under dedup, replay determinism) into property tests, and Chapter 27
assembles the whole system.

---

## 4. Worked Examples

**WE1 — The fold is the posterior; replay reproduces it.**
Prior precision $\Lambda_0 = 4$ (the Chapter 22 scale prior). Two logged events contribute precision
$25$ (a lift with $s = 0.2$) then $16$ (a lift with $s = 0.25$). Folding the log gives $4 \to 29 \to
45$, identical to the batch precision $4 + 25 + 16 = 45$ (posterior std $1/\sqrt{45} \approx 0.149$).
Replaying to version $1$ (first event only) reproduces $29$ exactly — a deterministic time-travel query
against the immutable log.

**WE2 — Out-of-order arrival changes nothing.**
The same two events arrive in the opposite order. The fold gives $4 \to 20 \to 45$ — a different
intermediate, the same final posterior $45$ — because the precision contributions add and addition
commutes (Rung 2). The store can therefore accept late-arriving events without recomputing history in a
fixed order.

**WE3 — At-least-once delivery, and why dedup is mandatory.**
The first event is delivered twice (a producer retry after a lost acknowledgement). Without
deduplication the fold gives $4 + 25 + 25 + 16 = 70$ — the evidence of one experiment counted twice, a
silently overconfident posterior (std $1/\sqrt{70} \approx 0.120$ instead of $0.149$). Deduplicating by
`event_id` folds the repeated id once and restores $45$. The model bought commutativity for free but not
idempotency; the dedup key is what makes ordinary retry safe.

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + Python standard library; `matplotlib` for the figure;
verified headless). No database/streaming library. It:
1. **The event log and the fold (WE1):** represent the log as a list of event records (`event_id`,
   precision contribution, timestamp); implement `fold(log)` summing prior precision plus contributions;
   `assert` the fold equals the batch precision $45$ and that replay-to-version-$1$ gives $29$.
2. **Order-independence (WE2):** fold the log and its reverse; `assert` both give $45$.
3. **Dedup / idempotency (WE3):** duplicate the first event; `assert` the no-dedup fold gives $70$ and
   the dedup-by-`event_id` fold gives $45$; re-deliver an already-folded id and `assert` the posterior is
   unchanged (idempotent).
4. **Figure:** posterior precision (and std) versus log version — the accumulation staircase — with the
   duplicated/no-dedup run overlaid to show the spurious jump.
Every numeric claim asserted; deterministic (no RNG needed).

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why is the posterior stored as *derived state* (a fold over the log) rather than as the source
  of truth? Why does the store get commutativity "for free" but must *engineer* idempotency? What goes
  wrong under at-least-once delivery without a dedup key?
- **B:** Given a small event log (prior precision and a few contributions), compute the fold, the
  replay-to-version-$k$ value, the out-of-order fold, and the duplicate (no-dedup vs dedup) folds.
- **P:** (i) Prove the fold is order-independent (commutativity/associativity of the additive precision
  update). (ii) Prove that deduplication by `event_id` makes the fold under at-least-once delivery equal
  to the exactly-once posterior, and hence idempotent under redelivery.
- **A:** Extend the Code Tie-in — implement an append-only log with a `replay(version)` query and a
  dedup set; add a late-arriving (out-of-order) event and a duplicate event, and show the posterior is
  invariant to both; then add a schema field to new events and show old events still fold (forward
  compatibility).

---

## 7. Appendix solutions

Append `## Chapter 25 — Data Engineering` to `appendix/solutions.qmd`, **in chapter order** (after the
Chapter 24 block), inside the existing `content-visible` gated div. Full C/B/P/A; the P-block carries
both proofs (fold order-independence; dedup ⇒ at-least-once equals exactly-once). Part IV airy proof
spacing.

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities/quantities:
the posterior as a fold $\Lambda_k = \Lambda_0 + \sum_{j\le k} w_j s_j^{-2} c_j c_j^\top$; order-
independence (commutativity/associativity of the additive update); non-idempotence of the raw Bayesian
update; the dedup theorem (dedup by `event_id` ⇒ at-least-once = exactly-once); reproducible
replay-to-version. Close tying back to Chapter 22 (sequential = batch as the correctness guarantee) and
Chapter 24 (replay the ledger), and forward to Chapters 26–27.

---

## 9. Cross-references

- **Chapter 21 (Advanced Calibration):** the calibration factor each log event records.
- **Chapter 22 (The Prior Store):** sequential = batch, reused as the fold's order-independence.
- **Chapter 24 (Software Architecture):** the store interface and "replay the ledger," realized here.
- **Chapter 18 (Budget Optimization):** the optimizer that consumes the served posterior.
- **Chapters 26–27:** property-testing the invariants; the full system.

---

## 10. Bibliography

Reuse `@kleppmann2017` (Designing Data-Intensive Applications — event logs, at-least-once vs exactly-once
delivery, schema evolution; already present). No additions expected; the plan confirms all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor verified (done: fold in-order $= 45$, reversed $= 45$; replay-to-version-$1 = 29$;
  duplicate no-dedup $= 70$; duplicate dedup $= 45$).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`); no database/streaming
  library.
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no `\begin{align}`,
  no `psmallmatrix`, valid citation keys, no banned library names; Part IV `### Rung` headings + airy
  proofs.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.

---

## 12. Build note

This worktree branched off `main` before Chapter 24 merged into the appendix. **Before building**,
rebase onto the latest `main` (which now contains the Chapter 24 appendix block) so the Chapter 25
solutions block appends in correct chapter order after Chapter 24 and avoids a `solutions.qmd` conflict.
