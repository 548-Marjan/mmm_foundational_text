# Chapter 25 — Data Engineering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author Chapter 25 — build the Chapter 22 prior store as a real versioned data product — with the keystone *a deduplicated append-only event log + a commutative/associative fold = a deterministic, reorder-safe, retry-safe versioned posterior*, where commutativity is free (Ch. 22 sequential=batch) but idempotency must be engineered (dedup by event id).

**Architecture:** Replace the stub body of `parts/07-swe-cs/03-data-engineering.qmd` with the six-heading template. Append a Chapter 25 solutions block to `appendix/solutions.qmd`. No BibTeX additions (`@kleppmann2017` already present). Verify via headless code run + CI `quarto render`, ship as a PR.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), `numpy` + Python standard library + `matplotlib` only (no DB/ORM/streaming library), `references.bib`.

---

## STEP 0 (before building): rebase onto current main

This branch was created before Chapter 24 merged. **First action of the build:** `git fetch origin && git rebase origin/main`, so the Chapter 25 appendix block appends in chapter order after the Chapter 24 block and `solutions.qmd` does not conflict. Confirm `appendix/solutions.qmd` contains `## Chapter 24 — Software Architecture` after the rebase.

## STYLING — match Part IV (`parts/04-optimization/01-convexity.qmd`)

- Rungs are **`### Rung N — Title` H3 headings**.
- Proofs are **airy, step-by-step**: `**Proof.**` then each step its own short paragraph leading into a `$$` display on its own lines; `\begin{aligned}` for multi-line; end `$\blacksquare$`.
- Blank line before every opening `$$` and after every closing `$$`.

## Other conventions (enforced every task)

- **HOUSE RULES:** `numpy` + Python standard library + `matplotlib` only — no database/ORM/streaming library; implement the log, fold, and dedup directly. NEVER name PyMC-Marketing or any MMM/PPL/sampler library.
- **KaTeX:** `aligned` inside `$$ … $$` only. `$$` on their own lines. **Even** `$$` count per file. NEVER `\begin{psmallmatrix}`.
- "Key identities" in the Summary must be a **bulleted list**.
- Keep H1 `# Data Engineering` and the anchors line `*Canonical anchors: Kleppmann (DDIA).*`. Remove the stub callout.
- Citation key: `@kleppmann2017` (already in `references.bib`).
- Commit identity `jlh530i` / `jlh530i@gmail.com`. Branch → PR → user merges; never self-merge; never commit to main.

## Verified anchors (confirmed before planning)

- Prior precision $\Lambda_0 = 4$; logged events contribute precision $25$ (lift $s=0.2$) then $16$ (lift $s=0.25$).
- **Fold = batch:** $4 \to 29 \to 45$; replay-to-version-1 $= 29$; posterior std $1/\sqrt{45} \approx 0.149$.
- **Out-of-order:** reversed gives $4 \to 20 \to 45$ — same final $45$.
- **At-least-once duplicate:** no dedup $\to 4+25+25+16 = 70$ (std $1/\sqrt{70}\approx 0.120$); dedup by `event_id` $\to 45$; redelivering a folded id leaves the posterior unchanged.

## Critical files

- `parts/07-swe-cs/03-data-engineering.qmd` — the chapter (currently a stub; replace body).
- `appendix/solutions.qmd` — append Chapter 25 block before the final `:::` (after the Chapter 24 block, post-rebase).
- Styling exemplar (read-only): `parts/04-optimization/01-convexity.qmd`. Predecessors (read-only): `parts/07-swe-cs/02-software-architecture.qmd` (Ch24), `parts/06-causal-calibration/04-prior-store.qmd` (Ch22, sequential=batch).
- Authoritative design: `docs/superpowers/specs/2026-06-18-chapter-25-data-engineering-design.md`.

---

### Task 1: Front matter + Motivation

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`

- [ ] **Step 1: Strip the stub.** Keep lines 1–3 (H1 + blank + anchors line). Delete the stub callout and placeholder italics. Keep the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises` with `### C/B/P/A`, `## Summary`).

- [ ] **Step 2: Write Motivation (3–4 paragraphs).** Driving question: *"Experiments arrive over years, out of order, sometimes twice, under a schema that keeps changing. How do you store them so the live posterior is always correct, reproducible, and safe to recompute?"* Cover: Ch24 defined the store's interface/position, this builds its implementation; the prior store as an append-only event log whose fold is the posterior; preview the keystone (commutativity free from Ch22, idempotency engineered via dedup); ground in DDIA [@kleppmann2017]. No library named.

- [ ] **Step 3: Verify** even `$$` count, headings present/ordered. Commit: `feat(ch25): motivation`.

---

### Task 2: Theory & Proofs — Rungs 1–2 (event log; order-independent fold)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`. Part IV styling. Open `## Theory & Proofs` with one short orienting paragraph listing the five rungs.

- [ ] **Step 1: `### Rung 1 — The store as an append-only event log`.** Event sourcing: immutable append-only log, one record per experiment (`event_id`, estimand class secant/tangent/validation, constraint $c$, variance $s^2$, weight $w$, timestamp). The live posterior is derived state (a fold), the log is source of truth (auditable, replayable, survives schema/code change; the posterior is a rebuildable cache).

- [ ] **Step 2: `### Rung 2 — The fold and its order-independence`.** Define the fold; display the additive precision update

  $$
  \Lambda_k = \Lambda_0 + \sum_{j=1}^{k} w_j\, s_j^{-2}\, c_j c_j^\top .
  $$

  **`**Proposition (order-independence).**`** + airy **Proof**: the sum is unchanged under permutation (commutativity/associativity of matrix addition) — the matrix form of Ch22 sequential=batch. End `$\blacksquare$`. Consequences: out-of-order ingest; replay-to-version-$k$ well-defined (Ch24 replay-the-ledger). Anchor: in-order = reversed = 45.

- [ ] **Step 3: Verify** even `$$` count, no bare `\begin{align}`/`psmallmatrix`, Rung 2 proof ends `$\blacksquare$`. Commit: `feat(ch25): theory rungs 1-2 (event log, order-independent fold)`.

---

### Task 3: Theory & Proofs — Rungs 3–5 (keystone dedup/idempotency; schema; synthesis)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`. Part IV styling.

- [ ] **Step 1: `### Rung 3 — Idempotency must be engineered; dedup makes at-least-once safe (keystone)`.** At-least-once delivery ⇒ duplicates. Bayesian update is NOT idempotent — display

  $$
  \Lambda_0 + s^{-2}cc^\top + s^{-2}cc^\top \neq \Lambda_0 + s^{-2}cc^\top ,
  $$

  so a duplicate double-counts (anchor $4+25+25+16=70$). **`**Theorem.**`** Dedup by `event_id` ⇒ at-least-once fold = exactly-once posterior. **Airy Proof:** dedup makes the fold a function of the *set* of distinct ids, not the delivered multiset; two deliveries of one id collapse to one application; hence equals exactly-once, and redelivery is idempotent. End `$\blacksquare$`. Thesis: **commutativity free (Rung 2), idempotency engineered (here).** Anchor: dup no-dedup 70, dedup 45.

- [ ] **Step 2: `### Rung 4 — Schema evolution and compatibility`.** Backward compatibility (new readers read old events) + forward compatibility (old readers tolerate new fields) ⇒ producers/consumers deploy independently. The secant/tangent/validation taxonomy = stable schema core (Ch20/21/24); evolution adds optional fields around it. Migration = pure function on the log/fold, versioned with the data, so historical versions stay reconstructible. Brief, DDIA-grounded `[@kleppmann2017]`.

- [ ] **Step 3: `### Rung 5 — Synthesis: the store as a data product`.** Append-only + deduplicated + schema-versioned log → fold → live posterior, served behind the Ch24 interface into the Ch18 optimizer. Three guarantees ↔ three rungs: reorder-safe (R2), retry-safe (R3), reproducible (replay, R2–4). Forward: Ch26 property-tests these invariants; Ch27 assembles the system.

- [ ] **Step 4: Verify** even `$$` count, Rung 3 proof ends `$\blacksquare$`, `@kleppmann2017` used, airy spacing. Commit: `feat(ch25): theory rungs 3-5 (keystone dedup/idempotency, schema, synthesis)`.

---

### Task 4: Worked Examples (WE1–WE3)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`. `### WE1/WE2/WE3 — Title` subheads, Part IV spacing. Use the EXACT anchor numbers.

- [ ] **Step 1: WE1 — The fold is the posterior; replay reproduces it.** $\Lambda_0=4$; contributions $25$ then $16$; fold $4\to29\to45$ = batch $4+25+16=45$ (std $\approx0.149$); replay-to-version-1 $=29$ (deterministic time-travel query).

- [ ] **Step 2: WE2 — Out-of-order arrival changes nothing.** Reversed order $4\to20\to45$ — different intermediate, same final $45$ (commutativity). Late-arriving events need no fixed-order recompute.

- [ ] **Step 3: WE3 — At-least-once delivery, and why dedup is mandatory.** First event delivered twice; no-dedup $\to70$ (std $\approx0.120$, silently overconfident); dedup by `event_id` $\to45$. Commutativity free, idempotency engineered.

- [ ] **Step 4: Verify** even `$$` count; numbers match anchors. Commit: `feat(ch25): worked examples WE1-WE3`.

---

### Task 5: Code Tie-in (single runnable cell)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`

- [ ] **Step 1: Write a single ```{python}``` cell** under `## Code Tie-in` with a short prose preface. `numpy` + Python standard library + `matplotlib` only — NO database/streaming library. Figure ends `plt.show()`. Four blocks:
  1. **Event log + fold (WE1):** log = list of records `{"event_id","prec","t"}`; `fold(log, prior=4.0)` sums prior + contributions; `assert fold == 45`; `replay(log, version=1) == 29`.
  2. **Order-independence (WE2):** `assert fold(log) == fold(list(reversed(log))) == 45`.
  3. **Dedup / idempotency (WE3):** duplicate the first event; `assert fold(dup, dedup=False) == 70`; `assert fold(dup, dedup=True) == 45`; redeliver a folded id and `assert` posterior unchanged.
  4. **Figure:** posterior precision (and/or std) vs log version — the accumulation staircase — with the no-dedup duplicate run overlaid showing the spurious jump to 70.

- [ ] **Step 2: Extract and run headless.** Copy to `/tmp/ch25_code.py`, run `MPLBACKEND=Agg python3`. Expected: all asserts pass; prints show fold 45 / replay 29 / reversed 45 / dup-no-dedup 70 / dup-dedup 45. Fix until clean (confirm Python exit 0).

- [ ] **Step 3: Verify** even `$$` count; no banned libraries (`grep -inwE 'sqlalchemy|psycopg|kafka|sqlite3|pandas|networkx|pymc|stan|numpyro|pyro'` — note `sqlite3` is stdlib but avoid for clarity; use plain lists/dicts). Commit: `feat(ch25): code tie-in (event log, fold, dedup, replay)`.

---

### Task 6: Exercises (C / B / P / A — self-contained)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`

- [ ] **Step 1: Write the four exercise blocks** under `### C/B/P/A`.
  - **C:** (C1) Why store the posterior as derived state (a fold) rather than as source of truth? (C2) Why is commutativity free but idempotency engineered? (C3) What goes wrong under at-least-once delivery without a dedup key?
  - **B:** (B1) Given prior precision and a few contributions, compute the fold and posterior std. (B2) Compute replay-to-version-$k$ and the out-of-order fold. (B3) Compute the duplicate fold with and without dedup.
  - **P:** (P1) Prove the fold is order-independent (commutativity/associativity of the additive update). (P2) Prove dedup by `event_id` makes the at-least-once fold equal to the exactly-once posterior, hence idempotent under redelivery.
  - **A:** (A1) Implement an append-only log with `replay(version)` and a dedup set; add a late-arriving and a duplicate event, show invariance. (A2) Add a schema field to new events and show old events still fold (forward compatibility).

- [ ] **Step 2: Verify** even `$$` count. Commit: `feat(ch25): exercises C/B/P/A`.

---

### Task 7: Summary (auto-included)

**Files:** Modify `parts/07-swe-cs/03-data-engineering.qmd`

- [ ] **Step 1: Write `## Summary`** with a one-paragraph wrap, then **bulleted** "Key concepts" and **bulleted** "Key identities" (inline math). Identities (each a bullet):
  - posterior as a fold $\Lambda_k = \Lambda_0 + \sum_{j\le k} w_j s_j^{-2} c_j c_j^\top$
  - order-independence: the additive update commutes and associates (Ch22 sequential=batch)
  - non-idempotence: $\Lambda_0 + s^{-2}cc^\top + s^{-2}cc^\top \neq \Lambda_0 + s^{-2}cc^\top$
  - dedup theorem: dedup by `event_id` ⇒ at-least-once $=$ exactly-once
  - reproducible replay-to-version-$k$
  Close tying back to Ch22 (sequential=batch as the correctness guarantee) and Ch24 (replay the ledger), forward to Ch26–27.

- [ ] **Step 2: Verify** Key identities bulleted; even `$$` count. Commit: `feat(ch25): summary`.

---

### Task 8: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Append `## Chapter 25 — Data Engineering`** immediately before the final `:::`, **after the Chapter 24 block** (confirm Ch24 block present from the Step 0 rebase). Part IV airy proof spacing. Full C/B/P/A:
  - **C1–C3:** the log is the auditable, replayable source of truth and the posterior a rebuildable cache; commutativity follows from the model (additive update) while idempotency does not (re-applying a factor double-counts), so it must be enforced with a dedup key; without dedup, a retry under at-least-once delivery silently double-counts evidence and overstates precision.
  - **B1–B3:** fold $= 45$, std $1/\sqrt{45}\approx0.149$; replay-v1 $=29$, out-of-order $=45$; duplicate no-dedup $=70$ (std $\approx0.120$), dedup $=45$.
  - **P1:** order-independence via commutativity/associativity of the additive precision sum (airy), end `$\blacksquare$`. **P2:** dedup ⇒ fold depends on the *set* of ids ⇒ at-least-once equals exactly-once and redelivery is idempotent (airy), end `$\blacksquare$`.
  - **A1, A2:** describe the log+replay+dedup implementation (invariant to a late and a duplicate event) and forward compatibility (old events lacking the new field still fold, e.g. via a default).

- [ ] **Step 2: Verify** even `$$` count in `appendix/solutions.qmd`; block before final `:::`; proofs end `$\blacksquare$`. Commit: `feat(ch25): appendix solutions (C/B/P/A)`.

---

### Task 9: Final review + PR (no bib additions)

**Files:** verify chapter + appendix

- [ ] **Step 1: Confirm `@kleppmann2017` resolves** in `references.bib` (already present) — no additions.

- [ ] **Step 2: Lint pass.** Even `$$` count in chapter and appendix; no bare `\begin{align}`; no `\begin{psmallmatrix}`; six headings in order; H1 + anchors line intact; stub callout gone; proofs end `$\blacksquare$`; Key identities bulleted; only `@kleppmann2017` cited and it resolves; no banned library names. **Styling check:** `### Rung` H3 headings, airy proofs, `\begin{aligned}` for multi-line.

- [ ] **Step 3: Re-run the extracted code cell headless** (`MPLBACKEND=Agg python3`); confirm Python exit 0, all asserts pass.

- [ ] **Step 4: Commit** any final lint fixes: `chore(ch25): final lint`.

- [ ] **Step 5: Push the branch and open a PR** against `main` (title `feat: Chapter 25 — Data Engineering`; body summarizing the log/fold keystone + commutativity-free / idempotency-engineered + the Ch22/Ch24 tie-in). Start the background merge-watcher.

- [ ] **Step 6: Watch CI `quarto render` (HTML + PDF) to a green conclusion** — the real gate. Report the green PR; the user merges. **Do not push follow-up commits after the PR is up** (keep each PR a single complete change).

---

## Self-Review (done at plan time)

- **Spec coverage:** Rungs 1–5 → Tasks 2–3; WE1–3 → Task 4; code 4 blocks → Task 5; exercises C/B/P/A → Task 6; summary → Task 7; appendix → Task 8. All mapped.
- **Anchor consistency:** fold 45 / replay 29 / out-of-order 45 / dup-no-dedup 70 / dup-dedup 45 identical across Theory, Worked Examples, Code, Appendix — all verified above.
- **Styling:** Part IV `### Rung` headings + airy proofs baked into Tasks 2, 3, 8 and the Task 9 styling check.
- **Rebase note:** Step 0 rebases onto Ch.24-bearing main before any appendix edit.
- **No placeholders; house rules** (no DB/streaming/library, KaTeX, even `$$`, bulleted Key identities, `$\blacksquare$`, branch→PR→user-merges) enforced per task.
