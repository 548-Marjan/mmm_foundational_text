# Chapter 24 — Software Architecture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author Chapter 24 — architect the MMM system — with the keystone *the module dependency graph is acyclic ⇔ a topological build/test order exists*, and Clean Architecture principles (instability metric, Stable Dependencies Principle, Dependency Inversion) applied to the MMM pipeline.

**Architecture:** Replace the stub body of `parts/07-swe-cs/02-software-architecture.qmd` with the six-heading template. Append a Chapter 24 solutions block to `appendix/solutions.qmd`. Add one BibTeX entry (`@parnas1972`). Verify via headless code run + CI `quarto render`, ship as a PR.

**Tech Stack:** Quarto (`.qmd`, KaTeX math), Python standard library + `matplotlib` only (NO `networkx`), `references.bib`.

---

## STYLING — match Part IV (`parts/04-optimization/01-convexity.qmd`) — USER REQUIREMENT

- Rungs are **`### Rung N — Title` H3 headings** (not inline `**Rung**` bold).
- Proofs are **airy, step-by-step**: `**Proof.**` then each logical step its own short paragraph ending in a colon/transition, blank line, `$$` display with each `$$` on its own line, blank line, next step. Multi-line math via `\begin{aligned}` inside `$$`.
- Theorems/Propositions lead a paragraph as `**Theorem.**` / `**Proposition.**`; proofs end `$\blacksquare$`.
- Blank line before every opening `$$` and after every closing `$$`.

## Other conventions (enforced every task)

- **HOUSE RULES (critical):** Python standard library + `matplotlib` only — **NO `networkx`** or any third-party graph/architecture library; implement topological sort and cycle detection directly. NEVER name PyMC-Marketing or any MMM/PPL/sampler library.
- **KaTeX:** `aligned` inside `$$ … $$` only. `$$` on their own lines. **Even** `$$` count per file. NEVER `\begin{psmallmatrix}`.
- "Key identities" in the Summary must be a **bulleted list**.
- Keep H1 `# Software Architecture` and the anchors line `*Canonical anchors: Martin (Clean Architecture).*`. Remove the stub callout.
- Citation keys: `@martin2017` (exists), `@clrs2009` (exists), `@parnas1972` (added in T9).
- Commit identity `jlh530i` / `jlh530i@gmail.com`. Branch → PR → user merges; never self-merge; never commit to main.

## Verified anchors (confirmed before planning)

- **MMM dependency graph** (edge `u→v` = "u depends on v"): `domain: []`, `ingestion:[domain]`, `transforms:[domain]`, `store:[domain]`, `fit:[domain,transforms]`, `calibration:[domain,store]`, `optimizer:[domain,store]`, `reporting:[fit,optimizer]`.
- **Instability** $I = \text{fan-out}/(\text{fan-in}+\text{fan-out})$: domain `0.00` (fan-in 6), store `0.33`, transforms `0.50`, fit `0.67`, optimizer `0.67`, ingestion `1.00`, calibration `1.00`, reporting `1.00`.
- **Topological build order:** `domain → ingestion → transforms → store → fit → calibration → optimizer → reporting` (all 8 placed).
- **Cycle:** adding `optimizer→calibration` and `calibration→optimizer` → Kahn places only **5 of 8** (cycle detected); routing both through `store` (no direct edges) → all 8 again.

## Critical files

- `parts/07-swe-cs/02-software-architecture.qmd` — the chapter (currently a stub; replace body).
- `appendix/solutions.qmd` — append Chapter 24 block before the final `:::` (one `content-visible` gated div).
- `references.bib` — add `@parnas1972`.
- Styling exemplar (read-only): `parts/04-optimization/01-convexity.qmd`. Predecessor (read-only): `parts/07-swe-cs/01-cs-foundations.qmd` (Ch23, Part VII ch.1).
- Authoritative design: `docs/superpowers/specs/2026-06-18-chapter-24-software-architecture-design.md`.

---

### Task 1: Front matter + Motivation

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`

- [ ] **Step 1: Strip the stub.** Keep lines 1–3 (H1 + blank + anchors line). Delete the stub callout and placeholder italics. Keep the six template headings in order (`## Motivation`, `## Theory & Proofs`, `## Worked Examples`, `## Code Tie-in`, `## Exercises` with `### C/B/P/A`, `## Summary`).

- [ ] **Step 2: Write Motivation (3–4 paragraphs).** Driving question: *"The MMM stack is many parts that evolve at different rates. How do you structure them so a component can be built, tested, and replaced without disturbing the rest — and so the calibration loop doesn't tangle everything together?"* Cover: Ch23 grounded computation, this grounds structure; the MMM modules (domain core, ingestion, transforms, fit, store, calibration, optimizer, reporting); preview the DAG keystone and the data-cycle-vs-dependency-cycle payoff (the Part VI loop is a runtime data cycle but an acyclic module graph). No library named.

- [ ] **Step 3: Verify** even `$$` count, headings present/ordered. Commit: `feat(ch24): motivation`.

---

### Task 2: Theory & Proofs — Rungs 1–2 (modules; dependency graph + instability)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`. Part IV styling. Open `## Theory & Proofs` with one short orienting paragraph listing the five rungs.

- [ ] **Step 1: `### Rung 1 — Modules, interfaces, coupling, and cohesion`.** Module = unit with a public interface hiding internal decisions (information hiding, Parnas `[@parnas1972]`). Cohesion (high, within) vs coupling (low, between). Define fan-out (# it depends on) and fan-in (# depending on it). Introduce the MMM modules: domain core ($S(x;\theta)$ + prior-store interface), ingestion, transforms, fit, store, calibration, optimizer, reporting.

- [ ] **Step 2: `### Rung 2 — The dependency graph and the instability metric`.** $G=(V,E)$, edge $u\to v$ = "$u$ depends on $v$." Display the instability metric

  $$
  I(v) = \frac{\text{fan-out}(v)}{\text{fan-in}(v) + \text{fan-out}(v)} \in [0,1],
  $$

  $I=0$ maximally stable, $I=1$ maximally unstable. State the Stable Dependencies Principle (depend toward decreasing instability). MMM: domain core $I=0$ (high fan-in), reporting/ingestion $I=1$, store deliberately stable.

- [ ] **Step 3: Verify** even `$$` count, no bare `\begin{align}`/`psmallmatrix`, `@parnas1972` used. Commit: `feat(ch24): theory rungs 1-2 (modules, dependency graph, instability)`.

---

### Task 3: Theory & Proofs — Rungs 3–5 (keystone DAG; dependency rule; layering)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`. Part IV styling.

- [ ] **Step 1: `### Rung 3 — Acyclicity ⇔ a topological order exists (keystone)`.** Define topological order (every edge points forward; dependencies before dependents). **`**Theorem.**`** A finite directed graph admits a topological order iff it is acyclic. **Airy Proof:**
  - ($\Rightarrow$) a topological order ⇒ every edge forward ⇒ a cycle would force a backward edge, contradiction ⇒ acyclic.
  - ($\Leftarrow$) acyclic ⇒ a vertex of out-degree 0 exists (else following dependencies forever forces a repeat = cycle); place it, delete, induct = Kahn's algorithm, $O(|V|+|E|)$ (Ch23).
  End `$\blacksquare$`. **Architectural meaning:** topological order = valid build/init/test order; cycle ⇒ no module in it is independently buildable/testable.

- [ ] **Step 2: `### Rung 4 — The Dependency Rule: data-flow cycles versus dependency cycles`.** The Part VI calibration loop *looks* cyclic (optimizer depends on calibration, calibration depends on optimizer = 2-cycle ⇒ untestable by Rung 3). **Dependency Inversion:** both depend on the stable prior-store interface (both edges point to the store), neither on the other; cycle broken, DAG restored, store = stable core ($I\approx0$). State plainly: a runtime data/control cycle does NOT require a source-code dependency cycle. `[@martin2017]`

- [ ] **Step 3: `### Rung 5 — Layering and synthesis: the MMM architecture`.** Layered DAG, dependencies inward: stable core (domain abstractions: $S(x;\theta)$ + store interface), middle (fit, calibration, optimizer), outer (ingestion, reporting, config) which the core never depends on. Payoff: swap sampler/optimizer/data-source without touching the core. Forward to Ch25 (store as built product) and Ch26 (testing via mockable boundaries).

- [ ] **Step 4: Verify** even `$$` count, Rung 3 proof ends `$\blacksquare$`, `@martin2017` used, airy spacing. Commit: `feat(ch24): theory rungs 3-5 (keystone DAG, dependency rule, layering)`.

---

### Task 4: Worked Examples (WE1–WE3)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`. `### WE1/WE2/WE3 — Title` subheads, Part IV spacing. Use the EXACT anchor numbers.

- [ ] **Step 1: WE1 — Instability metrics on the MMM module graph.** Tabulate fan-in/fan-out/$I$: domain `0.00` (fan-in 6), store `0.33`, transforms `0.50`, fit `0.67`, optimizer `0.67`, ingestion/calibration/reporting `1.00`. Verify the Stable Dependencies Principle (edges point toward lower-or-equal $I$).

- [ ] **Step 2: WE2 — The topological build order.** Run Kahn by hand: `domain → ingestion → transforms → store → fit → calibration → optimizer → reporting`; confirm every module built after its dependencies; all 8 placed ⇒ acyclic.

- [ ] **Step 3: WE3 — Breaking and fixing the calibration loop's cycle.** Add optimizer→calibration and calibration→optimizer; Kahn places only 5 of 8 (optimizer, calibration, reporting never freed) ⇒ cycle detected, none independently testable. Apply Dependency Inversion (route both through store, drop direct edges) ⇒ DAG, all 8 sort. Runtime loop unchanged; only source-code dependencies inverted.

- [ ] **Step 4: Verify** even `$$` count; numbers match anchors. Commit: `feat(ch24): worked examples WE1-WE3`.

---

### Task 5: Code Tie-in (single runnable cell)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`

- [ ] **Step 1: Write a single ```{python}``` cell** under `## Code Tie-in` with a short prose preface. **Python standard library + `matplotlib` only — NO networkx.** Figures end `plt.show()`. Four blocks:
  1. **Instability (WE1):** encode the MMM dependency dict; compute fan_in, fan_out, $I$; print table; `assert` domain $I==0$ and domain has max fan-in; `assert` every edge `u→v` has `I(v) <= I(u)` (Stable Dependencies Principle).
  2. **Topological sort (WE2):** implement Kahn's algorithm (pure Python, `collections.deque`); `assert len(order)==8` and the order is a valid build order (for each edge u→v, v precedes u).
  3. **Cycle detection (WE3):** add optimizer↔calibration edges; `assert len(kahn(bad)) < 8`; then store-mediated version; `assert len(kahn(fixed))==8`.
  4. **Figure:** draw the MMM DAG with matplotlib (manual node positions in layers, `annotate` arrows) — left panel DAG (blue), right panel the injected cycle (red edges). No networkx.

- [ ] **Step 2: Extract and run headless.** Copy to `/tmp/ch24_code.py`, run `MPLBACKEND=Agg python3`. Expected: all asserts pass; table prints (domain 0.00 … reporting 1.00); topo order length 8; cycle version < 8. Fix until clean.

- [ ] **Step 3: Verify** even `$$` count; no banned libraries (`grep -inwE 'networkx|pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm'`). Commit: `feat(ch24): code tie-in (instability, topological sort, cycle detection)`.

---

### Task 6: Exercises (C / B / P / A — self-contained)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`

- [ ] **Step 1: Write the four exercise blocks** under `### C/B/P/A`.
  - **C:** (C1) Why does an acyclic graph permit independent build/test while a cycle does not? (C2) Why should dependencies point toward more stable (lower-$I$) modules? (C3) In what sense is the calibration loop a cycle, and in what sense not?
  - **B:** (B1) For a small dependency list, compute fan-in/fan-out/$I$ and identify most/least stable. (B2) Produce a topological order by hand. (B3) State where (if anywhere) the Stable Dependencies Principle is violated.
  - **P:** (P1) Prove a finite digraph admits a topological order iff acyclic (both directions; constructive = Kahn). (P2) Prove that introducing one abstraction node both mutually-dependent modules depend on (Dependency Inversion) converts their 2-cycle into a DAG.
  - **A:** (A1) Write a function returning a topological order or the modules in a cycle; run on the MMM graph with/without the loop edges. (A2) Compute each module's instability and flag any edge violating the Stable Dependencies Principle.

- [ ] **Step 2: Verify** even `$$` count. Commit: `feat(ch24): exercises C/B/P/A`.

---

### Task 7: Summary (auto-included)

**Files:** Modify `parts/07-swe-cs/02-software-architecture.qmd`

- [ ] **Step 1: Write `## Summary`** with a one-paragraph wrap, then **bulleted** "Key concepts" and **bulleted** "Key identities" (inline math). Identities (each a bullet):
  - instability $I(v) = \text{fan-out}/(\text{fan-in}+\text{fan-out}) \in [0,1]$
  - Stable Dependencies Principle: depend toward lower $I$
  - keystone: a finite digraph has a topological order $\iff$ it is acyclic
  - build/test order cost $O(|V|+|E|)$
  - Dependency Rule: route runtime cycles through a stable interface to keep the module graph acyclic
  Close tying back to Part VI (calibration loop as acyclic architecture) and forward to Ch25–26.

- [ ] **Step 2: Verify** Key identities bulleted; even `$$` count. Commit: `feat(ch24): summary`.

---

### Task 8: Appendix solutions

**Files:** Modify `appendix/solutions.qmd`

- [ ] **Step 1: Append `## Chapter 24 — Software Architecture`** immediately before the final `:::`. Part IV airy proof spacing. Full C/B/P/A:
  - **C1–C3:** acyclic ⇒ a build order exists so each module compiles/tests against already-built deps; a cycle ⇒ mutual entanglement, none isolatable. Depend toward stable modules so volatile code does not force churn in stable code. The loop is a runtime *data* cycle (information flows around it over time) but not a source-code *dependency* cycle (mediated by the store interface).
  - **B1–B3:** worked fan-in/fan-out/$I$ on the given list; a valid topological order; SDP-violation check.
  - **P1:** both-directions proof (airy) ending `$\blacksquare$`. **P2:** with edges $A\to B$ and $B\to A$ (2-cycle), introduce $C$ with $A\to C$, $B\to C$ and remove the direct edges; no path returns to $A$ or $B$, so no cycle — a DAG. End `$\blacksquare$`.
  - **A1, A2:** describe the function (Kahn; if `len(order)<|V|` the unplaced nodes are exactly those on/after a cycle) and the SDP-flagging (report edges `u→v` with `I(v) > I(u)`).

- [ ] **Step 2: Verify** even `$$` count in `appendix/solutions.qmd`; block before final `:::`; proofs end `$\blacksquare$`. Commit: `feat(ch24): appendix solutions (C/B/P/A)`.

---

### Task 9: Bibliography + final review + PR

**Files:** Modify `references.bib`; verify chapter + appendix

- [ ] **Step 1: Add the BibTeX entry** to `references.bib`:

```bibtex
@article{parnas1972,
  author  = {Parnas, David L.},
  title   = {On the Criteria To Be Used in Decomposing Systems into Modules},
  journal = {Communications of the ACM},
  year    = {1972},
  volume  = {15},
  number  = {12},
  pages   = {1053--1058}
}
```

- [ ] **Step 2: Lint pass.** Even `$$` count in chapter and appendix; no bare `\begin{align}`; no `\begin{psmallmatrix}`; six headings in order; H1 + anchors line intact; stub callout gone; proofs end `$\blacksquare$`; Key identities bulleted; only `@martin2017`/`@clrs2009`/`@parnas1972` cited and all resolve; no banned library names (`grep -inwE 'networkx|pymc|stan|numpyro|pyro|orbit|causalimpact|lightweight_mmm'`). **Styling check:** `### Rung` H3 headings, airy proofs, `\begin{aligned}` for multi-line.

- [ ] **Step 3: Re-run the extracted code cell headless** (`MPLBACKEND=Agg python3`); all asserts pass.

- [ ] **Step 4: Commit** `feat(ch24): information-hiding reference (Parnas 1972)`.

- [ ] **Step 5: Push the branch and open a PR** against `main` (title `feat: Chapter 24 — Software Architecture`; body summarizing the DAG keystone + the data-cycle-vs-dependency-cycle payoff). Start the background merge-watcher.

- [ ] **Step 6: Watch CI `quarto render` (HTML + PDF) to a green conclusion** — the real gate. Report the green PR; the user merges. **Do not push follow-up commits after the PR is up** (the user merges fast — keep each PR a single complete change).

---

## Self-Review (done at plan time)

- **Spec coverage:** Rungs 1–5 → Tasks 2–3; WE1–3 → Task 4; code 4 blocks → Task 5; exercises C/B/P/A → Task 6; summary → Task 7; appendix → Task 8; bib `@parnas1972` → Task 9. All mapped.
- **Anchor consistency:** instability table, topological order (8 modules), cycle (5 of 8) identical across Theory, Worked Examples, Code, Appendix — all verified above.
- **Styling:** Part IV `### Rung` headings + airy proofs baked into Tasks 2, 3, 8 and the Task 9 styling check.
- **No placeholders; house rules** (no networkx / no library names, KaTeX, even `$$`, bulleted Key identities, `$\blacksquare$`, branch→PR→user-merges) enforced per task.
