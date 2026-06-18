# Chapter 24 — Software Architecture: Design Spec

**Date:** 2026-06-18
**Part:** VII — Software Engineering & Computer Science (second chapter)
**File:** `parts/07-swe-cs/02-software-architecture.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter architects the MMM system. Where Chapter 23 grounded the *computation* (arithmetic, cost),
this chapter grounds the *structure*: how the many moving parts — ingestion, transforms, model fitting,
the prior store, calibration, optimization, reporting — are organized so that any one can be built,
tested, and replaced without disturbing the rest. The center of gravity is the **module dependency
graph as a directed acyclic graph (DAG)**, with Clean Architecture's principles (Martin) as the
supporting framework, applied throughout to the MMM pipeline.

**Center of gravity (locked, via AskUserQuestion):** the DAG keystone — *the module dependency graph is
acyclic if and only if a topological build/test order exists* — with Clean Architecture principles
(the Dependency Rule, the instability metric, information hiding) as support. The proof-bearing template
is kept.

**Driving question (locked):** *"The MMM stack is many parts that evolve at different rates. How do you
structure them so a component can be built, tested, and replaced without disturbing the rest — and so
the calibration loop doesn't tangle everything together?"*

**Unifying payoff:** the calibration–optimization loop of Part VI (Chapters 18 ↔ 20 ↔ 21 ↔ 22) is a
runtime *data* cycle, yet the *module dependency* graph stays acyclic, because the optimizer and the
calibration step communicate through the stable **prior-store interface** rather than depending on each
other directly. **Data-flow cycles are not dependency cycles** — and that distinction *is* Clean
Architecture's Dependency Rule (the Dependency Inversion Principle).

**Styling (locked, user directive):** match Part IV (`parts/04-optimization/01-convexity.qmd`) — rungs
as `### Rung N — Title` H3 headings, airy step-by-step proofs (each step its own short paragraph leading
into a `$$` display surrounded by blank lines, `\begin{aligned}` for multi-line).

---

## 2. Scope discipline

- Architects the MMM system already built in Parts I–VI; introduces no new modeling.
- Does **not** build the prior store as a deployed data product (storage, schema migration, pipelines) —
  that is Chapter 25. This chapter defines its *interface* and *position* in the architecture.
- Does **not** cover testing mechanics — Chapter 26. It establishes the acyclic, mockable structure that
  makes testing possible.
- **House rule (critical):** `numpy`/`scipy`/`matplotlib` and the Python standard library only.
  **No `networkx`** or any third-party graph/architecture library — the topological sort and cycle
  detection are implemented directly. Never name PyMC-Marketing or any MMM/PPL/sampler library.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their own lines;
  even `$$` count.

---

## 3. Theory & Proofs — the rungs (5; keystone at Rung 3)

**Rung 1 — Modules, interfaces, coupling, and cohesion.**
A **module** is a unit of software with a public interface that hides its internal decisions
(*information hiding*, Parnas [@parnas1972]): consumers depend on what it does, not how. Two qualities
govern a decomposition: **cohesion**, the relatedness of what lives inside a module, which should be
high; and **coupling**, the dependence between modules, which should be low. Make these countable:
**fan-out** of a module is the number of modules it depends on, **fan-in** the number that depend on it.
Introduce the MMM modules used throughout the chapter: a stable **domain** core (the response curve
$S(x;\theta)$ and the prior-store interface), **ingestion**, **transforms** (adstock/saturation),
**fit**, **store**, **calibration**, **optimizer**, and **reporting**.

**Rung 2 — The dependency graph and the instability metric.**
Model the system as a directed graph $G = (V, E)$ with $V$ the modules and an edge $u \to v$ meaning
"$u$ depends on $v$." Martin's **instability** of a module is
$$
I(v) = \frac{\text{fan-out}(v)}{\text{fan-in}(v) + \text{fan-out}(v)} \in [0, 1],
$$
with $I = 0$ maximally stable (only depended upon) and $I = 1$ maximally unstable (only depends). State
the **Stable Dependencies Principle**: dependencies should point in the direction of decreasing
instability — toward stable modules. In the MMM graph the domain core has $I = 0$ (high fan-in, nothing
below it) and reporting/ingestion have $I = 1$ (leaves that only consume); the prior store is
deliberately stable ($I$ small, high fan-in) because many components rely on it.

**Rung 3 — Proof P (KEYSTONE): acyclicity ⇔ a topological order exists.**
A **topological order** of $G$ is a linear arrangement of $V$ in which every edge points forward (every
dependency is satisfied before the module that needs it). **Theorem:** a finite directed graph admits a
topological order if and only if it is acyclic. **Proof.** ($\Rightarrow$) If a topological order exists,
every edge goes forward in it; a cycle would force at least one edge to point backward, a contradiction —
so $G$ is acyclic. ($\Leftarrow$) If $G$ is acyclic, then it has a vertex of out-degree $0$ — a module
that depends on nothing: otherwise every vertex has an outgoing edge, and following dependencies forward
through finitely many vertices forces a repeat, hence a cycle, contradicting acyclicity. Place that
vertex first, delete it (the remaining graph is still acyclic), and induct; the construction is Kahn's
algorithm and runs in $O(|V| + |E|)$ (Chapter 23). $\blacksquare$
**Architectural meaning:** a topological order *is* a valid build, initialization, and test order —
each module is constructed and unit-tested using only modules already in place. A cycle means **no**
module in the cycle can be built or tested in isolation; they are mutually entangled. Acyclic
dependencies are therefore exactly the condition for independent buildability and testability.

**Rung 4 — The Dependency Rule: data-flow cycles versus dependency cycles.**
The calibration–optimization loop of Part VI appears to cycle: the optimizer flags a high-leverage spend
region, an experiment is run there, calibration folds the result in, the store updates, and the optimizer
re-runs. Written naively as module dependencies — optimizer depends on calibration (to read the
calibrated curve) and calibration depends on optimizer (to know where to aim) — this is a 2-cycle, and
by Rung 3 neither module is independently testable. The resolution is the **Dependency Inversion
Principle**: introduce a stable abstraction both depend on. Both the optimizer and calibration depend on
the **prior-store interface** (both edges point *to* the store); neither depends on the other. The cycle
is broken, the graph is again a DAG, and the store — now with high fan-in — is the stable core ($I
\approx 0$). The lesson stated plainly: a **runtime data/control cycle** (information flowing around the
loop over time) does **not** require a **source-code dependency cycle**; routing interaction through a
stable interface keeps the module graph acyclic. This is Clean Architecture's Dependency Rule
[@martin2017].

**Rung 5 — Layering and synthesis: the MMM architecture.**
Assemble the full system as a layered DAG with dependencies pointing inward toward stability:
a **stable core** (low $I$, high fan-in) of domain abstractions — the response curve $S(x;\theta)$ and
the prior-store schema/interface; a **middle layer** of model fitting, calibration, and optimization that
depends on the core; and an **outer layer** (high $I$) of ingestion adapters, reporting, and
configuration that the core never depends on. The payoff is concrete: because dependencies point inward,
the sampler, the optimizer, or the data source can each be swapped without touching the core or each
other. This is what makes the Part VI calibration–optimization loop a *maintainable* system rather than a
tangle, and it sets up Chapter 25 (the prior store built as a real, versioned data product behind the
interface defined here) and Chapter 26 (testing, made possible by the acyclic, mockable boundaries).

---

## 4. Worked Examples

**WE1 — Instability metrics on the MMM module graph.**
Given the dependency list (domain; ingestion→domain; transforms→domain; store→domain; fit→domain,
transforms; calibration→domain, store; optimizer→domain, store; reporting→fit, optimizer), tabulate
fan-in, fan-out, and $I$ for each module. Results: domain $I = 0.00$ (fan-in 6), store $I = 0.33$,
transforms $I = 0.50$, fit $I = 0.67$, optimizer $I = 0.67$, ingestion/calibration/reporting $I = 1.00$.
Verify the Stable Dependencies Principle — every edge points from a higher-$I$ module to a lower-or-equal
$I$ module, i.e. toward the stable core.

**WE2 — The topological build order.**
Run the constructive proof of Rung 3 (Kahn's algorithm) on the MMM DAG by hand: repeatedly remove a
module whose dependencies are all already placed. One valid order is domain → ingestion → transforms →
store → fit → calibration → optimizer → reporting. Confirm it is a legal build order: every module is
built only after the modules it depends on. All eight modules are placed, certifying acyclicity.

**WE3 — Breaking and fixing the calibration loop's cycle.**
Add the naive loop edges optimizer → calibration and calibration → optimizer. Kahn's algorithm now places
only five of the eight modules (the optimizer, calibration, and reporting that depends on the optimizer
are never freed) — the cycle is detected, and none of the three is independently testable. Apply
Dependency Inversion: drop the two direct edges and route both modules through the store (optimizer →
store, calibration → store, already present). The graph is a DAG again and all eight modules sort. The
runtime loop is unchanged; only the source-code dependencies were inverted.

---

## 5. Code Tie-in

A single runnable `{python}` cell (Python standard library + `matplotlib`; verified headless). **No
third-party graph library.** It:
1. **Instability metrics (WE1):** encode the MMM dependency dictionary; compute fan-in, fan-out, and $I$
   for each module; print the table; `assert` the domain has $I = 0$ and the maximum fan-in, and that
   every edge points toward a module of $I$ no greater than its source (Stable Dependencies Principle).
2. **Topological sort (WE2):** implement Kahn's algorithm; `assert` it returns all eight modules and that
   the order is a valid build order (every dependency precedes its dependent).
3. **Cycle detection (WE3):** add the optimizer ↔ calibration edges; `assert` Kahn now returns fewer than
   eight modules (cycle detected); then route through the store and `assert` it returns all eight again.
4. **Figure:** draw the MMM dependency DAG with `matplotlib` (manual node positions, arrow annotations —
   no `networkx`), and a second panel showing the injected cycle in red.
Every non-figure claim asserted; deterministic (no RNG needed).

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** Why does an acyclic dependency graph permit independent building and testing, while a cycle does
  not? Why should dependencies point toward more stable modules? In what sense is the calibration loop a
  cycle, and in what sense is it not?
- **B:** For a small given dependency list, compute fan-in, fan-out, and the instability $I$ of each
  module, and identify the most and least stable; produce a topological order by hand and state where the
  Stable Dependencies Principle is violated, if anywhere.
- **P:** Prove that a finite directed graph admits a topological order if and only if it is acyclic
  (both directions; the constructive direction is Kahn's algorithm). Then prove that introducing a single
  abstraction node that two mutually-dependent modules both depend on (Dependency Inversion) converts
  their 2-cycle into a DAG.
- **A:** Extend the Code Tie-in — write a function that, given any dependency dictionary, returns a
  topological order or reports the modules involved in a cycle; run it on the MMM graph with and without
  the loop edges, and compute each module's instability to flag any dependency that violates the Stable
  Dependencies Principle.

---

## 7. Appendix solutions

Append `## Chapter 24 — Software Architecture` to `appendix/solutions.qmd`, **in chapter order** (after
the Chapter 23 block), inside the existing `content-visible` gated div. Full C/B/P/A; the P-block carries
both proofs (acyclicity ⇔ topological order; Dependency Inversion breaks a 2-cycle). Part IV airy proof
spacing.

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities/quantities:
instability $I(v) = \text{fan-out}/(\text{fan-in} + \text{fan-out}) \in [0,1]$; the Stable Dependencies
Principle (depend toward lower $I$); the keystone (acyclic ⇔ topological order exists), with build/test
cost $O(|V| + |E|)$; the Dependency Rule (route runtime cycles through a stable interface to keep the
module graph acyclic). Close tying back to Part VI (the calibration loop as an acyclic architecture) and
forward to Chapters 25–26.

---

## 9. Cross-references

- **Chapter 7 (Markov Chains):** directed-graph structure reused for the dependency graph.
- **Chapter 23 (CS Foundations):** the $O(|V| + |E|)$ cost of the topological sort.
- **Chapters 18, 20, 21, 22 (Part VI):** the calibration–optimization loop architected here as a DAG.
- **Chapter 25 (Data Engineering):** the prior store built as a versioned data product behind this
  chapter's interface.
- **Chapter 26 (Testing & Reliability):** testing enabled by the acyclic, mockable boundaries.

---

## 10. Bibliography

Reuse `@martin2017` (Clean Architecture — already present) and `@clrs2009` (topological sort / graph
algorithms — already present). Add:
- `@parnas1972` — Parnas, "On the Criteria To Be Used in Decomposing Systems into Modules"
  (Communications of the ACM, 1972) — information hiding.
The plan confirms exact BibTeX and that all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor verified (done: instability table — domain $0.00$, store $0.33$, transforms
  $0.50$, fit/optimizer $0.67$, ingestion/calibration/reporting $1.00$; topological order over all 8
  modules; injected optimizer↔calibration cycle leaves Kahn placing only 5 of 8).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`); no third-party graph
  library.
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, no `\begin{align}`,
  no `psmallmatrix`, valid citation keys, no banned library names; Part IV `### Rung` headings + airy
  proofs.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.
