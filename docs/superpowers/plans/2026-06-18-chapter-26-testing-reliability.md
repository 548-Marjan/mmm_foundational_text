# Chapter 26 — Testing & Reliability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Write Chapter 26 (Testing & Reliability), closing Part VII: establish how to know the MMM stack is correct when there is no oracle, using the book's proved theorems as metamorphic/property-test oracles, with the keystone being the probabilistic detection power of randomized testing.

**Architecture:** Replace the stub body of `parts/07-swe-cs/04-testing-reliability.qmd` with the six-heading chapter template (Motivation → Theory & Proofs → Worked Examples → Code Tie-in → Exercises → Summary), keeping the H1 and adding the anchors line. Append a gated solutions block to `appendix/solutions.qmd`. Add two bib entries. Verify via CI `quarto render`.

**Tech Stack:** Quarto `.qmd` (KaTeX math), single `{python}` cell using `numpy` + Python standard library + `matplotlib` only (no third-party testing library), `references.bib`.

---

## Conventions (enforced every task)

- **Styling = Part IV.** Rungs are `### Rung N — Title` H3 headings (NOT inline `**Rung**` bold). Proofs are airy: each step its own short paragraph leading into a `$$` display on its own lines, blank line before and after each `$$`; `\begin{aligned}` for multi-line. Exemplar (read-only): `parts/04-optimization/01-convexity.qmd`. Proofs end with `$\blacksquare$`.
- **House rules (critical).** `numpy` + Python stdlib + `matplotlib` only. **No `hypothesis`, `pytest`, `networkx`, `pandas`, or any third-party testing/MMM/PPL library** — build the property-test harness directly. NEVER name PyMC-Marketing or any MMM/PPL/sampler/causal library. NEVER `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only, never bare `\begin{align}`; `$$` on their own lines; **even `$$` count per file**. Key identities must be **bulleted**, not a run-on paragraph.
- **Math anchors (NumPy-verified, locked):** $P_{\text{detect}}(p{=}0.05, N{=}100) = 0.9941$ (empirical $\approx 0.9944$, seed 26); bound sample size $p{=}0.01, \delta{=}0.01 \Rightarrow N \ge \lceil \ln(1/\delta)/p \rceil = 461$ (exact requirement $\lceil \ln\delta/\ln(1-p)\rceil = 459$ — state the bound is conservative); $P_{\text{detect}}(p{=}0.01, N{=}100) = 0.634$; WE1 buggy store $70$ vs correct $45$; WE3 at $\sigma{=}0.5, n{=}10^4$: standard error $0.005$, $3\sigma/\sqrt n = 0.015$ (tol $0.002$ flakes, $0.02$ safe).
- **Git identity:** jlh530i / jlh530i@gmail.com. NEVER commit to main; this is a worktree branch → PR → user merges. One complete commit-set per PR; do not push follow-up commits after opening.

---

## Task 0: Rebase onto latest main (MUST run before any writing)

**Files:** none (git only).

This worktree branched off `main` at `5fa44de`, **before** Chapter 25 merged. `appendix/solutions.qmd` here ends at the Chapter 24 block; the Chapter 25 block is on `origin/main` but not here. Rebasing first makes the Chapter 26 solutions block append in correct chapter order (after Ch25) and avoids a `solutions.qmd` conflict.

- [ ] **Step 1: Fetch and rebase**

Run:
```bash
git fetch origin
git rebase origin/main
```
Expected: the `spec(ch26)` commit replays cleanly on top of the Ch25-bearing main. If a conflict surfaces in `appendix/solutions.qmd` (it should not, since the spec commit only touches `docs/`), resolve by keeping both the incoming Ch25 block and our content, then `git rebase --continue`.

- [ ] **Step 2: Verify the Ch25 block is now present**

Run: `grep -n "^## Chapter 25" appendix/solutions.qmd`
Expected: one match (the Chapter 25 — Data Engineering block). This confirms the rebase landed and the appendix is in chapter order, so the Ch26 block (Task 8) appends after it.

---

## Task 1: Front matter + Motivation

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (replace stub; keep H1 `# Testing & Reliability`).

- [ ] **Step 1: Replace the stub callout with the anchors line and Motivation**

Keep the H1 `# Testing & Reliability`. Immediately under it, add the anchors line (italic), then the `## Motivation` section:

```
# Testing & Reliability

*Canonical anchors: property-based testing (Claessen & Hughes); metamorphic testing (Chen et al.).*

## Motivation
```

Motivation prose (2–4 short paragraphs) must establish:
- Parts I–VI produced the models; Part VII grounded computation (Ch.23), structured architecture (Ch.24), made data durable (Ch.25). This chapter answers: **how do we know any of it is correct?**
- The **oracle problem**: a sampler returns draws, an optimizer returns an allocation, the store returns a posterior — none has a closed-form "right answer" to assert against, so the standard `assert f(x) == y` example test has nothing to compare to.
- The driving question (state it): *"You can't assert the sampler returns the 'right' posterior — there is no oracle. So how do you test statistical code, and how many random cases does it take to trust it?"*
- The thesis: **the theorems proved across the book are the test oracles** — metamorphic/property tests check the invariant relations those theorems guarantee, and a probabilistic argument (the keystone) says how many random cases it takes to trust the result.

- [ ] **Step 2: Verify structure**

Run: `grep -nE "^# |^## |Canonical anchors" parts/07-swe-cs/04-testing-reliability.qmd`
Expected: H1 intact, anchors line present, `## Motivation` present.

- [ ] **Step 3: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): motivation + anchors line"
```

---

## Task 2: Theory & Proofs — Rungs 1–2 (oracle problem; metamorphic testing)

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add `## Theory & Proofs` then Rungs 1–2).

- [ ] **Step 1: Add `## Theory & Proofs` and Rung 1**

`### Rung 1 — The oracle problem: example tests versus property tests`

Prose must define:
- **Example test:** fixes one input, asserts one expected output (`assert f(x) == y`); certifies a function only at the enumerated points and presupposes a **test oracle** — a known correct output.
- The **oracle problem:** much of the MMM stack has no oracle (sampler draws, optimizer allocation, store posterior — no closed-form right answer).
- **Property test:** asserts a relation that must hold for *all* inputs — a universally quantified statement checked on many (typically random) cases. State the shift: from "is the output equal to this value?" to "does the output satisfy this invariant?"

- [ ] **Step 2: Add Rung 2 (metamorphic testing)**

`### Rung 2 — Metamorphic testing: the theorems are the oracles`

Prose must define a **metamorphic relation** (a relation between outputs on related inputs that must hold even when neither exact value is known; canonical example "permuting the input does not change the result"), then state the book's central observation — **its proved theorems are exactly such relations**, the strongest possible oracles. Enumerate the relations as a bulleted list:
- **Chapter 25 — permutation invariance:** `fold(log) == fold(permute(log))`.
- **Chapter 25 — idempotency under dedup:** `fold(log + duplicate) == fold(log)`.
- **Chapter 24 — acyclicity:** a topological sort places all modules iff the dependency graph is acyclic.
- **Chapter 22 — EVPI monotonicity:** adding a ridge-aimed calibration factor never increases the EVPI.
- **Chapter 23 — conditioning:** $\kappa(X^\top X) \approx \kappa(X)^2$.

Close with the principle (cite `[@chen1998]`): a proved theorem asserts a relation that *must* hold, so a test of it fails only if the implementation is wrong; a property test is the empirical sampling of that theorem over the input space.

- [ ] **Step 3: Verify**

Run: `grep -nE "### Rung [12]" parts/07-swe-cs/04-testing-reliability.qmd`
Expected: Rung 1 and Rung 2 headings present.

- [ ] **Step 4: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): theory rungs 1-2 (oracle problem, metamorphic testing)"
```

---

## Task 3: Theory & Proofs — Rung 3 KEYSTONE (probabilistic power of randomized testing)

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add Rung 3).

- [ ] **Step 1: Add Rung 3 with the keystone theorem, proof, and corollary**

`### Rung 3 — The probabilistic power of randomized testing`

Frame: random property testing turns "how much testing is enough?" into a probability question.

**Theorem.** A defect causes a property to fail on a fraction $p \in (0,1]$ of the input space; $N$ cases are drawn independently and uniformly. Then the detection probability is

$$
P_{\text{detect}}(N) = 1 - (1-p)^N \;\ge\; 1 - e^{-pN}.
$$

**Proof (airy, Part IV spacing — each step its own short paragraph leading into a display):**

A single random case lands in the failing region with probability $p$ and misses it with probability $1-p$.

By independence, all $N$ cases miss simultaneously with probability

$$
(1-p)^N,
$$

so at least one case hits the failing region with probability $1-(1-p)^N$.

The exponential bound follows from $1-p \le e^{-p}$ for every $p$, hence

$$
(1-p)^N \le e^{-pN}, \qquad\text{so}\qquad P_{\text{detect}}(N) \ge 1 - e^{-pN}. \qquad \blacksquare
$$

**Corollary (sample size).** To detect with confidence at least $1-\delta$ it suffices that $(1-p)^N \le \delta$, i.e.

$$
N \;\ge\; \frac{\ln(1/\delta)}{p}.
$$

Reading: rarer bugs (small $p$) need proportionally more cases; a random property test carries a **quantifiable** detection guarantee an example test cannot. State the anchors in prose: $p=0.05, N=100 \Rightarrow P_{\text{detect}} \approx 0.994$; and $p=0.01, \delta=0.01 \Rightarrow N \ge \ln(100)/0.01 \approx 461$ (note this bound is conservative — the exact requirement $(0.99)^N \le 0.01$ gives $N \ge 459$).

- [ ] **Step 2: Verify even `$$` count and the proof terminator**

Run: `grep -c '\$\$' parts/07-swe-cs/04-testing-reliability.qmd` (must be even) and `grep -n 'blacksquare' parts/07-swe-cs/04-testing-reliability.qmd`.
Expected: even `$$` count; `\blacksquare` present in the Rung 3 proof.

- [ ] **Step 3: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): keystone Rung 3 (detection probability, sample-size rule)"
```

---

## Task 4: Theory & Proofs — Rungs 4–5 (numerical reliability; synthesis)

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add Rungs 4–5).

- [ ] **Step 1: Add Rung 4 (numerical reliability)**

`### Rung 4 — Numerical reliability: tolerances, seeds, and flakiness`

Three hazards, each a short paragraph:
- **Tolerances:** finite precision (Ch.23) makes exact-equality assertions wrong; compare within $|a-b| \le \texttt{atol} + \texttt{rtol}\,|b|$, sized to conditioning (achievable accuracy $\approx \kappa u$).
- **Reproducibility:** a stochastic test must fix its RNG seed so a failure is replayable and a minimal counterexample can be isolated; an unseeded failing property test cannot be debugged.
- **Flakiness:** a Monte-Carlo estimate (e.g. the Ch.22 EVPI) has standard error $\approx \sigma/\sqrt{n}$, so a test asserting the estimate within $\texttt{tol}$ flakes whenever $\texttt{tol} \lesssim 3\sigma/\sqrt{n}$ — not because the code is wrong but because $n$ and $\texttt{tol}$ were sized inconsistently (Ch.3). Fix: size $n$ up or $\texttt{tol}$ to the estimator's noise.

Use a `$$` display for at least the tolerance relation and the standard-error relation (airy spacing). Keep `$$` count even.

- [ ] **Step 2: Add Rung 5 (reliability of the running system + synthesis)**

`### Rung 5 — Reliability of the running system, and synthesis`

Prose must:
- Define reliability as correctness under partial failure; note the Ch.25 prior store already has the two properties that provide it — **idempotent** (dedup makes a retry safe) and **replayable** (posterior rebuilds from the immutable log after a crash) — so at-least-once delivery under failure does not corrupt state.
- Place property tests in the testing landscape: many fast unit/example tests, a layer of property/metamorphic tests guarding invariants, few slow end-to-end tests.
- Close the synthesis: the relations proved across the book now serve double duty as **specifications and as tests**; Part VII ends with the system built correctly (Ch.23), structured soundly (Ch.24), stored durably (Ch.25), and now **verified** to honor its own theorems.
- Forward to **Chapter 27**, which assembles and exercises the whole calibration–optimization system end to end.

- [ ] **Step 3: Verify**

Run: `grep -nE "### Rung [45]" parts/07-swe-cs/04-testing-reliability.qmd` and `grep -c '\$\$' parts/07-swe-cs/04-testing-reliability.qmd` (even).

- [ ] **Step 4: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): theory rungs 4-5 (numerical reliability, synthesis)"
```

---

## Task 5: Worked Examples

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add `## Worked Examples`, WE1–WE3).

- [ ] **Step 1: Add `## Worked Examples` with three subsections**

**WE1 — A metamorphic test catches what an example test misses.** Two folds: correct (dedup by `event_id`) and buggy (no dedup). An **example test** — fold a fixed log of distinct events, assert posterior precision $= 45$ — passes for both (the bug only manifests on a duplicate). The **metamorphic idempotency test** of Ch.25 — fold a log, then fold the same log with one event duplicated, assert the two posteriors equal — passes for the correct store and **fails for the buggy one ($45$ vs $70$)**. (Use the Ch.25 anchors: prior precision $4$, two events $+25$, $+16$ → $45$; the duplicated $+25$ pushes the buggy fold to $70$.)

**WE2 — How many random cases?** The dedup bug triggers only when a randomly generated log contains a duplicate — a fraction $p=0.05$ of generated logs. An example test using distinct events never catches it; $N=100$ random property checks catch it with probability $1 - 0.95^{100} \approx 0.994$. For a rarer defect, $p=0.01$, reaching $99\%$ confidence ($\delta=0.01$) needs $N \ge \ln(100)/0.01 \approx 461$ cases. The sample-size rule converts desired confidence directly into a number of random cases.

**WE3 — A flaky Monte-Carlo test, and how to size it.** The Ch.22 EVPI is estimated by Monte Carlo with $n$ draws; standard error $\approx \sigma/\sqrt{n}$. A regression test asserting the estimate within $\texttt{tol}$ flakes whenever $\texttt{tol} < 3\sigma/\sqrt{n}$. With $\sigma \approx 0.5, n = 10^4$, the standard error is $0.005$, so $\texttt{tol} = 0.002$ flakes while $0.02$ is safe; alternatively raise $n$. Fixing the seed makes a single run deterministic, but the principled fix is to size $n$ and $\texttt{tol}$ to each other.

Use `$$` displays where a relation is stated; keep the count even.

- [ ] **Step 2: Verify**

Run: `grep -nE "## Worked Examples|WE1|WE2|WE3|metamorphic" parts/07-swe-cs/04-testing-reliability.qmd`
Expected: section + three worked examples present.

- [ ] **Step 3: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): worked examples (metamorphic catch, detection numbers, flaky MC)"
```

---

## Task 6: Code Tie-in (single `{python}` cell, headless-verified)

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add `## Code Tie-in` with one `{python}` cell).
- Scratch (not committed): write the cell body to `/tmp/ch26_code.py` and run headless before pasting.

- [ ] **Step 1: Write and verify the code cell headless**

Build the cell to do exactly this (numpy + stdlib + matplotlib only; **no third-party testing library**):
1. **Prior store folds** — `fold(events, dedup=True)` (precision = prior $4$ + sum of distinct `event_id` contributions) and the buggy `dedup=False`. Reuse the Ch.25 contributions ($+25, +16$).
2. **A property-test harness** — `make_log(rng, p)` generates a log of distinct events and, with probability $p$, inserts a duplicate of one event; metamorphic checks `mr_permutation(fold_fn, log)` and `mr_idempotent(fold_fn, log)` each returning pass/fail.
3. **WE1 example vs metamorphic** — assert the example test (fixed distinct log → precision $45$) passes for **both** folds; assert `mr_idempotent` passes for the correct fold and **fails** for the buggy fold (45 vs 70).
4. **WE2/keystone detection probability** — seed `np.random.default_rng(26)`; generate $N=100$ logs at $p=0.05$, count how often `mr_idempotent` catches the buggy store; assert the empirical rate is within tolerance (sized to sampling noise, e.g. `abs(emp - (1-(1-p)**N)) < 0.05`) of $1-(1-p)^N \approx 0.994$.
5. **Figure** — plot $1-(1-p)^N$ vs $N$ for $p \in \{0.01, 0.05, 0.1\}$ with empirical detection points overlaid; end the cell with `plt.show()`.

All asserts use tolerances sized to noise (Rung 4), not exact equality for Monte-Carlo quantities. Verify:
```bash
MPLBACKEND=Agg python3 /tmp/ch26_code.py
```
Expected: prints (correct fold 45, buggy duplicate 70, empirical detection ≈ 0.99), all asserts pass, no exception.

- [ ] **Step 2: Paste verified body into the `{python}` cell**

Wrap exactly the verified body in:
````
## Code Tie-in

```{python}
# ... verified body ...
```
````
Do not edit the body after verification except whitespace.

- [ ] **Step 3: Re-verify the extracted cell**

Re-extract the cell body to `/tmp/ch26_code.py` and re-run `MPLBACKEND=Agg python3 /tmp/ch26_code.py`. Expected: identical clean run.

- [ ] **Step 4: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): code tie-in (property-test harness, detection-rate experiment, figure)"
```

---

## Task 7: Exercises + Summary

**Files:**
- Modify: `parts/07-swe-cs/04-testing-reliability.qmd` (add `## Exercises` and `## Summary`).

- [ ] **Step 1: Add `## Exercises` (C / B / P / A — self-contained, no solution links)**

- **C — Conceptual:** What is the oracle problem, and why does it make a sampler/optimizer hard to test with example tests? Why is a proved theorem the strongest possible test oracle? Why does an unseeded stochastic test that fails resist debugging?
- **B — By Hand:** A defect fails on a fraction $p$ of inputs. Compute $P_{\text{detect}}$ for given $(p,N)$, and the number of random cases needed for target confidence $1-\delta$. Size a tolerance for a Monte-Carlo test of standard error $\sigma/\sqrt{n}$ so it is not flaky.
- **P — Prove It:** (i) Prove $N$ independent uniform cases detect a defect of failing-fraction $p$ with probability $1-(1-p)^N$, and derive $N \ge \ln(1/\delta)/p$. (ii) Show a metamorphic relation guaranteed by a proved theorem is a sound oracle — a test of it fails only if the implementation violates the theorem.
- **A — Applied:** Extend the Code Tie-in — add a third metamorphic relation (Ch.24 acyclicity or Ch.22 EVPI monotonicity) and a defect that violates it; estimate the empirical detection rate vs $N$ and compare to the keystone curve; then make a Monte-Carlo test flaky by tightening its tolerance and fix it by sizing $n$.

Use the four template subheadings (`### C — Conceptual / Reading Comprehension`, `### B — By Hand`, `### P — Prove It`, `### A — Applied / Code`) matching sibling chapters.

- [ ] **Step 2: Add `## Summary` (bulleted Key concepts + bulleted Key identities)**

`## Summary` with two bulleted subsections.
- **Key concepts (bulleted):** oracle problem; example vs property test; metamorphic relation; theorems as oracles; numerical reliability hazards (tolerances/seeds/flakiness); reliability via idempotency + replay; tests double as specifications.
- **Key identities (bulleted — NOT a run-on paragraph):**
  - example test certifies points; property test certifies an invariant over all inputs
  - metamorphic relation as oracle: `fold(log) == fold(permute(log))`, `fold(log+dup) == fold(log)`
  - detection probability $P_{\text{detect}} = 1-(1-p)^N \ge 1 - e^{-pN}$
  - sample-size rule $N \ge \ln(1/\delta)/p$
  - Monte-Carlo flakiness threshold $\texttt{tol} \gtrsim 3\sigma/\sqrt{n}$
  - idempotent + replayable $\Rightarrow$ reliable under retry

Close with one sentence tying back to Ch.22–25 (their theorems as oracles) and forward to Ch.27 (the capstone).

- [ ] **Step 3: Verify structure (six template headings in order) and even `$$`**

Run: `grep -nE "^## (Motivation|Theory|Worked|Code|Exercises|Summary)" parts/07-swe-cs/04-testing-reliability.qmd` and `grep -c '\$\$' parts/07-swe-cs/04-testing-reliability.qmd` (even).
Expected: the six headings appear once each in template order; even `$$` count.

- [ ] **Step 4: Commit**

```bash
git add parts/07-swe-cs/04-testing-reliability.qmd
git commit -m "feat(ch26): exercises (C/B/P/A) + summary"
```

---

## Task 8: Appendix solutions block

**Files:**
- Modify: `appendix/solutions.qmd` (append `## Chapter 26 — Testing & Reliability` after the Chapter 25 block, inside the existing gated div).

- [ ] **Step 1: Locate the insertion point**

Run: `grep -n "^## Chapter 25" appendix/solutions.qmd` (must exist — Task 0 ensured it). Find the end of the Ch.25 block (just before the next `## Chapter` or the closing `:::` of the gated div). The Ch.26 block goes there, in chapter order, **inside** the `content-visible when-meta="show-solutions"` div (do not open a new div).

- [ ] **Step 2: Append the Chapter 26 solutions block**

`## Chapter 26 — Testing & Reliability` with full C / B / P / A solutions:
- **C:** the oracle problem is that statistical/optimization outputs have no closed-form expected value to assert against; example tests need an oracle, so they cannot certify such code beyond trivially. A proved theorem is the strongest oracle because it is a relation that *must* hold for all inputs — a failing test is a genuine bug, never a wrong expectation. An unseeded failing stochastic test resists debugging because the failing input cannot be reproduced or minimized.
- **B:** worked numbers — e.g. $p=0.05, N=100 \Rightarrow 1-0.95^{100} \approx 0.994$; $p=0.01, \delta=0.01 \Rightarrow N \ge \ln(100)/0.01 \approx 461$; tolerance sizing $\texttt{tol} \ge 3\sigma/\sqrt{n}$, e.g. $\sigma=0.5, n=10^4 \Rightarrow \texttt{tol} \ge 0.015$.
- **P:** both proofs in full, Part IV airy spacing — (i) the detection-probability theorem + sample-size corollary (independence ⇒ $(1-p)^N$ miss; $1-p \le e^{-p}$ bound; $(1-p)^N \le \delta \Rightarrow N \ge \ln(1/\delta)/p$); (ii) soundness of a theorem-backed metamorphic relation (if the theorem holds for all valid inputs, an output pair violating the relation witnesses an implementation that does not compute the specified function — the test cannot raise a false alarm). End each proof with `$\blacksquare$`.
- **A:** sketch the extension — a third metamorphic relation (e.g. acyclicity: build a DAG, assert a topological order exists; inject a back-edge defect, assert detection), empirical detection rate vs $N$ tracking $1-(1-p)^N$, and the flaky→fixed Monte-Carlo tolerance demonstration.

- [ ] **Step 3: Verify even `$$` count and gating**

Run: `grep -c '\$\$' appendix/solutions.qmd` (even) and confirm the new block sits before the final `:::` that closes the gated div (`grep -n ':::' appendix/solutions.qmd | tail -3`).

- [ ] **Step 4: Commit**

```bash
git add appendix/solutions.qmd
git commit -m "feat(ch26): appendix solutions (C/B/P/A, both proofs, Part IV spacing)"
```

---

## Task 9: Bibliography + final review + PR

**Files:**
- Modify: `references.bib` (add `@claessen2000`, `@chen1998`).
- Read-only review: `parts/07-swe-cs/04-testing-reliability.qmd`, `appendix/solutions.qmd`.

- [ ] **Step 1: Add the two bib entries**

Append to `references.bib` (verify exact keys resolve against the chapter's `[@claessen2000]` / `[@chen1998]` citations):

```bibtex
@inproceedings{claessen2000,
  author    = {Claessen, Koen and Hughes, John},
  title     = {{QuickCheck}: A Lightweight Tool for Random Testing of {Haskell} Programs},
  booktitle = {Proceedings of the Fifth ACM SIGPLAN International Conference on Functional Programming (ICFP)},
  pages     = {268--279},
  year      = {2000},
  publisher = {ACM}
}

@techreport{chen1998,
  author      = {Chen, T. Y. and Cheung, S. C. and Yiu, S. M.},
  title       = {Metamorphic Testing: A New Approach for Generating Next Test Cases},
  number      = {HKUST-CS98-01},
  institution = {Department of Computer Science, Hong Kong University of Science and Technology},
  year        = {1998}
}
```

- [ ] **Step 2: Verify both keys are cited and defined**

Run:
```bash
grep -nE "claessen2000|chen1998" references.bib
grep -nE "@claessen2000|@chen1998" parts/07-swe-cs/04-testing-reliability.qmd
```
Expected: each key defined once in the bib and cited at least once in the chapter.

- [ ] **Step 3: Full lint pass**

Run the structure/KaTeX checks over the chapter and appendix:
```bash
grep -c '\$\$' parts/07-swe-cs/04-testing-reliability.qmd   # even
grep -c '\$\$' appendix/solutions.qmd                       # even
grep -nE 'begin\{align\}|psmallmatrix' parts/07-swe-cs/04-testing-reliability.qmd appendix/solutions.qmd  # none
grep -niE 'pymc|stan|hypothesis|pytest|networkx|pandas|numpyro|pyro|sqlalchemy|kafka' parts/07-swe-cs/04-testing-reliability.qmd  # none (banned libs)
grep -nE '^## ' parts/07-swe-cs/04-testing-reliability.qmd  # six template headings in order
```
Expected: even `$$` counts; no banned LaTeX; no banned library names; six headings in template order; H1 + anchors line intact.

- [ ] **Step 4: Final headless code re-run**

Re-extract the `{python}` cell to `/tmp/ch26_code.py` and run `MPLBACKEND=Agg python3 /tmp/ch26_code.py`. Expected: clean run, all asserts pass (correct fold 45, buggy 70, empirical detection ≈ 0.99).

- [ ] **Step 5: Commit bib + push branch + open PR**

```bash
git add references.bib
git commit -m "feat(ch26): bibliography (Claessen & Hughes; Chen et al.)"
git push -u origin worktree-ch26-testing-reliability
gh pr create --title "Chapter 26 — Testing & Reliability" --body "$(cat <<'EOF'
## Summary
- Closes Part VII: tests the MMM stack with the book's theorems as metamorphic/property-test oracles.
- Keystone: randomized-testing detection probability $P_{\text{detect}} = 1-(1-p)^N \ge 1-e^{-pN}$ and sample-size rule $N \ge \ln(1/\delta)/p$.
- Worked examples (metamorphic catch 45 vs 70; detection numbers; flaky Monte-Carlo sizing), runnable harness (numpy + stdlib only), C/B/P/A exercises + appendix solutions, two new references.

## Test Plan
- [ ] CI `quarto render` (HTML + PDF) green
- [ ] Code cell runs headless (`MPLBACKEND=Agg python3`)

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 6: Watch CI to green, then report**

Start a background watcher on the PR's render check; report the PR number and watch to a green render conclusion. Do **not** push follow-up commits after opening (the user merges fast — keep one complete commit-set). The user merges the PR.

---

## Self-Review (run after drafting, before execution)

1. **Spec coverage:** every spec section maps to a task — role/motivation (T1), 5 rungs incl. keystone (T2–T4), 3 worked examples (T5), code tie-in with harness + figure (T6), C/B/P/A + summary (T7), appendix (T8), bib + verification (T9), rebase build-note (T0). ✓
2. **Placeholder scan:** no TBD/TODO; every prose section enumerates required content; code cell is specified by behavior + verified anchors. ✓
3. **Anchor consistency:** $P_{\text{detect}}(0.05,100)=0.994$, $N\ge461$ (bound; exact 459), buggy 70 vs correct 45, $3\sigma/\sqrt n=0.015$ at $\sigma{=}0.5,n{=}10^4$ — identical across rungs, worked examples, code, appendix. ✓
4. **House rules:** no banned libraries, even `$$`, bulleted Key identities, Part IV `### Rung` headings + airy proofs, branch→PR→user-merges. ✓
