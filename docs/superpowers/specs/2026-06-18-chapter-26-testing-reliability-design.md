# Chapter 26 — Testing & Reliability: Design Spec

**Date:** 2026-06-18
**Part:** VII — Software Engineering & Computer Science (fourth chapter; closes Part VII)
**File:** `parts/07-swe-cs/04-testing-reliability.qmd`
**Status:** Approved design → ready for implementation plan

---

## 1. Role in the book

This chapter tests the system the book has built. Parts I–VI produced the models; Part VII grounded their
computation (Chapter 23), structured their architecture (Chapter 24), and made their data durable
(Chapter 25). This chapter establishes how to know any of it is correct — and the difficulty specific to
statistical and optimization code is the **oracle problem**: there is no closed-form expected output to
assert a sampler, an optimizer, or the prior-store fold against. The chapter's thesis is that **the
theorems proved across the book are the test oracles**: metamorphic and property-based testing check the
invariant relations those theorems guarantee, and a probabilistic argument says how many random cases it
takes to trust the result.

**Center of gravity (locked, via AskUserQuestion):** metamorphic / property-based testing, with the
book's proven theorems as the oracles, and the **keystone** being the probabilistic power of randomized
testing — a defect on a fraction $p$ of inputs is caught by $N$ random cases with probability
$1 - (1-p)^N$.

**Driving question (locked):** *"You can't assert the sampler returns the 'right' posterior — there is no
oracle. So how do you test statistical code, and how many random cases does it take to trust it?"*

**Anchors (the stub has no anchors line; add one):**
`*Canonical anchors: property-based testing (Claessen & Hughes); metamorphic testing (Chen et al.).*`

**Styling (locked, user directive):** match Part IV (`parts/04-optimization/01-convexity.qmd`) — rungs
as `### Rung N — Title` H3 headings, airy step-by-step proofs (each step its own short paragraph leading
into a `$$` display on its own lines; `\begin{aligned}` for multi-line).

---

## 2. Scope discipline

- Tests components already built (the prior store, the optimizer, the dependency graph, the conditioning
  results); introduces no new modeling.
- Does **not** re-prove the relations being tested (Chapters 22–25 own them); it *uses* them as oracles.
- **House rule (critical):** `numpy` + Python standard library + `matplotlib` only. **No `hypothesis`,
  `pytest`, or any third-party testing library** — the property-test harness (random input generation,
  metamorphic checks, detection-rate estimation) is implemented directly. Never name PyMC-Marketing or
  any MMM/PPL/sampler library.
- **Never** `\begin{psmallmatrix}`. KaTeX: `aligned` inside `$$ … $$` only; `$$` on their own lines;
  even `$$` count.

---

## 3. Theory & Proofs — the rungs (5; keystone at Rung 3)

**Rung 1 — The oracle problem: example tests versus property tests.**
An **example test** fixes one input and asserts one expected output — the standard `assert f(x) == y`. It
certifies a function only at the points enumerated, and it presupposes a **test oracle**: a known correct
output to compare against. For much of the MMM stack there is none — a sampler returns draws, an
optimizer returns an allocation, the store returns a posterior, and none has a closed-form "right answer"
to assert against (the **oracle problem**). A **property test** sidesteps the missing oracle by asserting
a relation that must hold for *all* inputs — a universally quantified statement checked on many
(typically random) cases rather than a single point. Define both, and state the shift: from "is the
output equal to this value?" to "does the output satisfy this invariant?"

**Rung 2 — Metamorphic testing: the theorems are the oracles.**
A **metamorphic relation** is a relation between the outputs of a program on related inputs that must hold
even when neither output's exact value is known — the canonical example being "permuting the input does
not change the result." Metamorphic testing checks such relations, and it is the natural fit for code
with no oracle. The crucial observation for this book: **its proved theorems are exactly such relations**,
and therefore the strongest possible oracles. Enumerate the metamorphic relations the book supplies:
- **Chapter 25 — permutation invariance:** `fold(log) == fold(permute(log))`.
- **Chapter 25 — idempotency under dedup:** `fold(log + duplicate) == fold(log)`.
- **Chapter 24 — acyclicity:** a topological sort places all modules iff the dependency graph is acyclic.
- **Chapter 22 — EVPI monotonicity:** adding a ridge-aimed calibration factor never increases the EVPI.
- **Chapter 23 — conditioning:** $\kappa(X^\top X) \approx \kappa(X)^2$.
State the principle: a proved theorem asserts a relation that *must* hold, so a test of it can only fail
if the implementation is wrong; a property test is the empirical sampling of that theorem over the input
space. [@chen1998]

**Rung 3 — Proof P (KEYSTONE): the probabilistic power of randomized testing.**
Random property testing turns "how much testing is enough?" into a probability question. **Theorem.**
Suppose a defect causes a property to fail on a fraction $p \in (0, 1]$ of the input space, and $N$ test
cases are drawn independently and uniformly. Then the probability the defect is detected is
$$
P_{\text{detect}}(N) = 1 - (1-p)^N \;\ge\; 1 - e^{-pN}.
$$
**Proof.** A single random case lands in the failing region with probability $p$ and misses it with
probability $1 - p$. By independence, all $N$ cases miss with probability $(1-p)^N$, so at least one hits
with probability $1 - (1-p)^N$. The bound follows from $1 - p \le e^{-p}$ for all $p$, hence
$(1-p)^N \le e^{-pN}$. $\blacksquare$
**Corollary (sample size).** To detect with confidence at least $1 - \delta$ it suffices that
$(1-p)^N \le \delta$, i.e.
$$
N \;\ge\; \frac{\ln(1/\delta)}{p},
$$
using the same bound. State the reading: rarer bugs (small $p$) need proportionally more cases, and a
property test with random generation carries a *quantifiable* detection guarantee an example test cannot.
Anchors: $p = 0.05, N = 100 \Rightarrow P_{\text{detect}} \approx 0.994$; $p = 0.01, \delta = 0.01
\Rightarrow N \ge 461$.

**Rung 4 — Numerical reliability: tolerances, seeds, and flakiness.**
Three reliability hazards specific to numerical and stochastic tests. **Tolerances:** finite precision
(Chapter 23) makes exact-equality assertions wrong; compare within a tolerance $|a - b| \le
\texttt{atol} + \texttt{rtol}\,|b|$, sized to the conditioning of the computation (achievable accuracy
$\approx \kappa\, u$). **Reproducibility:** a stochastic test must fix its RNG seed so a failure is
replayable and a minimal counterexample can be isolated; an unseeded property test that fails cannot be
debugged. **Flakiness:** a Monte-Carlo estimate of a quantity (e.g. the Chapter 22 EVPI) has standard
error $\approx \sigma/\sqrt{n}$, so a test asserting the estimate is within a tolerance $\texttt{tol}$
fails intermittently whenever $\texttt{tol}$ is set below roughly $3\sigma/\sqrt{n}$ — the test is flaky
not because the code is wrong but because $n$ and $\texttt{tol}$ were sized inconsistently (Chapter 3).
The fix is to size $n$ up or $\texttt{tol}$ to the estimator's own noise.

**Rung 5 — Reliability of the running system and synthesis.**
Reliability is correctness under partial failure, and the prior store of Chapter 25 already has the two
properties that provide it: operations are **idempotent** (dedup makes a retry safe) and **replayable**
(the posterior rebuilds from the immutable log after a crash), so at-least-once delivery under failure
does not corrupt state. Place property tests in the testing landscape — many fast unit/example tests, a
layer of property/metamorphic tests guarding the invariants, few slow end-to-end tests — and close the
synthesis: the relations proved across the book now serve double duty as **specifications and as tests**,
so Part VII ends with the system not only built correctly (Chapter 23), structured soundly (Chapter 24),
and stored durably (Chapter 25) but **verified** to honor its own theorems. Forward to Chapter 27, which
assembles and exercises the whole calibration–optimization system end to end.

---

## 4. Worked Examples

**WE1 — A metamorphic test catches what an example test misses.**
Two implementations of the prior-store fold: a correct one that deduplicates by `event_id`, and a buggy
one that does not. An **example test** — fold a fixed log of distinct events, assert the posterior
precision equals $45$ — **passes for both**, because the bug only manifests on a duplicate. The
**metamorphic test** of Chapter 25's idempotency relation — fold a log, then fold the same log with one
event duplicated, assert the two posteriors are equal — **passes for the correct store and fails for the
buggy one** ($45$ versus $70$). The theorem, used as an oracle, detects a defect the example test cannot
see.

**WE2 — How many random cases? the detection probability in numbers.**
Suppose the dedup bug is triggered only when a randomly generated log happens to contain a duplicate, a
fraction $p = 0.05$ of generated logs. A single example test that happens to use distinct events
($p$-region missed) never catches it; $N = 100$ random property checks catch it with probability
$1 - 0.95^{100} \approx 0.994$. For a rarer defect, $p = 0.01$, reaching $99\%$ confidence
($\delta = 0.01$) requires $N \ge \ln(100)/0.01 \approx 461$ cases. The sample-size rule converts a
desired confidence directly into a number of random cases.

**WE3 — A flaky Monte-Carlo test, and how to size it.**
The Chapter 22 EVPI is estimated by Monte Carlo with $n$ draws; the estimate has standard error
$\approx \sigma/\sqrt{n}$. A regression test asserting the estimate is within $\texttt{tol}$ of a target
EVPI will pass or fail at random — flaky — whenever $\texttt{tol} < 3\sigma/\sqrt{n}$, because the
estimator's own noise exceeds the tolerance. With (say) $\sigma \approx 0.5$ and $n = 10^4$, the standard
error is $0.005$, so a tolerance of $0.002$ flakes while $0.02$ is safe; alternatively raise $n$ to
shrink the standard error. Fixing the seed makes any single run deterministic, but the principled fix is
to size $n$ and $\texttt{tol}$ to each other.

---

## 5. Code Tie-in

A single runnable `{python}` cell (`numpy` + Python standard library; `matplotlib` for the figure;
verified headless). **No third-party testing library** (`hypothesis`, `pytest`) — the harness is built
directly. It:
1. **A property-test harness:** a random-log generator (`rng`, with a tunable probability of inserting a
   duplicate) and two metamorphic checks — permutation invariance and idempotency under dedup — each a
   function returning pass/fail on a generated log.
2. **Correct versus buggy store (WE1):** run an example test (fixed distinct log → precision $45$) on both
   the dedup and no-dedup folds and `assert` both pass; run the idempotency metamorphic test on both and
   `assert` the correct store passes and the buggy store fails.
3. **Detection probability (WE2 / keystone):** generate $N$ random logs at duplicate-fraction $p$, count
   how often the metamorphic test catches the buggy store, and `assert` the empirical detection rate is
   close to $1 - (1-p)^N$.
4. **Figure:** the detection-probability curve $1 - (1-p)^N$ versus $N$ for a few $p \in \{0.01, 0.05,
   0.1\}$, with the empirical detection points overlaid.
Every numeric claim asserted; seed `np.random.default_rng(26)`; Monte-Carlo asserts use tolerances sized
to the sampling noise (Rung 4), not exact equality.

---

## 6. Exercises (C / B / P / A — self-contained, no inline solution links)

- **C:** What is the oracle problem, and why does it make a sampler or optimizer hard to test with example
  tests? Why is a proved theorem the strongest possible test oracle? Why does an unseeded stochastic test
  that fails resist debugging?
- **B:** A defect fails on a fraction $p$ of inputs. Compute the detection probability for given $(p, N)$,
  and the number of random cases needed for a target confidence $1 - \delta$. Size a tolerance for a
  Monte-Carlo test of standard error $\sigma/\sqrt{n}$ so that it is not flaky.
- **P:** (i) Prove that $N$ independent uniform random cases detect a defect of failing-fraction $p$ with
  probability $1 - (1-p)^N$, and derive the sample-size rule $N \ge \ln(1/\delta)/p$. (ii) Show that a
  metamorphic relation guaranteed by a proved theorem is a sound oracle — a test of it fails only if the
  implementation violates the theorem.
- **A:** Extend the Code Tie-in — add a third metamorphic relation (e.g. Chapter 24 acyclicity or Chapter
  22 EVPI monotonicity) and a defect that violates it; estimate the empirical detection rate as a function
  of $N$ and compare it to the keystone curve; then make a Monte-Carlo test flaky by tightening its
  tolerance and fix it by sizing $n$.

---

## 7. Appendix solutions

Append `## Chapter 26 — Testing & Reliability` to `appendix/solutions.qmd`, **in chapter order** (after
the Chapter 25 block), inside the existing `content-visible` gated div. Full C/B/P/A; the P-block carries
both proofs (the detection-probability theorem + sample-size rule; metamorphic relation as a sound
oracle). Part IV airy proof spacing.

---

## 8. Summary (auto-included)

Bulleted **Key concepts** + bulleted **Key identities** (inline math, bulleted). Identities/quantities:
example vs property test; metamorphic relation as oracle; detection probability $P_{\text{detect}} =
1 - (1-p)^N \ge 1 - e^{-pN}$; sample-size rule $N \ge \ln(1/\delta)/p$; Monte-Carlo flakiness threshold
$\texttt{tol} \gtrsim 3\sigma/\sqrt{n}$; idempotent + replayable ⇒ reliable under retry. Close tying back
to Chapters 22–25 (their theorems as oracles) and forward to Chapter 27 (the capstone).

---

## 9. Cross-references

- **Chapter 3 (Probability & Statistics):** independence and the Monte-Carlo standard error.
- **Chapter 22 (Prior Store):** EVPI monotonicity as a metamorphic relation.
- **Chapter 23 (CS Foundations):** finite precision and conditioning, sizing test tolerances.
- **Chapter 24 (Software Architecture):** acyclicity as a metamorphic relation.
- **Chapter 25 (Data Engineering):** permutation invariance and idempotency as metamorphic relations;
  replay and dedup as reliability mechanisms.
- **Chapter 27 (Capstone):** the whole system, exercised end to end.

---

## 10. Bibliography

The chapter has no existing testing references; add two:
- `@claessen2000` — Claessen & Hughes, "QuickCheck: A Lightweight Tool for Random Testing of Haskell
  Programs" (ICFP 2000) — property-based testing.
- `@chen1998` — Chen, Cheung & Yiu, "Metamorphic Testing: A New Approach for Generating Next Test Cases"
  (technical report, 1998) — metamorphic testing.
The plan confirms exact BibTeX and that all keys resolve.

---

## 11. Verification (the real gate)

- Every numeric anchor verified (done: $P_{\text{detect}}(p{=}0.05, N{=}100) = 0.9941$, empirical
  $\approx 0.9944$; sample size $p{=}0.01, \delta{=}0.01 \Rightarrow N \ge 461$ (exact $459$); WE1 buggy
  store $70$ vs correct $45$; WE3 standard error $\sigma/\sqrt{n}$ sizing).
- The single code cell runs top-to-bottom headless (`MPLBACKEND=Agg python3`); no third-party testing
  library.
- KaTeX/structure lint: even `$$` count, six template headings in order, H1 intact, anchors line added,
  no `\begin{align}`, no `psmallmatrix`, valid citation keys, no banned library names; Part IV
  `### Rung` headings + airy proofs.
- **CI `quarto render` (HTML + PDF) green on the PR** — the gate. User merges.

---

## 12. Build note

This worktree branched off `main` before Chapter 25 merged into the appendix. **Before building**, rebase
onto the latest `main` (which now contains the Chapter 24 and Chapter 25 appendix blocks) so the Chapter
26 solutions block appends in correct chapter order and avoids a `solutions.qmd` conflict.
