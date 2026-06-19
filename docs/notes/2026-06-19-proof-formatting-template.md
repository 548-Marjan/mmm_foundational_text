# Proof formatting template (annotated-steps style)

**Decision:** standardize proofs on the *annotated-steps + justification-column* style, so the
mathematics is visually separated from the prose and each algebraic step carries its reason. This
codifies the style several later proofs already use (Ch. 1 normal equations, Ch. 19 OVB, Ch. 17
Kalman) and brings the cramped early-foundations proofs up to it.

Canonical worked example: the steepest-ascent proof in §2.2.1
(`parts/01-foundations/02-multivariable-calculus.qmd`). Before/after below.

---

## The template

A proof is built from these elements, in this order:

1. **Optional shorthand line.** If a long symbol recurs, abbreviate it once up front to keep the
   aligned lines short and the justification column from overflowing:
   `**Proof.** Write $g \equiv \nabla f(\beta)$.`

2. **Logical parts get a short italic sub-label.** When a proof has distinct movements (e.g. an
   existence claim and a uniqueness claim, or a maximizer and an orthogonality claim), introduce
   each with an italic lead-in like `*The maximizing direction.*`. A single-movement proof skips
   these.

3. **Each chain of (in)equalities goes in an `aligned` block with a justification column.** Use
   `&` before the relation symbol and `&&` before the reason. Reasons are wrapped in `\text{...}`
   and kept to ≈4 words:

   ```
   $$
   \begin{aligned}
   D_u f(\beta) &= g^\top u                         && \text{(directional derivative)} \\
                &\le \lVert g\rVert\,\lVert u\rVert  && \text{(Cauchy–Schwarz)} \\
                &= \lVert g\rVert                    && (\lVert u\rVert = 1),
   \end{aligned}
   $$
   ```

4. **Prose is short connective tissue only** — one sentence to set up the block, one to read off
   the conclusion. No key equation stays inline in a sentence; lift it into a display block.

5. **Close with `$\qquad\blacksquare$`** (or `$\blacksquare$`), as today.

### When the chain-of-steps form does not fit

Some proofs are not a chain of equalities (induction, "choose a curve…", contradiction,
measure-theoretic arguments). For those:

- Still split movements with italic sub-labels and keep prose short.
- Lift the **key relation(s)** into standalone display `$$ ... $$` (no justification column needed
  if there is only one line), and state the justification in the sentence that introduces it.
- Do not force an `aligned` block onto a step that has no algebraic siblings.

The justification column is for *derivations*; the display-equation lift is the floor every proof
meets.

### Rendering notes / guardrails

- `aligned` with `&&\text{...}` renders in MathJax (Quarto HTML default) and KaTeX. Keep
  annotations short so they don't overflow on narrow screens or in PDF; abbreviate the subject
  (the shorthand line) rather than writing long reasons.
- Use an en-dash in names (`Cauchy–Schwarz`) consistently.
- Put the trailing comma/period **inside** the aligned block on the last line, as shown.
- One blank line between the prose and the `$$` block.

---

## Before / after (§2.2.1, steepest ascent)

### Before
```
**Proof.** By Cauchy–Schwarz, $\nabla f(\beta)^\top u \le
\lVert\nabla f(\beta)\rVert\,\lVert u\rVert = \lVert\nabla f(\beta)\rVert$, with
equality iff $u$ is a nonnegative multiple of $\nabla f(\beta)$; the unit
maximizer is therefore $u^\star=\nabla f(\beta)/\lVert\nabla f(\beta)\rVert$
(assuming $\nabla f(\beta)\neq 0$). For the orthogonality claim, let
$\gamma(t)$ be any curve lying in the level set $\{x:f(x)=f(\beta)\}$ with
$\gamma(0)=\beta$. Then $f(\gamma(t))$ is constant, so differentiating and
using the chain rule, $0=\tfrac{d}{dt}f(\gamma(t))\big|_{0}
=\nabla f(\beta)^\top\gamma'(0)$. Hence $\nabla f(\beta)$ is orthogonal to every
tangent direction of the level set. $\qquad\blacksquare$
```

### After
```
**Proof.** Write $g \equiv \nabla f(\beta)$.

*The maximizing direction.* For any unit vector $u$,
$$
\begin{aligned}
D_u f(\beta) &= g^\top u                         && \text{(directional derivative)} \\
             &\le \lVert g\rVert\,\lVert u\rVert  && \text{(Cauchy–Schwarz)} \\
             &= \lVert g\rVert                    && (\lVert u\rVert = 1),
\end{aligned}
$$
with equality iff $u$ is a nonnegative multiple of $g$. The unit maximizer is
therefore $u^\star = g/\lVert g\rVert$ (assuming $g \neq 0$).

*Orthogonality to the level set.* Let $\gamma(t)$ be any curve in the level set
$\{x : f(x) = f(\beta)\}$ with $\gamma(0) = \beta$. Then
$$
\begin{aligned}
0 &= \tfrac{d}{dt}\, f(\gamma(t))\big|_{0}  && \text{($f\circ\gamma$ is constant)} \\
  &= g^\top \gamma'(0)                       && \text{(chain rule).}
\end{aligned}
$$
Hence $g$ is orthogonal to every tangent direction $\gamma'(0)$ of the level set. $\qquad\blacksquare$
```

---

## Rollout checklist (for the later whole-text pass)

- [ ] Foundations proofs first (Ch. 1–3) — the most cramped.
- [ ] Per proof: add shorthand if a symbol recurs; split movements with italic sub-labels; convert
      each equality chain to an annotated `aligned` block; pull any remaining inline key equations
      to display; verify it renders (no column overflow).
- [ ] Leave already-clean proofs (Ch. 1 normal equations, Ch. 19 OVB, Ch. 17 Kalman, Ch. 20 DiD/IV,
      Ch. 21/22 calibration) alone unless a justification column clearly improves them.
- [ ] Spot-render a build after each part to catch any `aligned` breakage early.
- [ ] No mathematical content changes — formatting only.
