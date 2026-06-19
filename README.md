# The Mathematics and Engineering of Marketing Mix Modeling

A self-contained textbook that takes a reader from **calculus, linear algebra,
and calculus-based statistics** to **advanced Marketing Mix Modeling (MMM)
practice** — covering both the mathematics and the software engineering /
computer science that a production MMM relies on.

📖 **Read it online:** <https://jlhf80.github.io/mmm_foundational_text/>

The book is built with [Quarto](https://quarto.org) and runs **foundations →
applications**: linear algebra and calculus, regression and Bayesian inference,
MCMC / HMC / NUTS, convex and constrained optimization (with duality and shadow
prices), the MMM itself and its causal calibration, and the engineering that
makes it reliable.

## Structure

| Part | Topic |
|------|-------|
| I | Mathematical Foundations (linear algebra, multivariable calculus, probability & statistics) |
| II | Regression & Bayesian Inference |
| III | Computation: Sampling (Markov chains, MCMC, HMC & NUTS, model checking) |
| IV | Optimization (convexity, linear programming, SLSQP; duality & shadow prices) |
| V | Marketing Mix Modeling: Modeling (adstock, saturation, fitting, DLM / state-space, budget optimization) |
| VI | Causal Grounding & Calibration (causal inference, quasi-experimental design, calibration, prior store) |
| VII | Software Engineering & Computer Science |
| VIII | Capstone (end-to-end worked engagement) |

Each chapter follows the same rhythm: intuition → theory & proofs → worked
examples → runnable code → four tiers of exercises (Conceptual, By hand, Prove
it, Applied/code). Solutions live in a single appendix at the back.

## Building locally

Requires [Quarto](https://quarto.org/docs/get-started/) and Python 3.11+.

```bash
pip install -r requirements.txt

# Reader edition (no solutions)
quarto render

# Instructor edition (with worked solutions)
quarto render --metadata show-solutions:true

# Live preview
quarto preview
```

Output is written to `_book/` (HTML site + PDF).

## Status

In progress, authored part by part. **Parts I–IV** (the prerequisite mathematics —
linear algebra and calculus through Bayesian inference, sampling, and
optimization) are written; **Parts V–VIII** (the MMM itself, causal grounding &
calibration, software engineering, and the capstone) are being authored, with
stub chapters following the standard template (see any `parts/**/**.qmd`).

## A note on code and data

All code is original and self-contained — minimal NumPy, SciPy, and Matplotlib
examples that teach the same patterns a production MMM uses. The book contains
**no proprietary source, firm priors, or client data**.
