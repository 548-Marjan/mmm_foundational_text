# The Mathematics and Engineering of Marketing Mix Modeling

A self-contained textbook that takes a reader from **Calculus 2, linear algebra,
and calculus-based statistics** to **advanced Marketing Mix Modeling (MMM)
practice** — covering both the mathematics and the software engineering /
computer science that a production MMM relies on.

📖 **Read it online:** <https://jlhf80.github.io/mmm_foundational_text/>

The book is built with [Quarto](https://quarto.org) and runs **foundations →
applications**: linear algebra and calculus, regression and Bayesian inference,
MCMC / HMC / NUTS, convex and constrained optimization (with duality and shadow
prices), the MMM synthesis itself, and the engineering that makes it reliable.

## Structure

| Part | Topic |
|------|-------|
| I | Mathematical Foundations (linear algebra, multivariable calculus, probability & statistics) |
| II | Regression & Bayesian Inference |
| III | Computation: Sampling (Markov chains, MCMC, HMC & NUTS, calibration) |
| IV | Optimization (convexity, linear programming, SLSQP; duality & shadow prices) |
| V | Marketing Mix Modeling: Synthesis (adstock, saturation, fitting, budget optimization) |
| VI | Software Engineering & Computer Science |
| VII | Capstone (end-to-end worked engagement) |

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

In progress, authored part-by-part. **Part I:** Chapter 1 (Linear Algebra) and
Chapter 2 (Multivariable Calculus & the Optimization Toolkit) are written; the
remaining chapters are stubs following the same template (see any
`parts/**/**.qmd` for the standard structure).

## A note on code and data

All code is original and self-contained — minimal PyMC / SciPy examples that
teach the same patterns a production MMM uses. The book contains **no
proprietary source, firm priors, or client data**.
