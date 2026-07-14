# The Selling Playbook

A private, personal Quarto book — a field guide to founder-led selling, synthesized *by sales stage* from six books:

1. **Fanatical Prospecting** (Jeb Blount) — build pipeline
2. **SPIN Selling** (Neil Rackham) — run discovery
3. **The Challenger Sale** (Dixon & Adamson) — create demand through insight
4. **Influence** (Robert Cialdini) — understand buying psychology
5. **Let's Get Real or Let's Not Play** (Khalsa & Illig) — sell as a trusted advisor
6. **Never Split the Difference** (Chris Voss) — close and negotiate

Plus a curated founder-led sales resource list and a weekly operating cadence.

> **A pre-rendered PDF is committed at [`The-Selling-Playbook.pdf`](The-Selling-Playbook.pdf)** so you can read it straight from GitHub without building anything. Re-render it (see below) after editing any chapter to keep it current.

## This is not published online

`sales/` is a **separate Quarto project** from the root MMM book (it has its own `_quarto.yml`). The root `quarto render` and the Render Book CI workflow do **not** descend into it, so this book is **never built in CI and never deployed to GitHub Pages**. The source is committed to the private repo; you render it locally when you want to read it.

## Build it locally

From this `sales/` directory:

```bash
# HTML book (opens as a navigable mini-site)
quarto render --to html
open _book/index.html          # macOS

# PDF (single file)
quarto render --to pdf
open _book/The-Selling-Playbook.pdf

# Both formats at once
quarto render

# Live-reloading preview while editing
quarto preview
```

The PDF build needs a LaTeX engine. If you don't have one, install TinyTeX once:

```bash
quarto install tinytex
```

## Layout

```
sales/
├── _quarto.yml   # book config (HTML + PDF); its presence makes this a separate project
├── index.qmd     # preface
└── chapters/
    ├── 00-cheatsheet.qmd          # the one-page system
    ├── 01-pipeline.qmd            # Fanatical Prospecting
    ├── 02-discovery.qmd           # SPIN
    ├── 03-demand.qmd              # Challenger
    ├── 04-psychology.qmd          # Influence (cross-cutting lens)
    ├── 05-trusted-advisor.qmd     # Let's Get Real
    ├── 06-closing-negotiation.qmd # Never Split the Difference
    ├── 07-founder-resources.qmd
    └── 08-operating-cadence.qmd
```

Build output (`_book/`, `.quarto/`) is gitignored.
