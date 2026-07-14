# Rendered PDF

This folder holds the **rendered PDF of the book**, *The Mathematics and Engineering
of Marketing Mix Modeling*.

It exists because Quarto's own output directory, `_book/`, is git-ignored (see
`.gitignore`) — anything Quarto renders there is local-only and never committed.
Dropping the PDF here instead keeps a version-controlled copy in the repo and on
GitHub.

## How to generate the PDF

From the repo root:

```bash
quarto render --to pdf
```

Quarto writes the PDF into `_book/` (its configured output dir). Copy or move it here, e.g.:

```bash
cp _book/*.pdf "pdf/"
```

The PDF format is configured in `_quarto.yml` under `format: pdf:` (`documentclass: scrreprt`).
Rendering requires a LaTeX engine (e.g. `quarto install tinytex`).

> Note: `.gitignore` ignores `*.tex` but **not** `*.pdf`, so PDFs placed in this
> folder are tracked normally.
