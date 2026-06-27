# Literature has moved to `../priorlit`

The reference literature for this project is no longer stored in this repo. It was
consolidated into the shared **`../priorlit`** collection (one level above this repo,
alongside the other zfish projects), so every project draws from a single literature
library.

## Where things went

| Was here | Now in `../priorlit` |
|----------|----------------------|
| `lit/*.pdf` (16 PDFs) | `priorlit/pdf/` — renamed to the `author_year_words` convention |
| `lit/converted/*.md` (PDF→markdown text) | `priorlit/curated/conversions/` (preserved verbatim) |
| `vault/*.md` (Obsidian notes, 30) | `priorlit/curated/vault_p01bruno/` (preserved verbatim) |
| — | `priorlit/parsed/` + `priorlit/vault/notes/` also hold auto-generated, unified text + notes for every PDF |

A copy of this repo's `references.bib` is mirrored at
`priorlit/curated/references_p01bruno.bib` for reference.

## What stayed here

- **`references.bib`** — the manuscript build (`manuscript/make`) cites it directly, so
  it remains in this repo. `build_vault.py` also reads it.

## Filename mapping (PDFs)

`AuthorYYYY_words.pdf` → `author_year_words.pdf`, e.g.
`Audira2020_zebrafish_strains_behavioral.pdf` → `audira_2020_zebrafish_strains_behavioral.pdf`.
Johnson 2025 and Philpott 2012 are shared with the P5 project and are stored once in
`priorlit/pdf/`.
