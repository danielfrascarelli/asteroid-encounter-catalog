# A&A submission package

LaTeX source for the dataset paper *"A systematic catalogue of real 3D close
encounters between numbered asteroids during the Gaia DR3 window, with a
measured completeness budget."* Content is a faithful conversion of
[`../dataset_paper_draft.md`](../dataset_paper_draft.md) into the Astronomy &
Astrophysics document class.

## Files

- `aa_encounters.tex` — main manuscript (A&A `aa.cls`).
- `references.bib` — BibTeX; volumes/pages verified against ADS/publisher pages.
- `aa_encounters.pdf` — the compiled manuscript (committed for convenience).
- `aa.cls`, `aa.bst`, `linenoaa.sty` — the official A&A macro package
  (v9.4, March 2026), from
  <https://www.aanda.org/doc_journal/instructions/macro/aa/macro-latex-aa.zip>.
  Bundled here for a reproducible build; they remain © EDP Sciences / A&A and are
  provided under the terms of the A&A LaTeX kit (`readme.txt` in that archive).
  `lineno.sty` is not bundled (it ships with TeX Live).
- Figures are pulled from `../figures/` via `\graphicspath`.

## Build

```bash
pdflatex aa_encounters
bibtex   aa_encounters
pdflatex aa_encounters
pdflatex aa_encounters
```

For the actual upload, copy the four `fig*.pdf` from `../figures/` into this
directory and delete the `\graphicspath` line (A&A expects a flat submission).

## Open TODOs before submission (author-owned)

Marked inline in the `.tex` with `TODO(author)`:

1. **Author list, affiliations, ORCID, contact e-mail** (`\author`,
   `\institute`).
2. **Acknowledgements** — funding, software stack (astropy, rebound, ASSIST,
   polars, scipy), Gaia/DPAC boilerplate.
3. **Data availability DOI** — register the catalogue (Zenodo/VizieR) and insert
   the DOI in Sect. "Data availability".
4. **`fuentesmunoz2024` LPSC abstract number** — confirm `55, 2388` against the
   LPSC 2024 (55th) programme on ADS. Now that the full AJ paper
   (`fuentesmunoz2025`) is out, consider whether the LPSC abstract is still
   needed as a separate citation.

## Known content notes

- **Figure 3 needs regeneration.** The on-disk `../figures/fig3_aei_map.pdf`
  predates the fix in `scripts/bench/make_paper_figures.py` that removed the
  duplicated "Inclination *i* (deg)" label (colorbar + marginal). Regenerate:

  ```bash
  docker compose run --rm pipeline bash -c \
    "pip install --quiet matplotlib && python -m scripts.bench.make_paper_figures"
  mv data/output/figures/fig3_aei_map.png data/output/figures/fig3_aei_map.pdf docs/figures/
  ```

- Two figure-only aggregate numbers — **93,010 encountering bodies** and the
  **4.25 km/s median relative velocity** — are computed live by
  `make_paper_figures.py` and are not stated in any prose companion doc. They are
  internally self-consistent across the manuscript; to be fully referee-proof,
  consider recording them in a provenance note alongside the figures.

- References were checked against ADS in 2026-07: Goffin 2014 (A&A 565, A56),
  Michalak 2000 (A&A 360, 363), Park 2016 (Nature 537, 515), Russell 2012
  (Science 336, 684), Tanga 2023 (A&A 674, A12), Vernazza 2020 (Nat. Astron. 4,
  136), Fuentes-Muñoz 2025 (AJ 170, 353 — authors Fuentes-Muñoz, Farnocchia,
  Giorgini & Park; the draft had erroneously copied Scheeres/Tanga, now fixed).
