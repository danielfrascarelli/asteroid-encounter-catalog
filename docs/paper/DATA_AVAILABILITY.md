# Data availability — deposit checklist

Steps to satisfy the A&A data policy and fill the `TODO(author)` DOI in
`aa_encounters.tex` (Sect. "Data availability").

> **⚠️ Post-B1 regeneration (2026-07):** the SHA-256 prefixes and row counts below
> refer to the pre-2026-07 freeze. The catalogue is being **regenerated** after the
> refinement-window fix (`refine.py`, tribunal B1), which changes epochs/distances of
> a large fraction of rows. Re-record all SHAs and the exact row count from the new
> `encounters_catalog_rebound_005au_b1fix.parquet` (and its characterised/hybrid
> derivatives) before depositing. The full column dictionary is Appendix A
> (`\ref{tab:schema}`) of the paper.

## What to deposit

| File | Rows | Size | SHA-256 prefix | Role |
|------|-----:|-----:|----------------|------|
| `encounters_catalog_rebound_005au.parquet` | 72,236,904 | ~2.8 GB | `b0272be7…` | Kepler-refined candidate catalogue (primary product) |
| `encounters_catalog_hybrid_stageb.parquet` | 72,236,904 | ~6.1 GB | (record at deposit) | Hybrid: per-row `refinement_method ∈ {kepler, nbody}` + both distances/epochs/velocities |
| `encounters_characterized_full.parquet` | 72,236,904 | ~5.8 GB | (record at deposit) | Adds observability, elongation, magnitudes, diameters, taxonomy |
| provenance sidecar (JSON) | — | small | — | MPCORB hash `3e44e7d3…`, frozen config, pipeline params |

Total ~15 GB. Zenodo accepts this (default 50 GB/record on request); the CDS/VizieR
route may prefer a reduced table plus a link to the full Zenodo record.

## A&A data policy (two tracks, do both if possible)

1. **VizieR / CDS** — A&A expects machine-readable tables at the CDS. For a 72 M-row
   catalogue, submit a described schema + a representative/reduced table (e.g. the
   observable subset, or a per-perturber target list) and reference the full Zenodo
   record for the bulk parquet. Contact CDS (`cds-question@unistra.fr`) with the
   ReadMe. Cite as `J/A+A/VVV/LXX` once assigned.
2. **Zenodo** — mint a DOI for the full parquet bundle. Metadata template:
   [`zenodo_data_deposit.json`](zenodo_data_deposit.json). Fill `creators`
   (name/affiliation/ORCID) and, after acceptance, the paper DOI in
   `related_identifiers`.

## Procedure

1. Freeze the three parquet files (already frozen; see `FROZEN_RUN.md`). Recompute and
   record the SHA-256 of all three:
   ```bash
   sha256sum data/output/encounters_catalog_rebound_005au.parquet \
             data/output/encounters_catalog_hybrid_stageb.parquet \
             data/output/encounters_characterized_full.parquet
   ```
2. Edit `zenodo_data_deposit.json`: creators, and (post-acceptance) the paper DOI.
3. Create the Zenodo record (web UI or API), upload the files, publish → get the DOI.
4. Paste the DOI into `aa_encounters.tex` where the `TODO(author)` in Sect. "Data
   availability" is, e.g. `The catalogue is archived at Zenodo (\doi{10.5281/zenodo.XXXXXXX}).`
5. If using CDS/VizieR, add the VizieR designation alongside the Zenodo DOI.

## Reproducibility note

MPCORB is a living file; the catalogue is bit-reproducible only against the frozen
snapshot `MPCORB_20160217.DAT` (SHA `3e44e7d3…`). The measured diameters/albedos in the
characterised catalogue are a second living input: the JPL SBDB physical-data snapshot
(`data/raw/sbdb_physical.parquet`, download date and SHA-256 in its JSON sidecar). Both
input hashes must be recorded at deposit. Anything derived downstream must cite the
catalogue SHA and the code commit at its generation time.
