# Figure / aggregate-number provenance

Every headline number in the manuscript traces to a companion doc or a script. Two
aggregate figures are computed live by `scripts/bench/make_paper_figures.py` and were
not previously stated in any prose companion doc; the numerical-consistency audit
flagged them. Recorded here so a referee has a citable source.

| Quantity | Value | How it is computed | Where cited |
|----------|------:|--------------------|-------------|
| Encountering bodies (Fig. 3 population) | **93,010** | `figure3_aei_map`: `data/raw/gaia_orbits.parquet`, drop nulls, filter to the main-belt frame `a ∈ (1.5, 4.5) AU`, `e ∈ [0, 0.5)`, `i < 40°`. Row count of that filtered set. | §2.4 text and Fig. 3 caption |
| Median relative velocity at closest approach | **4.25 km s⁻¹** (4.2548) | `figure2_relvel_hist`: streaming `rel_vel_km_s.median()` over `data/output/encounters_characterized_full.parquet` (72,236,904 rows), `polars` lazy + `engine="streaming"`. | §2.4 text and Fig. 2 caption |

Both were re-derived on 2026-07-04 against the frozen catalogue (SHA `b0272be7…`) and
agree with the values in the manuscript. Reproduce:

```bash
docker compose run --rm pipeline python -c "
import polars as pl
df = pl.read_parquet('data/raw/gaia_orbits.parquet', columns=['a_au','e','i_deg']).drop_nulls()
df = df.filter((pl.col('a_au')>1.5)&(pl.col('a_au')<4.5)&(pl.col('e')>=0.0)&(pl.col('e')<0.5)&(pl.col('i_deg')<40.0))
print('bodies', df.height)
print('median v_rel', (pl.scan_parquet('data/output/encounters_characterized_full.parquet')
    .select(pl.col('rel_vel_km_s').median()).collect(engine='streaming').item()))
"
```

All other figures: Fig. 1 (separation histogram) and Fig. 4 (threshold censoring) are
streaming/loaded aggregations described in `make_paper_figures.py`; Fig. 4's numbers
are also in `docs/kepler_threshold_bias_paper.md`.
