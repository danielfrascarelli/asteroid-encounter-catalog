"""Streamlit dashboard — Gaia DR3 Asteroid Close-Encounter Catalog.

Four tabs:
  1. Encounter Catalog  — full 119 k encounter catalog with filters
  2. Novel Encounters   — 379 science-grade candidates (filtered subset)
  3. Gaia Coverage      — per-encounter transit coverage audit
  4. Mass Candidates    — Cat B candidates + published mass-fit results
"""

from __future__ import annotations

from pathlib import Path

import plotly.express as px
import plotly.graph_objects as go
import polars as pl
import streamlit as st

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Gaia Asteroid Encounters",
    page_icon="🪨",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_CATALOG_PATH = Path("data/output/encounters_characterized.parquet")
_NOVEL_PATH = Path("data/output/relevant_novel_encounters.csv")
_COVERAGE_PATH = Path("data/output/gaia_coverage_audit.csv")
_CANDIDATES_PATH = Path("data/output/mass_candidates.csv")

# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


@st.cache_data(show_spinner="Loading encounter catalog…")
def _load_catalog() -> pl.DataFrame | None:
    if not _CATALOG_PATH.exists():
        return None
    return pl.read_parquet(_CATALOG_PATH)


@st.cache_data(show_spinner="Loading novel encounters…")
def _load_novel() -> pl.DataFrame | None:
    if not _NOVEL_PATH.exists():
        return None
    return pl.read_csv(_NOVEL_PATH)


@st.cache_data(show_spinner="Loading coverage audit…")
def _load_coverage() -> pl.DataFrame | None:
    if not _COVERAGE_PATH.exists():
        return None
    return pl.read_csv(_COVERAGE_PATH)


@st.cache_data(show_spinner="Loading mass candidates…")
def _load_candidates() -> pl.DataFrame | None:
    if not _CANDIDATES_PATH.exists():
        return None
    return pl.read_csv(_CANDIDATES_PATH)


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

st.title("🪨 Gaia DR3 Asteroid Close Encounters")
st.caption(
    "Systematic catalog of asteroid pair separations during the Gaia observation window "
    "(July 2014 – May 2017). "
    "N-body propagation (REBOUND) · KD-tree detection · AL-weighted mass fitting."
)

tab_catalog, tab_novel, tab_coverage, tab_mass = st.tabs([
    "📦 Encounter Catalog",
    "🔭 Novel Encounters",
    "📡 Gaia Coverage",
    "⚖️ Mass Candidates",
])

# ===========================================================================
# TAB 1 — ENCOUNTER CATALOG
# ===========================================================================

with tab_catalog:
    df_all = _load_catalog()
    if df_all is None:
        st.warning("Catalog not found. Run the pipeline first (`scripts/run_pipeline.py`).")
    else:
        # ── Metrics ────────────────────────────────────────────────────────
        n_total = len(df_all)
        n_gaia = int(df_all["gaia_observable"].sum())
        d_min = float(df_all["dist_au"].min())  # type: ignore[arg-type]

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total encounters", f"{n_total:,}")
        c2.metric("Gaia-observable", f"{n_gaia:,}", f"{n_gaia/n_total*100:.0f}%")
        c3.metric("Closest approach", f"{d_min:.5f} AU")
        c4.metric("Time span", "Jul 2014 – May 2017")

        st.divider()

        # ── Sidebar-style filters inside tab ───────────────────────────────
        with st.expander("🔧 Filters", expanded=False):
            fc1, fc2, fc3 = st.columns(3)
            dist_max = fc1.slider(
                "Max separation (AU)",
                min_value=0.001, max_value=0.05,
                value=0.05, step=0.001, format="%.3f",
                key="cat_dist",
            )
            gaia_only = fc2.checkbox("Gaia-observable only", key="cat_gaia")
            all_classes = sorted(
                set(df_all["class_1"].to_list()) | set(df_all["class_2"].to_list())
            )
            sel_classes = fc3.multiselect(
                "Orbit class", all_classes, default=all_classes, key="cat_class"
            )

        df = df_all.filter(pl.col("dist_au") <= dist_max)
        if gaia_only:
            df = df.filter(pl.col("gaia_observable"))
        if sel_classes:
            df = df.filter(
                pl.col("class_1").is_in(sel_classes) | pl.col("class_2").is_in(sel_classes)
            )

        st.caption(f"Showing **{len(df):,}** of {n_total:,} encounters")

        # ── Charts ─────────────────────────────────────────────────────────
        ch1, ch2 = st.columns(2)
        with ch1:
            fig = px.histogram(
                df.to_pandas(), x="dist_au", nbins=80,
                labels={"dist_au": "Min. separation (AU)"},
                color_discrete_sequence=["#4a90d9"],
                title="Distance distribution",
            )
            fig.update_layout(margin=dict(t=40, b=10), height=300, showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

        with ch2:
            fig = px.scatter(
                df.sample(min(10_000, len(df))).to_pandas(),
                x="dist_au", y="rel_vel_km_s",
                color="class_1", opacity=0.4,
                hover_data=["designation_1", "designation_2", "date_utc"],
                labels={
                    "dist_au": "Separation (AU)",
                    "rel_vel_km_s": "Rel. velocity (km/s)",
                    "class_1": "Class",
                },
                title="Velocity vs separation (sample 10 k)",
            )
            fig.update_layout(margin=dict(t=40, b=10), height=300)
            st.plotly_chart(fig, use_container_width=True)

        ch3, ch4 = st.columns(2)
        with ch3:
            # Timeline histogram
            fig = px.histogram(
                df.to_pandas(), x="date_utc", nbins=60,
                labels={"date_utc": "Date"},
                color_discrete_sequence=["#e67e22"],
                title="Encounters per month",
            )
            fig.update_layout(margin=dict(t=40, b=10), height=280, showlegend=False)
            st.plotly_chart(fig, use_container_width=True)

        with ch4:
            # Orbit class breakdown
            class_counts = (
                pl.concat([
                    df.select(pl.col("class_1").alias("class")),
                    df.select(pl.col("class_2").alias("class")),
                ])
                .group_by("class").agg(pl.len().alias("count"))
                .sort("count", descending=True)
            )
            fig = px.pie(
                class_counts.to_pandas(), values="count", names="class", hole=0.4,
                title="Encounters by orbit class",
            )
            fig.update_layout(margin=dict(t=40, b=10), height=280)
            st.plotly_chart(fig, use_container_width=True)

        # ── Table ──────────────────────────────────────────────────────────
        st.subheader("Closest encounters")
        show_cols = [c for c in [
            "number_1", "designation_1", "number_2", "designation_2",
            "date_utc", "dist_au", "dist_km", "rel_vel_km_s",
            "diameter_1_km", "diameter_2_km", "class_1", "class_2",
            "solar_elongation_deg", "gaia_observable",
        ] if c in df.columns]
        top_df = df.sort("dist_au").head(500).select(show_cols)
        st.dataframe(top_df.to_pandas(), use_container_width=True, height=380)

        st.download_button(
            "⬇ Download filtered catalog (CSV)",
            data=df.select(show_cols).to_pandas().to_csv(index=False).encode(),
            file_name="gaia_encounters_filtered.csv",
            mime="text/csv",
        )

# ===========================================================================
# TAB 2 — NOVEL ENCOUNTERS
# ===========================================================================

with tab_novel:
    df_nov = _load_novel()
    if df_nov is None:
        st.warning(
            "Novel encounters file not found: `data/output/relevant_novel_encounters.csv`. "
            "Run `scripts/analyze_mass_candidates.py`."
        )
    else:
        n_nov = len(df_nov)
        n_unk = int(df_nov["mass_unknown"].sum()) if "mass_unknown" in df_nov.columns else 0
        n_gaia_nov = int(df_nov["gaia_observable"].sum()) if "gaia_observable" in df_nov.columns else 0

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Novel candidates", f"{n_nov:,}")
        c2.metric("Unknown-mass perturbers", f"{n_unk:,}")
        c3.metric("Gaia-observable", f"{n_gaia_nov:,}")
        c4.metric("Unique perturbers", str(df_nov["number_1"].n_unique()))

        st.divider()

        # Filters
        with st.expander("🔧 Filters", expanded=True):
            fc1, fc2, fc3 = st.columns(3)
            nov_dist = fc1.slider(
                "Max separation (AU)", 0.001, 0.05, 0.05, 0.001,
                format="%.3f", key="nov_dist",
            )
            nov_vel = fc2.slider(
                "Max rel. velocity (km/s)", 0.5, 20.0, 10.0, 0.5,
                key="nov_vel",
            )
            nov_gaia = fc3.checkbox("Gaia-observable only", key="nov_gaia")

        df_n = df_nov.filter(pl.col("dist_au") <= nov_dist)
        df_n = df_n.filter(pl.col("rel_vel_km_s") <= nov_vel)
        if nov_gaia and "gaia_observable" in df_n.columns:
            df_n = df_n.filter(pl.col("gaia_observable"))

        # Deflection score chart
        ch1, ch2 = st.columns(2)
        with ch1:
            top20 = df_n.sort("deflection_score", descending=True).head(20)
            fig = px.bar(
                top20.to_pandas(),
                x="deflection_score",
                y=top20["designation_1"].to_list(),
                orientation="h",
                color="dist_au",
                color_continuous_scale="Blues_r",
                labels={"deflection_score": "Deflection score", "y": "Perturber"},
                title="Top 20 by deflection score",
            )
            fig.update_layout(
                yaxis=dict(autorange="reversed"),
                margin=dict(t=40, b=10), height=420, coloraxis_colorbar_title="AU",
            )
            st.plotly_chart(fig, use_container_width=True)

        with ch2:
            fig = px.scatter(
                df_n.to_pandas(),
                x="dist_au", y="deflection_score",
                color="rel_vel_km_s",
                size="diameter_1_km",
                hover_data=["designation_1", "designation_2", "date_utc"],
                color_continuous_scale="Viridis",
                log_y=True,
                labels={
                    "dist_au": "Separation (AU)",
                    "deflection_score": "Deflection score",
                    "rel_vel_km_s": "Vel. (km/s)",
                },
                title="Deflection score landscape",
            )
            fig.update_layout(margin=dict(t=40, b=10), height=420)
            st.plotly_chart(fig, use_container_width=True)

        st.subheader(f"All {len(df_n):,} filtered novel encounters")
        show = [c for c in [
            "number_1", "designation_1", "diameter_1_km",
            "number_2", "designation_2", "diameter_2_km",
            "date_utc", "dist_au", "rel_vel_km_s",
            "deflection_score", "gaia_observable",
        ] if c in df_n.columns]
        st.dataframe(
            df_n.sort("deflection_score", descending=True).select(show).to_pandas(),
            use_container_width=True, height=400,
        )

        st.download_button(
            "⬇ Download novel encounters (CSV)",
            data=df_n.select(show).to_pandas().to_csv(index=False).encode(),
            file_name="novel_encounters_filtered.csv", mime="text/csv",
        )

# ===========================================================================
# TAB 3 — GAIA COVERAGE
# ===========================================================================

with tab_coverage:
    df_cov = _load_coverage()
    if df_cov is None:
        st.info(
            "Coverage audit not yet run. Execute:\n\n"
            "```\ndocker compose run --rm pipeline python -m scripts.audit_gaia_coverage\n```"
        )
    else:
        n_viable = int(df_cov["viable_coverage"].sum())
        n_data = int(df_cov["has_gaia_data"].sum())
        n_total_cov = len(df_cov)

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Candidates audited", f"{n_total_cov:,}")
        c2.metric("Have Gaia data", f"{n_data:,}", f"{n_data/n_total_cov*100:.0f}%")
        c3.metric("Viable for LOO fit", f"{n_viable:,}", f"{n_viable/n_total_cov*100:.0f}%")
        c4.metric("Not viable", f"{n_total_cov - n_viable:,}")

        st.divider()

        # Coverage scatter: pre vs post transits, colored by viability
        ch1, ch2 = st.columns(2)
        with ch1:
            fig = px.scatter(
                df_cov.to_pandas(),
                x="n_pre_transits",
                y="n_post_transits",
                color="viable_coverage",
                hover_data=["perturber_name", "target_designation", "date_utc",
                            "nearest_post_days", "dist_au"],
                color_discrete_map={True: "#2ecc71", False: "#e74c3c"},
                labels={
                    "n_pre_transits": "Pre-encounter transits",
                    "n_post_transits": "Post-encounter transits",
                    "viable_coverage": "Viable",
                },
                title="Pre vs post-encounter Gaia transit counts",
            )
            # Viability thresholds
            fig.add_hline(y=3, line_dash="dot", line_color="gray",
                          annotation_text="min post=3")
            fig.add_vline(x=5, line_dash="dot", line_color="gray",
                          annotation_text="min pre=5")
            fig.update_layout(margin=dict(t=40, b=10), height=400)
            st.plotly_chart(fig, use_container_width=True)

        with ch2:
            # Gap to nearest post-encounter transit
            viable_only = df_cov.filter(pl.col("viable_coverage")).drop_nulls(["nearest_post_days"])
            fig = px.histogram(
                viable_only.to_pandas(),
                x="nearest_post_days",
                nbins=40,
                color_discrete_sequence=["#2ecc71"],
                labels={"nearest_post_days": "Days to first post-encounter transit"},
                title="Gap: encounter → first post-encounter transit (viable only)",
            )
            fig.add_vline(x=30, line_dash="dash", line_color="orange",
                          annotation_text="30 d")
            fig.update_layout(margin=dict(t=40, b=10), height=400)
            st.plotly_chart(fig, use_container_width=True)

        # Top viable candidates table
        st.subheader("Top viable candidates (sorted by deflection score)")

        # Join with novel encounters to get mass_unknown flag if available
        df_nov2 = _load_novel()
        if df_nov2 is not None and "mass_unknown" in df_nov2.columns:
            df_cov_aug = df_cov.join(
                df_nov2.select(["number_2", "mass_unknown"]).rename(
                    {"number_2": "target_number"}
                ),
                on="target_number", how="left",
            )
        else:
            df_cov_aug = df_cov.with_columns(pl.lit(True).alias("mass_unknown"))

        show_cov = [c for c in [
            "perturber_number", "perturber_name", "target_designation",
            "date_utc", "dist_au", "deflection_score",
            "n_pre_transits", "n_post_transits",
            "nearest_pre_days", "nearest_post_days",
            "n_window_transits", "viable_coverage", "note",
        ] if c in df_cov_aug.columns]

        viable_df = (
            df_cov_aug
            .filter(pl.col("viable_coverage"))
            .sort("deflection_score", descending=True)
            .select(show_cov)
        )
        st.dataframe(viable_df.to_pandas(), use_container_width=True, height=420)

        with st.expander("Show all (including non-viable)"):
            all_cov = df_cov_aug.sort("deflection_score", descending=True).select(show_cov)
            st.dataframe(all_cov.to_pandas(), use_container_width=True, height=400)

        st.download_button(
            "⬇ Download full coverage audit (CSV)",
            data=df_cov.to_pandas().to_csv(index=False).encode(),
            file_name="gaia_coverage_audit.csv", mime="text/csv",
        )

# ===========================================================================
# TAB 4 — MASS CANDIDATES & FIT RESULTS
# ===========================================================================

with tab_mass:
    df_cand = _load_candidates()

    # ── Published mass-fit results (hardcoded from our LOO runs) ──────────
    st.subheader("Published mass-fit results from this pipeline")

    fit_results = [
        {
            "Perturber": "(111) Ate",
            "Target": "2000 NT3",
            "Date": "2016-06-08",
            "dist (AU)": 0.000472,
            "M_fit (kg)": 5.43e17,
            "M_lit (kg)": 1.76e18,
            "ρ_fit (g/cm³)": 0.51,
            "ρ_lit (g/cm³)": 1.15,
            "χ²_red": 1.03,
            "Verdict": "⚠️ Signal < noise (137-day gap)",
        },
        {
            "Perturber": "(165) Loreley",
            "Target": "1996 TF50",
            "Date": "2014-12-08",
            "dist (AU)": 0.00254,
            "M_fit (kg)": 7.0e17,
            "M_lit (kg)": 7.25e17,
            "ρ_fit (g/cm³)": 2.4,
            "ρ_lit (g/cm³)": "—",
            "χ²_red": 0.82,
            "Verdict": "✅ Consistent with mass estimate",
        },
    ]

    import pandas as pd
    st.dataframe(pd.DataFrame(fit_results), use_container_width=True)

    st.caption(
        "Loreley: mass within 3% of our ρ=1.5 g/cm³ estimate. "
        "Ate: unphysical density — Gaia had a 137-day observation gap around closest approach, "
        "so the signal (expected 1.8 mas) is buried in 22 mas unmodelled perturbation noise."
    )

    st.divider()

    # ── Cat B candidates ───────────────────────────────────────────────────
    if df_cand is None:
        st.info("Mass candidates file not found: `data/output/mass_candidates.csv`.")
    else:
        st.subheader("Cat B: mass determination candidates")

        n_cat = len(df_cand)
        n_viable_cand = int(df_cand["viable"].sum()) if "viable" in df_cand.columns else 0

        c1, c2, c3 = st.columns(3)
        c1.metric("Candidates", n_cat)
        c2.metric("Viable (deflection > threshold)", n_viable_cand)
        c3.metric("Deflection threshold", "100 μas (Gaia single-transit precision)")

        # Deflection bar chart
        show_top = min(25, n_cat)
        df_top = df_cand.head(show_top) if isinstance(df_cand, pl.DataFrame) else df_cand
        if isinstance(df_top, pl.DataFrame):
            df_top_pd = df_top.to_pandas()
        else:
            df_top_pd = df_top

        if "deflection_muas" in df_top_pd.columns and "perturber_name" in df_top_pd.columns:
            fig = px.bar(
                df_top_pd,
                x="deflection_muas",
                y="perturber_name",
                orientation="h",
                color="dist_au",
                color_continuous_scale="RdYlGn_r",
                labels={
                    "deflection_muas": "Expected deflection (μas)",
                    "perturber_name": "Perturber",
                    "dist_au": "Dist (AU)",
                },
                title=f"Expected astrometric deflection — top {show_top} candidates",
            )
            fig.add_vline(x=100, line_dash="dash", line_color="red",
                          annotation_text="Gaia threshold (100 μas)")
            fig.update_layout(
                yaxis=dict(autorange="reversed"),
                margin=dict(t=40, b=10), height=500,
            )
            st.plotly_chart(fig, use_container_width=True)

        # Full table
        show_cols_c = [c for c in [
            "rank", "perturber_name", "perturber_diameter_km",
            "target_designation", "target_diameter_km",
            "date_utc", "dist_au", "rel_vel_km_s",
            "deflection_score", "mass_est_kg", "deflection_muas",
            "gaia_precision_muas", "viable",
        ] if c in (df_cand.columns if isinstance(df_cand, pl.DataFrame) else df_cand.columns)]

        st.subheader("All mass candidates")
        st.dataframe(
            (df_cand.to_pandas() if isinstance(df_cand, pl.DataFrame) else df_cand)[show_cols_c],
            use_container_width=True, height=380,
        )

    # ── Coverage context for top candidates ───────────────────────────────
    df_cov2 = _load_coverage()
    if df_cov2 is not None:
        st.divider()
        st.subheader("Gaia coverage for top 10 Cat B candidates")

        # Best candidates from coverage audit (unknown mass, viable, top deflection)
        df_nov3 = _load_novel()
        if df_nov3 is not None and "mass_unknown" in df_nov3.columns:
            merged = df_cov2.join(
                df_nov3.select(["number_2", "mass_unknown"]).rename({"number_2": "target_number"}),
                on="target_number", how="left",
            )
            top10 = (
                merged
                .filter(pl.col("viable_coverage"))
                .filter(pl.col("mass_unknown").fill_null(True))
                .sort("deflection_score", descending=True)
                .head(10)
            )
        else:
            top10 = df_cov2.filter(pl.col("viable_coverage")).head(10)

        show_top10 = [c for c in [
            "perturber_name", "target_designation", "date_utc", "dist_au",
            "deflection_score", "n_pre_transits", "n_post_transits",
            "nearest_post_days", "viable_coverage",
        ] if c in top10.columns]

        st.dataframe(top10.select(show_top10).to_pandas(), use_container_width=True)

        st.caption(
            "**Recommended next fits**: Alkeste (gap=30d, 76 post obs), "
            "Germania (gap=21d, 204 post obs), Sirona (gap=17d), "
            "Industria (gap=9d, 47 post obs), Urda (gap=11d, 55 post obs)."
        )
