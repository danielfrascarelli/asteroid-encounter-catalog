"""Streamlit dashboard for the Gaia asteroid close-encounter catalog."""

from __future__ import annotations

from pathlib import Path

import plotly.express as px
import polars as pl
import streamlit as st

from src.catalog.query import filter_encounters, top_encounters

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

_CATALOG_PATH = Path("data/output/encounters_characterized.parquet")
_PAGE_ROWS = 500

st.set_page_config(
    page_title="Gaia Asteroid Encounters",
    page_icon="🪨",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


@st.cache_data(show_spinner="Loading catalog…")
def _load_catalog() -> pl.DataFrame:
    if not _CATALOG_PATH.exists():
        st.error(f"Catalog not found: {_CATALOG_PATH}. Run the pipeline first.")
        st.stop()
    return pl.read_parquet(_CATALOG_PATH)


df_all = _load_catalog()

# Derived ranges for sliders
dist_min_data = float(df_all["dist_au"].min())  # type: ignore[arg-type]
dist_max_data = float(df_all["dist_au"].max())  # type: ignore[arg-type]
date_min = str(df_all["date_utc"].min())
date_max = str(df_all["date_utc"].max())

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _date(s: str):  # type: ignore[no-untyped-def]
    from datetime import date

    return date.fromisoformat(s)


# ---------------------------------------------------------------------------
# Sidebar — filters
# ---------------------------------------------------------------------------

st.sidebar.title("🔭 Filters")

dist_range = st.sidebar.slider(
    "Min. separation (AU)",
    min_value=dist_min_data,
    max_value=dist_max_data,
    value=(dist_min_data, dist_max_data),
    step=0.0001,
    format="%.4f",
)

all_classes = sorted(set(df_all["class_1"].to_list()) | set(df_all["class_2"].to_list()))
selected_classes = st.sidebar.multiselect(
    "Orbit class (body 1 or 2)",
    options=all_classes,
    default=all_classes,
)

date_range = st.sidebar.date_input(
    "Date range",
    value=(
        _date(date_min),
        _date(date_max),
    ),
    min_value=_date(date_min),
    max_value=_date(date_max),
)

gaia_only = st.sidebar.checkbox("Gaia-observable only", value=False)

st.sidebar.divider()
st.sidebar.markdown(
    f"**Total encounters:** {len(df_all):,}  \n"
    f"**Gaia-observable:** {int(df_all['gaia_observable'].sum()):,}"
)


# ---------------------------------------------------------------------------
# Apply filters
# ---------------------------------------------------------------------------

date_start_str = date_range[0].isoformat() if len(date_range) >= 1 else date_min
date_end_str = date_range[1].isoformat() if len(date_range) >= 2 else date_max

df = filter_encounters(
    df_all,
    min_dist_au=dist_range[0],
    max_dist_au=dist_range[1],
    date_start=date_start_str,
    date_end=date_end_str,
    gaia_observable_only=gaia_only,
)

# Class filter (body 1 OR body 2 in selected classes)
if selected_classes:
    df = df.filter(
        pl.col("class_1").is_in(selected_classes) | pl.col("class_2").is_in(selected_classes)
    )

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

st.title("🪨 Gaia DR3 Asteroid Close Encounters")
st.caption(
    "Systematic catalog of asteroid pair separations during the Gaia observation window "
    "(July 2014 – May 2017). Kepler 2-body propagation, 1-hour time step."
)

col1, col2, col3, col4 = st.columns(4)
col1.metric("Filtered encounters", f"{len(df):,}")
col2.metric("Gaia-observable", f"{int(df['gaia_observable'].sum()):,}")
col3.metric(
    "Closest approach",
    f"{float(df['dist_au'].min()):.5f} AU" if len(df) else "—",  # type: ignore[arg-type]
)
col4.metric(
    "Fastest encounter",
    f"{float(df['rel_vel_km_s'].max()):.2f} km/s" if len(df) else "—",  # type: ignore[arg-type]
)

st.divider()

# ---------------------------------------------------------------------------
# Plots — row 1: histogram + scatter
# ---------------------------------------------------------------------------

c1, c2 = st.columns(2)

with c1:
    st.subheader("Distance distribution")
    fig_hist = px.histogram(
        df.to_pandas(),
        x="dist_au",
        nbins=60,
        labels={"dist_au": "Min. separation (AU)"},
        color_discrete_sequence=["#4a90d9"],
    )
    fig_hist.update_layout(margin=dict(t=10, b=10), height=320)
    st.plotly_chart(fig_hist, use_container_width=True)

with c2:
    st.subheader("Velocity vs separation")
    fig_scatter = px.scatter(
        df.to_pandas(),
        x="dist_au",
        y="rel_vel_km_s",
        color="class_1",
        hover_data=["designation_1", "designation_2", "date_utc"],
        labels={
            "dist_au": "Min. separation (AU)",
            "rel_vel_km_s": "Relative velocity (km/s)",
            "class_1": "Class (body 1)",
        },
        opacity=0.5,
    )
    fig_scatter.update_layout(margin=dict(t=10, b=10), height=320)
    st.plotly_chart(fig_scatter, use_container_width=True)

# ---------------------------------------------------------------------------
# Plots — row 2: diameter scatter + observable pie
# ---------------------------------------------------------------------------

c3, c4 = st.columns(2)

with c3:
    st.subheader("Encounter size: body 1 vs body 2")
    fig_diam = px.scatter(
        df.drop_nulls(["diameter_1_km", "diameter_2_km"]).to_pandas(),
        x="diameter_1_km",
        y="diameter_2_km",
        color="gaia_observable",
        hover_data=["designation_1", "designation_2", "dist_au"],
        labels={
            "diameter_1_km": "Diameter body 1 (km)",
            "diameter_2_km": "Diameter body 2 (km)",
            "gaia_observable": "Gaia-observable",
        },
        log_x=True,
        log_y=True,
        opacity=0.5,
        color_discrete_map={True: "#2ecc71", False: "#95a5a6"},
    )
    fig_diam.update_layout(margin=dict(t=10, b=10), height=320)
    st.plotly_chart(fig_diam, use_container_width=True)

with c4:
    st.subheader("Encounters by orbit class")
    class_counts = (
        pl.concat(
            [
                df.select(pl.col("class_1").alias("class")),
                df.select(pl.col("class_2").alias("class")),
            ]
        )
        .group_by("class")
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )
    fig_pie = px.pie(
        class_counts.to_pandas(),
        values="count",
        names="class",
        hole=0.4,
    )
    fig_pie.update_layout(margin=dict(t=10, b=10), height=320)
    st.plotly_chart(fig_pie, use_container_width=True)

# ---------------------------------------------------------------------------
# Top encounters table
# ---------------------------------------------------------------------------

st.divider()
st.subheader("Top encounters (closest approach)")

display_cols = [
    "number_1",
    "number_2",
    "designation_1",
    "designation_2",
    "date_utc",
    "dist_au",
    "dist_km",
    "rel_vel_km_s",
    "diameter_1_km",
    "diameter_2_km",
    "class_1",
    "class_2",
    "solar_elongation_deg",
    "gaia_observable",
]
present_cols = [c for c in display_cols if c in df.columns]

top = top_encounters(df, n=_PAGE_ROWS, by="dist_au", ascending=True).select(present_cols)
st.dataframe(top.to_pandas(), use_container_width=True, height=400)

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------

st.download_button(
    label="⬇ Download filtered catalog (CSV)",
    data=df.select(present_cols).to_pandas().to_csv(index=False).encode(),
    file_name="gaia_encounters_filtered.csv",
    mime="text/csv",
)
