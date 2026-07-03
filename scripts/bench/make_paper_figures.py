"""Generate publication-quality figures for the dataset paper.

This script produces the four main figures referenced by
``docs/dataset_paper_draft.md`` describing the Gaia DR3 asteroid close-encounter
catalogue and the Kepler-vs-N-body threshold-bias appendix.

Figures
-------
1. ``fig1_separation_hist.png`` -- histogram of minimum encounter separations
   (``dist_au``) with a log Y axis and the 0.05 AU detection threshold marked.
2. ``fig2_relvel_hist.png`` -- distribution of relative velocity at closest
   approach (``rel_vel_km_s``), with the median annotated.
3. ``fig3_aei_map.png`` -- orbital-element map (a vs e coloured by inclination,
   plus an inclination marginal) of the population that participates in
   encounters, with the main-belt boundaries marked.
4. ``fig4_kepler_nbody_threshold.png`` -- Delta dist (N-body minus Kepler)
   versus the Kepler distance for the near-threshold validation sample, showing
   the downward crossings that quantify the (~0.70 %) false-negative censoring.

Data handling
-------------
The two headline catalogues have ~72 M rows each and must **not** be loaded into
memory whole (there was a prior OOM on this dataset with only ~4 GB free). All
per-row statistics for figures 1 and 2 are computed as **streaming histogram
aggregations** with polars lazy + ``engine="streaming"``: we bin ``dist_au`` /
``rel_vel_km_s`` inside the query and only pull back the (small) bin-count table.
The orbital-elements file (fig 3) is tiny (~130 k bodies) and is loaded directly,
then subsampled to <=100 k points for scatter rendering. Figure 4 uses a small
validation parquet loaded whole.

Environment notes
-----------------
* Runs inside the project Docker image
  (``docker compose run --rm pipeline python -m scripts.bench.make_paper_figures``).
* Inside the container **only** ``data``, ``logs``, ``src`` and ``scripts`` are
  mounted -- ``docs`` is NOT. Therefore figures are written to
  ``data/output/figures/`` from the container; move them to ``docs/figures/`` on
  the host afterwards (``mv data/output/figures/*.png docs/figures/``).
* ``matplotlib`` may not be baked into the image; install it in the same
  ``docker compose run`` invocation before calling this module, e.g.::

      docker compose run --rm pipeline bash -c \
          "pip install --quiet matplotlib && \
           python -m scripts.bench.make_paper_figures"
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless backend for containers
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.figure import Figure

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger("make_paper_figures")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

OUT_DIR = Path("data/output/figures")

ENCOUNTERS_PARQUET = Path("data/output/encounters_catalog_rebound_005au.parquet")
CHARACTERIZED_PARQUET = Path("data/output/encounters_characterized_full.parquet")
ORBITS_PARQUET = Path("data/raw/gaia_orbits.parquet")
KEPLER_NBODY_PARQUET = Path("data/output/kepler_false_negatives/band_refined.parquet")

THRESHOLD_AU = 0.05
DPI = 300
RNG_SEED = 42
MAX_SCATTER = 100_000

# Sober journal aesthetics -------------------------------------------------- #
plt.rcParams.update(
    {
        "font.size": 11,
        "font.family": "serif",
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "axes.edgecolor": "0.2",
        "axes.linewidth": 0.8,
        "figure.dpi": 110,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
    }
)

NEUTRAL = "0.35"
ACCENT = "#8c1515"  # muted dark red for reference lines/annotations


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #


def _save(fig: Figure, stem: str) -> list[Path]:
    """Save a figure as PNG (and PDF) into :data:`OUT_DIR`.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to write.
    stem : str
        File name without extension.

    Returns
    -------
    list of pathlib.Path
        Paths actually written.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for ext in ("png", "pdf"):
        path = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(path)
        written.append(path)
        logger.info("wrote %s (%d bytes)", path, path.stat().st_size)
    plt.close(fig)
    return written


def _streaming_histogram(
    parquet: Path,
    column: str,
    lo: float,
    hi: float,
    n_bins: int,
    *,
    predicate: pl.Expr | None = None,
    clip: bool = True,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Compute a histogram of one column via a streaming polars aggregation.

    The full column is never materialised in memory: rows are assigned to a bin
    index inside the lazy query and only the (n_bins-row) group counts are
    collected.

    Parameters
    ----------
    parquet : pathlib.Path
        Source parquet file.
    column : str
        Numeric column to histogram.
    lo, hi : float
        Histogram range (inclusive of ``lo``; ``hi`` is the right edge).
    n_bins : int
        Number of equal-width bins.
    predicate : polars.Expr, optional
        Optional row filter applied before binning.
    clip : bool, default True
        If True, values outside ``[lo, hi]`` are clipped into the edge bins
        (so ``total`` equals all filtered rows). If False, out-of-range values
        are dropped, avoiding an artificial pile-up spike in the edge bins.

    Returns
    -------
    counts : numpy.ndarray
        Per-bin counts, shape ``(n_bins,)``.
    edges : numpy.ndarray
        Bin edges, shape ``(n_bins + 1,)``.
    total : int
        Total number of rows counted (after the predicate and range filter).
    """
    edges = np.linspace(lo, hi, n_bins + 1)
    width = (hi - lo) / n_bins

    lf = pl.scan_parquet(parquet)
    if predicate is not None:
        lf = lf.filter(predicate)
    if not clip:
        lf = lf.filter((pl.col(column) >= lo) & (pl.col(column) <= hi))

    # Map to an integer bin id in [0, n_bins - 1].
    bin_expr = (
        ((pl.col(column).clip(lo, hi) - lo) / width)
        .floor()
        .clip(0, n_bins - 1)
        .cast(pl.Int32)
        .alias("bin")
    )
    grouped = (
        lf.select(bin_expr).group_by("bin").agg(pl.len().alias("count")).collect(engine="streaming")
    )

    counts = np.zeros(n_bins, dtype=np.int64)
    for row in grouped.iter_rows():
        counts[int(row[0])] = row[1]
    return counts, edges, int(counts.sum())


# --------------------------------------------------------------------------- #
# Figure 1 -- minimum-separation histogram
# --------------------------------------------------------------------------- #


def figure1_separation_hist() -> list[Path]:
    """Histogram of minimum encounter separations (``dist_au``), log Y."""
    logger.info("Figure 1: separation histogram from %s", ENCOUNTERS_PARQUET)
    n_bins = 100
    counts, edges, total = _streaming_histogram(
        ENCOUNTERS_PARQUET, "dist_au", 0.0, THRESHOLD_AU, n_bins
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.bar(
        centers,
        counts,
        width=width,
        color=NEUTRAL,
        edgecolor="none",
        align="center",
    )
    ax.set_yscale("log")
    ax.set_xlabel("Minimum separation $d_{\\min}$ (AU)")
    ax.set_ylabel("Number of encounters per bin")
    ax.set_xlim(0.0, THRESHOLD_AU * 1.02)

    ax.axvline(THRESHOLD_AU, color=ACCENT, linestyle="--", linewidth=1.2)
    ax.annotate(
        "detection threshold\n0.05 AU",
        xy=(THRESHOLD_AU, counts.max()),
        xytext=(THRESHOLD_AU * 0.72, counts.max() * 0.35),
        color=ACCENT,
        fontsize=9,
        ha="center",
        arrowprops=dict(arrowstyle="->", color=ACCENT, linewidth=0.9),
    )
    ax.text(
        0.03,
        0.94,
        f"$N = {total/1e6:.1f}\\times10^{{6}}$ encounters",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
    )
    ax.grid(True, which="major", axis="y", linestyle=":", linewidth=0.5, color="0.8")
    fig.tight_layout()
    return _save(fig, "fig1_separation_hist")


# --------------------------------------------------------------------------- #
# Figure 2 -- relative-velocity histogram
# --------------------------------------------------------------------------- #


def figure2_relvel_hist() -> list[Path]:
    """Distribution of relative velocity at encounter (``rel_vel_km_s``)."""
    logger.info("Figure 2: relative-velocity histogram from %s", CHARACTERIZED_PARQUET)

    # Robust range + median in a single streaming pass.
    stats = (
        pl.scan_parquet(CHARACTERIZED_PARQUET)
        .select(
            pl.col("rel_vel_km_s").quantile(0.001).alias("q_lo"),
            pl.col("rel_vel_km_s").quantile(0.999).alias("q_hi"),
            pl.col("rel_vel_km_s").median().alias("median"),
        )
        .collect(engine="streaming")
    )
    q_lo = float(stats["q_lo"][0])
    q_hi = float(stats["q_hi"][0])
    median = float(stats["median"][0])
    lo = max(0.0, np.floor(q_lo))
    hi = np.ceil(q_hi)
    logger.info("rel_vel range [%.2f, %.2f] km/s, median %.3f", lo, hi, median)

    n_bins = 80
    counts, edges, _ = _streaming_histogram(
        CHARACTERIZED_PARQUET, "rel_vel_km_s", lo, hi, n_bins, clip=False
    )
    total = (
        pl.scan_parquet(CHARACTERIZED_PARQUET).select(pl.len()).collect(engine="streaming").item()
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]

    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.bar(
        centers,
        counts,
        width=width * 0.95,
        color=NEUTRAL,
        edgecolor="none",
        align="center",
    )
    ax.set_xlabel("Relative velocity at closest approach $v_{\\rm rel}$ (km s$^{-1}$)")
    ax.set_ylabel("Number of encounters per bin")
    ax.set_xlim(lo, hi)

    ax.axvline(median, color=ACCENT, linestyle="--", linewidth=1.2)
    ax.annotate(
        f"median = {median:.2f} km s$^{{-1}}$",
        xy=(median, counts.max() * 0.9),
        xytext=(median + 0.12 * (hi - lo), counts.max() * 0.9),
        color=ACCENT,
        fontsize=9,
        va="center",
        arrowprops=dict(arrowstyle="->", color=ACCENT, linewidth=0.9),
    )
    ax.text(
        0.97,
        0.94,
        f"$N = {total/1e6:.1f}\\times10^{{6}}$ encounters",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="right",
    )
    ax.grid(True, which="major", axis="y", linestyle=":", linewidth=0.5, color="0.8")
    fig.tight_layout()
    return _save(fig, "fig2_relvel_hist")


# --------------------------------------------------------------------------- #
# Figure 3 -- (a, e, i) map
# --------------------------------------------------------------------------- #


def figure3_aei_map() -> list[Path]:
    """Orbital-element map (a vs e coloured by i) of the encounter population."""
    logger.info("Figure 3: (a, e, i) map from %s", ORBITS_PARQUET)
    df = pl.read_parquet(ORBITS_PARQUET, columns=["a_au", "e", "i_deg"]).drop_nulls()
    # Keep the main-belt-ish region for a legible frame; drop extreme outliers.
    df = df.filter(
        (pl.col("a_au") > 1.5)
        & (pl.col("a_au") < 4.5)
        & (pl.col("e") >= 0.0)
        & (pl.col("e") < 0.5)
        & (pl.col("i_deg") < 40.0)
    )
    n_total = df.height
    logger.info("%d bodies in plotted a-e-i window", n_total)

    if n_total > MAX_SCATTER:
        df = df.sample(n=MAX_SCATTER, seed=RNG_SEED)
        logger.info("subsampled to %d points for scatter", df.height)

    a = df["a_au"].to_numpy()
    e = df["e"].to_numpy()
    i = df["i_deg"].to_numpy()

    fig, (ax0, ax1) = plt.subplots(
        1,
        2,
        figsize=(9.0, 4.4),
        gridspec_kw={"width_ratios": [3.0, 1.0]},
        layout="constrained",
    )

    # Main panel: a vs e coloured by inclination.
    sc = ax0.scatter(
        a,
        e,
        c=i,
        s=2,
        alpha=0.35,
        cmap="viridis",
        vmin=0.0,
        vmax=np.percentile(i, 99),
        edgecolors="none",
        rasterized=True,
    )
    ax0.set_xlabel("Semi-major axis $a$ (AU)")
    ax0.set_ylabel("Eccentricity $e$")
    ax0.set_xlim(1.9, 3.6)
    ax0.set_ylim(0.0, 0.45)

    # Main-belt Kirkwood-gap boundaries (approx. mean-motion resonances).
    belt_lines = {
        "3:1 (2.50)": 2.50,
        "5:2 (2.82)": 2.82,
        "2:1 (3.28)": 3.28,
    }
    for label, a_res in belt_lines.items():
        ax0.axvline(a_res, color="0.5", linestyle=":", linewidth=0.8)
        ax0.text(
            a_res,
            0.44,
            label,
            rotation=90,
            va="top",
            ha="right",
            fontsize=7,
            color="0.4",
        )
    ax0.annotate(
        "inner belt",
        xy=(2.25, 0.02),
        fontsize=8,
        color="0.3",
        ha="center",
    )
    ax0.annotate(
        "outer belt",
        xy=(3.05, 0.02),
        fontsize=8,
        color="0.3",
        ha="center",
    )

    cbar = fig.colorbar(sc, ax=ax0, pad=0.02)
    cbar.set_label("Inclination $i$ (deg)")
    cbar.solids.set_alpha(1.0)

    # Marginal panel: inclination distribution.
    ax1.hist(
        i,
        bins=60,
        orientation="horizontal",
        color=NEUTRAL,
        edgecolor="none",
    )
    ax1.set_ylim(0.0, np.percentile(i, 99.5))
    ax1.set_xlabel("count")
    ax1.set_ylabel("Inclination $i$ (deg)")
    ax1.grid(True, axis="x", linestyle=":", linewidth=0.5, color="0.85")

    ax0.text(
        0.02,
        0.97,
        f"{n_total:,} bodies" + (f" ({df.height:,} shown)" if df.height < n_total else ""),
        transform=ax0.transAxes,
        fontsize=8,
        va="top",
    )
    return _save(fig, "fig3_aei_map")


# --------------------------------------------------------------------------- #
# Figure 4 -- Kepler vs N-body near the threshold
# --------------------------------------------------------------------------- #


def figure4_kepler_nbody_threshold() -> tuple[list[Path], bool]:
    """Delta dist (N-body - Kepler) vs Kepler distance near the threshold.

    Uses the real near-threshold validation sample
    (``kepler_false_negatives/band_refined.parquet``): Kepler-selected pairs in
    the [0.05, 0.06) AU band re-refined with N-body, which lets us *measure* the
    downward-crossing (false-negative) rate that the distance-thresholded
    catalogue censors.

    Returns
    -------
    written : list of pathlib.Path
        Paths written.
    real_data : bool
        True if real data were used (always True here unless the parquet is
        missing, in which case an illustrative schematic is drawn).
    """
    if KEPLER_NBODY_PARQUET.exists():
        return _figure4_real(), True
    logger.warning("Kepler/N-body parquet missing; drawing schematic figure")
    return _figure4_schematic(), False


def _figure4_real() -> list[Path]:
    """Real-data version of figure 4."""
    logger.info("Figure 4 (real data) from %s", KEPLER_NBODY_PARQUET)
    df = pl.read_parquet(KEPLER_NBODY_PARQUET)
    dk = df["dist_au_kepler"].to_numpy()
    delta = df["delta_dist_au"].to_numpy() * 1e3  # AU -> milli-AU for readability
    n = dk.size

    # Downward crossings: Kepler said >= threshold, N-body says < threshold.
    dn = df["dist_au_nbody"].to_numpy()
    crossing_down = (dk >= THRESHOLD_AU) & (dn < THRESHOLD_AU)
    n_cross = int(crossing_down.sum())
    rate = n_cross / n

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.scatter(
        dk[~crossing_down],
        delta[~crossing_down],
        s=6,
        color=NEUTRAL,
        alpha=0.35,
        edgecolors="none",
        rasterized=True,
        label="no crossing",
    )
    ax.scatter(
        dk[crossing_down],
        delta[crossing_down],
        s=16,
        color=ACCENT,
        alpha=0.9,
        edgecolors="none",
        label=f"censored below 0.05 AU ($N={n_cross}$)",
    )
    ax.axhline(0.0, color="0.5", linewidth=0.7)
    ax.axvline(THRESHOLD_AU, color="0.2", linestyle="--", linewidth=1.0)
    ax.text(
        THRESHOLD_AU,
        ax.get_ylim()[1],
        " 0.05 AU threshold",
        color="0.2",
        fontsize=8,
        va="top",
        ha="left",
    )

    ax.set_xlabel("Kepler minimum distance $d_{\\rm Kep}$ (AU)")
    ax.set_ylabel(r"$\Delta d = d_{\rm N\text{-}body} - d_{\rm Kep}$ (mAU)")
    ax.set_title("")

    med = float(np.median(delta))
    std = float(np.std(delta))
    ax.text(
        0.03,
        0.05,
        f"median $\\Delta d = {med:+.3f}$ mAU\n"
        f"$\\sigma_{{\\Delta d}} = {std:.3f}$ mAU\n"
        f"downward crossings: {n_cross}/{n} = {rate*100:.2f}%",
        transform=ax.transAxes,
        fontsize=8.5,
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.7"),
    )
    ax.legend(loc="upper left", frameon=False)
    ax.grid(True, linestyle=":", linewidth=0.5, color="0.85")
    fig.tight_layout()
    return _save(fig, "fig4_kepler_nbody_threshold")


def _figure4_schematic() -> list[Path]:
    """Illustrative schematic from published summary numbers (fallback only)."""
    logger.info("Figure 4 (schematic) from published summary numbers")
    rng = np.random.default_rng(RNG_SEED)
    # Published: median Delta ~ 12 uAU, sigma ~ 3.7e-4 AU, censoring ~0.70%.
    median_au = 12e-6
    sigma_au = 3.7e-4
    n = 17_000
    dk = rng.uniform(0.05, 0.06, n)
    delta = rng.normal(median_au, sigma_au, n)
    dn = dk + delta
    crossing_down = dn < THRESHOLD_AU

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.scatter(
        dk[~crossing_down],
        delta[~crossing_down] * 1e3,
        s=6,
        color=NEUTRAL,
        alpha=0.3,
        edgecolors="none",
        rasterized=True,
    )
    ax.scatter(
        dk[crossing_down],
        delta[crossing_down] * 1e3,
        s=14,
        color=ACCENT,
        alpha=0.9,
        edgecolors="none",
    )
    ax.axhline(0.0, color="0.5", linewidth=0.7)
    ax.axvline(THRESHOLD_AU, color="0.2", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Kepler minimum distance $d_{\\rm Kep}$ (AU)")
    ax.set_ylabel(r"$\Delta d = d_{\rm N\text{-}body} - d_{\rm Kep}$ (mAU)")
    ax.text(
        0.03,
        0.95,
        "ILLUSTRATIVE SCHEMATIC\n(published summary numbers,\nnot per-pair data)",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
        color=ACCENT,
        fontweight="bold",
    )
    ax.grid(True, linestyle=":", linewidth=0.5, color="0.85")
    fig.tight_layout()
    return _save(fig, "fig4_kepler_nbody_threshold")


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #


def main() -> None:
    """Generate all four paper figures into :data:`OUT_DIR`."""
    logger.info("Output directory: %s", OUT_DIR.resolve())
    figure1_separation_hist()
    figure2_relvel_hist()
    figure3_aei_map()
    _, real = figure4_kepler_nbody_threshold()
    logger.info("Figure 4 used %s data", "REAL" if real else "SCHEMATIC")
    logger.info("Done.")


if __name__ == "__main__":
    main()
