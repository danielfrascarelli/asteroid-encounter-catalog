"""Validate top encounters against JPL Horizons barycentric positions."""

from __future__ import annotations

import numpy as np
import polars as pl
from astropy.time import Time
from astroquery.jplhorizons import Horizons

cat = pl.read_parquet("data/output/encounters_catalog.parquet")
top5 = cat.head(5)

header = "{:<3} {:<44} {:<14} {:<14} {:<12}".format(
    "#", "Par", "Nuestro (AU)", "JPL (AU)", "Error (AU)"
)
print(header)
print("-" * 90)

for i, row in enumerate(top5.iter_rows(named=True), 1):
    t = Time(row["jd_tdb"], format="jd", scale="tdb")
    epochs = {
        "start": (t - 0.25).utc.iso[:19],
        "stop": (t + 0.25).utc.iso[:19],
        "step": "1h",
    }
    try:
        v1 = Horizons(id=str(row["number_1"]), location="@0", epochs=epochs).vectors()
        v2 = Horizons(id=str(row["number_2"]), location="@0", epochs=epochs).vectors()
        jd_target = row["jd_tdb"]
        idx1 = int(np.argmin(np.abs(np.array(v1["datetime_jd"]) - jd_target)))
        idx2 = int(np.argmin(np.abs(np.array(v2["datetime_jd"]) - jd_target)))
        p1 = np.array([float(v1["x"][idx1]), float(v1["y"][idx1]), float(v1["z"][idx1])])
        p2 = np.array([float(v2["x"][idx2]), float(v2["y"][idx2]), float(v2["z"][idx2])])
        jpl_dist = float(np.linalg.norm(p1 - p2))
        err = abs(jpl_dist - row["dist_au"])
        label = "({}) — ({})".format(row["number_1"], row["number_2"])
        print(
            "{:<3} {:<44} {:<14.6f} {:<14.6f} {:<12.6f}".format(
                i, label, row["dist_au"], jpl_dist, err
            )
        )
    except Exception as e:
        label = "({}) — ({})".format(row["number_1"], row["number_2"])
        print(f"{i:<3} {label:<44} ERROR: {e}")
