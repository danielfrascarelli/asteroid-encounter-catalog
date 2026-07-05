"""Pipeline configuration: loads config.yaml (+ optional config.local.yaml override)."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

# Gaia DR3 defaults, used to synthesise a "dr3" release block when a legacy
# config carries only the flat ``gaia_sso.table`` field (no ``releases`` map).
# Keeps pre-FPR configs valid without edits. These mirror the constants the
# DR3 ingest/mass layer used before the release selector existed.
_DR3_EPOCH_REF_JD_TCB = 2455197.5  # epoch = JD_TCB - this (days since J2010.0 TCB)
_DR3_WINDOW_START = "2014-07-25T00:00:00"
_DR3_WINDOW_END = "2017-05-28T00:00:00"
_DR3_MP_MAX = 160_000
_DR3_DEFAULT_TABLE = "gaiadr3.sso_observation"


# ---------------------------------------------------------------------------
# Leaf dataclasses
# ---------------------------------------------------------------------------


@dataclass
class RunConfig:
    name: str
    description: str
    seed: int


@dataclass
class PathsConfig:
    data_root: str
    raw: str
    cache: str
    output: str
    logs: str


@dataclass
class MpcorbSourceConfig:
    url: str
    local_filename: str
    refresh_days: int


@dataclass
class GaiaReleaseConfig:
    """Per-release metadata that differs between Gaia DR3 and FPR.

    Everything that is coupled to the Gaia data release lives here so a single
    ``release`` flag swaps all of it atomically (avoiding states like
    "FPR table + DR3 window"). See ``docs/gaia_fpr_data_model.md``.
    """

    table: str
    epoch_ref_jd_tcb: float
    window_start: str  # ISO 8601 UTC
    window_end: str  # ISO 8601 UTC
    mp_max: int
    # Columns to drop from the shared base set for this release. FPR's
    # sso_observation has no ``g_mag`` (no photometry), so a SELECT including it
    # would fail. Empty for DR3.
    columns_drop: list[str] = field(default_factory=list)
    # Boolean column flagging astrometric rejections, if the release exposes one
    # (FPR: ``is_rejected``). When set, the mass-layer fetch filters it out so
    # rejected transits don't corrupt the fit. None for DR3 (no such column).
    reject_flag_column: str | None = None


@dataclass
class GaiaSSOSourceConfig:
    archive_url: str
    columns: list[str]
    # Active release selector ("dr3" | "fpr") and the per-release metadata map.
    release: str = "dr3"
    releases: dict[str, GaiaReleaseConfig] = field(default_factory=dict)
    # Legacy/deprecated flat table name. Only consulted when ``releases`` is
    # empty (pre-FPR configs); a synthetic dr3 release is built from it.
    table: str | None = None
    batch_size: int = 5_000
    n_workers: int | str = "auto"
    max_retries: int = 3

    def active(self) -> GaiaReleaseConfig:
        """Resolve the currently selected release block.

        Falls back to a synthetic DR3 release when a legacy config provides only
        the flat ``table`` field and no ``releases`` map.
        """
        if self.releases:
            if self.release not in self.releases:
                raise ValueError(
                    f"gaia_sso.release = '{self.release}' not found in releases "
                    f"{sorted(self.releases)}"
                )
            return self.releases[self.release]
        # Legacy path: synthesise dr3 from flat fields + DR3 defaults.
        return GaiaReleaseConfig(
            table=self.table or _DR3_DEFAULT_TABLE,
            epoch_ref_jd_tcb=_DR3_EPOCH_REF_JD_TCB,
            window_start=_DR3_WINDOW_START,
            window_end=_DR3_WINDOW_END,
            mp_max=_DR3_MP_MAX,
            columns_drop=[],
        )

    def active_columns(self) -> list[str]:
        """Base columns minus the active release's ``columns_drop`` (order kept)."""
        drop = set(self.active().columns_drop)
        return [c for c in self.columns if c not in drop]


@dataclass
class GaiaOrbitsSourceConfig:
    archive_url: str
    batch_size: int = 5_000
    n_workers: int | str = "auto"
    max_retries: int = 3


@dataclass
class JplHorizonsSourceConfig:
    api_url: str
    rate_limit_seconds: float


@dataclass
class Fienga2003SourceConfig:
    vizier_catalog: str
    output_filename: str


@dataclass
class Goffin2014SourceConfig:
    vizier_catalog: str
    output_filename: str


@dataclass
class Galad2002SourceConfig:
    source_url: str
    output_filename: str


@dataclass
class SourcesConfig:
    mpcorb: MpcorbSourceConfig
    gaia_sso: GaiaSSOSourceConfig
    gaia_orbits: GaiaOrbitsSourceConfig
    jpl_horizons: JplHorizonsSourceConfig
    fienga_2003: Fienga2003SourceConfig
    goffin_2014: Goffin2014SourceConfig
    galad_2002: Galad2002SourceConfig


@dataclass
class SemimajorAxisRange:
    min: float
    max: float


@dataclass
class SubsetConfig:
    only_numbered: bool
    max_asteroids: int | None
    exclude_neas: bool
    semimajor_axis_au: SemimajorAxisRange


@dataclass
class TimeWindowConfig:
    start: str
    end: str
    scale: str


@dataclass
class ReboundConfig:
    integrator: str
    include_planets: list[str]
    include_major_asteroids: bool


@dataclass
class PropagationConfig:
    method: str
    time_step_hours: float
    reference_frame: str
    cache_results: bool
    rebound: ReboundConfig
    # Optional coarser step (hours) for the bulk-scan trajectory cache. When
    # set and larger than ``time_step_hours`` the bulk cache is built at this
    # step (12× smaller at 12 h vs 1 h) and the KD-tree scan widens its query
    # radius accordingly. Refinement falls back to Kepler-on-demand when the
    # cache is too coarse for quadratic interpolation.
    coarse_step_hours: float | None = None
    # On-disk format for the bulk trajectory cache. ``"zarr"`` (default) stores
    # the (T, N, 3) float32 array as a chunked, Blosc-zstd-bitshuffled Zarr v2
    # directory — typical ratio ~5× on smooth orbital data. ``"memmap"`` keeps
    # the original raw ``.npy`` layout. Cache keys are stable across formats,
    # so a stored memmap cache survives a config change to zarr (it just won't
    # be read; a new zarr cache will be computed alongside).
    cache_format: str = "zarr"


@dataclass
class PrefilterConfig:
    enabled: bool
    semimajor_diff_max_au: float
    inclination_diff_max_deg: float


@dataclass
class KdTreeConfig:
    leaf_size: int


@dataclass
class RefinementConfig:
    enabled: bool
    fine_time_step_seconds: float
    window_hours: float


@dataclass
class DetectionConfig:
    threshold_au: float
    prefilter: PrefilterConfig
    kdtree: KdTreeConfig
    refinement: RefinementConfig
    # Conservative upper bound on inter-asteroid relative velocity used to
    # widen the KD-tree query radius when the bulk trajectory is sampled
    # coarsely. ~25 km/s safely covers MBA-only subsets (max heliocentric
    # speeds ≲ 30 km/s at perihelion; mutual encounter velocities much lower).
    max_relative_velocity_km_s: float = 25.0
    # Out-of-core detection: stream scan candidates to disk shards, dedup with
    # DuckDB, and refine in batches. Bounds parent RAM for the full numbered
    # population, which OOMs the in-memory path in a 24 GB Docker Desktop VM.
    out_of_core: bool = False


@dataclass
class CharacterizeConfig:
    compute_relative_velocity: bool
    compute_phase_angle: bool
    estimate_diameters: bool
    default_albedo: float
    # Tabla de diámetros/albedos medidos (SBDB). None = solo albedo por zona/default.
    physical_data: str | None = None


@dataclass
class ParallelConfig:
    enabled: bool
    n_workers: int | str  # "auto" o entero
    backend: str
    chunk_size_days: int


@dataclass
class OutputConfig:
    format: str
    filename: str
    compression: str
    include_metadata: bool


@dataclass
class LoggingConfig:
    level: str
    format: str
    file_logging: bool
    show_progress_bars: bool


@dataclass
class KnownPair:
    perturber: int
    test: int


@dataclass
class ValidationConfig:
    known_pairs: list[KnownPair]
    compare_with_jpl: bool
    jpl_validation_top_n: int


@dataclass
class PipelineConfig:
    run: RunConfig
    paths: PathsConfig
    sources: SourcesConfig
    subset: SubsetConfig
    time_window: TimeWindowConfig
    propagation: PropagationConfig
    detection: DetectionConfig
    characterize: CharacterizeConfig
    parallel: ParallelConfig
    output: OutputConfig
    logging: LoggingConfig
    validation: ValidationConfig


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into *base*, returning a new dict."""
    result = dict(base)
    for key, val in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(val, dict):
            result[key] = _deep_merge(result[key], val)
        else:
            result[key] = val
    return result


def _require(d: dict[str, Any], *keys: str) -> None:
    """Raise ValueError if any key is missing from *d*."""
    for k in keys:
        if k not in d:
            raise ValueError(f"Missing required config key: '{k}'")


def _build_gaia_sso(d: dict[str, Any]) -> GaiaSSOSourceConfig:
    """Construct GaiaSSOSourceConfig, converting the nested ``releases`` map."""
    releases = {
        name: GaiaReleaseConfig(
            table=rel["table"],
            epoch_ref_jd_tcb=float(rel["epoch_ref_jd_tcb"]),
            window_start=rel["window_start"],
            window_end=rel["window_end"],
            mp_max=int(rel["mp_max"]),
            columns_drop=list(rel.get("columns_drop", [])),
            reject_flag_column=rel.get("reject_flag_column"),
        )
        for name, rel in (d.get("releases") or {}).items()
    }
    return GaiaSSOSourceConfig(
        archive_url=d["archive_url"],
        columns=d["columns"],
        release=d.get("release", "dr3"),
        releases=releases,
        table=d.get("table"),
        batch_size=d.get("batch_size", 5_000),
        n_workers=d.get("n_workers", "auto"),
        max_retries=d.get("max_retries", 3),
    )


def _build(raw: dict[str, Any]) -> PipelineConfig:
    """Construct a PipelineConfig from a raw merged dict, failing fast on missing keys."""
    _require(
        raw,
        "run",
        "paths",
        "sources",
        "subset",
        "time_window",
        "propagation",
        "detection",
        "characterize",
        "parallel",
        "output",
        "logging",
        "validation",
    )

    r = raw["run"]
    _require(r, "name", "description", "seed")

    p = raw["paths"]
    _require(p, "data_root", "raw", "cache", "output", "logs")

    s = raw["sources"]
    _require(
        s,
        "mpcorb",
        "gaia_sso",
        "gaia_orbits",
        "jpl_horizons",
        "fienga_2003",
        "goffin_2014",
        "galad_2002",
    )
    _require(s["mpcorb"], "url", "local_filename", "refresh_days")
    _require(s["gaia_sso"], "archive_url", "columns")
    # New release-selector schema: validate each release block. Legacy configs
    # without ``releases`` are accepted (synthetic dr3 built at .active()).
    gaia_releases = s["gaia_sso"].get("releases")
    if gaia_releases:
        active = s["gaia_sso"].get("release", "dr3")
        if active not in gaia_releases:
            raise ValueError(
                f"gaia_sso.release = '{active}' not in releases {sorted(gaia_releases)}"
            )
        for rel_name, rel in gaia_releases.items():
            _require(
                rel,
                "table",
                "epoch_ref_jd_tcb",
                "window_start",
                "window_end",
                "mp_max",
            )
    elif "table" not in s["gaia_sso"]:
        raise ValueError("gaia_sso needs either a 'releases' map or a legacy 'table' field")
    _require(s["gaia_orbits"], "archive_url")
    _require(s["jpl_horizons"], "api_url", "rate_limit_seconds")
    _require(s["fienga_2003"], "vizier_catalog", "output_filename")
    _require(s["goffin_2014"], "vizier_catalog", "output_filename")
    _require(s["galad_2002"], "source_url", "output_filename")

    sub = raw["subset"]
    _require(sub, "only_numbered", "max_asteroids", "exclude_neas", "semimajor_axis_au")
    _require(sub["semimajor_axis_au"], "min", "max")

    tw = raw["time_window"]
    _require(tw, "start", "end", "scale")

    prop = raw["propagation"]
    _require(prop, "method", "time_step_hours", "reference_frame", "cache_results", "rebound")
    _require(prop["rebound"], "integrator", "include_planets", "include_major_asteroids")
    # coarse_step_hours and max_relative_velocity_km_s are optional with sensible
    # defaults — older configs without these keys still validate.

    det = raw["detection"]
    _require(det, "threshold_au", "prefilter", "kdtree", "refinement")
    _require(det["prefilter"], "enabled", "semimajor_diff_max_au", "inclination_diff_max_deg")
    _require(det["kdtree"], "leaf_size")
    _require(det["refinement"], "enabled", "fine_time_step_seconds", "window_hours")

    char = raw["characterize"]
    _require(
        char,
        "compute_relative_velocity",
        "compute_phase_angle",
        "estimate_diameters",
        "default_albedo",
    )

    par = raw["parallel"]
    _require(par, "enabled", "n_workers", "backend", "chunk_size_days")

    out = raw["output"]
    _require(out, "format", "filename", "compression", "include_metadata")

    log = raw["logging"]
    _require(log, "level", "format", "file_logging", "show_progress_bars")

    val = raw["validation"]
    _require(val, "known_pairs", "compare_with_jpl", "jpl_validation_top_n")

    return PipelineConfig(
        run=RunConfig(**r),
        paths=PathsConfig(**p),
        sources=SourcesConfig(
            mpcorb=MpcorbSourceConfig(**s["mpcorb"]),
            gaia_sso=_build_gaia_sso(s["gaia_sso"]),
            gaia_orbits=GaiaOrbitsSourceConfig(**s["gaia_orbits"]),
            jpl_horizons=JplHorizonsSourceConfig(**s["jpl_horizons"]),
            fienga_2003=Fienga2003SourceConfig(**s["fienga_2003"]),
            goffin_2014=Goffin2014SourceConfig(**s["goffin_2014"]),
            galad_2002=Galad2002SourceConfig(**s["galad_2002"]),
        ),
        subset=SubsetConfig(
            only_numbered=sub["only_numbered"],
            max_asteroids=sub["max_asteroids"],
            exclude_neas=sub["exclude_neas"],
            semimajor_axis_au=SemimajorAxisRange(**sub["semimajor_axis_au"]),
        ),
        time_window=TimeWindowConfig(**tw),
        propagation=PropagationConfig(
            method=prop["method"],
            time_step_hours=prop["time_step_hours"],
            reference_frame=prop["reference_frame"],
            cache_results=prop["cache_results"],
            rebound=ReboundConfig(**prop["rebound"]),
            coarse_step_hours=prop.get("coarse_step_hours"),
            cache_format=str(prop.get("cache_format", "zarr")),
        ),
        detection=DetectionConfig(
            threshold_au=det["threshold_au"],
            prefilter=PrefilterConfig(**det["prefilter"]),
            kdtree=KdTreeConfig(**det["kdtree"]),
            refinement=RefinementConfig(**det["refinement"]),
            max_relative_velocity_km_s=float(det.get("max_relative_velocity_km_s", 25.0)),
            out_of_core=bool(det.get("out_of_core", False)),
        ),
        characterize=CharacterizeConfig(**char),
        parallel=ParallelConfig(**par),
        output=OutputConfig(**out),
        logging=LoggingConfig(**log),
        validation=ValidationConfig(
            known_pairs=[KnownPair(**kp) for kp in val["known_pairs"]],
            compare_with_jpl=val["compare_with_jpl"],
            jpl_validation_top_n=val["jpl_validation_top_n"],
        ),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def load_config(
    base_path: str | Path = "config.yaml",
    local_path: str | Path | None = "config.local.yaml",
) -> PipelineConfig:
    """Load pipeline configuration from *base_path*, optionally merging *local_path*.

    Parameters
    ----------
    base_path:
        Path to the main config file (typically ``config.yaml``).
    local_path:
        Path to the local override file (typically ``config.local.yaml``).
        If the file does not exist it is silently ignored.

    Returns
    -------
    PipelineConfig
        Fully validated, typed configuration object.

    Raises
    ------
    FileNotFoundError
        If *base_path* does not exist.
    ValueError
        If required keys are missing from the merged configuration.
    """
    base_path = Path(base_path)
    if not base_path.exists():
        raise FileNotFoundError(f"Config file not found: {base_path}")

    with base_path.open() as fh:
        raw: dict[str, Any] = yaml.safe_load(fh)

    if local_path is not None:
        local_path = Path(local_path)
        if local_path.is_file():
            with local_path.open() as fh:
                local_raw: dict[str, Any] = yaml.safe_load(fh) or {}
            raw = _deep_merge(raw, local_raw)
            logger.info("Loaded local config override from %s", local_path)

    return _build(raw)
