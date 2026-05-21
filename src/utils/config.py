"""Pipeline configuration: loads config.yaml (+ optional config.local.yaml override)."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)


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
class GaiaSSOSourceConfig:
    table: str
    archive_url: str
    columns: list[str]
    batch_size: int = 5_000
    n_workers: int | str = "auto"
    max_retries: int = 3


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


@dataclass
class CharacterizeConfig:
    compute_relative_velocity: bool
    compute_phase_angle: bool
    estimate_diameters: bool
    default_albedo: float


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
    _require(s["gaia_sso"], "table", "archive_url", "columns")
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
            gaia_sso=GaiaSSOSourceConfig(**s["gaia_sso"]),
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
