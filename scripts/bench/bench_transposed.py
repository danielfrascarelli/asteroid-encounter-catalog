"""Bench: transposed storage layout (3, N, T) where time is contiguous.

Hypothesis: Delta + Blosc only exploits time-smoothness when the time axis is
contiguous in memory. With shape (T, N, 3) (current), adjacent floats in flat
buffer order are unrelated asteroid positions → no compression. Storing as
(3, N, T) makes per-asteroid-per-coord time series contiguous → Delta produces
small diffs → high compression.

Trade-off: slab access (read positions[t] expecting (N, 3)) requires slicing
the LAST axis of (3, N, T) and transposing → decompresses the whole T_chunk
window. With LRU cache + sequential scan, this is amortised.
"""

from __future__ import annotations

import json
import logging
import shutil
import time
from pathlib import Path

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
log = logging.getLogger("trans")


def _on_disk(p: Path) -> int:
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


def _slab_read_via_axis2(z_ro, t_max: int = 200) -> tuple[float, float]:
    """Sequential slab reads from a (3, N, T)-shaped zarr.

    Each `positions[t]` would map to z[:, :, t] returning shape (3, N); we
    transpose to (N, 3) to match the KD-tree consumer's expectation.
    """
    cksum = 0.0
    t0 = time.monotonic()
    for k in range(t_max):
        sl = np.asarray(z_ro[:, :, k]).T  # (N, 3)
        cksum += float(sl[0, 0])
    return time.monotonic() - t0, cksum


def _bench(data_TN3, out_dir, label, *, t_chunk, filters, clevel=5):
    """Write data with shape (3, N, T) using given t_chunk along last axis."""
    import zarr
    from numcodecs import Blosc

    # Transpose to (3, N, T)
    data_3NT = np.ascontiguousarray(np.transpose(data_TN3, (2, 1, 0)))
    T = data_3NT.shape[2]
    N = data_3NT.shape[1]
    chunks = (3, N, t_chunk)

    if out_dir.exists():
        shutil.rmtree(out_dir)
    compressor = Blosc(cname="zstd", clevel=clevel, shuffle=Blosc.BITSHUFFLE)
    t0 = time.monotonic()
    z = zarr.open(
        str(out_dir),
        mode="w",
        shape=data_3NT.shape,
        chunks=chunks,
        dtype=data_3NT.dtype.str,
        compressor=compressor,
        filters=filters,
    )
    z[:] = data_3NT
    write_t = time.monotonic() - t0

    z_ro = zarr.open(str(out_dir), mode="r")
    slab_t, _ = _slab_read_via_axis2(z_ro, t_max=min(200, T))

    # Sanity: a couple of reconstructed slabs match the source
    max_err = 0.0
    for k in np.linspace(0, T - 1, 10).astype(int):
        recon = np.asarray(z_ro[:, :, k]).T  # (N, 3)
        diff = np.abs(recon - data_TN3[k])
        max_err = max(max_err, float(diff.max()))

    on_disk = _on_disk(out_dir)
    raw = data_TN3.nbytes
    ratio = raw / max(on_disk, 1)
    log.info(
        "%-60s ratio=%6.2fx  write=%5.1fs  slab[200]=%5.2fs  max|Δ|=%.2e AU",
        label,
        ratio,
        write_t,
        slab_t,
        max_err,
    )
    return ratio, slab_t, max_err


def main():
    from numcodecs import BitRound, Delta

    cache_dir = Path("data/cache/bench_zarr")
    src = sorted(cache_dir.glob("trajectory_*.npy"))[0]
    meta = json.loads(sorted(cache_dir.glob("trajectory_*.json"))[0].read_text())
    T, N, _ = meta["shape"]
    data = np.array(np.memmap(src, dtype=np.float32, mode="r", shape=(T, N, 3)))
    # Materialise to avoid memmap surprises during transpose
    log.info("Loaded shape (T=%d, N=%d, 3)  raw=%.2f GB", T, N, data.nbytes / 1e9)

    out_root = cache_dir / "transposed"
    out_root.mkdir(exist_ok=True)

    # Pure Delta (no precision loss)
    for t_chunk in (64, 128, 256, 512, 1024):
        _bench(
            data,
            out_root / f"T_delta_tc{t_chunk}",
            f"(3,N,{t_chunk}) Delta(t) zstd5 bitshuffle",
            t_chunk=t_chunk,
            filters=[Delta(dtype="float32")],
        )

    # Delta + BitRound
    for keepbits in (16, 12, 10, 8):
        _bench(
            data,
            out_root / f"T_br{keepbits}_delta_tc256",
            f"(3,N,256) BitRound({keepbits})+Delta(t) zstd5 bitshuffle",
            t_chunk=256,
            filters=[BitRound(keepbits=keepbits), Delta(dtype="float32")],
        )

    # Filter order matters! Try Delta-first (so BitRound rounds delta values)
    for keepbits in (16, 12, 10, 8):
        _bench(
            data,
            out_root / f"T_delta_br{keepbits}_tc256",
            f"(3,N,256) Delta(t)+BitRound({keepbits}) zstd5 bitshuffle  [delta first]",
            t_chunk=256,
            filters=[Delta(dtype="float32"), BitRound(keepbits=keepbits)],
        )


if __name__ == "__main__":
    main()
