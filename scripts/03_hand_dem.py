"""Step 3: Height Above Nearest Drainage (HAND) and slope from the 3DEP DEM.

The USGS 3DEP 1/3 arc-second DEM (averaged to the 20 m grid) is conditioned
with pysheds (fill pits, fill depressions, resolve flats), D8 flow directions
and flow accumulation are computed, and cells draining more than 1 km2 are
treated as the drainage network. HAND is the elevation of each cell above
the drainage cell it flows to. Cells with HAND above 15 m are later excluded
from flood mapping, a common rule for SAR flood maps (for example Twele and
others, 2016) that removes radar shadow and dark dry surfaces on high ground.
"""
import sys
from pathlib import Path

import numpy as np

if not hasattr(np, "in1d"):      # pysheds 0.5 still calls np.in1d (removed in recent NumPy)
    np.in1d = np.isin
from pysheds.grid import Grid  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

STREAM_KM2 = 1.0


def main():
    dem = C.read_tif(C.RAW / "ancillary" / "dem_3dep_13.tif").astype("float64")
    ok = np.isfinite(dem) & (dem > -100)
    print("DEM valid fraction", ok.mean(), "range", np.nanmin(dem[ok]), np.nanmax(dem[ok]),
          flush=True)
    # Ocean and bay cells without DEM values are set to -1 m so they act as outlets.
    filled = np.where(ok, dem, -1.0)
    tmp = C.WORK / "dem_for_pysheds.tif"
    C.write_tif(tmp, filled.astype("float32"), nodata=-9999.0)

    grid = Grid.from_raster(str(tmp))
    d = grid.read_raster(str(tmp))
    d = grid.fill_pits(d)
    print("pits filled", flush=True)
    d = grid.fill_depressions(d)
    print("depressions filled", flush=True)
    d = grid.resolve_flats(d)
    print("flats resolved", flush=True)
    fdir = grid.flowdir(d)
    acc = grid.accumulation(fdir)
    thr = STREAM_KM2 * 1e6 / (C.RES * C.RES)
    streams = acc > thr
    hand = grid.compute_hand(fdir, grid.read_raster(str(tmp)), streams)
    hand = np.asarray(hand, dtype="float32")
    hand[~ok] = np.nan
    hand[hand < 0] = 0

    gy, gx = np.gradient(np.where(ok, dem, np.nan), C.RES)
    slope = np.degrees(np.arctan(np.hypot(gx, gy))).astype("float32")

    C.write_tif(C.WORK / "hand_m.tif", hand)
    C.write_tif(C.WORK / "slope_deg.tif", slope)
    C.write_tif(C.WORK / "dem_m.tif", np.where(ok, dem, np.nan).astype("float32"))
    C.write_tif(C.WORK / "streams_1km2.tif", streams.astype("uint8"), dtype="uint8", nodata=255)
    print("HAND percentiles 5/50/95:", np.nanpercentile(hand, [5, 50, 95]), flush=True)


if __name__ == "__main__":
    main()
