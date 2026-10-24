"""Shared settings for the Hurricane Harvey Sentinel-1 flood mapping pipeline.

All scripts import this module. Large inputs and intermediate rasters go to a
cache folder that is NOT part of the repository. Set the environment variable
HARVEY_CACHE to choose where it lives (default: ./cache next to the repo root).
"""
import os
from pathlib import Path

import numpy as np
from pyproj import Transformer
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path(os.environ.get("HARVEY_CACHE", ROOT / "cache")).resolve()
RAW = CACHE / "raw"          # downloaded files (treated as untrusted data)
WORK = CACHE / "work"        # intermediate rasters
OUT = ROOT / "outputs"
FIG = OUT / "figures"
TAB = OUT / "tables"
GIS = OUT / "gis"
DATA = ROOT / "data"
for p in (RAW, WORK, FIG, TAB, GIS, DATA):
    p.mkdir(parents=True, exist_ok=True)

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

# Area of interest (lon/lat): west Houston (Addicks and Barker reservoirs,
# Buffalo Bayou), north Houston (Cypress and Spring creeks, Lake Houston) and
# the lower Brazos River lowlands from Simonton to the coast near Freeport.
AOI_LONLAT = (-96.25, 29.00, -95.00, 30.35)

# Analysis grid: UTM zone 15N, 20 m pixels (2 x 2 multilook of 10 m RTC).
CRS = "EPSG:32615"
RES = 20.0


def grid():
    """Return (transform, width, height) of the analysis grid."""
    tr = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    lon0, lat0, lon1, lat1 = AOI_LONLAT
    xs, ys = tr.transform([lon0, lon0, lon1, lon1], [lat0, lat1, lat0, lat1])
    xmin = np.floor(min(xs) / RES) * RES
    xmax = np.ceil(max(xs) / RES) * RES
    ymin = np.floor(min(ys) / RES) * RES
    ymax = np.ceil(max(ys) / RES) * RES
    width = int(round((xmax - xmin) / RES))
    height = int(round((ymax - ymin) / RES))
    return from_origin(xmin, ymax, RES, RES), width, height


# Sentinel-1 RTC scenes (Planetary Computer collection "sentinel-1-rtc").
# Primary pair: descending relative orbit 143, about 12:22 UTC (07:22 CDT).
# Secondary pair: ascending relative orbit 34, about 00:26 UTC (19:26 CDT the
# previous calendar day).
SCENES = {
    "desc143_pre": {
        "dates": ["2017-08-06", "2017-08-12", "2017-08-18"],
        "orbit": 143,
        "note": "pre-event median of three dates",
    },
    "desc143_post": {"dates": ["2017-08-30"], "orbit": 143, "note": "post-event"},
    "asc034_pre": {"dates": ["2017-08-05"], "orbit": 34, "note": "pre-event"},
    "asc034_post": {"dates": ["2017-08-29"], "orbit": 34, "note": "post-event"},
}

GTIFF_PROFILE = dict(
    driver="GTiff", tiled=True, blockxsize=512, blockysize=512,
    compress="deflate", predictor=2, BIGTIFF="IF_SAFER",
)


def profile(dtype="float32", nodata=np.nan, count=1, **kw):
    transform, width, height = grid()
    p = dict(GTIFF_PROFILE, dtype=dtype, nodata=nodata, count=count,
             crs=CRS, transform=transform, width=width, height=height)
    if dtype in ("uint8", "int16", "uint16", "int32") or nodata is None:
        p["predictor"] = 2 if dtype != "uint8" else 1
    else:
        p["predictor"] = 3
    p.update(kw)
    return p


def write_tif(path, arr, dtype=None, nodata=None, **kw):
    import rasterio
    arr = np.asarray(arr)
    dtype = dtype or str(arr.dtype)
    if nodata is None and np.issubdtype(np.dtype(dtype), np.floating):
        nodata = np.nan
    p = profile(dtype=dtype, nodata=nodata, **kw)
    with rasterio.open(path, "w", **p) as dst:
        dst.write(arr.astype(dtype), 1)


def read_tif(path):
    import rasterio
    with rasterio.open(path) as src:
        return src.read(1)
