"""Step 1: download every public input used by the pipeline.

Everything is anonymous and keyless:
  * Sentinel-1 RTC (gamma0, VV and VH) from the Microsoft Planetary Computer
    STAC API, collection "sentinel-1-rtc" (assets signed with the
    planetary-computer package; no account or key).
  * JRC Global Surface Water v1.3 occurrence ("jrc-gsw"), ESA WorldCover 2020
    v100 ("esa-worldcover") and the USGS 3DEP 1/3 arc-second DEM
    ("3dep-seamless") from the same STAC API.
  * USGS Hurricane Harvey flood-inundation polygons, mapped-area boundaries and
    high-water marks (Watson and others, 2018, doi:10.5066/F7VH5N3N) from
    ScienceBase.
  * Census 2017 ACS 5-year tract geodatabase for Texas (TIGER/Line with
    selected demographic tables) and 2017 cartographic ZCTA boundaries.
  * USGS NWIS instantaneous gage height for nine stream and reservoir gauges.

Rasters are cut to the analysis grid in common.py and stored in the cache.
Re-running skips files that already exist.
"""
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import planetary_computer as pc
import pystac_client
import rasterio
import requests
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import from_bounds

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works with python -I
import common as C  # noqa: E402

GDAL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    GDAL_HTTP_MULTIPLEX="YES",
    GDAL_HTTP_VERSION="2",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="3",
    GDAL_CACHEMAX=2048,
    VSI_CACHE="TRUE",
)
LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def catalog():
    return pystac_client.Client.open(C.STAC_URL, modifier=pc.sign_inplace)


def warp_to_grid(href, resampling, overview_level=None, src_nodata=None,
                 dst_dtype="float32"):
    """Read a remote raster onto the analysis grid (any source CRS)."""
    transform, width, height = C.grid()
    kw = {} if overview_level is None else {"overview_level": overview_level}
    with rasterio.open(href, **kw) as src:
        nod = src.nodata if src_nodata is None else src_nodata
        with WarpedVRT(src, crs=C.CRS, transform=transform, width=width,
                       height=height, resampling=resampling, src_nodata=nod,
                       nodata=nod) as vrt:
            arr = vrt.read(1).astype(dst_dtype)
            if nod is not None and np.issubdtype(arr.dtype, np.floating):
                arr[arr == nod] = np.nan
    return arr


def read_s1(href):
    """Read one Sentinel-1 RTC band onto the 20 m grid.

    The RTC COGs are already in UTM 15N at 10 m. The first internal overview
    is 20 m and is aligned with the analysis grid, so it is read directly with
    a single window request; this avoids pulling the full 10 m file across
    the network. If the scene is in another CRS, fall back to a warp.
    """
    transform, width, height = C.grid()
    bounds = (transform.c, transform.f - height * C.RES,
              transform.c + width * C.RES, transform.f)
    with rasterio.open(href) as full:
        same = (full.crs.to_string() == C.CRS)
    if not same:
        return warp_to_grid(href, Resampling.average, src_nodata=-32768)
    with rasterio.open(href, overview_level=0) as src:
        win = from_bounds(*bounds, src.transform).round_offsets().round_lengths()
        arr = src.read(1, window=win, boundless=True, fill_value=-32768,
                       out_shape=(height, width)).astype("float32")
    arr[(arr == -32768) | (arr <= 0)] = np.nan
    return arr


def download_s1(cat):
    out = C.RAW / "s1"
    out.mkdir(exist_ok=True)
    rows = []
    for key, cfg in C.SCENES.items():
        for date in cfg["dates"]:
            items = list(cat.search(collections=["sentinel-1-rtc"],
                                    bbox=C.AOI_LONLAT,
                                    datetime=f"{date}T00:00:00Z/{date}T23:59:59Z").items())
            items = [i for i in items
                     if i.properties["sat:relative_orbit"] == cfg["orbit"]]
            items.sort(key=lambda i: i.id)
            for pol in ("vv", "vh"):
                path = out / f"s1_{date}_o{cfg['orbit']:03d}_{pol}.tif"
                for it in items:
                    rows.append(dict(set=key, date=date, orbit=cfg["orbit"],
                                     orbit_state=it.properties["sat:orbit_state"],
                                     platform=it.properties.get("platform"),
                                     polarization=pol.upper(), item_id=it.id,
                                     datetime_utc=str(it.datetime)))
                if path.exists():
                    continue
                mosaic = None
                for it in items:
                    t0 = time.time()
                    arr = read_s1(pc.sign(it.assets[pol].href))
                    log(f"  read {it.id} {pol} in {time.time() - t0:.0f} s, "
                        f"valid {np.isfinite(arr).mean():.3f}")
                    mosaic = arr if mosaic is None else np.where(
                        np.isfinite(mosaic), mosaic, arr)
                C.write_tif(path, mosaic)
                log(f"wrote {path.name}")
    pd.DataFrame(rows).drop_duplicates().to_csv(C.TAB / "sentinel1_scenes_used.csv",
                                                index=False)


def download_ancillary(cat):
    out = C.RAW / "ancillary"
    out.mkdir(exist_ok=True)
    jobs = [
        # name, collection, query filter, asset, resampling, dtype
        ("jrc_occurrence.tif", "jrc-gsw", None, "occurrence", Resampling.nearest, "uint8"),
        ("worldcover2020.tif", "esa-worldcover", None, "map", Resampling.mode, "uint8"),
        ("dem_3dep_13.tif", "3dep-seamless", None, "data", Resampling.average, "float32"),
    ]
    for name, col, query, asset, rs, dtype in jobs:
        path = out / name
        if path.exists():
            continue
        items = list(cat.search(collections=[col], bbox=C.AOI_LONLAT).items())
        if col == "3dep-seamless":
            items = [i for i in items if i.id.endswith("-13")]
        if col == "esa-worldcover":
            items = [i for i in items if "2020" in i.id]
        log(f"{col}: {len(items)} items: {[i.id for i in items]}")
        mosaic = None
        for it in items:
            nod = 255 if col == "jrc-gsw" else (0 if col == "esa-worldcover" else None)
            arr = warp_to_grid(pc.sign(it.assets[asset].href), rs, src_nodata=nod,
                               dst_dtype="float32")
            mosaic = arr if mosaic is None else np.where(np.isfinite(mosaic), mosaic, arr)
        if dtype == "uint8":
            fill = 255 if col == "jrc-gsw" else 0
            m = np.where(np.isfinite(mosaic), mosaic, fill).astype("uint8")
            C.write_tif(path, m, dtype="uint8", nodata=fill)
        else:
            C.write_tif(path, mosaic.astype("float32"))
        log(f"wrote {name}")


def fetch(url, path, timeout=600):
    if path.exists():
        return
    log(f"GET {url}")
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)


SCIENCEBASE = {
    # basin: (item id, file name) from doi:10.5066/F7VH5N3N
    "Brazos_upper": ("5aa02374e4b0b1c392e68400", "Brazos_upper.zip"),
    "Brazos_lower": ("5aa02374e4b0b1c392e68400", "Brazos_lower.zip"),
    "San_Jacinto": ("5aa023ebe4b0b1c392e6881b", "San_Jacinto.zip"),
}


def download_usgs_harvey():
    out = C.RAW / "usgs_harvey"
    out.mkdir(exist_ok=True)
    for basin, (item, fname) in SCIENCEBASE.items():
        zpath = out / fname
        if not zpath.exists():
            meta = requests.get(f"https://www.sciencebase.gov/catalog/item/{item}",
                                params={"format": "json"}, timeout=120).json()
            url = [f["url"] for f in meta["files"] if f["name"] == fname][0]
            fetch(url, zpath)
        with zipfile.ZipFile(zpath) as z:
            z.extractall(out)


def download_census():
    out = C.RAW / "census"
    out.mkdir(exist_ok=True)
    fetch("https://www2.census.gov/geo/tiger/TIGER_DP/2017ACS/ACS_2017_5YR_TRACT_48.gdb.zip",
          out / "ACS_2017_5YR_TRACT_48.gdb.zip", timeout=1800)
    fetch("https://www2.census.gov/geo/tiger/GENZ2017/shp/cb_2017_us_zcta510_500k.zip",
          out / "cb_2017_us_zcta510_500k.zip", timeout=1800)


GAUGES = ["08114000", "08116650", "08072500", "08073000", "08073500",
          "08073600", "08074000", "08068500", "08069000"]


def download_nwis():
    """15-minute gage height (00065) or reservoir pool elevation (62614)."""
    path = C.DATA / "nwis_gauges_harvey_2017.csv"
    if path.exists():
        return
    url = "https://waterservices.usgs.gov/nwis/iv/"
    params = dict(format="json", sites=",".join(GAUGES), startDT="2017-08-20",
                  endDT="2017-09-15", parameterCd="00065,62614,62615", siteStatus="all")
    js = requests.get(url, params=params, timeout=300).json()
    rows = []
    for ts in js["value"]["timeSeries"]:
        si = ts["sourceInfo"]
        site = si["siteCode"][0]["value"]
        name = si["siteName"]
        lat = si["geoLocation"]["geogLocation"]["latitude"]
        lon = si["geoLocation"]["geogLocation"]["longitude"]
        pcode = ts["variable"]["variableCode"][0]["value"]
        unit = ts["variable"]["unit"]["unitCode"]
        nodata = ts["variable"]["noDataValue"]
        for block in ts["values"]:
            for v in block["value"]:
                val = float(v["value"])
                if val == nodata:
                    continue
                rows.append(dict(site_no=site, site_name=name, lat=lat, lon=lon,
                                 parameter=pcode, unit=unit, datetime=v["dateTime"],
                                 value=val, qualifiers=";".join(v["qualifiers"])))
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)
    log(f"NWIS: {len(df)} values, {df.site_no.nunique()} sites")


def main():
    t0 = time.time()
    with rasterio.Env(**GDAL_ENV):
        cat = catalog()
        download_usgs_harvey()
        download_census()
        download_nwis()
        download_ancillary(cat)
        download_s1(cat)
    stamp = dict(access_date_utc=time.strftime("%Y-%m-%d", time.gmtime()),
                 runtime_min=round((time.time() - t0) / 60, 1), log=LOG)
    (C.TAB / "download_log.json").write_text(json.dumps(stamp, indent=1))


if __name__ == "__main__":
    main()
