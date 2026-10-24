"""Step 2: build pre- and post-event backscatter images and filter speckle.

* The pre-event image for the descending track is the per-pixel median of
  three dates (6, 12 and 18 August 2017), which lowers speckle and removes
  short-lived wet-soil effects. The ascending track has one pre-event date.
* A Lee filter (7 x 7 window) is applied to linear gamma0. The equivalent
  number of looks (ENL) used by the filter is estimated from the image itself
  as the median local ENL of the most homogeneous 10 percent of windows.
* Output: filtered gamma0 in decibels for VV and VH.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

WIN = 7


def local_stats(img, size):
    valid = np.isfinite(img)
    x = np.where(valid, img, 0.0).astype("float64")
    w = ndimage.uniform_filter(valid.astype("float64"), size)
    m = ndimage.uniform_filter(x, size)
    m2 = ndimage.uniform_filter(x * x, size)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = m / w
        var = np.maximum(m2 / w - mean ** 2, 0)
    return mean, var, w


def estimate_enl(img, size=WIN):
    mean, var, w = local_stats(img, size)
    ok = (w > 0.99) & (mean > 0) & (var > 0)
    cv = np.sqrt(var[ok]) / mean[ok]
    enl = (mean[ok] ** 2) / var[ok]
    cut = np.quantile(cv, 0.10)
    return float(np.median(enl[cv <= cut]))


def lee_filter(img, enl, size=WIN):
    """Classic Lee (1980) filter for multiplicative speckle."""
    mean, var, _ = local_stats(img, size)
    cu2 = 1.0 / enl                       # squared speckle coefficient of variation
    with np.errstate(invalid="ignore", divide="ignore"):
        ci2 = var / mean ** 2
        k = np.clip((ci2 - cu2) / (ci2 * (1 + cu2)), 0, 1)
    out = mean + k * (np.where(np.isfinite(img), img, mean) - mean)
    out[~np.isfinite(img)] = np.nan
    return out.astype("float32")


def to_db(x):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (10 * np.log10(np.clip(x, 1e-5, None))).astype("float32")


def main():
    info = {}
    for key, cfg in C.SCENES.items():
        for pol in ("vv", "vh"):
            stack = [C.read_tif(C.RAW / "s1" / f"s1_{d}_o{cfg['orbit']:03d}_{pol}.tif")
                     for d in cfg["dates"]]
            raw = stack[0] if len(stack) == 1 else np.nanmedian(np.stack(stack), axis=0)
            enl = estimate_enl(raw)
            filt = lee_filter(raw, enl)
            C.write_tif(C.WORK / f"{key}_{pol}_db.tif", to_db(filt))
            C.write_tif(C.WORK / f"{key}_{pol}_db_unfiltered.tif", to_db(raw))
            info[f"{key}_{pol}"] = dict(dates=cfg["dates"], enl_estimate=round(enl, 2),
                                        valid_fraction=round(float(np.isfinite(raw).mean()), 4))
            print(key, pol, info[f"{key}_{pol}"], flush=True)
    (C.TAB / "speckle_filter_enl.json").write_text(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
