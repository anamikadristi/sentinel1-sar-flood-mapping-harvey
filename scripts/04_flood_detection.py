"""Step 4: flood detection by Otsu thresholding and log-ratio change detection.

For each image pair (descending 30 Aug 2017 vs. 6-18 Aug median; ascending
29 Aug 2017 vs. 5 Aug):

1. Log-ratio change image: d = 10*log10(post/pre) = post_dB - pre_dB.
2. Split-based Otsu thresholds. The image is cut into 4 x 4 km tiles. Tiles
   whose post-event backscatter is clearly bimodal (Ashman D > 2, dark class
   between 10 and 90 percent of the tile) are pooled, and Otsu's method is run
   on the pooled pixels. This follows the tile-selection idea of Martinis and
   others (2009) and Chini and others (2017) because a single global
   histogram is dominated by dry land. The same is done for the change
   image d (tiles with a clear decrease mode).
3. Open-water flood = post_dB < T_post AND d < T_change, for VV and for VH.
   The primary map uses VV; VH is kept as a sensitivity check.
4. Masks: permanent water (JRC GSW occurrence >= 75 percent or ESA
   WorldCover class 80), HAND > 15 m, and no-data pixels.
5. Clean-up: connected flood patches smaller than 10 pixels (0.4 ha) are
   removed.
6. Urban double-bounce flag: built-up pixels (WorldCover class 50) whose VV
   backscatter rose by more than 3 dB are flagged as "possible flooded urban".
   They are not counted as flood because a backscatter increase also comes
   from other causes (moving vehicles, debris, construction, wet walls).

Class codes in the output GeoTIFF:
  0 not flooded, 1 flooded (open water), 2 permanent water,
  3 excluded (HAND > 15 m), 4 possible flooded urban (VV increase > 3 dB,
  not counted), 255 no data.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from skimage.filters import threshold_otsu

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

TILE = 200            # pixels (4 km at 20 m)
HAND_MAX = 15.0
MIN_PATCH = 10
URBAN_INC_DB = 3.0
PAIRS = {
    "desc143_20170830": ("desc143_pre", "desc143_post"),
    "asc034_20170829": ("asc034_pre", "asc034_post"),
}


def ashman_d(x, t):
    a, b = x[x < t], x[x >= t]
    if a.size < 50 or b.size < 50:
        return 0.0, a.size / max(x.size, 1)
    d = np.sqrt(2) * abs(a.mean() - b.mean()) / np.sqrt(a.var() + b.var())
    return float(d), a.size / x.size


def split_otsu(img, valid, max_tiles=30, seed=0):
    """Pool pixels from the most bimodal tiles and return an Otsu threshold."""
    h, w = img.shape
    cands = []
    for r in range(0, h - TILE + 1, TILE):
        for c in range(0, w - TILE + 1, TILE):
            v = valid[r:r + TILE, c:c + TILE]
            if v.mean() < 0.8:
                continue
            x = img[r:r + TILE, c:c + TILE][v]
            if x.std() < 1e-3:
                continue
            t = threshold_otsu(x)
            d, frac = ashman_d(x, t)
            if d > 2 and 0.1 <= frac <= 0.9:
                cands.append((d, r, c))
    cands.sort(reverse=True)
    chosen = cands[:max_tiles]
    if len(chosen) < 3:
        raise RuntimeError("too few bimodal tiles for split-based Otsu")
    pool = np.concatenate([img[r:r + TILE, c:c + TILE][valid[r:r + TILE, c:c + TILE]]
                           for _, r, c in chosen])
    return float(threshold_otsu(pool)), len(cands), [(int(r), int(c)) for _, r, c in chosen], pool


def remove_small(mask, n):
    lab, k = ndimage.label(mask, structure=np.ones((3, 3)))
    if k == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    keep = sizes >= n
    keep[0] = False
    return keep[lab]


def main():
    jrc = C.read_tif(C.RAW / "ancillary" / "jrc_occurrence.tif")
    wc = C.read_tif(C.RAW / "ancillary" / "worldcover2020.tif")
    hand = C.read_tif(C.WORK / "hand_m.tif")
    permanent = ((jrc >= 75) & (jrc <= 100)) | (wc == 80)
    high = ~(hand <= HAND_MAX)          # also true where HAND is NaN
    built = wc == 50
    summary = {}
    for name, (pre_k, post_k) in PAIRS.items():
        res = {}
        maps = {}
        for pol in ("vv", "vh"):
            pre = C.read_tif(C.WORK / f"{pre_k}_{pol}_db.tif")
            post = C.read_tif(C.WORK / f"{post_k}_{pol}_db.tif")
            valid = np.isfinite(pre) & np.isfinite(post)
            d = post - pre
            land = valid & ~permanent & ~high
            t_post, n_post, tiles_post, pool_post = split_otsu(post, valid)
            t_chg, n_chg, tiles_chg, pool_chg = split_otsu(np.where(land, d, np.nan), land)
            flood = land & (post < t_post) & (d < t_chg)
            flood = remove_small(flood, MIN_PATCH)
            maps[pol] = dict(flood=flood, valid=valid, d=d)
            if pol == "vv":
                # single-criterion variants kept for the method comparison in step 5
                otsu_only = remove_small(land & (post < t_post), MIN_PATCH)
                change_only = remove_small(land & (d < t_chg), MIN_PATCH)
                C.write_tif(C.WORK / f"flood_vv_otsu_only_{name}.tif", otsu_only.astype("uint8"),
                            dtype="uint8", nodata=255)
                C.write_tif(C.WORK / f"flood_vv_change_only_{name}.tif",
                            change_only.astype("uint8"), dtype="uint8", nodata=255)
                res["vv_otsu_only_km2"] = round(otsu_only.sum() * C.RES ** 2 / 1e6, 1)
                res["vv_change_only_km2"] = round(change_only.sum() * C.RES ** 2 / 1e6, 1)
            res[pol] = dict(T_post_db=round(t_post, 2), T_change_db=round(t_chg, 2),
                            bimodal_tiles_post=n_post, bimodal_tiles_change=n_chg,
                            tiles_used=len(tiles_post),
                            flood_km2=round(flood.sum() * C.RES ** 2 / 1e6, 1))
            np.save(C.WORK / f"hist_{name}_{pol}_post.npy",
                    np.random.default_rng(0).choice(pool_post, min(200000, pool_post.size), replace=False))
            np.save(C.WORK / f"hist_{name}_{pol}_change.npy",
                    np.random.default_rng(0).choice(pool_chg, min(200000, pool_chg.size), replace=False))
            C.write_tif(C.WORK / f"logratio_{name}_{pol}.tif", d.astype("float32"))
        vv, vh = maps["vv"], maps["vh"]
        valid = vv["valid"]
        cls = np.zeros(valid.shape, dtype="uint8")
        urban_flag = valid & built & ~permanent & ~high & (vv["d"] > URBAN_INC_DB)
        cls[urban_flag] = 4
        cls[valid & high] = 3
        cls[valid & permanent] = 2
        cls[vv["flood"]] = 1
        cls[~valid] = 255
        C.write_tif(C.WORK / f"flood_classes_{name}.tif", cls, dtype="uint8", nodata=255)
        C.write_tif(C.WORK / f"flood_vh_{name}.tif", vh["flood"].astype("uint8"), dtype="uint8",
                    nodata=255)
        both = vv["flood"] & vh["flood"]
        either = vv["flood"] | vh["flood"]
        px = C.RES ** 2 / 1e6
        res["classes_km2"] = {k: round(float((cls == v).sum() * px), 1) for k, v in
                              [("not_flooded", 0), ("flooded", 1), ("permanent_water", 2),
                               ("excluded_hand", 3), ("possible_flooded_urban", 4)]}
        res["vv_vh_agreement"] = dict(both_km2=round(both.sum() * px, 1),
                                      either_km2=round(either.sum() * px, 1),
                                      jaccard=round(both.sum() / max(either.sum(), 1), 3))
        res["valid_km2"] = round(valid.sum() * px, 1)
        summary[name] = res
        print(name, json.dumps(res, indent=1), flush=True)
    (C.TAB / "flood_detection_summary.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
