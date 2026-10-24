"""Step 5: compare the SAR flood maps with independent USGS Harvey data.

Reference: USGS flood-inundation polygons and mapped-area boundaries for the
upper and lower Brazos River and the San Jacinto River (Watson and others,
2018, doi:10.5066/F7VH5N3N). These were interpolated from surveyed
high-water marks (HWMs), so they show PEAK inundation, which on many reaches
came one to four days after the radar passes. They are an independent
estimate, not ground truth.

Pixel comparison is limited to each mapped-area boundary, to pixels with
valid SAR data, and outside permanent water and HAND > 15 m (the same
pixels the SAR map could call flooded). Metrics: overall accuracy (OA),
precision (user's accuracy for flood), recall (producer's accuracy, hit
rate), F1 and Cohen's kappa, also split by ESA WorldCover 2020 land cover:
built-up (class 50), tree cover (class 10) and open land (all other land
classes). For reference-flood pixels the script also reports how often VV
backscatter rose by more than 3 dB, the signature of double bounce between
water and walls or tree trunks.

Point check: share of USGS HWM sites with SAR flood within 60 m (3 x 3
pixel window), by height of water above ground.
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from rasterio.features import rasterize
from scipy import ndimage
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

BASINS = {
    "Brazos upper (Simonton to Rosharon)": ("Brazos_upper", "brazos_upper"),
    "Brazos lower (Rosharon to Freeport)": ("Brazos_lower", "brazos_lower"),
    "San Jacinto (Cypress, Spring, Lake Houston)": ("San_Jacinto", "sj"),
}
PAIRS = ["desc143_20170830", "asc034_20170829"]


def metrics(pred, ref):
    tp = int((pred & ref).sum()); fp = int((pred & ~ref).sum())
    fn = int((~pred & ref).sum()); tn = int((~pred & ~ref).sum())
    n = tp + fp + fn + tn
    oa = (tp + tn) / n if n else np.nan
    prec = tp / (tp + fp) if tp + fp else np.nan
    rec = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * prec * rec / (prec + rec) if prec and rec else np.nan
    pe = ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n ** 2 if n else np.nan
    kappa = (oa - pe) / (1 - pe) if n else np.nan
    px = C.RES ** 2 / 1e6
    return dict(n_pixels=n, ref_flood_km2=round((tp + fn) * px, 1),
                sar_flood_km2=round((tp + fp) * px, 1), TP=tp, FP=fp, FN=fn, TN=tn,
                OA=round(oa, 3), precision=round(prec, 3), recall=round(rec, 3),
                F1=round(f1, 3), kappa=round(kappa, 3))


def main():
    transform, width, height = C.grid()
    shape = (height, width)
    wc = C.read_tif(C.RAW / "ancillary" / "worldcover2020.tif")
    built = wc == 50
    tree = wc == 10
    openl = ~built & ~tree
    LANDS = lambda m: [("all land", m), ("built-up", m & built), ("tree cover", m & tree),
                       ("open land", m & openl)]
    aoi = gpd.GeoSeries([box(*C.AOI_LONLAT)], crs=4326).to_crs(C.CRS).iloc[0]
    rows, hwm_rows = [], []
    ref_union = np.zeros(shape, bool)
    bnd_union = np.zeros(shape, bool)
    for label, (folder, stem) in BASINS.items():
        d = C.RAW / "usgs_harvey" / folder
        bnd = gpd.read_file(d / f"{stem}_bndy.shp").to_crs(C.CRS)
        ref = gpd.read_file(d / f"{stem}_final.shp").to_crs(C.CRS)
        hwm = gpd.read_file(d / f"{stem}_hwm.shp").to_crs(C.CRS)
        bmask = rasterize([(g, 1) for g in bnd.geometry], out_shape=shape,
                          transform=transform, fill=0, dtype="uint8").astype(bool)
        rmask = rasterize([(g, 1) for g in ref.geometry], out_shape=shape,
                          transform=transform, fill=0, dtype="uint8").astype(bool)
        ref_union |= rmask & bmask
        bnd_union |= bmask
        frac_in_aoi = bnd.geometry.intersection(aoi).area.sum() / bnd.geometry.area.sum()
        for pair in PAIRS:
            cls = C.read_tif(C.WORK / f"flood_classes_{pair}.tif")
            vh = C.read_tif(C.WORK / f"flood_vh_{pair}.tif") == 1
            oo = C.read_tif(C.WORK / f"flood_vv_otsu_only_{pair}.tif") == 1
            co = C.read_tif(C.WORK / f"flood_vv_change_only_{pair}.tif") == 1
            evalm = bmask & np.isin(cls, [0, 1, 4])
            for variant, pred in [("VV", cls == 1), ("VH", vh), ("VV or VH", (cls == 1) | vh),
                                  ("VV Otsu only", oo), ("VV log-ratio only", co)]:
                for lc, lcm in LANDS(evalm):
                    m = metrics(pred[lcm], rmask[lcm])
                    rows.append(dict(basin=label, pair=pair, map=variant, land=lc,
                                     boundary_share_in_aoi=round(frac_in_aoi, 3), **m))
            # High-water marks
            flood = cls == 1
            near = ndimage.maximum_filter(flood.astype("uint8"), size=3).astype(bool)
            for _, h in hwm.iterrows():
                col, row = ~transform * (h.geometry.x, h.geometry.y)
                r, c = int(row), int(col)
                if not (0 <= r < height and 0 <= c < width):
                    continue
                hwm_rows.append(dict(basin=label, pair=pair, site_no=h["site_no"],
                                     waterbody=h.get("waterbody"),
                                     peak_date=h.get("peak_date"),
                                     height_above_ground_ft=h.get("height_abo"),
                                     sar_class=int(cls[r, c]), built_up=bool(built[r, c]),
                                     sar_flood_within_60m=bool(near[r, c])))
    df = pd.DataFrame(rows)
    df.to_csv(C.TAB / "validation_usgs_inundation.csv", index=False)
    hw = pd.DataFrame(hwm_rows)
    hw.to_csv(C.TAB / "validation_hwm_points.csv", index=False)

    # Pooled metrics over all three basins
    pooled = []
    for pair in PAIRS:
        cls = C.read_tif(C.WORK / f"flood_classes_{pair}.tif")
        evalm = bnd_union & np.isin(cls, [0, 1, 4])
        d = C.read_tif(C.WORK / f"logratio_{pair}_vv.tif")
        for lc, lcm in LANDS(evalm):
            r = ref_union & lcm
            pooled.append(dict(pair=pair, land=lc, **metrics((cls == 1)[lcm], ref_union[lcm]),
                               ref_share_vv_increase_gt3db=round(float((d[r] > 3).mean()), 3),
                               dry_share_vv_increase_gt3db=round(float((d[lcm & ~ref_union] > 3).mean()), 3),
                               ref_median_logratio_db=round(float(np.nanmedian(d[r])), 2)))
    pd.DataFrame(pooled).to_csv(C.TAB / "validation_pooled.csv", index=False)
    var = []
    for pair in PAIRS:
        cls = C.read_tif(C.WORK / f"flood_classes_{pair}.tif")
        evalm = bnd_union & np.isin(cls, [0, 1, 4])
        for variant, f in [("VV Otsu + log-ratio (primary)", None), ("VH Otsu + log-ratio", "flood_vh"),
                           ("VV Otsu only", "flood_vv_otsu_only"),
                           ("VV log-ratio only", "flood_vv_change_only")]:
            pred = (cls == 1) if f is None else (C.read_tif(C.WORK / f"{f}_{pair}.tif") == 1)
            var.append(dict(pair=pair, map=variant, **metrics(pred[evalm], ref_union[evalm])))
    var = pd.DataFrame(var)
    var.to_csv(C.TAB / "validation_method_variants.csv", index=False)
    print(var.drop(columns=["TP", "FP", "FN", "TN"]).to_string(), flush=True)

    hw["depth_bin_ft"] = pd.cut(hw["height_above_ground_ft"].astype(float),
                                [-0.01, 0.5, 2, 4, 100],
                                labels=["0-0.5", "0.5-2", "2-4", ">4"])
    summ = (hw.groupby(["pair", "depth_bin_ft"], observed=False)["sar_flood_within_60m"]
            .agg(["count", "mean"]).round(3).reset_index())
    summ.to_csv(C.TAB / "validation_hwm_by_depth.csv", index=False)
    C.write_tif(C.WORK / "usgs_ref_flood.tif", ref_union.astype("uint8"), dtype="uint8", nodata=255)
    C.write_tif(C.WORK / "usgs_ref_boundary.tif", bnd_union.astype("uint8"), dtype="uint8",
                nodata=255)
    print(pd.DataFrame(pooled).to_string(), flush=True)
    print(df[(df["map"] == "VV") & (df.land == "all land")].drop(columns=["TP", "FP", "FN", "TN"]).to_string(), flush=True)
    print(summ.to_string(), flush=True)
    print(hw.groupby(["pair", "built_up"])["sar_flood_within_60m"].agg(["count", "mean"]))


if __name__ == "__main__":
    main()
