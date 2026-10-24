"""Step 6: flooded area and exposed population by census tract and ZCTA.

Population, housing units and people below poverty come from the 2017 ACS
5-year estimates in the Census TIGER/Line geodatabase for Texas tracts
(tables B01001, B25001, B17001). Two exposure estimates are reported:

* area-weighted: tract population times the flooded share of the tract's
  land pixels (all land cover classes);
* built-up dasymetric: tract population is spread evenly over the tract's
  ESA WorldCover built-up pixels (class 50) and counted where those pixels
  are mapped as flooded. Tracts without built-up pixels fall back to the
  area-weighted rule.

Because C-band SAR misses most flooding between buildings (see README), the
built-up estimate is a lower bound for people living in flooded areas. The
ZIP-level table (2010 ZCTA boundaries, 2017 cartographic file) sums the same
pixel-level population surface, so it adds up to the tract totals.

Outputs: CSV tables, a GeoPackage with tract and ZCTA layers plus flood
polygons, and a shapefile copy of the tract layer for ArcGIS users.
"""
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio
from rasterio.features import rasterize, shapes
from scipy import ndimage
from shapely.geometry import box, shape

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

PAIR = "desc143_20170830"     # primary map
GDB = "ACS_2017_5YR_TRACT_48.gdb.zip"


def load_tracts():
    p = f"/vsizip/{C.RAW / 'census' / GDB}"
    g = pyogrio.read_dataframe(p, layer="ACS_2017_5YR_TRACT_48_TEXAS",
                               columns=["GEOID", "NAMELSAD", "COUNTYFP", "ALAND", "GEOID_Data"],
                               bbox=C.AOI_LONLAT)
    t1 = pyogrio.read_dataframe(p, layer="X01_AGE_AND_SEX", columns=["GEOID", "B01001e1"],
                                read_geometry=False)
    t2 = pyogrio.read_dataframe(p, layer="X25_HOUSING_CHARACTERISTICS",
                                columns=["GEOID", "B25001e1"], read_geometry=False)
    t3 = pyogrio.read_dataframe(p, layer="X17_POVERTY", columns=["GEOID", "B17001e1", "B17001e2"],
                                read_geometry=False)
    tab = t1.merge(t2, on="GEOID").merge(t3, on="GEOID").rename(columns={"GEOID": "GEOID_Data"})
    g = g.merge(tab, on="GEOID_Data", how="left")
    g = g.rename(columns={"B01001e1": "pop_2017", "B25001e1": "housing_units",
                          "B17001e1": "pov_universe", "B17001e2": "pop_below_poverty"})
    return g.to_crs(C.CRS)


def main():
    transform, width, height = C.grid()
    shp = (height, width)
    px_km2 = C.RES ** 2 / 1e6
    cls = C.read_tif(C.WORK / f"flood_classes_{PAIR}.tif")
    wc = C.read_tif(C.RAW / "ancillary" / "worldcover2020.tif")
    valid_land = np.isin(cls, [0, 1, 3, 4])
    flood = cls == 1
    built = (wc == 50) & valid_land
    aoi = gpd.GeoSeries([box(*C.AOI_LONLAT)], crs=4326).to_crs(C.CRS).iloc[0]

    tr = load_tracts()
    tr = tr[tr.geometry.centroid.within(aoi)].reset_index(drop=True)
    tr["tid"] = np.arange(1, len(tr) + 1)
    tid = rasterize(zip(tr.geometry, tr.tid), out_shape=shp, transform=transform, fill=0,
                    dtype="int32")
    n = len(tr) + 1
    land_px = np.bincount(tid[valid_land], minlength=n)
    flood_px = np.bincount(tid[flood], minlength=n)
    built_px = np.bincount(tid[built], minlength=n)
    bflood_px = np.bincount(tid[built & flood], minlength=n)
    urb_px = np.bincount(tid[cls == 4], minlength=n)
    nod_px = np.bincount(tid[cls == 255], minlength=n)
    all_px = np.bincount(tid.ravel(), minlength=n)
    i = tr.tid.values
    tr["land_km2_mapped"] = land_px[i] * px_km2
    tr["sar_coverage"] = 1 - nod_px[i] / np.maximum(all_px[i], 1)
    tr["flood_km2"] = flood_px[i] * px_km2
    tr["flood_pct_land"] = 100 * flood_px[i] / np.maximum(land_px[i], 1)
    tr["builtup_km2"] = built_px[i] * px_km2
    tr["builtup_flood_pct"] = 100 * bflood_px[i] / np.maximum(built_px[i], 1)
    tr["urban_flag_km2"] = urb_px[i] * px_km2
    f_area = flood_px[i] / np.maximum(land_px[i], 1)
    f_built = np.where(built_px[i] > 0, bflood_px[i] / np.maximum(built_px[i], 1), f_area)
    for col in ["pop_2017", "housing_units", "pop_below_poverty"]:
        tr[f"{col}_exposed_areaw"] = tr[col] * f_area
        tr[f"{col}_exposed_builtup"] = tr[col] * f_built

    # Pixel population surface (built-up dasymetric) for ZCTA sums
    w = np.where(built, 1.0, 0.0)
    no_built = built_px == 0
    w = np.where(no_built[tid] & valid_land, 1.0, w)
    w[tid == 0] = 0
    wsum = np.bincount(tid.ravel(), weights=w.ravel(), minlength=n)
    pop_t = np.zeros(n); pop_t[i] = tr.pop_2017.fillna(0).values
    pop_px = (pop_t / np.maximum(wsum, 1e-9))[tid] * w
    pop_px[tid == 0] = 0

    z = gpd.read_file(f"zip://{C.RAW / 'census' / 'cb_2017_us_zcta510_500k.zip'}",
                      bbox=C.AOI_LONLAT).to_crs(C.CRS)
    z = z[z.geometry.centroid.within(aoi)].reset_index(drop=True)
    z["zid"] = np.arange(1, len(z) + 1)
    zid = rasterize(zip(z.geometry, z.zid), out_shape=shp, transform=transform, fill=0,
                    dtype="int32")
    m = len(z) + 1
    j = z.zid.values
    zl = np.bincount(zid[valid_land], minlength=m)
    zf = np.bincount(zid[flood], minlength=m)
    zp = np.bincount(zid.ravel(), weights=pop_px.ravel(), minlength=m)
    zpf = np.bincount(zid[flood], weights=pop_px[flood], minlength=m)
    z["land_km2_mapped"] = zl[j] * px_km2
    z["flood_km2"] = zf[j] * px_km2
    z["flood_pct_land"] = 100 * zf[j] / np.maximum(zl[j], 1)
    z["pop_2017_est"] = zp[j]
    z["pop_exposed_builtup"] = zpf[j]
    z = z.rename(columns={"ZCTA5CE10": "zcta"})[["zcta", "land_km2_mapped", "flood_km2",
                                                  "flood_pct_land", "pop_2017_est",
                                                  "pop_exposed_builtup", "geometry"]]

    rcols = ["flood_km2", "flood_pct_land", "builtup_flood_pct", "land_km2_mapped",
             "builtup_km2", "urban_flag_km2", "sar_coverage"] + \
            [c for c in tr.columns if "exposed" in c]
    tr[rcols] = tr[rcols].round(3)
    z[["land_km2_mapped", "flood_km2", "flood_pct_land", "pop_2017_est",
       "pop_exposed_builtup"]] = z[["land_km2_mapped", "flood_km2", "flood_pct_land",
                                    "pop_2017_est", "pop_exposed_builtup"]].round(2)
    keep = ["GEOID", "NAMELSAD", "COUNTYFP", "pop_2017", "housing_units", "pop_below_poverty"] + rcols
    tr[keep].sort_values("flood_km2", ascending=False).to_csv(
        C.TAB / "flood_by_census_tract.csv", index=False)
    z.drop(columns="geometry").sort_values("flood_km2", ascending=False).to_csv(
        C.TAB / "flood_by_zcta.csv", index=False)

    counties = {"201": "Harris", "157": "Fort Bend", "039": "Brazoria", "473": "Waller",
                "339": "Montgomery", "167": "Galveston", "291": "Liberty", "015": "Austin",
                "481": "Wharton", "071": "Chambers", "407": "San Jacinto", "321": "Matagorda",
                "477": "Washington", "185": "Grimes"}
    tr["county"] = tr.COUNTYFP.map(counties).fillna(tr.COUNTYFP)
    cty = tr.groupby("county")[["pop_2017", "flood_km2", "land_km2_mapped",
                                "pop_2017_exposed_areaw", "pop_2017_exposed_builtup",
                                "housing_units_exposed_builtup",
                                "pop_below_poverty_exposed_builtup"]].sum()
    cty["flood_pct_land"] = 100 * cty.flood_km2 / cty.land_km2_mapped
    cty = cty.round(1).sort_values("flood_km2", ascending=False)
    cty.loc["All tracts in AOI"] = cty.sum(numeric_only=True)
    cty.loc["All tracts in AOI", "flood_pct_land"] = round(
        100 * cty.loc["All tracts in AOI", "flood_km2"] / cty.loc["All tracts in AOI", "land_km2_mapped"], 1)
    cty = cty.round(1)
    cty.to_csv(C.TAB / "flood_exposure_by_county.csv")
    print(cty.to_string(), flush=True)

    # Flood area by land cover
    names = {10: "Tree cover", 20: "Shrubland", 30: "Grassland", 40: "Cropland", 50: "Built-up",
             60: "Bare", 80: "Permanent water", 90: "Herbaceous wetland", 95: "Mangroves"}
    lc = pd.DataFrame({"class": list(names.values()),
                       "land_km2": [((wc == k) & valid_land).sum() * px_km2 for k in names],
                       "flood_km2": [((wc == k) & flood).sum() * px_km2 for k in names]})
    lc = lc[lc.land_km2 > 0]
    lc["flood_pct"] = 100 * lc.flood_km2 / lc.land_km2
    lc.round(2).to_csv(C.TAB / "flood_by_landcover.csv", index=False)
    print(lc.round(2).to_string(), flush=True)

    # Vector outputs
    gpkg = C.GIS / "harvey_sar_flood_2017.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    out_tr = tr[keep + ["geometry"]]
    out_tr.to_file(gpkg, layer="tract_flood_exposure", driver="GPKG")
    z.to_file(gpkg, layer="zcta_flood", driver="GPKG")
    # Flood polygons: patches of 2 ha or more, simplified to 20 m
    lab, k = ndimage.label(flood, structure=np.ones((3, 3)))
    sizes = np.bincount(lab.ravel())
    big = (sizes >= 50)[lab] & flood
    polys = [shape(g) for g, v in shapes(big.astype("uint8"), mask=big, transform=transform)
             if v == 1]
    fp = gpd.GeoDataFrame({"source": ["S1 VV 2017-08-30 vs Aug median"] * len(polys)},
                          geometry=polys, crs=C.CRS)
    fp["geometry"] = fp.geometry.simplify(20)
    fp["area_ha"] = (fp.area / 1e4).round(2)
    fp.to_file(gpkg, layer="flood_extent_20170830", driver="GPKG")
    shp_dir = C.GIS / "tract_flood_exposure_shp"
    shp_dir.mkdir(exist_ok=True)
    short = {"NAMELSAD": "NAME", "pop_2017": "POP17", "housing_units": "HU17",
             "pop_below_poverty": "POV17", "flood_km2": "FLD_KM2", "flood_pct_land": "FLD_PCT",
             "builtup_flood_pct": "BLT_FLDPCT", "land_km2_mapped": "LAND_KM2",
             "builtup_km2": "BLT_KM2", "urban_flag_km2": "URBFLG_KM2", "sar_coverage": "SAR_COV",
             "pop_2017_exposed_areaw": "POP_EXP_AW", "pop_2017_exposed_builtup": "POP_EXP_BU",
             "housing_units_exposed_areaw": "HU_EXP_AW",
             "housing_units_exposed_builtup": "HU_EXP_BU",
             "pop_below_poverty_exposed_areaw": "POV_EXP_AW",
             "pop_below_poverty_exposed_builtup": "POV_EXP_BU"}
    out_tr.rename(columns=short).to_file(shp_dir / "tract_flood_exposure.shp")
    print("tracts", len(tr), "zctas", len(z), "flood polygons", len(fp), flush=True)


if __name__ == "__main__":
    main()
