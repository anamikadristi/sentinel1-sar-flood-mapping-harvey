"""Step 7: figures for the README (PNG, 150 dpi) and the small GeoTIFF export."""
import json
import sys
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm, LinearSegmentedColormap
from matplotlib.patches import Patch
from pyproj import Transformer
from rasterio.enums import Resampling

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

PAIR = "desc143_20170830"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
NAVY, GRAY = "#0d366b", "#c9c8c3"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titleweight": "bold",
                     "axes.titlesize": 10, "figure.dpi": 150, "savefig.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False})
TR, W, H = C.grid()
EXTENT = (TR.c, TR.c + W * C.RES, TR.f - H * C.RES, TR.f)
LL2UTM = Transformer.from_crs("EPSG:4326", C.CRS, always_xy=True)
ZOOMS = {
    "Addicks and Barker reservoirs (west Houston)": (-95.72, 29.70, -95.55, 29.87),
    "Brazos River at Richmond and Rosenberg": (-95.90, 29.48, -95.70, 29.65),
}
PLACES = {"Addicks Res.": (-95.63, 29.83), "Barker Res.": (-95.66, 29.76),
          "Richmond": (-95.76, 29.58), "Lake Houston": (-95.24, 29.95),
          "Downtown Houston": (-95.37, 29.76), "Freeport": (-95.36, 29.05)}


def win(lonlat):
    x0, y0 = LL2UTM.transform(lonlat[0], lonlat[1])
    x1, y1 = LL2UTM.transform(lonlat[2], lonlat[3])
    c0, r1 = ~TR * (x0, y0)
    c1, r0 = ~TR * (x1, y1)
    r0, r1, c0, c1 = int(r0), int(r1), int(c0), int(c1)
    return (slice(max(r0, 0), r1), slice(max(c0, 0), c1)), (x0, x1, y0, y1)


def sub(a, s=4):
    return a[::s, ::s]


def km_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(""); ax.set_ylabel("")
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(GRAY)


def scalebar(ax, x0, y0, km=10):
    ax.plot([x0, x0 + km * 1000], [y0, y0], color=INK, lw=2)
    ax.text(x0 + km * 500, y0 + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.012, f"{km} km",
            ha="center", va="bottom", fontsize=8, color=INK)


def load_ref():
    out = {}
    for folder, stem in [("Brazos_upper", "brazos_upper"), ("Brazos_lower", "brazos_lower"),
                         ("San_Jacinto", "sj")]:
        d = C.RAW / "usgs_harvey" / folder
        out[folder] = dict(bnd=gpd.read_file(d / f"{stem}_bndy.shp").to_crs(C.CRS),
                           ref=gpd.read_file(d / f"{stem}_final.shp").to_crs(C.CRS),
                           hwm=gpd.read_file(d / f"{stem}_hwm.shp").to_crs(C.CRS))
    return out


def fig_study_area(ref, gauges):
    pre = C.read_tif(C.WORK / "desc143_pre_vv_db.tif")
    fig, ax = plt.subplots(figsize=(7.2, 8.2))
    ax.imshow(sub(pre), extent=EXTENT, cmap="gray", vmin=-22, vmax=0)
    for i, (k, v) in enumerate(ref.items()):
        v["bnd"].boundary.plot(ax=ax, color=AQUA, lw=1.1)
        v["hwm"].plot(ax=ax, color=ORANGE, markersize=9, marker="o", edgecolor="white",
                      linewidth=0.4)
    g = gpd.GeoDataFrame(gauges, geometry=gpd.points_from_xy(gauges.lon, gauges.lat),
                         crs=4326).to_crs(C.CRS)
    g.plot(ax=ax, color=BLUE, marker="^", markersize=40, edgecolor="white", linewidth=0.6)
    for name, (lo, la) in PLACES.items():
        x, y = LL2UTM.transform(lo, la)
        ax.text(x, y, name, color="white", fontsize=8, fontweight="bold",
                bbox=dict(facecolor=INK, alpha=0.55, lw=0, pad=1.5))
    for name, b in ZOOMS.items():
        _, (x0, x1, y0, y1) = win(b)
        ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="white",
                                   lw=1.2, ls="--"))
    ax.set_xlim(EXTENT[0], EXTENT[1]); ax.set_ylim(EXTENT[2], EXTENT[3])
    km_axes(ax)
    scalebar(ax, EXTENT[0] + 5000, EXTENT[2] + 5000, 20)
    ax.legend(handles=[Patch(fc="none", ec=AQUA, label="USGS mapped-area boundary"),
                       plt.Line2D([], [], color=ORANGE, marker="o", ls="", label="USGS high-water mark"),
                       plt.Line2D([], [], color=BLUE, marker="^", ls="", label="USGS stream or reservoir gauge"),
                       Patch(fc="none", ec="white", ls="--", label="Zoom windows (Fig. 2)")],
              loc="upper right", fontsize=7.5, framealpha=0.9)
    ax.set_title("Study area: Sentinel-1 VV backscatter, pre-event median (6 to 18 Aug 2017)",
                 loc="left")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig1_study_area.png"); plt.close(fig)


def fig_pre_post(pair_keys=("desc143_pre", "desc143_post")):
    pre = C.read_tif(C.WORK / f"{pair_keys[0]}_vv_db.tif")
    post = C.read_tif(C.WORK / f"{pair_keys[1]}_vv_db.tif")
    d = post - pre
    cmap_div = LinearSegmentedColormap.from_list("div", ["#104281", "#5598e7", "#f0efec",
                                                         "#e66767", "#a61b1b"])
    fig, axes = plt.subplots(2, 3, figsize=(10, 7))
    for row, (name, b) in enumerate(ZOOMS.items()):
        s, (x0, x1, y0, y1) = win(b)
        ext = (x0, x1, y0, y1)
        for col, (img, title) in enumerate([(pre, "Pre-event VV (dB)"), (post, "30 Aug 2017 VV (dB)"),
                                            (d, "Log-ratio, post minus pre (dB)")]):
            ax = axes[row, col]
            if col < 2:
                im = ax.imshow(img[s], extent=ext, cmap="gray", vmin=-22, vmax=0)
            else:
                im = ax.imshow(img[s], extent=ext, cmap=cmap_div,
                               norm=TwoSlopeNorm(vcenter=0, vmin=-10, vmax=10))
            km_axes(ax)
            if row == 0:
                ax.set_title(title, loc="left")
            if col == 0:
                ax.set_ylabel(name, fontsize=8.5)
            if row == 1:
                fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03)
        scalebar(axes[row, 0], x0 + 800, y0 + 800, 5)
    fig.suptitle("Flooded fields and reservoir pools turn dark (low backscatter); urban flooding often brightens",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig2_pre_post_logratio_zooms.png"); plt.close(fig)


CLS_COLORS = [GRAY, BLUE, NAVY, "#e7e6e1", ORANGE]
CLS_LABELS = ["Not flooded", "Flooded (open water)", "Permanent water", "Excluded, HAND > 15 m",
              "Possible flooded urban (VV +3 dB, not counted)"]


def class_rgba(cls):
    rgba = np.zeros(cls.shape + (4,), dtype="float32")
    for v, colr in [(1, BLUE), (2, NAVY), (4, ORANGE)]:
        rgba[cls == v] = matplotlib.colors.to_rgba(colr)
    rgba[cls == 3] = matplotlib.colors.to_rgba("#ffffff", 0.45)
    return rgba


def fig_flood_map(summary):
    cls = C.read_tif(C.WORK / f"flood_classes_{PAIR}.tif")
    post = C.read_tif(C.WORK / "desc143_post_vv_db.tif")
    fig, ax = plt.subplots(figsize=(7.2, 8.2))
    ax.imshow(sub(post, 3), extent=EXTENT, cmap="gray", vmin=-22, vmax=2, alpha=0.9)
    # use max-pooling so thin flooded strips stay visible at overview scale
    c3 = cls[: H - H % 3, : W - W % 3].reshape(H // 3, 3, W // 3, 3)
    pooled = np.where((c3 == 1).any(axis=(1, 3)), 1, c3[:, 1, :, 1])
    ax.imshow(class_rgba(pooled), extent=(EXTENT[0], EXTENT[0] + (W - W % 3) * C.RES,
                                          EXTENT[3] - (H - H % 3) * C.RES, EXTENT[3]))
    km_axes(ax)
    ax.set_xlim(EXTENT[0], EXTENT[1]); ax.set_ylim(EXTENT[2], EXTENT[3])
    scalebar(ax, EXTENT[0] + 5000, EXTENT[2] + 5000, 20)
    for name, (lo, la) in PLACES.items():
        x, y = LL2UTM.transform(lo, la)
        ax.text(x, y, name, color="white", fontsize=8, fontweight="bold",
                bbox=dict(facecolor=INK, alpha=0.55, lw=0, pad=1.5))
    s = summary[PAIR]["classes_km2"]
    ax.legend(handles=[Patch(fc=BLUE, label=f"Flooded, open water ({s['flooded']:,.0f} km²)"),
                       Patch(fc=NAVY, label="Permanent water (JRC, WorldCover)"),
                       Patch(fc=ORANGE, label=f"Possible flooded urban, not counted ({s['possible_flooded_urban']:,.0f} km²)"),
                       Patch(fc="#ffffff", ec=GRAY, label="Excluded, HAND > 15 m")],
              loc="upper right", fontsize=7.5, framealpha=0.92)
    ax.set_title("Sentinel-1 flood extent, 30 Aug 2017 07:22 CDT (VV, Otsu + log-ratio)", loc="left")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig3_flood_extent_20170830.png"); plt.close(fig)


def fig_hist(summary):
    fig, axes = plt.subplots(2, 2, figsize=(9, 5.6))
    for i, pol in enumerate(["vv", "vh"]):
        for j, kind in enumerate(["post", "change"]):
            x = np.load(C.WORK / f"hist_{PAIR}_{pol}_{kind}.npy")
            t = summary[PAIR][pol]["T_post_db" if kind == "post" else "T_change_db"]
            ax = axes[i, j]
            ax.hist(x, bins=120, color=BLUE, alpha=0.85, edgecolor="none")
            ax.axvline(t, color=ORANGE, lw=2)
            ax.text(t, ax.get_ylim()[1] * 0.92, f"  Otsu {t:.1f} dB", color=INK, fontsize=8)
            ax.set_title(f"{pol.upper()} {'post-event backscatter' if kind == 'post' else 'log-ratio (post minus pre)'}",
                         loc="left")
            ax.set_xlabel("dB"); ax.set_ylabel("pixels in selected tiles")
            ax.grid(axis="y", color="#e7e6e1", lw=0.6)
    fig.suptitle("Pooled histograms from bimodal 4 km tiles, descending pass 30 Aug 2017",
                 x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig4_otsu_histograms.png"); plt.close(fig)


def fig_validation(ref):
    cls = C.read_tif(C.WORK / f"flood_classes_{PAIR}.tif")
    rf = C.read_tif(C.WORK / "usgs_ref_flood.tif") == 1
    bd = C.read_tif(C.WORK / "usgs_ref_boundary.tif") == 1
    post = C.read_tif(C.WORK / "desc143_post_vv_db.tif")
    ev = bd & np.isin(cls, [0, 1, 4])
    sar = cls == 1
    rgba = np.zeros(cls.shape + (4,), "float32")
    rgba[ev & sar & rf] = matplotlib.colors.to_rgba(BLUE)
    rgba[ev & sar & ~rf] = matplotlib.colors.to_rgba(ORANGE)
    rgba[ev & ~sar & rf] = matplotlib.colors.to_rgba(AQUA)
    b = ref["Brazos_upper"]["bnd"].total_bounds
    pad = 1500
    c0, r0 = ~TR * (b[0] - pad, b[3] + pad)
    c1, r1 = ~TR * (b[2] + pad, b[1] - pad)
    s = (slice(max(int(r0), 0), int(r1)), slice(max(int(c0), 0), int(c1)))
    ext = (b[0] - pad, b[2] + pad, b[1] - pad, b[3] + pad)
    v = pd.read_csv(C.TAB / "validation_usgs_inundation.csv")
    fig = plt.figure(figsize=(11, 5.4))
    ax = fig.add_axes([0.01, 0.06, 0.50, 0.86])
    ax.imshow(post[s], extent=ext, cmap="gray", vmin=-22, vmax=2)
    ax.imshow(rgba[s], extent=ext)
    ref["Brazos_upper"]["bnd"].boundary.plot(ax=ax, color="white", lw=0.8)
    km_axes(ax)
    ax.set_xlim(ext[0], ext[1]); ax.set_ylim(ext[2], ext[3])
    scalebar(ax, ext[0] + 2000, ext[2] + 2000, 10)
    ax.legend(handles=[Patch(fc=BLUE, label="Both flooded"),
                       Patch(fc=ORANGE, label="SAR only"),
                       Patch(fc=AQUA, label="USGS peak map only")],
              loc="lower right", fontsize=8, framealpha=0.92)
    ax.set_title("Brazos River (upper reach): SAR 30 Aug vs. USGS peak inundation", loc="left")
    ax2 = fig.add_axes([0.72, 0.12, 0.26, 0.76])
    sel = v[(v.pair == PAIR) & (v["map"] == "VV") & (v.land == "all land")]
    pool = pd.read_csv(C.TAB / "validation_pooled.csv")
    pool = pool[(pool.pair == PAIR) & (pool.land != "all land")]
    labels, f1s, recs, precs = [], [], [], []
    for _, r in sel.iterrows():
        labels.append(r.basin.split(" (")[0] + ", all land")
        f1s.append(r.F1); recs.append(r.recall); precs.append(r.precision)
    for _, r in pool.iterrows():
        labels.append(f"All basins, {r.land}")
        f1s.append(r.F1); recs.append(r.recall); precs.append(r.precision)
    y = np.arange(len(labels))
    hgt = 0.26
    ax2.barh(y - hgt, precs, height=hgt - 0.03, color=BLUE, label="Precision")
    ax2.barh(y, recs, height=hgt - 0.03, color=ORANGE, label="Recall")
    ax2.barh(y + hgt, f1s, height=hgt - 0.03, color=AQUA, label="F1")
    ax2.set_yticks(y, labels, fontsize=7.5)
    ax2.invert_yaxis()
    ax2.set_xlim(0, 1)
    ax2.grid(axis="x", color="#e7e6e1", lw=0.6)
    ax2.legend(fontsize=7.5, loc="lower right", frameon=False)
    ax2.set_title("Agreement with USGS (VV, 30 Aug)", loc="left")
    fig.savefig(C.FIG / "fig5_validation_usgs.png"); plt.close(fig)


def fig_gauges(gauges):
    g = gauges.copy()
    g["t"] = pd.to_datetime(g["datetime"], utc=True)
    sites = g.drop_duplicates("site_no")[["site_no", "site_name"]].values
    n = len(sites)
    ncol = 3
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(10, 2.3 * nrow), sharex=True)
    sar_times = [pd.Timestamp("2017-08-29T00:26Z"), pd.Timestamp("2017-08-30T12:22Z")]
    for ax, (site, name) in zip(axes.ravel(), sites):
        d = g[g.site_no == site]
        par = "00065" if (d.parameter == "00065").any() else d.parameter.iloc[0]
        d = d[d.parameter == par].sort_values("t")
        ax.plot(d.t, d.value, color=BLUE, lw=1.6)
        pk = d.loc[d.value.idxmax()]
        ax.plot(pk.t, pk.value, "o", color=BLUE, ms=5, mec="white")
        for st in sar_times:
            ax.axvline(st, color=ORANGE, lw=1.2, ls="--")
        lab = "stage (ft)" if par == "00065" else "pool elev. (ft)"
        ax.set_title(f"{name[:34]}", loc="left", fontsize=8.5)
        ax.set_ylabel(lab, fontsize=7.5)
        ax.grid(color="#e7e6e1", lw=0.6)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.xaxis.set_major_locator(matplotlib.dates.DayLocator(bymonthday=[22, 29, 5, 12]))
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%d %b"))
    fig.suptitle("USGS gauges: dashed lines mark the two Sentinel-1 passes (29 Aug 00:26 and 30 Aug 12:22 UTC)",
                 x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig6_gauge_hydrographs.png"); plt.close(fig)


def fig_tracts():
    tr = gpd.read_file(C.GIS / "harvey_sar_flood_2017.gpkg", layer="tract_flood_exposure")
    fig, ax = plt.subplots(figsize=(7.2, 7.6))
    cm = ListedColormap(["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
    bins = [0, 1, 5, 10, 20, 40, 100]
    norm = BoundaryNorm(bins, cm.N)
    tr.plot(column="flood_pct_land", ax=ax, cmap=cm, norm=norm, edgecolor="white", linewidth=0.25)
    km_axes(ax)
    ax.set_xlim(EXTENT[0], EXTENT[1]); ax.set_ylim(EXTENT[2], EXTENT[3])
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.01)
    cb.set_label("Share of tract land mapped as flooded (%)")
    for name, (lo, la) in PLACES.items():
        x, y = LL2UTM.transform(lo, la)
        ax.text(x, y, name, color=INK, fontsize=7.5,
                bbox=dict(facecolor="white", alpha=0.7, lw=0, pad=1.2))
    ax.set_title("SAR-mapped flood share by census tract, 30 Aug 2017", loc="left")
    fig.tight_layout()
    fig.savefig(C.FIG / "fig7_tract_flood_share.png"); plt.close(fig)


def export_small_geotiff():
    """Copy the 20 m flood class raster to outputs/gis with a color table."""
    src_p = C.WORK / f"flood_classes_{PAIR}.tif"
    for name, factor in [("flood_classes_20170830_20m.tif", 1)]:
        with rasterio.open(src_p) as src:
            data = src.read(1, out_shape=(src.height // factor, src.width // factor),
                            resampling=Resampling.mode if factor > 1 else Resampling.nearest)
            tr = src.transform * src.transform.scale(src.width / data.shape[1],
                                                     src.height / data.shape[0])
            prof = src.profile
            prof.update(width=data.shape[1], height=data.shape[0], transform=tr,
                        compress="deflate", predictor=1, zlevel=9)
            with rasterio.open(C.GIS / name, "w", **prof) as dst:
                dst.write(data, 1)
                dst.update_tags(1, classes="0 not flooded; 1 flooded open water; 2 permanent water; "
                                           "3 excluded HAND>15m; 4 possible flooded urban (not counted); 255 no data",
                                source="Sentinel-1 RTC (Planetary Computer), S1B 2017-08-30 vs 2017-08-06/12/18 median")
                cmap = {0: (230, 230, 225, 255), 1: (42, 120, 214, 255), 2: (13, 54, 107, 255),
                        3: (255, 255, 255, 255), 4: (235, 104, 52, 255), 255: (0, 0, 0, 0)}
                dst.write_colormap(1, cmap)


def main():
    summary = json.loads((C.TAB / "flood_detection_summary.json").read_text())
    gauges = pd.read_csv(C.DATA / "nwis_gauges_harvey_2017.csv", dtype={"site_no": str, "parameter": str})
    ref = load_ref()
    fig_study_area(ref, gauges.drop_duplicates("site_no"))
    fig_pre_post()
    fig_flood_map(summary)
    fig_hist(summary)
    fig_validation(ref)
    fig_gauges(gauges)
    fig_tracts()
    export_small_geotiff()


if __name__ == "__main__":
    main()
