#!/usr/bin/env python3
#  그림 규칙: 패널 문자 대문자·모두 같은 크기(PANEL_LABEL_PT), 컬러, N 화살표 크게, 범례에 색 이름 없음, A·B는 같은 색 척도
import os, unicodedata
import numpy as np, pandas as pd, geopandas as gpd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
NFC = lambda s: unicodedata.normalize("NFC", str(s))
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); AN = os.path.join(ROOT, "04_분석"); RES = os.path.join(AN, "FLU", "results_std3")
SHP = os.environ.get("ST_SHP", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "boundary", "final.shp"))
OUTD = os.path.join(AN, "FLU", "figs"); os.makedirs(OUTD, exist_ok=True)
PANEL_LABEL_PT = 60; plt.rcParams.update({"font.family": "Arial", "font.size": 16})
G = gpd.read_file(SHP); G["region"] = G.region.map(NFC)
FO = pd.read_csv(os.path.join(ROOT, "02_패널", "FLU_fh_outcomes.csv"), encoding="utf-8-sig"); FO["fbd_region"] = FO.fbd_region.map(NFC)
a = FO[FO.year <= 2015].groupby("fbd_region").SR.mean().rename("v").reset_index(); b = FO[FO.year >= 2020].groupby("fbd_region").SR.mean().rename("v").reset_index()
def ld(f, col="OR"):
    r = pd.read_csv(os.path.join(RES, f)); r["fbd_region"] = r.fbd_region.map(NFC); return r.rename(columns={col: "v"})
sp = ld("S3_spatial_OR.csv"); sf = ld("S0_forced_only_spatialOR.csv")
lo, hi = np.quantile(pd.concat([a.v, b.v]), [0.02, 0.98])
P = [("A", a, "seq", "Influenza vaccination, age-sex\nstandardised %,\n2011–2015 mean"), ("B", b, "seq", "Influenza vaccination, age-sex\nstandardised %,\n2020–2024 mean"),
     ("C", sf, "div", "Spatial odds ratio,\nstructural controls only"), ("D", sp, "div", "Spatial odds ratio,\nmain model")]
fig, axes = plt.subplots(2, 2, figsize=(16, 15))
for ax, (lab, r, kind, title) in zip(axes.flat, P):
    cols = ["fbd_region", "v"] + (["class"] if "class" in r else [])
    g = G.merge(r[cols], left_on="region", right_on="fbd_region", how="left", validate="one_to_one"); assert g.v.notna().all()
    if kind == "div":
        x = np.abs(np.log(g.v)).max(); norm = TwoSlopeNorm(vmin=np.exp(-x), vcenter=1.0, vmax=np.exp(x)); cmap = "RdBu_r"
    else:
        norm = plt.Normalize(lo, hi); cmap = "YlGnBu"
    g.plot(column="v", cmap=cmap, norm=norm, ax=ax, edgecolor="#9a9a9a", linewidth=0.15)
    if "class" in g:
        g[g["class"] == "hotspot"].boundary.plot(ax=ax, color="black", linewidth=0.9)
        g[g["class"] == "coldspot"].boundary.plot(ax=ax, color="black", linewidth=0.9, linestyle=(0, (2, 1.5)))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); cb = fig.colorbar(sm, ax=ax, shrink=0.45, pad=0.01, location="right"); cb.set_label(title, fontsize=15); cb.ax.tick_params(labelsize=13)
    ax.set_axis_off()
    ax.annotate("", xy=(0.86, 0.36), xytext=(0.86, 0.10), xycoords="axes fraction", arrowprops=dict(facecolor="black", width=27, headwidth=81, headlength=70))
    ax.text(0.86, 0.40, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=48, fontweight="bold")
plt.subplots_adjust(left=0.05, right=0.97, top=0.95, bottom=0.05, wspace=0.05, hspace=0.12); fig.canvas.draw()
for (lab, *_), ax in zip(P, axes.flat):
    bb = ax.get_position(); fig.text(bb.x0 - 0.03, bb.y1 + 0.01, lab, fontsize=PANEL_LABEL_PT, fontweight="bold", ha="left", va="bottom")
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([0], [0], color="black", lw=1.5, label="Hotspot (C, D): P(spatial OR > 1) > 0.95"), Line2D([0], [0], color="black", lw=1.5, ls=(0, (2, 1.5)), label="Coldspot (C, D): P(spatial OR > 1) < 0.05")],
           loc="lower center", ncol=2, fontsize=14, frameon=False, bbox_to_anchor=(0.5, 0.0))
for ext in ["png", "tiff"]: fig.savefig(os.path.join(OUTD, f"Fig_maps_FLU.{ext}"), dpi=300 if ext == "png" else 600, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None)
print("[99FL4] 저장", OUTD)
