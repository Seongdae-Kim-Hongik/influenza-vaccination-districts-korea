#  그림 규칙: 컬러, 패널 문자 대문자·같은 크기, 600 dpi TIFF, 그림 안 설명문 금지
import os, shutil
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
import content_flu as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(OUT, exist_ok=True)
AN = os.path.normpath(os.path.join(HERE, "..", "..")); N = C.N
im = Image.open(os.path.join(AN, "figs", "Fig_maps_FLU.tiff")); im.save(os.path.join(OUT, "Figure1.tiff"), compression="tiff_lzw", dpi=(600, 600))
shutil.copy(os.path.join(AN, "figs", "Fig_maps_FLU.png"), os.path.join(OUT, "Figure1.png"))
plt.rcParams.update({"font.family": "Arial", "font.size": 11})
def save(fig, k):
    for ext, dpi in (("png", 300), ("tiff", 600)): fig.savefig(os.path.join(OUT, f"Figure{k}.{ext}"), dpi=dpi, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None, bbox_inches="tight")
VARS = [v for g, vs in C.GROUPS[:2] for v in vs]
NA = dict(OR=np.nan, lo=np.nan, hi=np.nan, cred=False)
SER = [("Main", lambda v: N["main"][v], "#1f4e79", "o"), ("Knorr-Held IV", lambda v: N["kh"][f"KH_IV|{v}"], "#2e8b57", "s"),
       ("Spatial+", lambda v: N["splus"][v], "#6a3d9a", "D"), ("Without survey covariates", lambda v: N["same"].get(f"non_survey_covariates_only|{v}", NA), "#d95f02", "P"),
       ("Between", lambda v: N["bw"][f"{v}|between"], "#b8860b", "^"), ("Within", lambda v: N["bw"][f"{v}|within"], "#8b1a1a", "v")]
fig, ax = plt.subplots(figsize=(8.5, 7.5)); y0 = np.arange(len(VARS))[::-1]; off = np.linspace(-0.34, 0.34, len(SER))
for (lab, f, col, mk), dy in zip(SER, off):
    first = True
    for v, y in zip(VARS, y0):
        d = f(v)
        if np.isnan(d["OR"]): continue
        ax.plot([d["lo"], d["hi"]], [y + dy, y + dy], color=col, lw=1.4)
        ax.plot(d["OR"], y + dy, marker=mk, color=col, mfc=col if d["cred"] else "white", ms=6.5, label=lab if first else None); first = False
ax.axvline(1, color="#555", lw=0.8, ls="--")
ax.set_yticks(y0); ax.set_yticklabels([C.LABEL[v] for v in VARS]); ax.set_xscale("log")
ticks = [0.92, 0.96, 1.0, 1.04, 1.08]; ax.set_xticks(ticks); ax.set_xticklabels([f"{t:.2f}" for t in ticks]); ax.minorticks_off()
vals = [f(v) for _, f, _, _ in SER for v in VARS]; lo = np.nanmin([d["lo"] for d in vals]); hi = np.nanmax([d["hi"] for d in vals]); assert 0.91 <= lo and hi <= 1.10, (lo, hi)
ax.set_xlim(0.91, 1.10); ax.set_xlabel("Odds ratio of influenza vaccination per standard deviation (log scale)"); ax.set_ylim(-0.7, len(VARS) - 0.3)
for s in ["top", "right"]: ax.spines[s].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.40, -0.22), ncol=3, frameon=False, fontsize=10)
plt.tight_layout(); save(fig, 2); plt.close(fig)
print("figs ok", round(lo, 3), round(hi, 3))
