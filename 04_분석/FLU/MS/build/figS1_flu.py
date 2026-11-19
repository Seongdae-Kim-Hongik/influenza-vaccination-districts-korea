import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import content_flu as C
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(os.path.dirname(HERE), "out", "Figures"); os.makedirs(OUT, exist_ok=True)
N = C.N; plt.rcParams.update({"font.family": "Arial", "font.size": 11})
yr = sorted(int(y) for y in N["moran_year"]); mi = [N["moran_year"][y] if y in N["moran_year"] else N["moran_year"][str(y)] for y in yr]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.plot(yr, mi, color="#1f4e79", lw=1.8, marker="o", ms=6)
ax.axhline(N["moran"], color="#8b1a1a", lw=1.0, ls="--", label=f"14-year district mean, pooled ({N['moran']:.2f})")
ax.axhline(0, color="#555", lw=0.8)
ax.set_xticks(yr); ax.set_xticklabels([str(y) for y in yr], rotation=45); ax.set_ylim(0, 0.5)
ax.set_xlabel("Year"); ax.set_ylabel("Moran's I"); ax.legend(frameon=False, loc="upper left")
for s in ("top", "right"): ax.spines[s].set_visible(False)
for ext, dpi in (("png", 300), ("tiff", 600)): fig.savefig(os.path.join(OUT, f"FigureS1.{ext}"), dpi=dpi, pil_kwargs={"compression": "tiff_lzw"} if ext == "tiff" else None, bbox_inches="tight")
print("FigureS1 saved", list(zip(yr, mi)))
