#!/usr/bin/env python3
#  입력: 02_패널/MASTER_시군구_연도_wide.parquet
#   주 결과 kchs__DT_INFLUENZA(최근 1년 인플루엔자 예방접종률, 19세 이상) SR·SR_SE·N·CR·CR_SE, 2011-2024 전 칸
#   조사 지역 효과(S21) kchs__DT_DIA_KIDNEY(DMKID), kchs__DT_117075_H_BELT_BACK(뒷좌석 안전벨트), 영양표시 활용률 → logit SR
#  산출: 02_패널/FLU_fh_outcomes.csv
import os, unicodedata
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.normpath(os.path.join(HERE, ".."))
M = pd.read_parquet(os.path.join(ROOT, "02_패널", "MASTER_시군구_연도_wide.parquet")); M.columns = [NFC(c) for c in M.columns]
M = M[(M.year >= 2011) & (M.year <= 2024)].copy(); M["fbd_region"] = M.fbd_region.map(NFC)
tb = "kchs__DT_INFLUENZA__"
NUT = [c for c in M.columns if "NUT_LABEL" in c and c.endswith("__SR")]; assert len(NUT) == 1, NUT
d = pd.DataFrame({"year": M.year.values, "fbd_region": M.fbd_region.values, "population": M.population.astype(float).values})
for it in ["SR", "SR_SE", "N", "CR", "CR_SE"]: d[it] = M[tb + it].astype(float).values
assert len(d) == 229 * 14 and d[["SR", "SR_SE", "N", "CR", "CR_SE"]].notna().all().all() and d.SR.between(0, 100, inclusive="neither").all()
n0 = int((d.SR_SE <= 0).sum()); d["SE_used"] = d.SR_SE.clip(lower=0.05); d["CRSE_used"] = d.CR_SE.clip(lower=0.05)
lg = lambda p: np.log(p / (1 - p))
p = d.SR / 100; d["y_logit"] = lg(p); d["v_logit"] = (d.SE_used / 100) ** 2 / (p * (1 - p)) ** 2
pc = d.CR / 100; d["y_logit_cr"] = lg(pc); d["v_logit_cr"] = (d.CRSE_used / 100) ** 2 / (pc * (1 - pc)) ** 2
d["n_eff"] = np.maximum(np.round(p * (1 - p) / (d.SE_used / 100) ** 2), 10); d["y_eff"] = np.round(p * d.n_eff)
for k, c in {"dmkid": "kchs__DT_DIA_KIDNEY__SR", "belt": "kchs__DT_117075_H_BELT_BACK__SR", "nut": NUT[0]}.items():
    x = M[c].astype(float).values / 100; d[f"{k}_SR"] = x * 100
    d[f"{k}_logit"] = np.where((x > 0) & (x < 1), lg(np.clip(x, 1e-6, 1 - 1e-6)), np.nan)
print(f"[99FL] rows {len(d)} · SE≤0 칸 {n0} · SR {d.SR.min():.1f}-{d.SR.max():.1f} · SE 중앙 {d.SR_SE.median():.2f} · n_eff 중앙 {d.n_eff.median():.0f} · N 중앙 {d.N.median():.0f} · N 합 {d.N.sum():,.0f}")
for k in ["dmkid", "belt", "nut"]: print(f"[99FL] {k} 관측 연도별", d.groupby("year")[f"{k}_logit"].count().tolist())
print("[99FL] 연도별 인구가중 SR:", d.groupby("year").apply(lambda g: round(np.average(g.SR, weights=g.population), 2)).to_dict())
print("[99FL] 연도별 SR 단순평균:", d.groupby("year").SR.mean().round(2).to_dict())
fp = os.path.join(ROOT, "02_패널", "FLU_fh_outcomes.csv"); d.to_csv(fp, index=False, encoding="utf-8-sig"); print("[99FL] 저장", fp, d.shape)
