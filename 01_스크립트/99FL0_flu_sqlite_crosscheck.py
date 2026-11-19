#!/usr/bin/env python3
import os, sqlite3, unicodedata, json, hashlib
import numpy as np, pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DB = os.path.expandvars("$ST_GPT_DATASET/02_체육건강_시공간_GPT버전/02_체육건강_시공간_GPT버전.sqlite")
con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
cols_db = [r[1] for r in con.execute("PRAGMA table_info(region_year)")]; cols_db_n = {NFC(c): c for c in cols_db}
S = json.load(open(os.path.join(ROOT, "04_분석/FLU/config_FLU_std3.json"), encoding="utf-8"))
raw = {v["key"]: (v.get("raw_key") or v.get("column")) for v in S["vars"].values()}
M = pd.read_parquet(os.path.join(ROOT, "02_패널/MASTER_시군구_연도_wide.parquet")); M.columns = [NFC(c) for c in M.columns]; M["fbd_region"] = M.fbd_region.map(NFC)
M = M[(M.year >= 2011) & (M.year <= 2024)]
want = ["population"] + [f"kchs__DT_INFLUENZA__{i}" for i in ["SR", "SR_SE", "N", "CR", "CR_SE"]] + [NFC(r) for r in raw.values() if r and not NFC(r).startswith("ext__")]
want = list(dict.fromkeys(want)); rows = []
q = ", ".join(['"year"', '"fbd_region"'] + [f'"{cols_db_n[c]}"' for c in want if c in cols_db_n])
DBd = pd.read_sql(f"SELECT {q} FROM region_year WHERE year BETWEEN 2011 AND 2024", con); DBd.columns = [NFC(c) for c in DBd.columns]; DBd["fbd_region"] = DBd.fbd_region.map(NFC)
J = M[["year", "fbd_region"] + [c for c in want if c in M]].merge(DBd, on=["year", "fbd_region"], how="outer", suffixes=("_ours", "_db"), indicator=True)
print("키 결합", J._merge.value_counts().to_dict())
for c in want:
    if c not in cols_db_n: rows.append(dict(variable=c, status="SQLite에 열 없음")); continue
    a, b = J[c + "_ours"].astype(float), J[c + "_db"].astype(float)
    both = a.notna() & b.notna(); d = (a - b).abs()[both]
    rows.append(dict(variable=c, n_ours=int(a.notna().sum()), n_db=int(b.notna().sum()), n_both=int(both.sum()), na_mismatch=int((a.isna() != b.isna()).sum()),
                     n_diff_gt_1e9=int((d > 1e-9 * np.maximum(1, b[both].abs())).sum()), max_abs_diff=float(d.max()) if len(d) else np.nan,
                     status="동일" if (a.isna() == b.isna()).all() and (d <= 1e-9 * np.maximum(1, b[both].abs())).all() else "차이"))
out = pd.DataFrame(rows)
# 창원 2023 인구
cw = J[(J.year == 2023) & (J.fbd_region == NFC("경상남도창원시"))][["population_ours", "population_db"]]
print("창원 2023 인구 ours(마스터)/db", cw.values.tolist(), "· 패널 population_override = 1,009,038")
P = pd.read_csv(os.path.join(ROOT, "04_분석/FLU/panel_FLU_std3.csv"), encoding="utf-8-sig"); P["fbd_region"] = P.fbd_region.map(NFC)
print("FLU 패널 창원 2023 population", P.loc[(P.year == 2023) & (P.fbd_region == NFC("경상남도창원시")), "population"].tolist())
# ext 변수: 정본 CSV 기반 → SQLite region_year 에 대응 열 유무
ext = [c for c in cols_db if any(k in NFC(c) for k in ["livelihood", "tertiary", "restaurant", "instant", "drinking", "entertain", "singing", "수급", "대졸", "음식점"])]
print("SQLite region_year ext 후보 열", ext[:20])
dp = [r[1] for r in con.execute("PRAGMA table_info(disease_panel)")]; print("disease_panel 열 중 ext 후보", [c for c in dp if any(k in c for k in ["livelihood", "tertiary", "restaurant", "instant"])])
fp = os.path.join(ROOT, "04_분석/FLU/00_기준SQLite_대조_261119.csv"); out.to_csv(fp, index=False, encoding="utf-8-sig"); print(out.to_string()); print("저장", fp)
h = hashlib.sha256(open(DB, "rb").read()).hexdigest(); print("SQLite SHA-256", h)
