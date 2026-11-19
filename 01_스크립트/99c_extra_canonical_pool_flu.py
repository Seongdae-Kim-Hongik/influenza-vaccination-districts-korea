#!/usr/bin/env python3
#  산출: 02_패널/EXT_canonical_pool_FLU_261119.csv
import os, unicodedata
import pandas as pd
NFC = lambda s: unicodedata.normalize("NFC", str(s))
PAN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "02_패널")
d = pd.read_csv(os.path.join(PAN, "EXT_canonical_pool_STRS_261113.csv"), encoding="utf-8-sig"); d["fbd_region"] = d.fbd_region.map(NFC)
d = d[(d.year >= 2011) & (d.year <= 2024)][["year", "fbd_region", "basic_livelihood_rate", "tertiary_educated_share"]].sort_values(["fbd_region", "year"])
assert d.fbd_region.nunique() == 229 and d.groupby("fbd_region").year.nunique().eq(14).all() and not d.duplicated(["year", "fbd_region"]).any() and d.notna().all().all()
cw = d[(d.fbd_region == NFC("경상남도창원시")) & d.year.isin([2021, 2022, 2023])].basic_livelihood_rate; assert cw.between(2.3, 2.4).all(), cw.tolist()
fp = os.path.join(PAN, "EXT_canonical_pool_FLU_261119.csv"); d.to_csv(fp, index=False, encoding="utf-8-sig"); print("[99c-flu] 저장", fp, d.shape); print(d.describe().T[["min", "50%", "max"]].round(2))
