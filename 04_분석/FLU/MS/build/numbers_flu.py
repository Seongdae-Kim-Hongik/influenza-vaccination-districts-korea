import os, json, unicodedata
import pandas as pd, numpy as np
NFC = lambda s: unicodedata.normalize("NFC", str(s))
HERE = os.path.dirname(os.path.abspath(__file__)); AN = os.path.normpath(os.path.join(HERE, "..", "..", "..")); ROOT = os.path.dirname(AN)
R = os.path.join(AN, "FLU", "results_std3")
rd = lambda f: pd.read_csv(os.path.join(R, f))
D = lambda r: dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), cred=bool(r.cred))
HYP = ["walking_practice_rate", "moderate_pa_rate", "public_sports_facilities_per_100k", "village_sports_facilities_per_100k", "tertiary_educated_share", "basic_livelihood_rate", "single_person_household_share"]
N = {"hyp": HYP}
FO = pd.read_csv(os.path.join(ROOT, "02_패널", "FLU_fh_outcomes.csv"), encoding="utf-8-sig")
assert len(FO) == 3206 and FO.SR.notna().all()
N["n_cells"] = len(FO); N["n_resp_total"] = int(FO.N.sum()); N["n_resp_median"] = int(FO.N.median()); N["n_resp_min"] = int(FO.N.min())
ym = FO.groupby("year").SR.mean(); N["sr"] = {int(y): round(v, 1) for y, v in ym.items()}
p = FO.groupby("fbd_region").SR.mean(); N["p10"] = round(p.quantile(.1), 1); N["p90"] = round(p.quantile(.9), 1); N["sr_dist_min"] = round(p.min(), 1); N["sr_dist_max"] = round(p.max(), 1)
_a = FO[FO.year <= 2015].groupby("fbd_region").SR.mean(); _b = FO[FO.year >= 2020].groupby("fbd_region").SR.mean(); _x = _b - _a; N["n_decline"] = int((_x < 0).sum()); N["chg_median"] = round(float(_x.median()), 1); N["chg_max"] = round(float(_x.max()), 1)
N["se_median"] = round(FO.SR_SE.median(), 1); N["sr_mean_all"] = round(float(FO.SR.mean()), 1)
s0 = rd("S0_screen.csv").iloc[0]; N["moran"] = round(s0.moran_I_pooled_logitSR, 2); N["moran_ymin"] = round(s0.moran_I_year_min, 2); N["moran_ymax"] = round(s0.moran_I_year_max, 2)
N["persist"] = round(s0.persistence_rank_rho, 2); N["phi_forced"] = round(s0.phi_forced_only, 2); N["p9010"] = round(s0.P90_P10_SR, 2)
N["kendallW"] = round(rd("S6_kendallW.csv").kendall_W[0], 2)
my = rd("S6_moran_year.csv").set_index("year"); N["moran_year"] = {int(y): round(v, 2) for y, v in my.moran_I.items()}; N["sr_popw"] = {int(y): round(v, 1) for y, v in my.SR_popweighted.items()}
f3 = rd("S3_main_fit.csv").iloc[0]; N["phi_main"] = round(f3.phi, 2); N["n_hot"] = int(f3.n_hotspot); N["n_cold"] = int(f3.n_coldspot); N["resid_moran"] = round(f3.resid_moran_I, 2); N["resid_moran_p"] = round(f3.resid_moran_p, 2)
N["largest_hot"] = int(f3.largest_hot); N["largest_cold"] = int(f3.largest_cold); N["or_sp_min"] = round(f3.OR_min, 2); N["or_sp_max"] = round(f3.OR_max, 2); N["waic_main"] = float(f3.WAIC)
sp = rd("S3_spatial_OR.csv"); sp["sido"] = sp.sido.map(NFC); sp["fbd_region"] = sp.fbd_region.map(NFC)
N["hot_sido"] = sp[sp["class"] == "hotspot"].sido.value_counts().to_dict(); N["cold_sido"] = sp[sp["class"] == "coldspot"].sido.value_counts().to_dict()
iscounty = sp.fbd_region.str.endswith("군"); N["hot_county"] = int((iscounty & (sp["class"] == "hotspot")).sum()); N["cold_county"] = int((iscounty & (sp["class"] == "coldspot")).sum())
sf = rd("S0_forced_only_spatialOR.csv"); N["hot_forced0"] = int((sf["class"] == "hotspot").sum()); N["cold_forced0"] = int((sf["class"] == "coldspot").sum())
vf = rd("S3_vs_forced_only.csv").iloc[0]; N["sp_rho_vs_forced"] = round(vf.spearman_spatial_OR, 2); N["sd_main"] = round(vf.sd_logOR_main, 3); N["sd_forced"] = round(vf.sd_logOR_forced, 3)
N["hot_forced"] = int(vf.hot_forced); N["cold_forced"] = int(vf.cold_forced); N["waic_forced"] = float(vf.WAIC_forced)
vp = rd("S3_varpart.csv").set_index("component").share; N["vp"] = {k: round(100 * v, 0) for k, v in vp.items()}
m = rd("S3_main_OR.csv").set_index("cov"); N["main"] = {k: D(r) for k, r in m.iterrows()}
nc = rd("S3_negcontrol.csv").set_index("cov"); N["neg"] = {k: D(r) for k, r in nc.iterrows()}
s2 = rd("S2_setcompare.csv").set_index("set"); N["sets_waic"] = {k: float(v) for k, v in s2.WAIC.items()}
if os.path.exists(os.path.join(R, "S7_multiverse_summary.csv")): mv = rd("S7_multiverse_summary.csv").set_index("cov"); N["mv_provisional"] = False
else:   # 체인 도중 시험 빌드용(최종 빌드 전 assert로 막음)
    import glob as _g; _p = pd.concat([pd.read_csv(f) for f in _g.glob(os.path.join(R, "S7_multiverse_results_part_*.csv"))]); _s = {k: np.sign(np.log(v["OR"])) for k, v in N["main"].items()}
    mv = pd.DataFrame([dict(cov=v, cred_frac=_p[_p["cov"] == v].cred.mean(), sign_agree=(np.sign(np.log(_p[_p["cov"] == v].OR)) == _s[v]).mean(), n_specs=int((_p["cov"] == v).sum()), OR_median=_p[_p["cov"] == v].OR.median(), OR_min=_p[_p["cov"] == v].OR.min(), OR_max=_p[_p["cov"] == v].OR.max()) for v in HYP]).set_index("cov")
    mv["robust"] = (mv.cred_frac >= .75) & (mv.sign_agree >= .95); N["mv_provisional"] = True
mv = mv; N["mv"] = {k: dict(cred=float(r.cred_frac), sign=float(r.sign_agree), robust=bool(r.robust), n=int(r.n_specs), med=float(r.OR_median), mn=float(r.OR_min), mx=float(r.OR_max)) for k, r in mv.iterrows()}
N["n_specs"] = int(pd.read_csv(os.path.join(R, "S7_multiverse_specs.csv")).shape[0])
kh = rd("S8_KH.csv"); N["kh"] = {f"{r.model}|{r['cov']}": D(r) for _, r in kh.iterrows()}; N["kh_waic"] = kh.drop_duplicates("model").set_index("model").WAIC.to_dict()
bw = rd("S9_between_within.csv"); N["bw"] = {f"{r['cov']}|{r.component}": D(r) for _, r in bw.iterrows()}
w4 = rd("S4_window_OR.csv"); N["win"] = {f"{r.window}|{r['cov']}": D(r) for _, r in w4.iterrows()}
c5 = rd("S5_covid.csv"); N["covid"] = {r["cov"]: D(r) for _, r in c5.iterrows()}
lg = rd("S10_lag.csv"); N["lag"] = {f"{int(r.lag)}|{r['cov']}": D(r) for _, r in lg.iterrows()}
s11 = rd("S11_leave_one_sido.csv"); N["loo"] = {v: dict(n=int((s11["cov"] == v).sum()), n_cred=int(s11[s11["cov"] == v].cred.sum()), n_sign=int((np.sign(np.log(s11[s11["cov"] == v].OR)) == np.sign(np.log(N["main"][v]["OR"]))).sum()), mn=float(s11[s11["cov"] == v].OR.min()), mx=float(s11[s11["cov"] == v].OR.max())) for v in HYP}
s12 = rd("S12_likelihood.csv"); N["lik"] = {f"{r.likelihood}|{r['cov']}": D(r) for _, r in s12.iterrows()}
s13 = rd("S13_graph_prior.csv"); N["gp"] = {f"{r.variant}|{r['cov']}": D(r) for _, r in s13.iterrows()}
s14 = rd("S14_excl2020_21.csv"); N["ex2021"] = {r["cov"]: D(r) for _, r in s14.iterrows()}; ms = rd("S14_excl2020_21_summary.csv").iloc[0]; N["ex2021_sum"] = {k: (float(v) if not isinstance(v, str) else v) for k, v in ms.items()}
s16 = rd("S16_RSR.csv"); N["rsr"] = {r["cov"]: D(r) for _, r in s16.iterrows()}
s18 = rd("S18_influence.csv"); N["infl"] = {r["cov"]: D(r) for _, r in s18.iterrows()}; N["infl_variant"] = s18.variant.iloc[0]
s21 = rd("S21_survey_area.csv"); N["svy"] = {f"{r.variant}|{r['cov']}": D(r) for _, r in s21.iterrows()}; N["svy_sum"] = rd("S21_survey_area_summary.csv").set_index("variant").to_dict(orient="index")
s22 = rd("S22_spatialplus.csv"); N["splus"] = {r["cov"]: D(r) for _, r in s22.iterrows()}
s23 = rd("S23_same_sample.csv"); N["same"] = {f"{r.variant}|{r['cov']}": D(r) for _, r in s23.iterrows()}; N["same_sum"] = rd("S23_same_sample_summary.csv").iloc[0].to_dict()
s24 = rd("S24_manufacturing.csv"); N["mfg"] = {f"{r.variant}|{r['cov']}": D(r) for _, r in s24.iterrows()}
if os.path.exists(os.path.join(R, "S17_region_slopes.csv")):
    s17 = rd("S17_region_slopes.csv")
    if "covariate" in s17: N["region"] = {f"{r.covariate}|{r.region}": D(r) for _, r in s17.iterrows()}
# 시군구 변화(2011-15 → 2020-24)·지역군 집계·44 사후
a = FO[FO.year <= 2015].groupby("fbd_region").SR.mean(); b = FO[FO.year >= 2020].groupby("fbd_region").SR.mean(); ch = b - a
N["n_increase"] = int((ch > 0).sum()); N["chg_median"] = round(float(ch.median()), 1); N["chg_min"] = round(float(ch.min()), 1); N["chg_max"] = round(float(ch.max()), 1)
WEST = ["충청남도", "전라북도", "전라남도", "광주시", "세종시"]; SE = ["경상남도", "경상북도", "부산시", "대구시", "울산시"]
N["hot_west"] = int(sp[(sp["class"] == "hotspot") & sp.sido.isin(WEST)].shape[0]); N["cold_se"] = int(sp[(sp["class"] == "coldspot") & sp.sido.isin(SE)].shape[0])
N["n_west"] = int(sp.sido.isin(WEST).sum()); N["n_se"] = int(sp.sido.isin(SE).sum())
pm = FO.groupby("fbd_region").SR.mean(); gW = pm[pm.index.map(lambda r: any(r.startswith(x) for x in WEST))]; gS = pm[pm.index.map(lambda r: any(r.startswith(x) for x in SE))]
N["sr_west"] = round(float(gW.mean()), 1); N["sr_se"] = round(float(gS.mean()), 1)
p44 = rd("S44_reverse_posthoc.csv"); N["posthoc"] = [D(r) for _, r in p44.iterrows()]
N["svy_sido_waic"] = float(N["svy_sum"]["sido_year_iid"]["WAIC"])
# 강건 판정: 다중우주 credible ≥75%·부호 ≥95% + KH II·IV credible 같은 방향
def same(d, v): return d["cred"] and np.sign(np.log(d["OR"])) == np.sign(np.log(N["main"][v]["OR"]))
N["robust"] = [v for v in HYP if N["mv"][v]["robust"] and same(N["kh"][f"KH_II|{v}"], v) and same(N["kh"][f"KH_IV|{v}"], v)]
# 연간 응답·2011→2024 감소
N["sr_drop_pp"] = round(N["sr"][2011] - N["sr"][2024], 1)
json.dump(N, open(os.path.join(HERE, "numbers_flu.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("robust:", N["robust"]); print("mv:", {k: (round(v["cred"], 2), round(v["sign"], 2)) for k, v in N["mv"].items()})
