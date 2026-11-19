#  문체: 단독 저자 무인칭(we 금지), 미국식 철자, Results = 숫자, Discussion = 통계 없음. 다른 호 원고는 인용·언급하지 않음.
import os, json, math
from decimal import Decimal, ROUND_HALF_UP
import content_intro_methods_flu as IM
HERE = os.path.dirname(os.path.abspath(__file__))
N = json.load(open(os.path.join(HERE, "numbers_flu.json"), encoding="utf-8"))
M = N["main"]
def r3(x):
    s = str(Decimal(repr(float(x))).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))
    return str(Decimal(repr(float(x))).quantize(Decimal("0.0001"), rounding=ROUND_HALF_UP)) if s == "1.000" and abs(float(x) - 1) >= 0.00005 else s
def o(v, d=None):
    d = d or M[v]; return f"OR {r3(d['OR'])}, 95% CrI {r3(d['lo'])}-{r3(d['hi'])}"
def oa(d): return f"{r3(d['OR'])}; {r3(d['lo'])}-{r3(d['hi'])}"
def pct(x): return f"{100*x:.0f}%"
def pp(v, d=None):
    p = N["sr_mean_all"] / 100; return abs(100 * p * (1 - p) * math.log((d or M[v])["OR"]))
WK, PA, FAC, VIL, ED, LIV, ONE = N["hyp"]
LABEL = {WK: "Walking", PA: "Moderate-or-higher physical activity", FAC: "Public sports facilities", VIL: "Village sports facilities", ED: "University-educated adults",
         LIV: "Basic livelihood recipients", ONE: "Single-person households",
         "urban_population_share": "Urban-zoned population share", "clinics_per_100k": "Clinics", "sex_ratio": "Sex ratio", "fiscal_independence": "Fiscal autonomy", "pop80_share": "Population aged 80+ years"}
UNIT = {WK: "%", PA: "%", FAC: "per 100,000", VIL: "per 100,000", ED: "%", LIV: "%", ONE: "%", "urban_population_share": "%", "clinics_per_100k": "per 100,000", "sex_ratio": "men per 100 women",
        "fiscal_independence": "%", "pop80_share": "%"}
EXPECT = {WK: "+", PA: "+", FAC: "+", VIL: "+", ED: "+", LIV: "+", ONE: "−"}
GROUPS = [("Physical activity and sports facilities", [WK, PA, FAC, VIL]), ("Social composition", [ED, LIV, ONE]),
          ("Controls", ["urban_population_share", "clinics_per_100k", "sex_ratio", "fiscal_independence", "pop80_share"])]
KH = lambda m, v: N["kh"][f"KH_{m}|{v}"]
BW = lambda v, c: N["bw"][f"{v}|{c}"]
SAME = lambda var, v: N["same"][f"{var}|{v}"]
SVY = lambda var, v: N["svy"][f"{var}|{v}"]
MV = N["mv"]; SR = N["sr"]; NSPEC = N["n_specs"]; ROB = N["robust"]
assert ROB == [WK], ROB
assert not KH("IV", ED)["cred"] and KH("IV", WK)["cred"] and not KH("IV", PA)["cred"]
SAMEN = "non_survey_covariates_only"
SIDO = {"강원도": "Gangwon", "경기도": "Gyeonggi", "서울시": "Seoul", "전라남도": "South Jeolla", "인천시": "Incheon", "경상남도": "South Gyeongsang", "경상북도": "North Gyeongsang", "대구시": "Daegu", "제주도": "Jeju",
        "충청남도": "South Chungcheong", "충청북도": "North Chungcheong", "전라북도": "North Jeolla", "대전시": "Daejeon", "울산시": "Ulsan", "광주시": "Gwangju", "부산시": "Busan", "세종시": "Sejong"}
def sido_list(d, k): return ", ".join(f"{SIDO[a]} ({b})" for a, b in list(d.items())[:k])
NEGC = [k for k, d in N["neg"].items() if d["cred"]]

TITLE = "Influenza vaccination coverage is slightly higher in Korean districts where more adults walk, but a regional divide remains unexplained"
RUNNING = "District influenza vaccination in Korea"
KEYWORDS = ["influenza vaccination", "walking", "health behavior", "geographic variation", "small-area estimation", "Korea"]

ABSTRACT = {
 "Background": "Influenza vaccination is a preventive health behavior that varies between places, but whether its local geography follows other health behaviors, such as physical activity, is unclear.",
 "Purpose": "To describe the space-time pattern of adult influenza vaccination across Korean districts and to test whether district physical activity, sports facilities, and social composition are associated with it.",
 "Methods": (f"Age-sex-standardized influenza vaccination rates for adults from the Korea Community Health Survey ({N['n_resp_total']/1e6:.1f} million responses) were modeled for 229 districts, 2011-2024, "
             f"with a spatiotemporal Fay-Herriot model with spatial and temporal random effects. Robustness was assessed across all {NSPEC} covariate subsets, space-time interaction models, and spatial-confounding adjustments."),
 "Results": (f"Mean coverage rose in every district between 2011-2015 and 2020-2024. Coverage clustered (Moran's I {N['moran']:.2f}), with hotspots concentrated in the west and coldspots in the southeast, and the covariates barely changed the spatial pattern. "
             f"Districts with more adults who walked regularly had higher coverage (odds ratio per standard deviation {oa(M[WK])}), the only association that met all prespecified robustness criteria; it was not present with lagged covariates. "
             f"Moderate-or-higher physical activity and a higher university-educated share (inversely, opposite to the hypothesis) were associated in the main model but attenuated in the type IV space-time interaction model. Two of three negative controls were also weakly associated."),
 "Conclusions": "Districts where walking is common have slightly higher vaccination coverage, consistent with clustering of health behaviors, but residual confounding cannot be excluded. A regional divide in coverage was not explained by the district characteristics studied.",
}

INTRO = IM.intro(N)
METHODS = IM.methods(N, NSPEC)

RESULTS = [
 ("Vaccination coverage and its geography",
  (f"The unweighted mean of district rates rose from {SR['2011']}% in 2011 to {SR['2021']}% in 2021, fell to {SR['2022']}% in 2022, and was {SR['2024']}% in 2024. Between the 2011-2015 and 2020-2024 window means, coverage increased in all 229 districts (median increase {N['chg_median']} percentage points; range {N['chg_min']}-{N['chg_max']}) (Table 1). "
   f"Fourteen-year district means ranged from {N['sr_dist_min']}% to {N['sr_dist_max']}%, and the 90th/10th percentile ratio was {N['p9010']:.2f}. "
   f"Coverage clustered spatially (pooled Moran's I {N['moran']:.2f}; annual values {N['moran_ymin']:.2f}-{N['moran_ymax']:.2f}, all p = 0.001), and district ranks were moderately persistent (Spearman {N['persist']:.2f} between 2011-2015 and 2020-2024; Kendall's W {N['kendallW']:.2f})."),
  (f"In the main model, {N['n_hot']} districts were hotspots, mainly in {sido_list(N['hot_sido'], 4)}, and {N['n_cold']} were coldspots, mainly in {sido_list(N['cold_sido'], 3)} (Fig. 1). "
   f"In a descriptive grouping, {N['hot_west']} of the {N['n_west']} districts in the western provinces (South Chungcheong, North and South Jeolla, Gwangju, and Sejong) were hotspots and {N['cold_se']} of the {N['n_se']} districts in the southeastern provinces (North and South Gyeongsang, Busan, Daegu, and Ulsan) were coldspots; the 14-year mean coverage was {N['sr_west']}% and {N['sr_se']}%, respectively. "
   f"The spatial effects of the main model correlated {N['sp_rho_vs_forced']:.2f} with those of the model with structural controls only, and their spread was almost unchanged (SD of log spatial odds ratios {N['sd_main']:.3f} in the main model vs {N['sd_forced']:.3f} in the controls-only model). "
   f"In the descriptive decomposition, the covariates accounted for about {N['vp']['hypothesis']:.0f}% of the summed variance, the spatial effect for {N['vp']['spatial']:.0f}%, and the year effects for {N['vp']['temporal']:.0f}%. The BYM2 mixing parameter was {N['phi_main']:.2f}, and residual Moran's I was {N['resid_moran']:.2f} (two-sided p = {N['resid_moran_p']:.2f})."),
 ),
 ("District covariates",
  (f"In the main model, walking ({o(WK)}) and moderate-or-higher physical activity ({o(PA)}) were associated with higher coverage, and the university-educated share with lower coverage ({o(ED)}), opposite to the hypothesized direction (Table 2, Fig. 2). "
   f"Public sports facilities ({o(FAC)}), village sports facilities, basic livelihood receipt, and single-person households were not credibly associated. "
   f"At the overall mean district rate ({N['sr_mean_all']}%), 1 SD of walking corresponded to about {pp(WK):.1f} percentage points higher coverage and 1 SD of the university-educated share to about {pp(ED):.1f} percentage points lower coverage. "
   f"The full covariate set improved WAIC over the controls alone ({N['waic_main']:.1f} vs {N['waic_forced']:.1f}); a four-covariate core set had a slightly lower WAIC ({N['sets_waic']['S5_core']:.1f})."),
  (f"Walking, physical activity, and education were credible in all of the {MV[WK]['n']} multiverse models containing each of them (Table S7). Only walking also remained credible in the Knorr-Held type II and type IV models (type IV, which had the lowest WAIC: {oa(KH('IV', WK))}); "
   f"physical activity ({oa(KH('IV', PA))}) and education ({oa(KH('IV', ED))}) attenuated toward 1 in the type IV model, whereas education remained credible in the type II and III models and physical activity in the type III model. "
   f"Public sports facilities ({oa(KH('IV', FAC))}) were credible in the type IV model and single-person households in the type III and IV models ({oa(KH('IV', ONE))} in type IV), but neither was credible in any multiverse model (Table S2). "
   f"Walking remained credible with restricted spatial regression ({oa(N['rsr'][WK])}) and Spatial+ ({oa(N['splus'][WK])}), and had credible between-district ({oa(BW(WK, 'between'))}) and within-district ({oa(BW(WK, 'within'))}) components. "
   f"Public sports facilities ({oa(BW(FAC, 'between'))}) and education ({oa(BW(ED, 'between'))}) had credible between-district components only."),
  (f"Without the two survey-based physical activity covariates, the other estimates were nearly unchanged and the spatial effects correlated {N['same_sum']['spearman_spatialOR_vs_main']:.3f} with those of the main model; added back alone, walking was {oa(SAME('non_survey_plus_walking_practice_rate', WK))}. "
   f"Adding province-by-year effects improved WAIC ({N['svy_sum']['sido_year_iid']['WAIC']:.1f} vs {N['waic_main']:.1f}) and left walking credible ({oa(SVY('sido_year_iid', WK))}); with an index of other KCHS items, walking was {oa(SVY('plus_kchs_index', WK))} and the index itself was associated with coverage ({oa(SVY('plus_kchs_index', 'kchs_index'))}). "
   f"Walking was not credible with 1-year lagged covariates ({oa(N['lag']['1|'+WK])}), in the 2016-2019 window ({oa(N['win']['W2_2016-19|'+WK])}), or in the 2017-2024 subset ({oa(N['mfg']['main_2017_2024|'+WK])}), and its region-specific slopes were credible in the capital region ({oa(N['region']['walking_practice_rate|Capital'])}) and Honam ({oa(N['region']['walking_practice_rate|Honam'])}) but not elsewhere (Table S9). "
   f"Results were similar when 2020-2021 were excluded (walking {oa(N['ex2021'][WK])}), across provinces left out (walking credible in {N['loo'][WK]['n_cred']} of {N['loo'][WK]['n']} refits), likelihoods, graphs, priors, and influential observations (Tables S1-S6); the manufacturing share was not associated ({oa(N['mfg']['plus_manufacturing_2017_2024|mfg'])}). "
   f"In post hoc models for the inverse education association (education with the controls only, without the other facility and social covariates, and without the urban-zoned share), estimates were {r3(min(d['OR'] for d in N['posthoc'][:4]))}-{r3(max(d['OR'] for d in N['posthoc'][:4]))}, all credible; in the prespecified crude-rate model with an age-65 control, it was {oa(N['lik']['crude_rate_FH_plus_old65|'+ED])}."),
  (f"Two of three negative controls had small credible associations: rear-seat belt use, a KCHS item ({oa(N['neg']['rear_seatbelt_wearing_rate'])}), and fires ({oa(N['neg']['fires_per_10000'])}); multicultural marriages were not associated ({oa(N['neg']['multicultural_marriage_share'])}) (Table S8)."),
 ),
]

DISCUSSION = [
 ("Main findings",
  ("Across 229 Korean districts and 14 years, adult influenza vaccination coverage rose everywhere, while hotspots of higher-than-expected coverage were concentrated in the western provinces and coldspots in the southeast. "
   "District physical activity, sports facilities, and social composition explained little of this geography. "
   "Walking was the only covariate that met the prespecified robustness criteria, but its association was small, contemporaneous rather than lagged, and absent in some periods and regions. "
   "Physical activity in general and the university-educated share were associated in the main model, the latter opposite to the hypothesis, but both attenuated in the best-fitting space-time interaction model."),
 ),
 ("Interpretation",
  ("The co-occurrence of walking and vaccination at the district level is consistent with a social-ecological view in which preventive behaviors share community and service-level determinants, and with a KCHS analysis in which regular walkers were more likely to be vaccinated [@yang2014influenza; @tian2019characteristics]. "
   "Shared local conditions, such as health promotion by public health centers or community programs that reach older adults, could produce such co-occurrence. "
   "Because walking and vaccination were reported by the same respondents and the association was contemporaneous, part of it may also reflect shared reporting tendencies or the same residents being more engaged with health services, a pattern related to the healthy-user effect in vaccine studies [@jackson2006evidence]. "
   "Another KCHS item index and rear-seat belt use, a negative control from the same survey, were also associated with coverage, so a shared survey-level component is plausible; the walking association persisted after adjustment for these items, but the design cannot separate these explanations, and the effect size was small."),
  ("The inverse association with the university-educated share was opposite to the hypothesis and was examined post hoc. "
   "It resembles individual-level KCHS findings in which lower education was associated with higher vaccination [@yang2014influenza] and an area-level privilege paradox described for childhood vaccination in Australia [@bryden2019privilege], and it contrasts with studies from England and the United States, where deprived areas had lower coverage [@hungerford2018influenza; @strully2022county; @green2026learning]. "
   "Because the outcome was age-standardized, differences in the share of older adults do not explain it. Social gradients in vaccination depend on how vaccination is financed and delivered [@nagata2013social; @jemna2024socio], but this study did not measure delivery, and the attenuation of the association in the type IV model shows that it is sensitive to space-time structure; as an ecological association, it does not support inference about individuals' education."),
  ("The regional pattern was little changed by the covariates, by exclusion of the pandemic years, or by adjustment for spatial confounding, and adding province-by-year effects improved model fit. "
   "District ranks were only moderately persistent, however, and hotspots summarize 14 years, so the pattern describes the period as a whole rather than a fixed divide. "
   "Possible explanations include provincial differences in vaccination campaigns, outreach by public health centers, or the organization of private clinics that deliver most vaccinations, as well as district characteristics not measured here and measurement differences between survey areas. "
   "Administrative vaccination records could test these explanations directly."),
 ),
 ("Strengths and limitations",
  ("Strengths include complete coverage of all districts for 14 years, a model that uses the published sampling variances, and robustness checks covering model choice, space-time interaction, spatial confounding, and survey effects. "
   "Several limitations apply. Vaccination was self-reported; self-report has high sensitivity but can overstate vaccination [@macdonald1999validation; @daley2024influenza]. "
   "The outcome combines adults who are eligible for free vaccination with those who are not, and age-specific district rates were not available, so the analysis cannot show which age groups drive the geography. "
   "Walking and physical activity came from the same respondents as the outcome, and two negative controls, including a KCHS item, were weakly associated, so residual confounding or shared measurement error could account for associations of the size observed. "
   "The university-educated share was interpolated between censuses and carried forward after 2020, and single-person household shares were filled with district means for 2011-2014, which limits the within-district comparisons for these covariates. "
   "Individual characteristics such as race, ethnicity, and individual income were not available; the analysis used district aggregates from a frozen analytic copy of the district database (November 2026) checked against the published tables. "
   "Finally, the associations are ecological and cannot be interpreted as individual effects [@greenland1989ecological]."),
 ),
]
CONCLUSION = ("Adult influenza vaccination coverage in Korea rose in every district between 2011-2015 and 2020-2024, with higher-than-expected coverage concentrated in the west and lower in the southeast. Districts where more adults walked regularly had slightly higher coverage, consistent with clustering of health behaviors, whereas the regional pattern was not explained by the district characteristics studied. "
              "How vaccination is organized and delivered across provinces deserves direct study.")
FIGCAP = {
 1: ("District influenza vaccination coverage and spatial effects, Korea, 2011-2024. (A) Age-sex-standardized coverage among adults, 2011-2015 mean. (B) 2020-2024 mean, same color scale. "
     "(C) Spatial odds ratios from the model with structural controls only. (D) Spatial odds ratios from the main model with district covariates. Solid outlines, hotspots (posterior probability that the spatial odds ratio exceeds 1 above 0.95); dashed outlines, coldspots (below 0.05).",
     "Four maps of South Korea's 229 districts. Panels A and B show vaccination coverage rising from about 30-40% in 2011-2015 to about 40-52% in 2020-2024. Panels C and D show nearly identical spatial patterns, with higher-than-expected coverage concentrated in the western provinces and lower-than-expected coverage in the southeastern provinces."),
 2: ("Odds ratios of influenza vaccination per 1 standard deviation of district covariates across models. Filled symbols, 95% credible interval excludes 1. "
     "Main, main Fay-Herriot model; Knorr-Held IV, type IV space-time interaction; Spatial+, covariates residualized on spatial eigenvectors; Without survey covariates, model excluding walking and moderate-or-higher physical activity; Between and Within, district-mean and annual-deviation components.",
     "Forest plot of odds ratios for seven district covariates under six model variants. Walking has odds ratios above 1 in all variants. Moderate-or-higher physical activity is above 1 except in the Knorr-Held IV model. University education is below 1 in most variants, with the interval crossing 1 in the Knorr-Held IV and within-district models. Public sports facilities and single-person households are credible in some variants (Knorr-Held IV and between-district), and village sports facilities and welfare receipt have intervals crossing 1."),
}
