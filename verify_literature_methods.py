"""Do the methods of the literature review (Literature_metanalysis/review.md)
beat the estimates the tool already ships?

Methods the physics brief has already scored (inverse Arya-Paris, inverse
Saxton, the Fredlund back-analyses, Peters et al. 2023, Brutsaert 1967) are
not repeated. What is left that a single retention curve, with bulk density
where known, can drive:

Ks (against the shipped capillary bundle, Marshall 1958 capped at 50 cm,
divided by 1.09, and against the tool's kNN):
  cpa        critical-path analysis as Ghanbarian & Patrignani (2026, eq. 1-2)
             use it: Ks = (rho g / mu) dc^2 / (C F), C = 53.5 (Skaggs 2011),
             dc the pore-throat diameter at the inflection point of the vG
             curve, F from porosity by Ghanbarian et al. (2014) with the
             critical porosity at theta_r and the crossover at 0.75. The
             inflection is taken on log h (h = m^(-1/n) / alpha, the
             Dexter S-index point) and, as a variant, on linear h
             (h = m^(1/n) / alpha).
  jks        Johnson et al. (1986): Ks = (rho g / mu) Lambda^2 / (8 F),
             Lambda = 0.21 dc (their eq. 5).
  guarracino Guarracino (2007): Ks = 4.65e4 phi alpha^2 cm/d, alpha in 1/cm.
  mishra     Mishra & Parker (1990): Ks = 108 (theta_s - theta_r)^2.5
             alpha^2 cm/s, alpha in 1/cm.
  ahuja_kc   effective-porosity Kozeny-Carman, Ahuja et al. (1989):
             Ks = 764.5 phi_e^3.29 cm/h, phi_e = phi - theta(33 kPa).
  Porosity phi is theta_s; where bulk density is known, 1 - bd / 2.65 is
  scored as a variant. Each form is also scored matched: divided by the
  median over the other sources of each source's median ratio to the
  measurement (one vote per source), as verify_ks_physical.py does.

Texture (against the tool's fractions, study held out, depth and bulk
density supplied; figures/cache/predictive.csv from verify_predictive.py):
  dry-end clay  clay from the water held at 1500 kPa, the adsorption route
             (Tuller & Or 2005; Wuddivira et al. 2012; Arthur et al.). The
             fixed form is the NRCS laboratory rule, clay = 2.5 x gravimetric
             1500 kPa water (%), with and without organic carbon removed;
             the fitted form is a straight line of clay on that water,
             fitted on the other studies.

Usage:  python verify_literature_methods.py
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

import ks_physical as kp
import swcc_texture as st
from verify_common import GROUP

RHO_G_MU = kp.DENSITY * kp.GRAVITY / kp.VISCOSITY         # 1/(m s)
C_CPA = 53.5
PHI_X = 0.75
LOG2 = np.log10(2.0)


# ---------------------------------------------------------------------------
# Ks forms, all in cm/h; alpha in 1/kPa as stored
# ---------------------------------------------------------------------------

def inflection_h_cm(alpha_kpa, n, scale="log"):
    m = 1.0 - 1.0 / n
    a_cm = alpha_kpa * 0.0980665                         # 1/kPa -> 1/cm
    return (m ** (-1.0 / n) if scale == "log" else m ** (1.0 / n)) / a_cm


def formation_factor(phi, phi_c):
    """Ghanbarian et al. (2014), as eq. 2 of Ghanbarian & Patrignani (2026);
    reduces to Archie's F = phi^-2 for phi_c = 0, phi_x = 1."""
    phi_c = np.minimum(phi_c, 0.9 * phi)
    lo = (1 - phi_c) * (PHI_X - phi_c) / (phi - phi_c) ** 2
    hi = (1 - phi_c) / (phi - phi_c)
    return np.where(phi < PHI_X, lo, hi)


def cpa_ks(thr, phi, alpha, n, scale="log", jks=False):
    h_m = inflection_h_cm(alpha, n, scale) / 100.0
    d = 4.0 * kp.SURFACE_TENSION / (kp.DENSITY * kp.GRAVITY * h_m)   # m
    f = formation_factor(phi, thr)
    k = (0.21 * d) ** 2 / (8.0 * f) if jks else d ** 2 / (C_CPA * f)
    return k * RHO_G_MU * kp.M_S_TO_CM_H


def guarracino_ks(phi, alpha):
    return 4.65e4 * phi * (alpha * 0.0980665) ** 2 / 24.0


def mishra_ks(thr, ths, alpha):
    return 108.0 * (ths - thr) ** 2.5 * (alpha * 0.0980665) ** 2 * 3600.0


def ahuja_ks(thr, ths, alpha, n, phi):
    phi_e = np.clip(phi - st.vg_theta(33.0, thr, ths, alpha, n), 1e-4, None)
    return 764.5 * phi_e ** 3.29


def shipped_ks(thr, ths, alpha, n):
    return kp.marshall_ks(thr, ths, alpha, n) / kp.MATCHING_FACTOR


def forms(d):
    thr, ths, a, n = (d[c].to_numpy(float) for c in ("thetar", "thetas", "alpha_kpa", "n"))
    phi_bd = (1.0 - d.bd / 2.65).where(d.bd.between(0.1, 2.3)).to_numpy(float)
    return {
        "shipped: Marshall, 50 cm cap": shipped_ks(thr, ths, a, n),
        "CPA, log-h inflection": cpa_ks(thr, ths, a, n),
        "CPA, linear-h inflection": cpa_ks(thr, ths, a, n, "linear"),
        "CPA, log-h, phi from bd": cpa_ks(thr, phi_bd, a, n),
        "JKS (Johnson 1986)": cpa_ks(thr, ths, a, n, jks=True),
        "Guarracino 2007": guarracino_ks(ths, a),
        "Mishra & Parker 1990": mishra_ks(thr, ths, a),
        "KC eff. porosity (Ahuja)": ahuja_ks(thr, ths, a, n, ths),
        "KC eff. porosity, phi from bd": ahuja_ks(thr, ths, a, n, phi_bd),
    }


def matched(raw, obs, src):
    """Divide by the median over the other sources of each source's median
    log ratio to the measurement (leave the scored source out)."""
    e = np.log10(raw) - np.log10(obs)
    ok = np.isfinite(e)
    per = pd.Series(e[ok]).groupby(src[ok]).median()
    out = np.full(len(raw), np.nan)
    for s in np.unique(src):
        m = src == s
        out[m] = raw[m] / 10.0 ** np.median(per.drop(s, errors="ignore"))
    return out


def ks_table(d, preds, base_key, title):
    obs, src = d.ksat_cmh.to_numpy(float), d.source.to_numpy()
    print(f"\n{title}")
    print(f"{'form':<44s} {'n':>5s} {'med|e|':>6s} {'RMSLE':>6s} {'<2x':>5s} "
          f"{'<10x':>5s} {'bias':>6s} {'rho':>5s} {'src med':>7s} "
          f"{'vs base':>7s} {'p':>8s} {'srcs better':>11s}")
    print("-" * 124)
    eb = np.abs(np.log10(preds[base_key]) - np.log10(obs))
    for k, p in preds.items():
        ok = np.isfinite(p) & (p > 0) & np.isfinite(eb)
        lr = np.log10(p[ok]) - np.log10(obs[ok])
        e = np.abs(lr)
        per = pd.Series(e).groupby(src[ok]).median()
        per_b = pd.Series(eb[ok]).groupby(src[ok]).median()
        better = int((per < per_b).sum())
        pv = wilcoxon(e, eb[ok]).pvalue if k != base_key else np.nan
        print(f"{k:<44s} {ok.sum():5d} {np.median(e):6.2f} {np.sqrt(np.mean(lr**2)):6.2f} "
              f"{np.mean(e <= LOG2)*100:4.0f}% {np.mean(e <= 1)*100:4.0f}% "
              f"{np.mean(lr):+6.2f} {spearmanr(p[ok], obs[ok]).statistic:5.2f} "
              f"{per.median():7.2f} {np.median(e) - np.median(eb[ok]):+7.2f} "
              f"{pv:8.1e} {better:>5d}/{len(per):<5d}")


def ks_part():
    df = st.load_reference_df()
    df["source"] = df.source_db.fillna("KSSL")
    d = df[(df.ksat_cmh > 0) & df.alpha_kpa.notna()].reset_index(drop=True)

    # 1. On the tool's own source-blocked Ks targets, against its kNN.
    kn = pd.read_csv("figures/cache/ks_neighbors.csv")
    kn = kn[kn.design == "source"][["layer_id", "k150"]]
    kn["layer_id"] = kn.layer_id.astype(str)
    t = d.assign(layer_id=d.layer_id.astype(str)).merge(kn, on="layer_id")
    t = t.drop_duplicates("layer_id").reset_index(drop=True)
    obs, src = t.ksat_cmh.to_numpy(float), t.source.to_numpy()
    raw = forms(t)
    preds = {"kNN (the tool, source hidden)": t.k150.to_numpy(float)}
    for k, p in raw.items():
        preds[k + ", raw"] = p
        if not k.startswith("shipped"):
            preds[k + ", matched"] = matched(p, obs, src)
    sh = preds["shipped: Marshall, 50 cm cap, raw"]
    for k in ("CPA, log-h inflection, matched", "KC eff. porosity (Ahuja), raw",
              "KC eff. porosity, phi from bd, raw"):
        preds[f"half shipped, half {k}"] = np.sqrt(sh * preds[k])
    ks_table(t, preds, "shipped: Marshall, 50 cm cap, raw",
             f"KS ON THE TOOL'S {len(t)} SOURCE-BLOCKED TARGETS "
             f"({t.source.nunique()} sources); base = shipped physics")

    # 2. Every layer with a measured Ks (physics needs nothing held out).
    obs, src = d.ksat_cmh.to_numpy(float), d.source.to_numpy()
    raw = forms(d)
    preds = {}
    for k, p in raw.items():
        preds[k + ", raw"] = p
        if not k.startswith("shipped"):
            preds[k + ", matched"] = matched(p, obs, src)
    ks_table(d, preds, "shipped: Marshall, 50 cm cap, raw",
             f"KS ON ALL {len(d)} LAYERS WITH A MEASURED KS "
             f"({d.source.nunique()} sources); base = shipped physics")

    # By texture group, the two best candidates against the shipped physics.
    g = d.texture_class.map(GROUP).to_numpy()
    print(f"\nBY TEXTURE GROUP, all layers (median |log10 error|)")
    keys = ["shipped: Marshall, 50 cm cap, raw", "CPA, log-h inflection, raw",
            "CPA, log-h inflection, matched", "KC eff. porosity (Ahuja), raw"]
    print(f"{'group':<8s} {'n':>5s} " + " ".join(f"{k[:22]:>23s}" for k in keys))
    for gp in ("Sandy", "Loamy", "Silty", "Clayey"):
        m = g == gp
        print(f"{gp:<8s} {m.sum():5d} " + " ".join(
            f"{np.nanmedian(np.abs(np.log10(preds[k][m] / obs[m]))):23.2f}" for k in keys))


# ---------------------------------------------------------------------------
# Texture: clay from the dry end
# ---------------------------------------------------------------------------

def texture_part():
    from verify_network import load
    g = load()
    pool = g[g.depth_cm.notna() & g.bd.notna()]
    tg = pd.concat([x.sample(min(len(x), 150), random_state=0)
                    for _, x in pool.groupby("texture_class")]).reset_index(drop=True)
    est = pd.read_csv("figures/cache/predictive.csv")
    est = est[est.design == "source"].set_index("j").sort_index()
    if len(est) != len(tg):
        raise SystemExit(f"cache holds {len(est)} targets, rebuild gives {len(tg)}")
    t = tg.loc[est.index].reset_index(drop=True)
    tool_sand, tool_clay = est.sand.to_numpy(), est.clay.to_numpy()
    tool_clay = np.where(np.isfinite(tool_clay), tool_clay, np.nan)
    true_sand, true_clay = t.sand.to_numpy(float), t.clay.to_numpy(float)

    def w1500(x):
        th = st.vg_theta(1500.0, x.thetar, x.thetas, x.alpha_kpa, x.n)
        return (th / x.bd).to_numpy(float) * 100.0             # g/g, %

    w = w1500(t)
    oc = t.oc.to_numpy(float)
    rule = np.clip(2.5 * w, 0, 100)
    rule_oc = np.where(np.isfinite(oc), np.clip(2.5 * (w - oc), 0, 100), np.nan)

    # Straight line of clay on w1500 (and on theta1500), fitted on the
    # other studies' layers that carry bulk density.
    pool_w = w1500(pool)
    pool_th = st.vg_theta(1500.0, pool.thetar, pool.thetas, pool.alpha_kpa,
                          pool.n).to_numpy(float)
    th = st.vg_theta(1500.0, t.thetar, t.thetas, t.alpha_kpa, t.n).to_numpy(float)
    fit_w, fit_th = np.full(len(t), np.nan), np.full(len(t), np.nan)
    for s in t.study.unique():
        m = (t.study == s).to_numpy()
        tr = ((pool.study != s).to_numpy() & np.isfinite(pool_w)
              & np.isfinite(pool_th) & pool.clay.notna().to_numpy())
        b = np.polyfit(pool_w[tr], pool.clay.to_numpy(float)[tr], 1)
        fit_w[m] = np.clip(np.polyval(b, w[m]), 0, 100)
        b = np.polyfit(pool_th[tr], pool.clay.to_numpy(float)[tr], 1)
        fit_th[m] = np.clip(np.polyval(b, th[m]), 0, 100)

    print(f"\nCLAY FROM THE DRY END, {len(t)} targets with bulk density "
          f"({t.study.nunique()} studies), study held out from the tool")
    print(f"{'estimate':<46s} {'n':>5s} {'MAE':>6s} {'bias':>6s} {'rho':>5s} "
          f"{'vs tool':>8s} {'p':>8s} {'studies better':>14s}")
    print("-" * 104)
    e_tool = np.abs(tool_clay - true_clay)
    for lab, p in (("tool (kNN fractions, depth + bd)", tool_clay),
                   ("NRCS rule 2.5 x w1500", rule),
                   ("NRCS rule 2.5 x (w1500 - OC), OC known", rule_oc),
                   ("line on w1500, other studies", fit_w),
                   ("line on theta1500, other studies", fit_th),
                   ("half tool, half line on w1500", 0.5 * (tool_clay + fit_w))):
        ok = np.isfinite(p) & np.isfinite(tool_clay)
        e = np.abs(p[ok] - true_clay[ok])
        per = pd.Series(e).groupby(t.study.to_numpy()[ok]).mean()
        per_t = pd.Series(e_tool[ok]).groupby(t.study.to_numpy()[ok]).mean()
        pv = wilcoxon(e, e_tool[ok]).pvalue if lab[:4] != "tool" else np.nan
        print(f"{lab:<46s} {ok.sum():5d} {e.mean():6.2f} "
              f"{np.mean(p[ok] - true_clay[ok]):+6.2f} "
              f"{spearmanr(p[ok], true_clay[ok]).statistic:5.2f} "
              f"{e.mean() - e_tool[ok].mean():+8.2f} {pv:8.1e} "
              f"{int((per < per_t).sum()):>6d}/{len(per):<6d}")

    # Does swapping in the dry-end clay improve the class? Tool sand kept.
    def cls(sand, clay):
        sand = np.clip(sand, 0, 100)
        clay = np.clip(clay, 0, 100 - sand)
        return np.array([st.usda_class(a, 100 - a - c, c) for a, c in zip(sand, clay)])
    truth = t.texture_class.to_numpy()
    base = cls(tool_sand, tool_clay)
    print(f"\nCLASS FROM FRACTIONS (tool's sand kept, clay replaced)")
    print(f"{'fractions':<46s} {'exact':>6s} {'group':>6s} {'vs tool':>8s}")
    for lab, c in (("tool sand + tool clay", tool_clay),
                   ("tool sand + line on w1500", fit_w),
                   ("tool sand + half/half clay", 0.5 * (tool_clay + fit_w))):
        p = cls(tool_sand, c)
        ex = np.mean(p == truth)
        gr = np.mean([GROUP[a] == GROUP[b] for a, b in zip(p, truth)])
        print(f"{lab:<46s} {ex*100:5.1f}% {gr*100:5.1f}% "
              f"{(ex - np.mean(base == truth))*100:+7.1f}")
    print(f"(the tool's own class, from the classifier: "
          f"{np.mean(est.cls.to_numpy() == truth)*100:.1f}% exact)")

    print(f"\nBY MEASURED TEXTURE GROUP, clay MAE (points)")
    gg = np.array([GROUP[c] for c in truth])
    for gp in ("Sandy", "Loamy", "Silty", "Clayey"):
        m = gg == gp
        print(f"  {gp:<8s} n={m.sum():4d}  tool {np.nanmean(e_tool[m]):5.2f}  "
              f"line on w1500 {np.nanmean(np.abs(fit_w[m] - true_clay[m])):5.2f}  "
              f"NRCS rule {np.nanmean(np.abs(rule[m] - true_clay[m])):5.2f}")


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ks_part()
    texture_part()
