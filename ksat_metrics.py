"""Scoring for the Ksat side of the tool.

Texture is a classification problem; Ksat is a regression problem spanning
five orders of magnitude (GSHP: 5th pct 0.004, median 0.63, 95th pct 6.0
cm/h), so it is scored in log10 space. A "factor of N" hit means the estimate
is within N-fold of the measurement, which is the accuracy language the
hydrology literature uses for Ks.

Coverage is the fraction of soils whose measured Ksat falls inside the
reported p5-p95 band. The band is nominally 90 %, so coverage well below 0.90
means the tool is overconfident and well above it means the band is wider
than it needs to be.
"""

import numpy as np
from scipy.stats import spearmanr

FACTORS = (2.0, 5.0, 10.0)


def ksat_scores(obs, med, lo, hi):
    """Compare measured and estimated Ksat (all cm/h).

    Non-finite and non-positive values are dropped: Ksat is logged, and GSHP
    stores some layers as exactly 0, which carries no usable information.
    """
    obs = np.asarray(obs, float)
    med = np.asarray(med, float)
    lo = np.asarray(lo, float)
    hi = np.asarray(hi, float)

    ok = (np.isfinite(obs) & np.isfinite(med) & (obs > 0) & (med > 0))
    n_drop = int((~ok).sum())
    o, m = obs[ok], med[ok]
    e = np.log10(m) - np.log10(o)

    out = {
        "n": int(ok.sum()),
        "n_dropped": n_drop,
        "bias_log10": float(np.mean(e)) if len(e) else float("nan"),
        "rmse_log10": float(np.sqrt(np.mean(e ** 2))) if len(e) else float("nan"),
        "mae_log10": float(np.mean(np.abs(e))) if len(e) else float("nan"),
    }
    for f in FACTORS:
        out[f"within_{f:g}x"] = (float(np.mean(np.abs(e) <= np.log10(f)))
                                 if len(e) else float("nan"))

    cov_ok = ok & np.isfinite(lo) & np.isfinite(hi)
    out["coverage_p5_p95"] = (
        float(np.mean((obs[cov_ok] >= lo[cov_ok]) & (obs[cov_ok] <= hi[cov_ok])))
        if cov_ok.any() else float("nan"))
    out["n_coverage"] = int(cov_ok.sum())

    # Rank agreement: does the tool at least order soils correctly, even where
    # the absolute level is off?
    out["spearman"] = (float(spearmanr(o, m).statistic) if len(e) > 2
                       else float("nan"))
    return out


def header():
    return (f"{'n':>5s} {'bias':>6s} {'RMSE':>6s} "
            + " ".join(f"{'<'+f'{f:g}'+'x':>6s}" for f in FACTORS)
            + f" {'cover':>6s} {'rho':>6s}")


def row(s):
    return (f"{s['n']:5d} {s['bias_log10']:+6.2f} {s['rmse_log10']:6.2f} "
            + " ".join(f"{s[f'within_{f:g}x']*100:5.0f}%" for f in FACTORS)
            + f" {s['coverage_p5_p95']*100:5.0f}% {s['spearman']:6.2f}")


def report(title, named_scores):
    """named_scores: list of (label, scores dict)."""
    print(f"\n{title}")
    print(f"{'':<18s} " + header())
    print("-" * (19 + len(header())))
    for label, s in named_scores:
        print(f"{label:<18s} " + row(s))
    print("  bias/RMSE in log10 cm/h; <Nx = within a factor of N; "
          "cover = measured inside p5-p95 (nominal 90%)")
