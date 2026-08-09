"""Build a reference table from sDB (Vereecken, Van Looy, Weynants & Javaux
2017), the soil retention and conductivity curve database.

    doi:10.1594/PANGAEA.879233, licence CC-BY-3.0.

Unlike HYPRES and EU-HYDI, sDB is openly licensed and therefore *may* be
redistributed with attribution, so the derived table is committed.

Contents: 38 Belgian profiles / 182 horizons. Each horizon carries a measured
retention curve (typically ~24 (theta, h) pairs), sand/silt/clay, bulk density,
organic carbon and a topsoil flag. Van Genuchten parameters are fitted here
with the same routine the tool uses (Mualem, m = 1 - 1/n).

Units in the .mat file: h in cm of water, theta in cm3/cm3, S_L_A in percent
(sand, silt "leem", clay), conductivity in cm/day. Verified in main() by
checking that the fitted alpha values sit on the same scale as GSHP's.

Depth is not stored explicitly; the topsoil flag is mapped to representative
mid-depths so the horizons can still be used by the --depth variant.
"""

import numpy as np
import pandas as pd
import scipy.io as sio

import swcc_texture as st

MAT = "data/sdb_raw/sDB.mat"
OUT = "data/sdb_reference.csv"
CM_TO_KPA = 0.0980665
TOPSOIL_DEPTH_CM, SUBSOIL_DEPTH_CM = 15.0, 60.0
MIN_POINTS = 5
MAX_RMSE = 0.03


def main():
    m = sio.loadmat(MAT, squeeze_me=False, struct_as_record=False)
    profiles = m["sDB"][0]

    rows, n_skip = [], 0
    for p in profiles:
        pid = str(np.asarray(p.IDProfile).ravel()[0])
        for hi, hor in enumerate(np.asarray(p.Hor).ravel()):
            th_h = np.asarray(hor.th_h, dtype=float)
            sla = np.asarray(hor.S_L_A, dtype=float).ravel()
            if th_h.ndim != 2 or th_h.shape[0] < MIN_POINTS or sla.size < 3:
                n_skip += 1
                continue

            theta, h_cm = th_h[:, 0], th_h[:, 1]
            good = np.isfinite(theta) & np.isfinite(h_cm) & (h_cm >= 0) & (theta > 0)
            if good.sum() < MIN_POINTS:
                n_skip += 1
                continue
            h_kpa = h_cm[good] * CM_TO_KPA
            theta = theta[good]

            try:
                popt, _ = st.fit_vg(h_kpa, theta)
            except Exception:
                n_skip += 1
                continue
            tr, ts, la, ln1 = popt
            alpha, n = 10.0 ** la, 1.0 + 10.0 ** ln1
            rmse = float(np.sqrt(np.mean(
                (st.vg_theta(h_kpa, tr, ts, alpha, n) - theta) ** 2)))
            if not np.isfinite(rmse) or rmse > MAX_RMSE or n <= 1.0:
                n_skip += 1
                continue

            sand, silt, clay = sla[:3]
            tot = sand + silt + clay
            if not (95 <= tot <= 105):
                n_skip += 1
                continue
            sand, silt, clay = (100 * sand / tot, 100 * silt / tot,
                                100 * clay / tot)

            # Saturated conductivity: h_K_th holds (h, K, theta); take K at the
            # smallest pressure head as Ksat.
            ksat_cmh = np.nan
            hk = np.asarray(getattr(hor, "h_K_th", np.empty((0, 3))), dtype=float)
            if hk.ndim == 2 and hk.shape[0] and hk.shape[1] >= 2:
                ok = np.isfinite(hk[:, 0]) & np.isfinite(hk[:, 1])
                if ok.any():
                    j = np.argmin(hk[ok, 0])
                    k_day = hk[ok][j, 1]
                    if k_day > 0:
                        ksat_cmh = k_day / 24.0

            topsoil = bool(np.asarray(hor.topsoil).ravel()[0])
            rows.append(dict(
                layer_id=f"sDB_{pid}_{hi}", profile_id=f"sDB_{pid}",
                texture_class=st.usda_class(sand, silt, clay),
                alpha_kpa=alpha, n=n, thetar=tr, thetas=ts,
                sand=sand, silt=silt, clay=clay, ksat_cmh=ksat_cmh,
                depth_cm=TOPSOIL_DEPTH_CM if topsoil else SUBSOIL_DEPTH_CM,
                rmse=rmse, n_points=int(good.sum())))

    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(f"sDB horizons usable: {len(out)}  (skipped {n_skip})")
    print(f"  with ksat: {out.ksat_cmh.notna().sum()}")
    print(f"  median retention points per curve: {out.n_points.median():.0f}")
    print(f"  median fit RMSE: {out.rmse.median():.4f}")
    print(out.texture_class.value_counts())

    gshp = pd.read_csv(st.REFERENCE_CSV)
    print("\nunit check -- class-median alpha [kPa^-1] / n, sDB vs GSHP:")
    for cls in ["sandy loam", "loam", "silt loam", "clay loam"]:
        a, b = out[out.texture_class == cls], gshp[gshp.texture_class == cls]
        if len(a) and len(b):
            print(f"  {cls:<12s} {a.alpha_kpa.median():7.3f} vs "
                  f"{b.alpha_kpa.median():7.3f}   n {a.n.median():5.2f} vs "
                  f"{b.n.median():5.2f}   ({len(a)} vs {len(b)} soils)")


if __name__ == "__main__":
    main()
