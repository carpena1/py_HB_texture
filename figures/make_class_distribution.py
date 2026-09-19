"""Layers per USDA texture class and group in the default reference.

    python figures/make_class_distribution.py

Each class bar is split into the distributed tables and the restricted ones
(EU-HYDI, Laikipia), as the coverage map splits them by location. Only class
counts are drawn, so the figure can be shared even when a table is
restricted. Writes figures/fig5_texture_classes.png/.svg and prints the counts.
"""

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import swcc_texture as st                                     # noqa: E402
from verify_common import GROUP, GROUP_ORDER                  # noqa: E402
import validation_plots as vp                                 # noqa: E402
from matplotlib.patches import Rectangle                      # noqa: E402

OUT = "figures"


def counts():
    d = st.load_reference_df()
    public = set(st.load_reference_df("public").layer_id.astype(str))
    d["public"] = d.layer_id.astype(str).isin(public)
    t = pd.DataFrame({
        "public": d[d.public].texture_class.value_counts(),
        "restricted": d[~d.public].texture_class.value_counts(),
        "ks": d[d.ksat_cmh.notna()].texture_class.value_counts(),
    }).reindex(vp.GROUPED).fillna(0).astype(int)
    t["total"] = t.public + t.restricted
    t["group"] = [GROUP[c] for c in t.index]
    return t


def main():
    t = counts()
    n = t.total.sum()
    print(t.assign(pct=(t.total / n * 100).round(1)))
    print(t.groupby("group", sort=False)[["public", "restricted", "total", "ks"]].sum()
          .assign(pct=lambda g: (g.total / n * 100).round(1)))

    vp.style()
    fig, ax = vp.plt.subplots(figsize=(8.6, 5.6))
    fig.subplots_adjust(top=0.82, bottom=0.14, left=0.18, right=0.8)
    # Classes top to bottom in group order, with a gap between groups.
    ys, y = [], 0.0
    for i, c in enumerate(t.index):
        if i and t.group.iloc[i] != t.group.iloc[i - 1]:
            y += 0.6
        ys.append(y)
        y += 1.0
    ys = -np.array(ys)
    h = 0.72
    for yy, (c, r) in zip(ys, t.iterrows()):
        ax.barh(yy, r.public, height=h, color=vp.S1, edgecolor=vp.SURFACE, lw=1, zorder=2)
        ax.barh(yy, r.restricted, left=r.public, height=h, color=vp.S2,
                edgecolor=vp.SURFACE, lw=1, zorder=2)
        ax.text(r.total + 60, yy, f"{r.total:,}  ({r.total / n * 100:.1f} %)",
                va="center", fontsize=8.5, color=vp.INK)
    ax.set_yticks(ys, t.index, fontsize=9.5)
    xmax = 5900
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_formatter(lambda v, _: f"{v:,.0f}")
    vp.hgrid(ax, [0, 1000, 2000, 3000, 4000, 5000])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("layers", fontsize=9, color=vp.INK2)
    # Group totals in the right margin, level with each group's middle.
    g = t.assign(y=ys).groupby("group", sort=False)
    for name in GROUP_ORDER:
        rows = g.get_group(name)
        tot = rows.total.sum()
        ym = rows.y.mean()
        ax.text(1.03, ym, f"{name}\n{tot:,} ({tot / n * 100:.0f} %)",
                transform=ax.get_yaxis_transform(), va="center", fontsize=9,
                color=vp.INK, fontweight="bold", linespacing=1.4)
        ax.plot([1.015, 1.015], [rows.y.max() + h / 2, rows.y.min() - h / 2],
                transform=ax.get_yaxis_transform(), color=vp.AXIS, lw=1,
                clip_on=False)
    hd = [Rectangle((0, 0), 1, 1, color=c) for c in (vp.S1, vp.S2)]
    ax.legend(hd, ["distributed tables", "restricted tables (EU-HYDI, Laikipia)"],
              loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=9)
    big, small = t.total.idxmax(), t.total.idxmin()
    vp.title(fig, "Texture classes in the default reference",
             f"{n:,} layers. {big.capitalize()} outnumbers {small} "
             f"{t.total.max() / t.total.min():.0f} to 1. The neighbour vote and the "
             "classifier weight classes evenly,\nbut a weight adds no information: "
             f"{small}'s {t.total.min()} layers stay a weak point.")
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(OUT, f"fig5_texture_classes.{ext}"), dpi=200,
                    bbox_inches="tight", pad_inches=0.15)
    print(f"  wrote {OUT}/fig5_texture_classes.png/.svg")


if __name__ == "__main__":
    main()
