"""E2 plotter including the two extra WCF-involving configs
(none/GU and GU/none) alongside the original seven. Mirror of
plot_e2.py — same data source, same axes, same styling convention —
but adds two more series.

Run:
    python \\
        -m simulation.plot_e2_with_wcf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Repo root goes on sys.path so `simulation.*` imports resolve.
ROOT = Path(__file__).resolve().parents[1]
# Simulation outputs and figures live inside this package.
DATA = Path(__file__).resolve().parents[0]
CSV = DATA / "results" / "e2_noncanon" / "grid_summary.csv"
FIG_DIR = DATA / "figures" / "e2_noncanon"

# (config label, style). Same 7 as plot_e2.py, plus 2 WCF-involving rows
# with distinct colors so they read as a separate family.
SERIES = [
    # matched
    ("GU/GU",        dict(color="#1a3a6c", ls="-",  marker="o")),
    ("GA/GA",        dict(color="#2c7fb8", ls="-",  marker="s")),
    ("GU,GA/GU,GA",  dict(color="#41ab5d", ls="-",  marker="^")),
    # partial overlap (existing)
    ("GU/GU,GA",     dict(color="#d95f0e", ls="--", marker="v")),
    ("GU,GA/GU",     dict(color="#cc4c02", ls="--", marker="D")),
    # disjoint mismatch (existing)
    ("GA/GU",        dict(color="#999999", ls=":",  marker="x")),
    ("GU/GA",        dict(color="#bdbdbd", ls=":",  marker="+")),
    # WCF-on-one-side configs
    ("none/GU",      dict(color="#6a3d9a", ls="-",  marker="p")),
    ("GU/none",      dict(color="#67000d", ls="-",  marker="h")),
    # NEW: matched WCF-only baseline
    ("none/none",    dict(color="#000000", ls="-",  marker="*")),
]


def _load():
    d = pd.read_csv(CSV)
    d["cfg"] = d.data_noncanon + "/" + d.test_noncanon
    # Marginal power = n_alt_flagged / N_alt_total (not / n_alt_tested).
    # See §2.2 audit 2026-05-25: cross-N cells drop planted anchors whose
    # stems are not detectable under V(N_test); marginal counts those as
    # missed detections, conditional silently excludes them.
    N_ALT = 2500
    d["power_marginal"] = d["n_alt_flagged"].astype(int) / N_ALT
    return d


def plot(d, metric, ylabel, ymax, out_name):
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for cfg, st in SERIES:
        sub = d[d.cfg == cfg].sort_values("v")
        if sub.empty:
            continue
        ax.plot(sub.v, sub[metric], label=cfg, linewidth=1.8,
                markersize=6, **st)
    if metric == "empirical_fdr":
        ax.axhline(0.05, color="red", ls=":", lw=1.0, zorder=0)
    else:
        ax.axhline(0.8, color="0.6", ls=":", lw=0.8, zorder=0)
    ax.set_xlabel("Mutations per target (v)")
    ax.set_ylabel(ylabel)
    ax.set_xticks([3, 4, 5, 6, 7])
    ax.set_ylim(-0.005 if metric == "empirical_fdr" else 0, ymax)
    ax.grid(True, alpha=0.3)
    ax.legend(title="data / test", fontsize=8, frameon=False, ncol=3,
              loc="lower right" if metric == "power" else "upper right")
    fig.tight_layout()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / out_name
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    return out


def main():
    d = _load()
    p = plot(d, "power_marginal", "Power (BH-adj p < 0.05)", 1.05,
             "fig_e2_power_with_wcf.pdf")
    f = plot(d, "empirical_fdr", "Empirical FDR", 0.08,
             "fig_e2_fdr_with_wcf.pdf")
    print(f"Wrote: {p.with_suffix('.png')}")
    print(f"Wrote: {f.with_suffix('.png')}")


if __name__ == "__main__":
    main()
