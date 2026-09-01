"""Combined §2.2 figure (7 panels A-G) for the manuscript.

Row I (analytic PoC) — closed-form per-target p-value behavior:
    A. SVE branch invariant across (R, N)
    B. SVP branch depends on (R, N) at a representative observation
    C. H0 ECDF of indicator-switch per-target p (v=8): SVP tail stays
       on/under the diagonal across all tested (R, N).

Row II (Monte-Carlo, matched conditions) — operating envelope at
N=GU, R=0.5, uniform abundance, n_alt=n_null=1000, 25 replicates/cell:
    D. FDR (mean realized FDP) vs v across K ∈ {4, 5, 6, 8, 10}, with
       exact-binomial bands and the pi_0 * alpha = 0.025 bound
    E. Power vs v across K ∈ {4, 5, 6, 8, 10}, same bands

Row III (Monte-Carlo, misspecified test) — (data-R × test-R) heatmaps
at fixed N=GU, K=5, v=5, n_alt=n_null=2500, 10 replicates/cell:
    F. Power across the 7×7 (data_R, test_R) grid
    G. FDR across the same grid

Inputs:
    poc panels: closed-form p-value functions (no on-disk inputs)
    grid panels: results/resim/titv-data0.5-test0.5__nc-dataGU-testGU__n1000/
                 grid_summary.csv  (25 replicates per cell)
    heatmap panels: results/resim/e3-RxR-square__K5-v5__nc-dataGU-testGU__n2500/
                    grid_summary.csv  (10 replicates per cell)

Output:
    figures/component2/fig_sim_sec22_combined.{pdf,png}

Run:
    python \\
        -m simulation.plot_fig4_calibration_power
"""
from __future__ import annotations

from math import comb
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from struct_rna.src.non_wcf import build_valid_set
from struct_rna.src.get_pval import target_p_ext, target_p_marginal_ext

# Repo root goes on sys.path so `simulation.*` imports resolve.
ROOT = Path(__file__).resolve().parents[1]
# Simulation outputs and figures live inside this package.
DATA = Path(__file__).resolve().parents[0]
FIG_DIR = DATA / "figures" / "component2"
# Row II is now a replicated grid: R independent BH runs per (K, v) cell,
# so panel D reports an FDR estimate (mean realized FDP) rather than the
# single realized FDP the earlier single-run grid could support.
GRID_CSV = (DATA / "results" / "resim"
            / "titv-data0.5-test0.5__nc-dataGU-testGU__n1000"
            / "grid_summary.csv")
# Row III is likewise replicated: 10 independent BH runs per
# (data_R, test_R) cell at K=5, v=5, n_alt=n_null=2500.
E3_CSV = (DATA / "results" / "resim"
          / "e3-RxR-square__K5-v5__nc-dataGU-testGU__n2500"
          / "grid_summary.csv")

# =====================================================================
# Row I — PoC constants (lifted from poc_sve_svp.py)
# =====================================================================
K_POC = 27
B_POC = [("A", "T"), ("G", "C"), ("C", "G"),
         ("A", "T"), ("G", "C"), ("T", "A")]
L_POC = len(B_POC)

# Two N-families, three R values each. Blue family = GU; orange = GU,GA.
COMBOS_GU = [(0.5, "GU"), (1.5, "GU"), (3.0, "GU")]
COMBOS_GU_GA = [(0.5, "GU,GA"), (1.5, "GU,GA"), (3.0, "GU,GA")]
COLORS_GU = ["#9ecae1", "#4292c6", "#08306b"]      # R = 0.5, 1.5, 3.0
COLORS_GU_GA = ["#fdd0a2", "#fd8d3c", "#7f2704"]   # R = 0.5, 1.5, 3.0

V_AB = list(range(4, 13))
V_ECDF = 8
SVP_STEMMUT, SVP_E = 4, 2

# =====================================================================
# Row II — K × v sweep constants
# =====================================================================
K_VALUES = [4, 5, 6, 8, 10]
K_COLORS = plt.cm.viridis(np.linspace(0.05, 0.85, len(K_VALUES)))

# =====================================================================
# Row III — Ti/Tv heatmap constants
# =====================================================================
R_SQUARE = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0]
V_SLICE = 5


def _make_shifted_cmap(base_name, vcenter, n=256):
    """Return a colormap where the base colormap's midpoint (its
    0.5 position) is shifted to `vcenter` in the new linear scale.

    Used so that the power heatmap's green band sits at power = 0.7
    while keeping the colorbar evenly spaced (vmin=0, vmax=1 linear)
    — instead of using TwoSlopeNorm, which would give an evenly
    sampled colormap but an unevenly spaced colorbar.
    """
    base = plt.colormaps[base_name]
    new_colors = []
    for i in range(n):
        x = i / (n - 1)
        if x <= vcenter:
            y = (x / vcenter) * 0.5 if vcenter > 0 else 0.0
        else:
            y = 0.5 + ((x - vcenter) / (1.0 - vcenter)) * 0.5
        new_colors.append(base(y))
    return LinearSegmentedColormap.from_list(
        f"{base_name}_shift{int(vcenter*100)}", new_colors, N=n
    )


POWER_CMAP = _make_shifted_cmap("viridis_r", vcenter=0.7)

# =====================================================================
# Typography
# =====================================================================
LEGEND_FS = 8
TITLE_FS = 10
LABEL_FS = 9
ROW_TITLE_FS = 12


# =====================================================================
# Panel-label helper
# =====================================================================
def add_panel_label(ax, label, *, dx=-0.18, dy=1.04, fontsize=14):
    ax.text(dx, dy, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight="bold",
            va="bottom", ha="left")


def _r_label(r):
    """LaTeX-formatted R label (calligraphic R)."""
    return rf"$\mathcal{{R}} = {r}$"


def _n_title(text):
    """LaTeX-formatted N family title (calligraphic N)."""
    return rf"$\mathcal{{N}} = \{{\mathrm{{{text}}}\}}$"


def _grouped_legends(ax, lines_gu, lines_gu_ga, loc="lower right"):
    """Two side-by-side legends grouped by N family, anchored at ``loc``."""
    # First legend (GU)
    if loc == "lower right":
        leg1_anchor = (0.55, 0.02)
        leg2_anchor = (1.00, 0.02)
        leg_loc = "lower right"
    elif loc == "upper right":
        leg1_anchor = (0.55, 0.98)
        leg2_anchor = (1.00, 0.98)
        leg_loc = "upper right"
    else:
        leg1_anchor = (0.02, 0.98)
        leg2_anchor = (0.36, 0.98)
        leg_loc = "upper left"

    leg1 = ax.legend(
        handles=lines_gu,
        title=_n_title("GU"),
        loc=leg_loc, bbox_to_anchor=leg1_anchor,
        fontsize=LEGEND_FS, title_fontsize=LEGEND_FS,
        frameon=False, labelspacing=0.3,
    )
    ax.add_artist(leg1)
    leg2 = ax.legend(
        handles=lines_gu_ga,
        title=_n_title("GU,GA"),
        loc=leg_loc, bbox_to_anchor=leg2_anchor,
        fontsize=LEGEND_FS, title_fontsize=LEGEND_FS,
        frameon=False, labelspacing=0.3,
    )
    return leg1, leg2


# =====================================================================
# Row I — PoC panels
# =====================================================================
def draw_poc_sve(ax):
    """A: SVE p vs v — analytic invariant across (R, N)."""
    R0, N0 = 1.5, "GU,GA"
    V0 = build_valid_set(N0)
    ys = [target_p_ext(K_POC, L_POC, v, 0, 0, B_POC, R0, V0)
          for v in V_AB]
    ref = [comb(K_POC - 2 * L_POC, v) / comb(K_POC, v) for v in V_AB]
    ax.plot(V_AB, ref, "k--", lw=1.5, zorder=0,
            label=r"$\binom{k-2L}{v}/\binom{k}{v}$")
    ax.plot(V_AB, ys, marker="o", ms=5, color="#2c7fb8",
            label=rf"SVE p ($\mathcal{{R}}={R0}$, $\mathcal{{N}}=\{{\mathrm{{{N0}}}\}}$)")
    ax.set_xlabel(r"$v$ (mutations / target)", fontsize=LABEL_FS)
    ax.set_ylabel("SVE target $p$-value", fontsize=LABEL_FS)
    ax.set_title("SVE branch ($s=0$): identity-free", fontsize=TITLE_FS)
    ax.annotate(r"invariant to $\mathcal{R}$, $\mathcal{N}$, $b$",
                xy=(0.50, 0.18), xycoords="axes fraction",
                fontsize=LEGEND_FS, ha="left", style="italic",
                color="0.35")
    ax.legend(fontsize=LEGEND_FS, frameon=False)
    ax.grid(alpha=0.3)
    add_panel_label(ax, "A")


def draw_poc_svp(ax):
    """B: SVP p vs v — depends on (R, N). Grouped legends by N."""
    lines_gu = []
    lines_gu_ga = []
    for (R, N), c in zip(COMBOS_GU, COLORS_GU):
        V = build_valid_set(N)
        ys = [target_p_ext(K_POC, L_POC, v, SVP_STEMMUT, SVP_E,
                           B_POC, R, V) for v in V_AB]
        line, = ax.plot(V_AB, ys, marker="s", ms=4, color=c, alpha=0.9,
                         label=_r_label(R))
        lines_gu.append(line)
    for (R, N), c in zip(COMBOS_GU_GA, COLORS_GU_GA):
        V = build_valid_set(N)
        ys = [target_p_ext(K_POC, L_POC, v, SVP_STEMMUT, SVP_E,
                           B_POC, R, V) for v in V_AB]
        line, = ax.plot(V_AB, ys, marker="s", ms=4, color=c, alpha=0.9,
                         label=_r_label(R))
        lines_gu_ga.append(line)
    ax.set_xlabel(r"$v$ (mutations / target)", fontsize=LABEL_FS)
    ax.set_ylabel("SVP target $p$-value", fontsize=LABEL_FS)
    ax.set_ylim(top=0.8)
    ax.set_title(rf"SVP branch ($s={SVP_STEMMUT}$, $e={SVP_E}$)",
                 fontsize=TITLE_FS)
    _grouped_legends(ax, lines_gu, lines_gu_ga, loc="upper left")
    ax.grid(alpha=0.3)
    add_panel_label(ax, "B")


def draw_poc_ecdf(ax):
    """C: H0 ECDF at v=V_ECDF. q0 atom indicator removed per request."""
    # Diagonal reference (uniform)
    ax.plot([0, 1], [0, 1], color="0.55", ls=":", lw=1.0)
    ax.text(0.62, 0.55, "uniform", fontsize=LEGEND_FS, color="0.45",
            rotation=38, ha="center", va="center", style="italic")

    lines_gu = []
    lines_gu_ga = []
    for (R, N), c in zip(COMBOS_GU, COLORS_GU):
        V = build_valid_set(N)
        support, pmf = target_p_marginal_ext(K_POC, L_POC, V_ECDF,
                                              B_POC, R, V)
        xs, ys, cum = [], [], 0.0
        for p, m in zip(support, pmf):
            cum += m
            xs.append(p); ys.append(cum)
        line, = ax.step(xs, ys, where="post", color=c, alpha=0.9,
                         label=_r_label(R))
        lines_gu.append(line)
    for (R, N), c in zip(COMBOS_GU_GA, COLORS_GU_GA):
        V = build_valid_set(N)
        support, pmf = target_p_marginal_ext(K_POC, L_POC, V_ECDF,
                                              B_POC, R, V)
        xs, ys, cum = [], [], 0.0
        for p, m in zip(support, pmf):
            cum += m
            xs.append(p); ys.append(cum)
        line, = ax.step(xs, ys, where="post", color=c, alpha=0.9,
                         label=_r_label(R))
        lines_gu_ga.append(line)
    ax.set_xlabel("target $p$-value", fontsize=LABEL_FS)
    ax.set_ylabel(r"$H_0$ CDF", fontsize=LABEL_FS)
    ax.set_title(rf"SVP-branch $H_0$ CDF at $v={V_ECDF}$", fontsize=TITLE_FS)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    _grouped_legends(ax, lines_gu, lines_gu_ga, loc="upper left")
    ax.grid(alpha=0.3)
    add_panel_label(ax, "C")


# =====================================================================
# Row II — FDR / power across (K, v), averaged over replicates
# =====================================================================
def _clopper_pearson(k, n, lvl=0.95):
    """Exact binomial interval. Used instead of mean +- 1.96 SE because at
    K = 8, 10 the false-discovery count per replicate is near zero and the
    normal approximation returns a negative lower bound."""
    from scipy.stats import beta
    a = (1.0 - lvl) / 2.0
    lo = 0.0 if k == 0 else beta.ppf(a, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a, k + 1, n - k)
    return lo, hi


def aggregate_replicates(df):
    """Per (K, v, abundance): mean realized FDP and mean power over the
    replicates, each with an exact-binomial band on the pooled counts.

    The point estimate is the mean of the per-replicate ratios (the FDR
    estimator); the band is computed from the pooled event counts, which
    for these cells agrees with the mean to <5e-4 while staying inside
    [0, 1] at the boundary cells.
    """
    rows = []
    for (K, v, ab), g in df.groupby(["K", "v", "abundance"]):
        V = int(g.n_null_flagged.sum())               # pooled false discoveries
        A = int(g.n_alt_flagged.sum())                # pooled true discoveries
        n0 = int(g.n_null_tested.sum())
        n1 = int(g.n_alt_tested.sum())
        R = V + A
        f_lo, f_hi = _clopper_pearson(V, n0)
        p_lo, p_hi = _clopper_pearson(A, n1)
        rows.append({
            "K": K, "v": v, "abundance": ab, "reps": len(g),
            "empirical_fdr": g.empirical_fdp.mean(),
            "fdr_lo": f_lo * n0 / R if R else 0.0,
            "fdr_hi": f_hi * n0 / R if R else 0.0,
            "power": g.power.mean(),
            "power_lo": p_lo, "power_hi": p_hi,
        })
    return pd.DataFrame(rows)


def aggregate_rxr(df):
    """Per (data_R, test_R): mean realized FDP and mean power over the
    replicates.  Column named `empirical_fdr` because averaging over
    independent BH runs is what makes the FDR label legitimate."""
    g = (df.groupby(["data_R", "test_R"])
           .agg(empirical_fdr=("empirical_fdp", "mean"),
                power=("power", "mean"),
                v=("v", "first"))
           .reset_index())
    return g


def draw_grid_metric(ax, df, metric, ylabel, ymax, ref_y,
                      ref_color, ref_style, legend_loc, panel_label):
    lo_col = "fdr_lo" if metric == "empirical_fdr" else "power_lo"
    hi_col = "fdr_hi" if metric == "empirical_fdr" else "power_hi"
    for K, color in zip(K_VALUES, K_COLORS):
        sub = (df[(df.abundance == "uniform") & (df.K == K)]
               .sort_values("v"))
        ax.fill_between(sub.v, sub[lo_col], sub[hi_col],
                        color=color, alpha=0.18, linewidth=0)
        ax.plot(sub.v, sub[metric], marker="o", linewidth=1.8,
                markersize=5, color=color, label=f"$K = {K}$")
    if ref_y is not None:
        ax.axhline(ref_y, color=ref_color, linestyle=ref_style,
                   linewidth=1.0, zorder=0)
        # The line is pi_0 * alpha = 0.025 (pi_0 = 0.5 by design), which is
        # the FDR the BH procedure controls here. Stated in the caption
        # rather than annotated on the axes.
    ax.set_xlabel(r"$v$ (mutations / target)", fontsize=LABEL_FS)
    ax.set_xticks([2, 3, 4, 5, 6, 7, 8])
    ax.set_xlim(1.7, 8.3)
    ax.set_ylim(-0.005 if metric == "empirical_fdr" else 0, ymax)
    ax.set_ylabel(ylabel, fontsize=LABEL_FS)
    ax.legend(loc=legend_loc, frameon=False, fontsize=LEGEND_FS)
    ax.grid(True, alpha=0.3)
    add_panel_label(ax, panel_label)


# =====================================================================
# Row III — Ti/Tv (data-R × test-R) heatmaps at v=V_SLICE
# =====================================================================
def draw_heatmap(ax, d, metric, vmin, vmax, cmap, panel_label,
                 cbar_label):
    """Plain linear colormap; the colorbar is evenly spaced. To shift
    the colormap's midpoint (e.g. green at 0.7 instead of 0.5), pass
    in a pre-shifted colormap built by `_make_shifted_cmap`."""
    sub = d[(d.v == V_SLICE) & d.data_R.isin(R_SQUARE)
            & d.test_R.isin(R_SQUARE)]
    piv = sub.pivot_table(index="data_R", columns="test_R",
                           values=metric)
    piv = piv.reindex(index=R_SQUARE, columns=R_SQUARE)
    im = ax.imshow(piv.values, origin="lower", cmap=cmap,
                   vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(R_SQUARE)))
    ax.set_xticklabels(R_SQUARE)
    ax.set_yticks(range(len(R_SQUARE)))
    ax.set_yticklabels(R_SQUARE)
    ax.set_xlabel(r"test $\mathcal{R}$ (assumed)", fontsize=LABEL_FS)
    ax.set_ylabel(r"data $\mathcal{R}$ (true)", fontsize=LABEL_FS)
    for iy in range(len(R_SQUARE)):
        for ix in range(len(R_SQUARE)):
            val = piv.values[iy, ix]
            if not np.isnan(val):
                # Heatmap cell text — fontsize bumped for readability.
                # Power uses viridis_r with TwoSlopeNorm centered at
                # 0.7 (green), so the dark-teal/blue range starts
                # around val ≈ 0.78 and that's where white text wins.
                # Empirical FDR uses Reds (white=low, dark=high).
                if metric == "power":
                    text_color = "white" if val > 0.78 else "black"
                else:  # empirical_fdr (Reds): dark at high FDR
                    text_color = "white" if val > 0.12 else "black"
                ax.text(ix, iy, f"{val:.2f}", ha="center",
                        va="center", fontsize=9,
                        color=text_color)
    ax.plot(range(len(R_SQUARE)), range(len(R_SQUARE)),
            color="red", lw=0.8, ls=":")
    cb = plt.colorbar(im, ax=ax, shrink=0.85, pad=0.04)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(cbar_label, fontsize=LABEL_FS)
    add_panel_label(ax, panel_label)


# =====================================================================
# Row titles (placed above each row after layout)
# =====================================================================
ROW_TITLES = [
    "I. Closed-form per-target $p$-value (analytic)",
    "II. Empirical operating envelope, matched conditions",
    r"III. Sensitivity to misspecification of the Ti/Tv ratio $\mathcal{R}$",
]


def add_row_titles(fig, row_first_axes, row_last_axes):
    """Place a bold centered title above each row, lifted clear of the
    external panel labels in the left margin."""
    for first, last, title in zip(row_first_axes, row_last_axes, ROW_TITLES):
        pos_first = first.get_position()
        pos_last = last.get_position()
        x_center = (pos_first.x0 + pos_last.x1) / 2
        y = min(pos_first.y1 + 0.045, 0.99)
        fig.text(x_center, y, title,
                 fontsize=ROW_TITLE_FS, fontweight="bold",
                 ha="center", va="bottom", color="#222222")


# =====================================================================
# Assembly
# =====================================================================
def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df_grid = aggregate_replicates(pd.read_csv(GRID_CSV))
    df_e3 = aggregate_rxr(pd.read_csv(E3_CSV))

    fig = plt.figure(figsize=(11, 10.5))
    gs = gridspec.GridSpec(
        nrows=3, ncols=6,
        figure=fig,
        height_ratios=[1.0, 0.95, 1.10],
        hspace=0.55,
        wspace=0.85,
        left=0.07, right=0.97,
        top=0.93, bottom=0.06,
    )

    # Row I (PoC) — 3 panels
    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_c = fig.add_subplot(gs[0, 4:6])

    # Row II (empirical sweep) — 2 panels
    ax_d = fig.add_subplot(gs[1, 0:3])
    ax_e = fig.add_subplot(gs[1, 3:6])

    # Row III (heatmaps) — 2 panels
    ax_f = fig.add_subplot(gs[2, 0:3])
    ax_g = fig.add_subplot(gs[2, 3:6])

    draw_poc_sve(ax_a)
    draw_poc_svp(ax_b)
    draw_poc_ecdf(ax_c)
    draw_grid_metric(
        ax_d, df_grid,
        metric="empirical_fdr",
        # pi_0 = n_null / (n_null + n_alt) = 0.5 by design, and BH controls
        # the FDR at pi_0 * alpha, so 0.025 -- not alpha = 0.05 -- is the
        # bound these curves can approach.
        ylabel="Empirical FDR",
        ymax=0.06,
        ref_y=0.025,
        ref_color="red", ref_style=":",
        legend_loc="upper right",
        panel_label="D",
    )
    draw_grid_metric(
        ax_e, df_grid,
        metric="power",
        ylabel=r"Power (BH-adjusted $p < 0.05$)",
        ymax=1.05,
        ref_y=None,
        ref_color="0.6", ref_style=":",
        legend_loc="lower right",
        panel_label="E",
    )
    draw_heatmap(
        ax_f, df_e3,
        metric="power",
        vmin=0.0, vmax=1.0,
        cmap=POWER_CMAP,            # shifted so green sits at power = 0.7
        panel_label="F",
        cbar_label="Mean power",
    )
    draw_heatmap(
        ax_g, df_e3,
        metric="empirical_fdr",
        vmin=0.0, vmax=0.20,
        cmap="Reds",
        panel_label="G",
        cbar_label="Empirical FDR",
    )

    # Row titles after layout is fixed (need first and last axes per row
    # so the title can be centered above each row span)
    add_row_titles(fig, [ax_a, ax_d, ax_f], [ax_c, ax_e, ax_g])

    out = FIG_DIR / "fig_sim_sec22_combined.pdf"
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Wrote: {out}")
    print(f"Wrote: {out.with_suffix('.png')}")


if __name__ == "__main__":
    main()
