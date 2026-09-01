"""Build a synthetic SPLASH input for one parameter cell, run STRUCT target mode,
read the output, and return power / FDR metrics.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from pandarallel import pandarallel

from struct_rna.structure_target_mode import SS_target
from struct_rna.src.non_wcf import build_valid_set, V_EXT

from .core import (
    abundance_counts,
    build_anchor_row,
    make_base_target,
    make_null_target,
    make_sve_target,
    make_svp_target,
    random_anchor_seq,
)

SUB_MODES = ("sve", "svp", "mixed")


def _gen_alt_targets(rng, base, geom, K, v, sub_mode, ti_tv_ratio,
                     data_valid=V_EXT):
    """Return K-1 less-abundant targets following the per-anchor sub-mode
    rule. ``data_valid`` is the data-side V(N) used by the SVP path
    (which pairs count as structure-preserving when planting)."""
    out = []
    for _ in range(K - 1):
        if sub_mode == "sve":
            out.append(make_sve_target(rng, base, geom, v, ti_tv_ratio=ti_tv_ratio))
        elif sub_mode == "svp":
            out.append(make_svp_target(rng, base, geom, v, ti_tv_ratio=ti_tv_ratio,
                                       valid=data_valid))
        else:  # mixed
            if rng.random() < 0.5:
                out.append(make_sve_target(rng, base, geom, v, ti_tv_ratio=ti_tv_ratio))
            else:
                out.append(make_svp_target(rng, base, geom, v, ti_tv_ratio=ti_tv_ratio,
                                           valid=data_valid))
    return out


def build_cell_dataframe(
    K: int,
    v: int,
    abundance_kind: str,
    n_alt: int,
    n_null: int,
    rng: np.random.Generator,
    ti_tv_ratio: float = 0.5,
    data_valid=V_EXT,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (splash_input_df, meta_df).

    splash_input_df has the columns process_df expects. meta_df keeps per-anchor
    bookkeeping (is_alt, sub_mode, n_wobble, L) for joining onto STRUCT output.

    ``data_valid`` is the data-side V(N): planted stems draw pairs from it
    and SVP-preservation is judged against it. Default V_EXT (G·U).
    """
    rows = []
    meta = []

    for n in range(n_alt):
        sub_mode = SUB_MODES[int(rng.integers(0, len(SUB_MODES)))]
        # Resample geometry until the sub-mode is feasible at v.
        for _ in range(50):
            geom = make_base_target(rng, valid=data_valid)
            base = geom["sequence"]
            non_stem_count = 27 - 2 * geom["L"]
            if sub_mode == "sve" and v > non_stem_count:
                continue  # SVE can't fit v mutations outside stem
            try:
                less = _gen_alt_targets(rng, base, geom, K, v, sub_mode,
                                        ti_tv_ratio, data_valid)
                break
            except RuntimeError:
                continue
        else:
            raise RuntimeError(f"alt anchor n={n} infeasible at K={K}, v={v}, mode={sub_mode}")

        counts = abundance_counts(K, abundance_kind)
        anchor_seq = random_anchor_seq(rng)
        anchor_id = f"alt_{n:06d}_{anchor_seq}"
        rows.append(build_anchor_row(anchor_id, base, less, counts))
        meta.append({
            "anchor": anchor_id,
            "is_alt": True,
            "sub_mode": sub_mode,
            "L": geom["L"],
            "loop_len": geom["loop_len"],
            "n_wobble": geom["n_wobble"],
        })

    for n in range(n_null):
        geom = make_base_target(rng, valid=data_valid)
        base = geom["sequence"]
        less = [make_null_target(rng, base, v, ti_tv_ratio=ti_tv_ratio) for _ in range(K - 1)]
        counts = abundance_counts(K, abundance_kind)
        anchor_seq = random_anchor_seq(rng)
        anchor_id = f"null_{n:06d}_{anchor_seq}"
        rows.append(build_anchor_row(anchor_id, base, less, counts))
        meta.append({
            "anchor": anchor_id,
            "is_alt": False,
            "sub_mode": "null",
            "L": geom["L"],
            "loop_len": geom["loop_len"],
            "n_wobble": geom["n_wobble"],
        })

    return pd.DataFrame(rows), pd.DataFrame(meta)


def run_cell(
    K: int,
    v: int,
    abundance_kind: str,
    n_alt: int,
    n_null: int,
    out_dir: Path,
    seed: int = 42,
    nb_workers: int = 1,
    ti_tv_ratio: float = 0.5,
    test_ti_tv_ratio: float = 0.5,
    data_noncanon: str = "GU",
    test_noncanon: str = "GU",
) -> dict:
    """Run one parameter cell end-to-end and return summary metrics.

    Writes the synthesized SPLASH input and the STRUCT output under out_dir.

    ``ti_tv_ratio`` controls the simulator's mutation identity weighting
    (the data side). ``test_ti_tv_ratio`` is forwarded to SS_target's
    ``titv`` argument and controls the null assumed by the test. Matched
    settings (test = data) restore nominal FDR; mismatched (test = 0.5,
    data = 2.0) reproduce the anti-conservative drift quantified in
    Component 3.

    ``data_noncanon`` is the non-canonical set planted in the data;
    ``test_noncanon`` is the set the test admits (SS_target ``noncanon``).
    Their 2×2 (data uses N? × test admits N?) is the E2 design. Defaults
    "GU"/"GU" reproduce the prior G·U behaviour.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    pandarallel.initialize(nb_workers=nb_workers, progress_bar=False, verbose=0)

    data_valid = build_valid_set(data_noncanon)
    splash_df, meta_df = build_cell_dataframe(
        K, v, abundance_kind, n_alt, n_null, rng,
        ti_tv_ratio=ti_tv_ratio, data_valid=data_valid,
    )
    input_tsv = out_dir / "splash_input.tsv"
    splash_df.to_csv(input_tsv, sep="\t", index=False)
    meta_df.to_csv(out_dir / "meta.tsv", sep="\t", index=False)

    output_prefix = str(out_dir / "struct")
    # Non-canonical behaviour is carried entirely by `noncanon`; the legacy
    # `wobble` bool was removed upstream (splash-structure 158f157), where
    # wobble=True was equivalent to noncanon="GU".
    SS_target(
        output_prefix=output_prefix,
        splash_output_file=str(input_tsv),
        element_annotation=False,
        titv=test_ti_tv_ratio,
        noncanon=test_noncanon,
    )

    # SS_target writes to <prefix>_results/structure_on_targets.tsv
    output_tsv = Path(f"{output_prefix}_results/structure_on_targets.tsv")
    if not output_tsv.exists():
        return {
            "K": K, "v": v, "abundance": abundance_kind,
            "n_alt": n_alt, "n_null": n_null,
            "anchors_with_stem": 0,
            "power": float("nan"),
            "fdr": float("nan"),
            "note": "no STRUCT output (no anchors had a detectable stem)",
        }
    struct_df = pd.read_csv(output_tsv, sep="\t")
    # struct_df has one row per (anchor, target). Reduce to per-anchor.
    per_anchor = struct_df[["anchor", "anchor_p", "anchor_p_BH"]].drop_duplicates()
    merged = per_anchor.merge(meta_df, on="anchor", how="inner")

    flagged = merged["anchor_p_BH"] < 0.05
    n_alt_tested = int(merged["is_alt"].sum())
    n_null_tested = int((~merged["is_alt"]).sum())
    n_alt_flagged = int((flagged & merged["is_alt"]).sum())
    n_null_flagged = int((flagged & ~merged["is_alt"]).sum())
    total_flagged = int(flagged.sum())
    power = n_alt_flagged / n_alt_tested if n_alt_tested else float("nan")
    fdr = n_null_flagged / total_flagged if total_flagged else 0.0

    alts = merged[merged["is_alt"]].copy()
    alts["flagged"] = alts["anchor_p_BH"] < 0.05

    by_mode = (
        alts.groupby("sub_mode")["flagged"].mean().round(3).to_dict()
    )
    by_L = (
        alts.groupby("L")["flagged"]
        .agg(["mean", "count"])
        .round(3)
        .to_dict("index")
    )
    by_L = {int(k): {"power": float(v["mean"]), "n": int(v["count"])} for k, v in by_L.items()}
    # Cross-tab L × sub_mode (counts may be small at 100 alts)
    by_L_mode = {}
    for (L_val, mode), grp in alts.groupby(["L", "sub_mode"]):
        by_L_mode[f"L{int(L_val)}_{mode}"] = {
            "power": round(float(grp["flagged"].mean()), 3),
            "n": int(len(grp)),
        }

    # --- E2 stratification by realized non-canonical count ---
    # meta `n_wobble` is the count of planted non-WCF pairs (= realized
    # n_nc under any N, since core._plant_once counts pairs ∉ V_WCF).
    # Binary stratum: stem contains ≥1 non-canonical pair vs all-WCF.
    def _stratum(df_):
        if len(df_) == 0:
            return {"power": float("nan"), "n": 0}
        return {"power": round(float(df_["flagged"].mean()), 3),
                "n": int(len(df_))}

    has_nc = alts[alts["n_wobble"] >= 1]
    all_wcf = alts[alts["n_wobble"] == 0]
    power_by_nc_stratum = {
        "has_nc": _stratum(has_nc),   # planted structure uses ≥1 N pair
        "all_wcf": _stratum(all_wcf), # planted stem is all-WCF (thin in large-N arms)
    }
    # Realized n_nc histogram over alt anchors (the planted population,
    # post cap-truncation) — reported so the N-dependent truncation is
    # visible rather than hidden.
    nc_hist = {
        int(k): int(c)
        for k, c in alts["n_wobble"].value_counts().sort_index().items()
    }

    return {
        "K": K, "v": v, "abundance": abundance_kind,
        "ti_tv_ratio": ti_tv_ratio,
        "test_ti_tv_ratio": test_ti_tv_ratio,
        "data_noncanon": data_noncanon,
        "test_noncanon": test_noncanon,
        "n_alt_tested": n_alt_tested, "n_null_tested": n_null_tested,
        "n_alt_flagged": n_alt_flagged, "n_null_flagged": n_null_flagged,
        "power": round(power, 3),
        "empirical_fdr": round(fdr, 3),
        "power_by_submode": by_mode,
        "power_by_L": by_L,
        "power_by_L_submode": by_L_mode,
        "power_by_nc_stratum": power_by_nc_stratum,
        "nc_hist": nc_hist,
    }
