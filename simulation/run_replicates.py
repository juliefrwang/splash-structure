"""Replicated simulation grids: R independent replicates per (K, v) cell.

Each *simulation setting* (a fixed data/test configuration) gets its own
folder under ``results/resim/`` holding everything produced under that
setting -- provenance, the long-form summary, and the per-replicate anchor
p-values:

    results/resim/<setting_id>/
        config.json                       full config + provenance
        grid_summary.csv                  one row per (K, v, abundance, rep)
        pvals/rep000/K5_v5_uniform.tsv    anchor, is_alt, anchor_p

One `run_cell` call = one BH run over n_alt + n_null anchors = one realized
FDP.  Averaging the realized FDP over the R replicates of a cell estimates
that cell's FDR, and the spread across replicates gives its uncertainty.
A single replicate supports only the word "FDP".

Per-replicate artifacts are deliberately slim.  Keeping each cell's full
`splash_input.tsv` + STRUCT output would cost ~3.5 MB x cells x reps (a few
GB at R=25); the anchor p-values are all the downstream analysis needs.

Run from repo root:
    PYTHONHASHSEED=0 python \\
        -m simulation.run_replicates --setting titv0.5 --reps 25
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from itertools import product
from pathlib import Path

import pandas as pd

# Repo root goes on sys.path so `simulation.*` imports resolve.
ROOT = Path(__file__).resolve().parents[1]
# Simulation outputs and figures live inside this package.
DATA = Path(__file__).resolve().parents[0]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from simulation.runner import run_cell  # noqa: E402

RESIM_DIR = DATA / "results" / "resim"

# ---------------------------------------------------------------------
# Simulation settings.  One entry = one folder = one fixed configuration.
# `id` is the folder name and is meant to be readable at a glance.
# ---------------------------------------------------------------------
SETTINGS = {
    # --- grid "Kv": K x v sweep at one fixed (data_R, test_R) ------------
    # Panels D/E of the sec-2.2 combined figure: matched uniform
    # substitution model, G*U non-canonical set on both data and test.
    "titv0.5": {
        "grid": "Kv",
        "id": "titv-data0.5-test0.5__nc-dataGU-testGU__n1000",
        "description": "Matched Ti/Tv = 0.5 (uniform identities), G*U both sides. "
                       "Regenerates the panel D/E operating envelope.",
        "ti_tv_ratio": 0.5,
        "test_ti_tv_ratio": 0.5,
        "data_noncanon": "GU",
        "test_noncanon": "GU",
        "n_alt": 1000,
        "n_null": 1000,
    },
    # Matched at the biologically realistic Ti/Tv.
    "titv2": {
        "grid": "Kv",
        "id": "titv-data2.0-test2.0__nc-dataGU-testGU__n1000",
        "description": "Matched Ti/Tv = 2.0, G*U both sides.",
        "ti_tv_ratio": 2.0,
        "test_ti_tv_ratio": 2.0,
        "data_noncanon": "GU",
        "test_noncanon": "GU",
        "n_alt": 1000,
        "n_null": 1000,
    },
    # --- grid "RxR": (data_R x test_R) square at fixed K, v --------------
    # Panels F/G: the misspecified-test heatmaps.  Same 7x7 square, K, v,
    # and n as the existing single-run results/e3_titv_grid, so the
    # replicated version is directly comparable to it.
    "e3square": {
        "grid": "RxR",
        "id": "e3-RxR-square__K5-v5__nc-dataGU-testGU__n2500",
        "description": "7x7 (data_R x test_R) square at K=5, v=5, uniform, G*U "
                       "both sides, n=2500. Replicates panels F/G.",
        "K": 5,
        "v": 5,
        "abundance": "uniform",
        "data_noncanon": "GU",
        "test_noncanon": "GU",
        "n_alt": 2500,
        "n_null": 2500,
    },
}

# (data_R, test_R) square for the "RxR" grid -- matches R_SQUARE in
# e3_grid.py and the R_SQUARE actually plotted by plot_combined_sec22.py.
R_SQUARE = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0]

GRID_K = (4, 5, 6, 8, 10)
GRID_V = (2, 3, 4, 5, 6, 7, 8)
NB_WORKERS = 8
BASE_SEED = 1000
BASE_SEED_RXR = 3030  # distinct block from the Kv grid (matches e3_grid.py)
REP_STRIDE = 10007  # prime, >> the within-grid seed span, so replicate
                    # seed blocks cannot overlap

_COMMON_FIELDS = [
    "rep", "seed",
    "n_alt_tested", "n_null_tested", "n_alt_flagged", "n_null_flagged",
    "power", "empirical_fdp",
    "power_sve", "power_svp", "power_mixed",
    "elapsed_sec",
]
FIELDNAMES = {
    "Kv": ["K", "v", "abundance"] + _COMMON_FIELDS,
    "RxR": ["data_R", "test_R", "K", "v", "abundance"] + _COMMON_FIELDS,
}


def build_cells(cfg: dict, abundance_arg: str) -> list[dict]:
    """One dict per cell: the run_cell kwargs that vary across the grid,
    plus the key fields used for resume and for the summary row."""
    if cfg["grid"] == "Kv":
        abundances = (["uniform", "skewed"] if abundance_arg == "both"
                      else [abundance_arg])
        return [
            {"K": K, "v": v, "abundance": ab,
             "ti_tv_ratio": cfg["ti_tv_ratio"],
             "test_ti_tv_ratio": cfg["test_ti_tv_ratio"]}
            for K, v, ab in product(GRID_K, GRID_V, abundances)
        ]
    if cfg["grid"] == "RxR":
        return [
            {"K": cfg["K"], "v": cfg["v"], "abundance": cfg["abundance"],
             "data_R": dr, "test_R": tr,
             "ti_tv_ratio": dr, "test_ti_tv_ratio": tr}
            for dr in R_SQUARE for tr in R_SQUARE
        ]
    raise ValueError(f"unknown grid type {cfg['grid']!r}")


def cell_seed(cfg: dict, cell: dict, rep: int) -> int:
    """Deterministic per (cell, replicate).  REP_STRIDE is far larger than
    the within-grid span, so replicate seed blocks cannot overlap."""
    if cfg["grid"] == "Kv":
        return (BASE_SEED + REP_STRIDE * rep + 100 * cell["K"] + 10 * cell["v"]
                + (0 if cell["abundance"] == "uniform" else 1))
    # RxR: index the R values rather than hashing floats
    di = R_SQUARE.index(cell["data_R"])
    ti = R_SQUARE.index(cell["test_R"])
    return BASE_SEED_RXR + REP_STRIDE * rep + 100 * di + 10 * ti


def cell_key(cfg: dict, cell: dict, rep: int) -> tuple:
    if cfg["grid"] == "Kv":
        return (cell["K"], cell["v"], cell["abundance"], rep)
    return (cell["data_R"], cell["test_R"], rep)


def cell_tag(cfg: dict, cell: dict) -> str:
    """Filename stem for the per-replicate p-value file."""
    if cfg["grid"] == "Kv":
        return f"K{cell['K']}_v{cell['v']}_{cell['abundance']}"
    return f"dR{cell['data_R']}_tR{cell['test_R']}"


def _git_commit(repo: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def write_config(out_dir: Path, setting_key: str, cfg: dict, reps: int,
                 cells: list) -> None:
    import numpy, pandas as _pd
    import os
    provenance = {
        "setting_key": setting_key,
        "setting": cfg,
        "reps": reps,
        "n_cells": len(cells),
        "grid_K": list(GRID_K),
        "grid_v": list(GRID_V),
        "base_seed": BASE_SEED,
        "rep_stride": REP_STRIDE,
        "seed_formula": "BASE_SEED + REP_STRIDE*rep + 100*K + 10*v + (0 uniform / 1 skewed)",
        "nb_workers": NB_WORKERS,
        "struct_commit": _git_commit(ROOT),
        "python": sys.version.split()[0],
        "numpy": numpy.__version__,
        "pandas": _pd.__version__,
        "PYTHONHASHSEED": os.environ.get("PYTHONHASHSEED", "<unset>"),
    }
    (out_dir / "config.json").write_text(json.dumps(provenance, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--setting", required=True, choices=sorted(SETTINGS))
    ap.add_argument("--reps", type=int, default=25)
    ap.add_argument("--abundance", default="uniform",
                    choices=["uniform", "skewed", "both"])
    args = ap.parse_args()

    cfg = SETTINGS[args.setting]
    out_dir = RESIM_DIR / cfg["id"]
    (out_dir / "pvals").mkdir(parents=True, exist_ok=True)

    cells = build_cells(cfg, args.abundance)
    fieldnames = FIELDNAMES[cfg["grid"]]
    write_config(out_dir, args.setting, cfg, args.reps, cells)

    summary_csv = out_dir / "grid_summary.csv"
    done = set()
    if summary_csv.exists():
        with open(summary_csv) as f:
            for row in csv.DictReader(f):
                if cfg["grid"] == "Kv":
                    done.add((int(row["K"]), int(row["v"]),
                              row["abundance"], int(row["rep"])))
                else:
                    done.add((float(row["data_R"]), float(row["test_R"]),
                              int(row["rep"])))
        print(f"Resuming: {len(done)} (cell, rep) runs already complete")

    write_header = not summary_csv.exists()
    f_out = open(summary_csv, "a", newline="")
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    if write_header:
        writer.writeheader(); f_out.flush()

    jobs = [(r, c) for r in range(args.reps) for c in cells]
    todo = [j for j in jobs if cell_key(cfg, j[1], j[0]) not in done]
    print(f"setting  : {cfg['id']}  (grid {cfg['grid']})")
    print(f"cells    : {len(cells)}   reps: {args.reps}   runs to do: {len(todo)}")
    t0 = time.time()

    scratch = out_dir / "_scratch"
    for i, (rep, cell) in enumerate(todo, start=1):
        seed = cell_seed(cfg, cell, rep)
        if scratch.exists():
            shutil.rmtree(scratch)
        c0 = time.time()
        res = run_cell(
            K=cell["K"], v=cell["v"], abundance_kind=cell["abundance"],
            n_alt=cfg["n_alt"], n_null=cfg["n_null"],
            out_dir=scratch, seed=seed, nb_workers=NB_WORKERS,
            ti_tv_ratio=cell["ti_tv_ratio"],
            test_ti_tv_ratio=cell["test_ti_tv_ratio"],
            data_noncanon=cfg["data_noncanon"],
            test_noncanon=cfg["test_noncanon"],
        )
        elapsed = round(time.time() - c0, 1)

        # Slim per-replicate artifact: anchor-level p-values only.
        sdf = pd.read_csv(scratch / "struct_results/structure_on_targets.tsv",
                          sep="\t", usecols=["anchor", "anchor_p"]).drop_duplicates("anchor")
        mdf = pd.read_csv(scratch / "meta.tsv", sep="\t", usecols=["anchor", "is_alt"])
        pv_dir = out_dir / "pvals" / f"rep{rep:03d}"
        pv_dir.mkdir(parents=True, exist_ok=True)
        mdf.merge(sdf, on="anchor", how="inner").to_csv(
            pv_dir / f"{cell_tag(cfg, cell)}.tsv", sep="\t", index=False)

        by_mode = res.get("power_by_submode", {})
        row = {
            "K": cell["K"], "v": cell["v"], "abundance": cell["abundance"],
            "rep": rep, "seed": seed,
            "n_alt_tested": res.get("n_alt_tested", 0),
            "n_null_tested": res.get("n_null_tested", 0),
            "n_alt_flagged": res.get("n_alt_flagged", 0),
            "n_null_flagged": res.get("n_null_flagged", 0),
            "power": res.get("power", float("nan")),
            # one BH run = one realized FDP, not an FDR
            "empirical_fdp": res.get("empirical_fdr", float("nan")),
            "power_sve": by_mode.get("sve", float("nan")),
            "power_svp": by_mode.get("svp", float("nan")),
            "power_mixed": by_mode.get("mixed", float("nan")),
            "elapsed_sec": elapsed,
        }
        if cfg["grid"] == "RxR":
            row["data_R"] = cell["data_R"]
            row["test_R"] = cell["test_R"]
        writer.writerow(row)
        f_out.flush()

        tot = time.time() - t0
        eta = (tot / i) * (len(todo) - i)
        print(f"[{i}/{len(todo)}] rep{rep:03d} {cell_tag(cfg, cell):<22s} "
              f"power={res.get('power', float('nan')):.3f} "
              f"fdp={res.get('empirical_fdr', float('nan')):.3f} "
              f"({elapsed:.1f}s) | {tot/60:.1f}m elapsed, ETA {eta/60:.1f}m",
              flush=True)

    if scratch.exists():
        shutil.rmtree(scratch)
    f_out.close()
    print(f"done: {out_dir}")


if __name__ == "__main__":
    main()
