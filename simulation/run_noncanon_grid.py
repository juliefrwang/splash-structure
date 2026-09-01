"""Non-canonical pair set grid: power and FDR across (data-N, test-N).

7 configs (data_noncanon / test_noncanon) — GU is the shared baseline
anchor; two 2×2 blocks (GA-vs-GU, GU,GA-vs-GU):

    GU/GU      GA/GA   GA/GU   GU/GA
               GU,GA/GU,GA   GU,GA/GU   GU/GU,GA

× v ∈ {3,4,5,6,7}, K=5, uniform abundance, titv=0.5,
n_alt=n_null=2500/cell.  35 cells.

Per cell we record aggregate power/FDR AND the binary n_nc
stratification (planted stem has ≥1 non-canonical pair vs all-WCF)
plus the realized n_nc histogram, so the N-dependent cap-truncation
is visible and power is comparable at matched composition.

Run from repo root:
    python \\
        -m simulation.run_noncanon_grid
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

# Repo root goes on sys.path so `simulation.*` imports resolve.
ROOT = Path(__file__).resolve().parents[1]
# Simulation outputs and figures live inside this package.
DATA = Path(__file__).resolve().parents[0]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from simulation.runner import run_cell  # noqa: E402

# (data_noncanon, test_noncanon)
CONFIGS = [
    ("GU", "GU"),          # baseline anchor
    ("GA", "GA"),          # matched single non-default
    ("GA", "GU"),          # test-narrower (data GA, test GU)
    ("GU", "GA"),          # test-wider (data GU, test GA)
    ("GU,GA", "GU,GA"),    # matched multi-pair
    ("GU,GA", "GU"),       # test-narrower (data GU,GA, test GU)
    ("GU", "GU,GA"),       # test-wider (data GU, test GU,GA)
]
V_VALUES = (3, 4, 5, 6, 7)
K = 5
ABUNDANCE = "uniform"
N_ALT = N_NULL = 2500
TITV = 0.5
NB_WORKERS = 8
BASE_SEED = 2025
GRID_DIR = "e2_noncanon"


def _tag(dn, tn):
    return f"{dn}__{tn}".replace(",", "")


def main():
    out_root = DATA / "results" / GRID_DIR
    out_root.mkdir(parents=True, exist_ok=True)
    summary_csv = out_root / "grid_summary.csv"

    done = set()
    if summary_csv.exists():
        with open(summary_csv) as f:
            for row in csv.DictReader(f):
                done.add((row["data_noncanon"], row["test_noncanon"],
                          int(row["v"])))
        print(f"Resuming: {len(done)} cells already complete")

    fieldnames = [
        "data_noncanon", "test_noncanon", "K", "v", "abundance",
        "n_alt_tested", "n_null_tested", "n_alt_flagged", "n_null_flagged",
        "power", "empirical_fdr",
        "power_has_nc", "n_has_nc", "power_all_wcf", "n_all_wcf",
        "nc_hist_json", "elapsed_sec",
    ]
    write_header = not summary_csv.exists()
    f_out = open(summary_csv, "a", newline="")
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    if write_header:
        writer.writeheader()
        f_out.flush()

    cells = [(dn, tn, v) for (dn, tn) in CONFIGS for v in V_VALUES]
    total = len(cells)
    t0 = time.time()

    for i, (dn, tn, v) in enumerate(cells, start=1):
        if (dn, tn, v) in done:
            continue
        cell_dir = out_root / f"{_tag(dn, tn)}_v{v}"
        seed = BASE_SEED + 1000 * V_VALUES.index(v) + 13 * CONFIGS.index((dn, tn))

        c_t0 = time.time()
        result = run_cell(
            K=K, v=v, abundance_kind=ABUNDANCE,
            n_alt=N_ALT, n_null=N_NULL,
            out_dir=cell_dir, seed=seed, nb_workers=NB_WORKERS,
            ti_tv_ratio=TITV, test_ti_tv_ratio=TITV,
            data_noncanon=dn, test_noncanon=tn,
        )
        result["elapsed_sec"] = round(time.time() - c_t0, 1)

        with open(cell_dir / "summary.json", "w") as fj:
            json.dump(result, fj, indent=2)

        strat = result.get("power_by_nc_stratum", {})
        has_nc = strat.get("has_nc", {})
        all_wcf = strat.get("all_wcf", {})
        row = {
            "data_noncanon": dn, "test_noncanon": tn,
            "K": K, "v": v, "abundance": ABUNDANCE,
            "n_alt_tested": result.get("n_alt_tested", 0),
            "n_null_tested": result.get("n_null_tested", 0),
            "n_alt_flagged": result.get("n_alt_flagged", 0),
            "n_null_flagged": result.get("n_null_flagged", 0),
            "power": result.get("power", float("nan")),
            "empirical_fdr": result.get("empirical_fdr", float("nan")),
            "power_has_nc": has_nc.get("power", float("nan")),
            "n_has_nc": has_nc.get("n", 0),
            "power_all_wcf": all_wcf.get("power", float("nan")),
            "n_all_wcf": all_wcf.get("n", 0),
            "nc_hist_json": json.dumps(result.get("nc_hist", {})),
            "elapsed_sec": result["elapsed_sec"],
        }
        writer.writerow(row)
        f_out.flush()

        el = time.time() - t0
        eta = (el / i) * (total - i)
        print(f"[{i}/{total}] {dn:>6}/{tn:<6} v={v} "
              f"power={row['power']:.3f} fdr={row['empirical_fdr']:.3f} "
              f"has_nc={row['power_has_nc']} (n={row['n_has_nc']}) "
              f"wcf={row['power_all_wcf']} (n={row['n_all_wcf']}) "
              f"({result['elapsed_sec']:.0f}s) ETA {eta/60:.1f}m")

    f_out.close()
    print(f"\nDone. Summary at {summary_csv}")


if __name__ == "__main__":
    main()
