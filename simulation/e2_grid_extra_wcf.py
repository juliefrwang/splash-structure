"""E2 — two extra configurations bringing N=empty into the
(data N, test N) sweep.

The existing E2 grid (e2_grid.py / results/e2_noncanon/) covers seven
(data, test) configs that all involve GU and/or GA. Per discussion
2026-05-25, we add two more configs that put N=empty on one side:

    (data_noncanon, test_noncanon):
      ("none", "GU")   - data uses strict WCF; test admits G·U.
      ("GU",   "none") - data uses GU; test admits strict WCF only.

Same K, v, abundance, R, and n_alt=n_null=2500 as the existing E2
(matching the existing CSV's row schema). Appends to the same
grid_summary.csv via the resume logic.

BASE_SEED bumped to avoid colliding with seeds used by the original
seven configs.

Run:
    python \\
        -m simulation.e2_grid_extra_wcf
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

CONFIGS = [
    ("none", "GU"),
    ("GU",   "none"),
    ("none", "none"),   # WCF-only matched baseline; added 2026-05-25
]
V_VALUES = (3, 4, 5, 6, 7)
K = 5
ABUNDANCE = "uniform"
N_ALT = N_NULL = 2500
NB_WORKERS = 8
BASE_SEED = 9090                       # distinct from e2_grid.py's BASE_SEED (2025)
GRID_DIR = "e2_noncanon"               # share dir; append to the same CSV


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
        print(f"Resuming: {len(done)} cells already complete in {summary_csv}")

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
    print(f"{len(CONFIGS)} new configs × {len(V_VALUES)} v = {total} cells")
    t0 = time.time()

    for i, (dn, tn, v) in enumerate(cells, start=1):
        if (dn, tn, v) in done:
            print(f"[{i}/{total}] skip {dn}/{tn} v={v} (already done)")
            continue

        tag = _tag(dn, tn) + f"_v{v}"
        cell_dir = out_root / tag
        seed = BASE_SEED + 1000 * V_VALUES.index(v) + 13 * CONFIGS.index((dn, tn))

        c_t0 = time.time()
        result = run_cell(
            K=K, v=v, abundance_kind=ABUNDANCE,
            n_alt=N_ALT, n_null=N_NULL,
            out_dir=cell_dir, seed=seed, nb_workers=NB_WORKERS,
            data_noncanon=dn, test_noncanon=tn,
        )
        result["elapsed_sec"] = round(time.time() - c_t0, 1)
        with open(cell_dir / "summary.json", "w") as fj:
            json.dump(result, fj, indent=2)

        nc = result.get("power_by_n_nc", {})
        nc_hist = nc.get("hist", {})
        row = {
            "data_noncanon": dn, "test_noncanon": tn,
            "K": K, "v": v, "abundance": ABUNDANCE,
            "n_alt_tested": result.get("n_alt_tested", 0),
            "n_null_tested": result.get("n_null_tested", 0),
            "n_alt_flagged": result.get("n_alt_flagged", 0),
            "n_null_flagged": result.get("n_null_flagged", 0),
            "power": result.get("power", float("nan")),
            "empirical_fdr": result.get("empirical_fdr", float("nan")),
            "power_has_nc": nc.get("has_nc", {}).get("power", float("nan")),
            "n_has_nc": nc.get("has_nc", {}).get("n", 0),
            "power_all_wcf": nc.get("all_wcf", {}).get("power", float("nan")),
            "n_all_wcf": nc.get("all_wcf", {}).get("n", 0),
            "nc_hist_json": json.dumps(nc_hist),
            "elapsed_sec": result["elapsed_sec"],
        }
        writer.writerow(row)
        f_out.flush()

        el = time.time() - t0
        eta = (el / i) * (total - i)
        print(f"[{i}/{total}] {dn}/{tn} v={v} "
              f"power={row['power']:.3f} fdr={row['empirical_fdr']:.3f} "
              f"({result['elapsed_sec']:.0f}s) ETA {eta/60:.1f}m")

    f_out.close()
    print(f"\nDone. Appended to {summary_csv}")


if __name__ == "__main__":
    main()
