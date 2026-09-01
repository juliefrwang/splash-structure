# Reproducing the simulation figures

Every Monte-Carlo result in the manuscript is produced by the scripts in this
directory, running against the STRUCT implementation in `struct_rna/` of this
same repository. Each figure below lists the exact commands, seeds, and
runtime.


## Figure 4 — calibration and power

Panels A–C are closed-form and take no simulation input. They are evaluated
directly from `struct_rna.src.get_pval` when the figure is drawn.

Panels D–E, the (K, v) operating envelope. 35 cells x 25 replicates = 875 runs,
1000 alternative and 1000 null anchors each. Measured runtime 109 min at
`NB_WORKERS = 8`.

```bash
PYTHONHASHSEED=0 python -m simulation.run_replicates --setting titv0.5 --reps 25
```

Panels F–G, the 7x7 (data-R, test-R) square at K=5, v=5. 49 cells x 10
replicates = 490 runs, 2500 anchors of each class. Measured runtime 139 min.

```bash
PYTHONHASHSEED=0 python -m simulation.run_replicates --setting e3square --reps 10
```

Then draw the figure:

```bash
python -m simulation.plot_fig4_calibration_power
```

Seeds are deterministic per (cell, replicate) and recorded in the `seed` column
of every `grid_summary.csv`, so an individual point can be re-run on its own:

    (K, v) envelope   seed = 1000 + 10007 * rep + 100 * K + 10 * v + (0 uniform, 1 skewed)
    (R, R) square     seed = 3030 + 10007 * rep + 100 * i(data_R) + 10 * i(test_R)

where `i(.)` indexes `R_SQUARE = [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0]`. The
stride of 10007 exceeds the within-grid seed span, so replicate seed blocks
cannot overlap.

## Table 2 and Supplementary Figure 3 — non-canonical pair set

Ten (data-N, test-N) configurations x v in {3,4,5,6,7}, K=5, uniform abundance,
R=1/2, 2500 anchors of each class, one run per cell. Both drivers append to the
same summary file and must be run in this order. Measured runtime 13 min.

```bash
PYTHONHASHSEED=0 python -m simulation.run_noncanon_grid             # seven configurations
PYTHONHASHSEED=0 python -m simulation.run_noncanon_grid_wcf_configs   # three WCF-involving configurations
python -m simulation.plot_suppfig3_noncanon_power
```

    seed = BASE_SEED + 1000 * index(v) + 13 * index(config)

with `BASE_SEED = 2025` in `run_noncanon_grid.py` and `9090` in `run_noncanon_grid_wcf_configs.py`,
chosen so the two driver's seeds cannot collide.

Table 2 reports power as flagged alternative anchors divided by the 2500
planted, not by the number actually tested. Anchors in which the stem-finder
locates no stem under the test-side pair set are therefore counted as misses.
The `power` column of `grid_summary.csv` uses the tested denominator instead, so
the two differ in the configurations where the test admits fewer or different
non-canonical pairs than the data contains.

## Committed results

`results/` holds the summary file for each grid plus the `config.json` recording
the commit, seeds, and package versions of the run. These are enough to redraw
every figure without re-running the simulation:

```bash
python -m simulation.plot_fig4_calibration_power
python -m simulation.plot_suppfig3_noncanon_power
```

Per-replicate anchor p-values and the synthetic SPLASH inputs are not committed.
They come to several GB and are fully regenerable from the seeds above.
