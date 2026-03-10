# MPM Uncertainty Explorer

This Shiny app provides a matrix-population-model uncertainty explorer with fixed-dimension tabs and single-matrix (`A`) data entry.

## What it does

- Supports four fixed matrix sizes via tabs: `2x2`, `3x3`, `4x4`, `5x5`.
- Left panel: cell-level sampling distributions for `A` matrix entries, color-coded by internal `U`/`F` assignment.
- Right panel: sampling distributions for derived quantities.
- Bottom table: point estimate, posterior median, and 95% interval summary for derived quantities.

## Inputs and assumptions

- Stage sample sizes are entered as one value per source stage (`n[j]`).
- Data entry is through a single `A` matrix.
- Internal mapping is fixed: top row except `[1,1]` is treated as `F`; all other cells are treated as `U`.
- `U` is normalised by `n[j]`; excess `U` counts are snapped down to respect stage totals.
- Optional structural-zero mode fixes entered zeros at zero during posterior sampling.
- Posterior draws are selectable (`300`, `500`, `1000`; default `500`).

Posterior model:

- `U`: Dirichlet posterior over living transitions + death residual per source stage.
- `F`: Gamma posterior for Poisson count/rate formulation.

Derived quantities shown:

- Population growth rate (`lambda`)
- Mature life expectancy (`L`)
- Generation time (`T`)
- Damping ratio (`rho`)

## Run

From the repository root in R:

```r
shiny::runApp("apps/MPM_Uncertainty_Explorer")
```
