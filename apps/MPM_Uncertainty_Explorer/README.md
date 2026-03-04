# MPM Uncertainty Explorer

This Shiny app provides a matrix-population-model uncertainty explorer using separate `U` and `F` matrices with fixed-dimension tabs.

## What it does

- Supports three fixed matrix sizes via tabs: `2x2`, `3x3`, `4x4`.
- Left panel: cell-level sampling distributions for `U` and `F` matrix entries.
- Right panel: sampling distributions for derived quantities.
- Bottom table: point estimate, posterior median, and 95% interval summary for derived quantities.

## Inputs and assumptions

- Stage sample sizes are entered as one value per source stage (`n[j]`).
- `U` inputs are transition counts among living stages.
- `F` inputs are recruit counts from each source stage.
- `U` is normalised by `n[j]`; excess counts are snapped down to respect stage totals.
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
