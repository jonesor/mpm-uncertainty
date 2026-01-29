models/
-------
Stan model files used by the analysis scripts.

Contents
--------
- `regress*.stan`               Regression models (including error variants).
- `regress_hier*.stan`          Hierarchical regression models.
- `categorical_logit_hier*.stan` Categorical logit hierarchical models.
- `movbeta_gprc*.stan`          Moving beta GPRC models.
- `varcomp.stan`                Variance component model.
- `null*.stan`                  Null models.

Usage
-----
Compiled in scripts via `rstan::stan_model("models/<file>.stan")`.
