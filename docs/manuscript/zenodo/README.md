docs/manuscript/zenodo/
-----------------------
Workspace for preparing Zenodo deposition bundles.

Contents
--------
- `staging/`  Timestamped candidate bundles created by `scripts/prepare_zenodo_bundle.R`.
- `archive/`  Finalised bundles retained after internal checks and before/after deposition.

Workflow
--------
1. Run `Rscript --vanilla scripts/prepare_zenodo_bundle.R`.
2. Review the generated bundle in `docs/manuscript/zenodo/staging/`.
3. Validate `MANIFEST.csv` / `MANIFEST.md` and run `RUN_REPRODUCTION.R` inside the bundle.
4. Move approved bundle to `docs/manuscript/zenodo/archive/` and upload to Zenodo.
