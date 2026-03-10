# Manuscript files

This folder contains the paper source files, rendered outputs, style template, and Zenodo packaging workspace.

## Core files

- `manuscript_main.Rmd` main manuscript source.
- `manuscript_supplement.Rmd` supplementary materials source.
- `render_manuscripts.R` render script for both documents and Zenodo bundle staging.

## Rendered outputs

- `sampling_uncertainty_mpm_main.docx` rendered main manuscript.
- `sampling_uncertainty_mpm_supplement.docx` rendered supplement.

## Styling

- `styles/reference.docx` Word reference style template used by both Rmd files.

## Reproducibility packaging

- `zenodo/` Zenodo deposition workspace:
  - `staging/` temporary bundle builds (git-ignored except `README.md`).
  - `archive/` zipped bundles for deposition (git-ignored except `README.md`).

## Submission bundle

- `submission/` curated submission-ready files:
  - main and supplement `.docx`,
  - combined manuscript+supplement `.docx`,
  - latest Zenodo `.zip`,
  - plain-text Zenodo description,
  - draft cover letters,
  - main manuscript figure files.

## Usage

From the project root in R:

```r
source("docs/manuscript/render_manuscripts.R")
```

This renders both `.docx` files and then prepares a Zenodo staging bundle.
