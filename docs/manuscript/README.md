# Manuscript files

Sources and build tools for the paper. Rendered outputs are **not** stored here.

## Core files

- `manuscript_main.Rmd` main manuscript source.
- `manuscript_supplement.Rmd` supplementary materials source.
- `render_manuscripts.R` renders both documents, copies the main figures and builds the Zenodo bundle.
- `styles/reference.docx` Word reference style template used by both Rmd files.

## Outputs

All rendered files go to `final_version/` at the project root (git-tracked, accepted versions):

- `sampling_uncertainty_mpm_main_final.docx`
- `sampling_uncertainty_mpm_supplement_final.docx`
- `figures/Figure_1.png` to `Figure_6.png` (300 dpi, journal file names)

Figures themselves are created by the analysis scripts in `figures/` (the only copy of each).

## Zenodo

- `zenodo/` bundle workspace (`staging/` and `archive/` are git-ignored).
- `submission/` holds the deposited Zenodo zip, its description text and the cover letters.

## Usage

From the project root in R:

```r
source("docs/manuscript/render_manuscripts.R")
```

Rendering requires a UTF-8 locale and pandoc (RStudio bundles one).
Peer-review history and earlier versions are in `peer_review/`.
