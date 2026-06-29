# Oral Mucositis Score Analysis - code only

This standalone directory contains only the R code used for the manuscript:

**Oral microbial signatures of head and neck cancer patients highlight diverse longitudinal patterns of oral mucositis severity**

No participant-level, microbiome, demographic, processed, cached or generated data are included.

## Data availability

The study data are restricted and are not publicly available because they contain potentially identifiable participant information. Data may be requested from the corresponding authors at MD Anderson. Access is subject to institutional, ethical and privacy review, approval of a Material Transfer Agreement and any applicable regulatory or governance requirements.

Do not commit study data or derived data to this repository. The included `.gitignore` excludes the expected data, cache and output locations and common data-file extensions.

## Directory layout

All scripts are intentionally kept in this directory. Run commands from this directory.

After approved access, place the supplied restricted inputs in:

```text
data/
```

The workflow creates these local directories, which are excluded from version control:

```text
data_filtered/
output/
figs/
figure_panels/
```

## Software

The revision was run with R 4.6. Required packages include:

- CRAN: `cluster`, `circlize`, `devtools`, `dplyr`, `ggplot2`, `ggpubr`, `magrittr`, `matrixStats`, `nlme`, `nsprcomp`, `patchwork`, `pheatmap`, `qs`, `readr`, `rmarkdown`, `seriation`, `stringr`, `tibble`, `tidyr`
- Bioconductor: `ComplexHeatmap`, `SummarizedExperiment`
- Method and visualisation packages: `PhiSpace`, `vizOmics`, `PLSKO`

The GitHub packages used by the analysis can be installed with:

```r
devtools::install_github("jiadongm/PhiSpace", subdir = "pkg")
devtools::install_github("ByronSyun/vizOmics")
```

Install `mixOmics` through Bioconductor before installing PLSKO if it is not already available. PLSKO 0.2.0 can then be installed with:

```r
devtools::install_github("guannan-yang/PLSKO/PLSKO", quiet = TRUE)
```

## Reproduction workflow

1. Preprocess the restricted symptom and microbiome inputs:

   ```r
   source("1.Preprocessing.R")
   source("2.prepareData.R")
   ```

2. Render the main analysis. This generates the primary score, trajectory, clustering, PLSKO and mixed-model results:

   ```r
   rmarkdown::render("Analysis.Rmd")
   ```

   Alternatively:

   ```bash
   quarto render Analysis.Rmd
   ```

3. Run the Round 4 scripts from this directory in the following order:

   ```r
   source("round4_cluster_validation_metrics.R")
   source("round4_demographic_pairwise_effects.R")
   source("round4_repeated_measure_lmm_taxa_screen.R")
   source("round4_age_weight_adjusted_lmm_taxa_screen.R")
   source("round4_filter_threshold_sensitivity.R")
   source("round4_prepare_patient_level_plsko_inputs.R")
   source("round4_run_patient_level_plsko_clean.R")
   source("round4_selected_taxa_supplementary_table.R")
   source("round4_regenerate_readable_figures.R")
   ```

`0.utils.R` contains helper functions sourced by the notebook and Round 4 scripts. The figure and supplementary-table scripts in this upload-only package write exclusively inside this standalone directory; manuscript-specific copy paths have been removed.

PLSKO computations are stochastic but use fixed seeds in the supplied code. Depending on the machine, the PLSKO steps may be computationally intensive.
