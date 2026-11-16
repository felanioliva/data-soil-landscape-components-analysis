# Reproducibility files — soil landscape analyses

This folder contains the source data and five R scripts used to reconstruct the numerical analyses and selected figures accompanying the manuscript.

## Input data

- `data/données paysages corrigées 2.xlsx`, worksheet `Feuil1` (156 observations). The source workbook is read but never modified.
- `01_data_preparation.R` checks the source columns and creates the in-memory `df_clean` dataset. The source data use landscape codes P1–P3; the PCA graphics display the manuscript labels L1–L3.

## Requirements

R with packages `readxl`, `ggplot2`, `patchwork`, `vegan`, `lme4`, `performance`, `dplyr`, `tidyr`, `multcompView`, and `car`. Install missing packages using `install.packages()`.

Set the working directory to this folder **before** running the scripts, so relative data and output paths resolve correctly. In RStudio, open the folder as a project or use `setwd()` to point to it.

## Execution order

```r
source("01_data_preparation.R")
source("02_correlations_PCA.R")
source("03_PERMANOVA.R")
source("04_Partitionnement_variables.R")
source("05_Figure6.R")
```

Run all scripts in the same R session, because scripts 02–05 use the `df_clean` object created by script 01.

| Script | Main purpose | Manuscript outputs |
| --- | --- | --- |
| `01_data_preparation.R` | Data validation, renaming and derived variables | Input to subsequent analyses |
| `02_correlations_PCA.R` | Pearson correlation matrix and standardized PCA | Supplementary Table S1; Figures S2A, S2B, 4B and 4C |
| `03_PERMANOVA.R` | Euclidean-distance PERMANOVA and multivariate dispersion tests | PERMANOVA results associated with Table 2 |
| `04_Partitionnement_variables.R` | Sequential contributions based on marginal R² from mixed models | Numerical results for Figure 5 |
| `05_Figure6.R` | One-way ANOVA, diagnostics, Tukey comparisons and graphical panels | Figure 6 |

Script 02 writes CSV and graphic files to `outputs/`. Script 05 writes figure files and supporting statistical CSVs to `Figure6_outputs/`. Scripts 03 and 04 print their results to the R console; they do not automatically export a final figure.

## Interpretation notes

- The PCA uses 12 standardized original variables, whereas the Pearson correlation matrix includes 14 variables. PCA grouping ellipses are normal ellipses at the 75% level.
- PERMANOVA uses 999 unrestricted permutations with a fixed random seed (`set.seed(123)`). Landscape dispersion is heterogeneous; significant PERMANOVA results should therefore not be interpreted as location effects alone. Sequential sums of squares depend on factor order. Repeated depths within plots are not explicitly accommodated in the unrestricted permutations.
- Figure 5 contributions are sequential differences in marginal R². Some mixed models yield singular-fit warnings, which have been retained rather than suppressed. The unexplained fraction is `1 - marginal R²`, not a pure conditional residual variance.
- Figure 6 reproduces the original one-way ANOVA workflow using all 156 depth observations. Potential within-plot dependence is not modeled. Diagnostic tests indicate departures from some ANOVA assumptions; the reconstruction is not an independent validation of these assumptions.

The scripts were executed successfully in an RStudio session, but exact statistical validity and agreement with every manuscript table should be evaluated separately from successful execution. The published Table 3 analysis is not included in this set of scripts.
