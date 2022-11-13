# The Presence and Impact of Multimorbidity Clusters on Adverse Outcome Across the Spectrum of Kidney Function

This repository contains statistical code written in [R](https://cran.r-project.org) that can be used to replicate the analyses reported in paper titled ["The Presence and Impact of Multimorbidity Clusters on Adverse Outcome Across the Spectrum of Kidney Function"](https://doi.org/10.1186/s12916-022-02628-2), published in BMC Medicine and available Open Access.

The structure of the code published in this repository is described in each subsection below, by topic.
Please note that:
- File names of actual data are excluded from each script, while intermediate objects (e.g., model fits) that are imported at later stages are retained;
- All scripts exports output (data, tables, figures) to certain directories that are not included in this repository.

## eGFR Interpolation

This topic concerns the following files:
* `01-fit-lmm.R`
* `02-interpolate.R`

The file `01-fit-lmm.R` can be used to fit a linear mixed effects model, which is then used to interpolate eGFR trajectories in file `02-interpolate.R`.
The interpolation procedure is described in more detail in the manuscript, please refer to that for all methodological considerations.

## Defining Multimorbidity Conditions

This topic concerns the following files:
* `03-multimorbidity.R`

This script is used to calculate each multimorbidity condition at every index date (i.e., eGFR threshold) considered in our analysis.
For more information of how each multimorbidity condition is defined refer to the manuscript.

## Clustering Analysis

This topic concerns the following files:
* `04-fit-clustering.R`
* `05-summarise-clustering.R`

The `04-fit-clustering.R` script can be used to fit the k-modes algorithm at each eGFR threshold and for a variety of possible number of clusters.
Elbow plots are produced to select the _optimal_ number of clusters; remember that this is specified by the analyst and cannot be fully automated.

The `05-summarise-clustering.R` script summarises the clustering analysis by creating plots and tables, and runs the algorithm to identify significant clusters.

## Survival Analysis and Predictive Performance

This topic concerns the following files:
* `06-outcomes.R`
* `30-predictive-performance.R`

The `06-outcomes.R` script runs all the outcomes analyses (MACE3, MACE4, mortality) and summarises them with a variety of tables and plots.
The `30-predictive-performance.R` script calculates predictive performance metrics for all the models we considered and compared.

## Sensitivity Analyses

This topic concerns the following files:
* `40-age-specific-below.R`
* `41-age-specific-above.R`

These two files can be used to run the whole analysis for subset of data defined by age (below and above 65 years of age, in our study).
Outcome tables and plots are exported as well.

# References

* Sullivan, M.K., Carrero, JJ., Jani, B.D. et al. The presence and impact of multimorbidity clusters on adverse outcomes across the spectrum of kidney function. BMC Med 20, 420 (2022). https://doi.org/10.1186/s12916-022-02628-2
