# Pb210Modeler
An R script for performing classical Pb-210 age-depth modeling (CF/CRS, CFCS, and CA/CIC).
[![DOI](https://zenodo.org/badge/1058859255.svg)](https://doi.org/10.5281/zenodo.17155594)

This script is an R implementation of the popular spreadsheet-based Pb-210 modeling solutions standardized in Sanchez-Cabeza & Ruiz-Fernández, 2012. Calculations are performed automatically, guided by a TRUE/FALSE user input dialogue.

Pb210Modeler accepts both Alpha and Gamma activity data (dpm/g), and uses dry bulk density (g/cm^3) for mass-related calculations. Activity and mass data should be provided as separate Excel (.xlsx) files. See alpha.xlsx, gamma.xlsx, and dbd.xlsx under the activity data and mass data folders for formatting conventions.

Any missing dry bulk density and supported Pb-210 activities are calculated via linear interpolation, while missing total Pb-210 activity is calculated via exponential interpolation. Missing uncertainties are calculated via error propagation. 

Pb210Modeler performs automatic background determination for Alpha data with the changepoints package and optional user-guided background determination. Background activity is calculated as a weighted mean and uncertainty is calculated as the standard error of the mean via REML modeling with a KNA test using the metafor package. The user may instead opt to use a simple arithmetic mean and standard error for background activity as well. 

Also included are user-guided surface active zone (SAZ) determination, automatic CFCS model fitting with the segmented package, and optional manual CFCS model fitting. Initial activity C(0) and initial dry bulk density DBD(0) at the surface are extrapolated by the best-fit linear models (as determined by adjusted R-square) calculated from points 1-3 through 1-10. Data tables (.csv) and accompanying plots (.pdf) are saved to model-specific folders. After performing age-depth modeling, Pb210Modeler includes an option for calculating the accumulation rates of other materials within the core, either as a weight fraction (g/g, i.e. carbon, PFAS, etc.) or as a particle concentration (p/g, microplastic, fly ash, etc.).

Pb210Modeler checks for model assumption violations due to a SAZ, missing inventory, age inversions, or negative SARs, and calculates Pb-210 atmospheric flux to aid CF model validation.

Required libraries include readxl, writexl, ggplot2, scales, segmented, dplyr, zoo, changepoint, metafor, and rstudioapi.

[Example data](https://app.geosamples.org/sample/igsn/10.58052/IEJDH0002) included here were kindly provided by [Josh Himmelstein](https://github.com/joshimmel).
