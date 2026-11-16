# Pb210Modeler
An R script for performing classical Pb-210 age-depth modeling (CF/CRS, CFCS, and CA/CIC).
[![DOI](https://zenodo.org/badge/1058859255.svg)](https://doi.org/10.5281/zenodo.17155594)

Pb210Modeler is compatible with Windows, MacOS, and Linux operating systems, and is run within the RStudio IDE. All data read/write operations are performed with Excel (.xlsx) files.

This script is a revised R implementation of the popular spreadsheet-based Pb-210 modeling solutions standardized in [Sanchez-Cabeza & Ruiz-Fernández, 2012](https://www.sciencedirect.com/science/article/pii/S0016703711001402?casa_token=Mh8tjFdxeIUAAAAA:BH9R_Yq80HbTUD866QAw_jI5fU-AH69YkyaaA6gTNu0wk0Elya17x1n1_ZBHEntuLI_EHmUGdA). Calculations are performed automatically, guided by a TRUE/FALSE user input dialogue.

Pb210Modeler accepts both Alpha and Gamma activity data (dpm/g), and accepts dry bulk density (g/cm^3) for mass-related calculations. Activity and mass data should be provided as separate Excel (.xlsx) files. See alpha.xlsx, gamma.xlsx, and dbd.xlsx under the activity data and mass data folders for formatting conventions.

Any missing dry bulk density values are calculated with linear interpolation via linear depth, while missing excess Pb-210 activity is calculated via exponential interpolation via mass depth. Missing uncertainties are calculated via error propagation of boundary samples. Interpolated activities are excluded from all subsequent modeling calculations with the exception of the CF inventory calculations.  

Pb210Modeler performs automatic background determination for Alpha data and optional user-guided background determination. Background activity is calculated as a mean of all measurements below equilibrium depth and uncertainty is calculated as the standard error of the mean. Pb210Modeler can optionally correct Gamma total Pb-210 activity for self absorption effects when provided a .xlsx file of attenuation coefficients for each sediment interval, and can decay correct excess Pb-210 activity to the date of core collection when provided with a .xlsx file of sample measurement dates.

Pb210Modeler includes an optional dialogue for importing porosity data to correct linear depths for compaction effects by the methods described in [Matsumoto 1975](https://www.jstage.jst.go.jp/article/geochemj1966/9/3/9_3_167/_article/-char/ja/), which are passed along to CFCS and CA models to establish compaction-corrected SARs. It also includes a optional dialogue for estimating local atmospheric Pb-210 flux via global curves established in [Zhang et al. 2021](https://essd.copernicus.org/articles/13/2963/2021/) to investigate sediment focusing effects.

Also included are automatic and user-guided options for surface active zone (SAZ) determination, automatic CFCS model fitting with the segmented package, and optional manual CFCS model fitting. Initial activity C(0) and initial dry bulk density DBD(0) at the surface may be extrapolated by user-selected expanding window linear models (log scale activity for C(0)) calculated from points 1-10, or provided as direct measurements. Data tables (.xlsx) and accompanying plots (.pdf) are saved to model-specific folders. After performing age-depth modeling, Pb210Modeler includes a dialogue for calculating the accumulation rates of other materials within the core, either as a weight fraction (g/g, i.e. carbon, PFAS, etc.) or as a particle concentration (p/g, microplastics, fly ash, etc.).

Pb210Modeler checks for model assumption violations due to a SAZ, missing inventory, age inversions, or negative SARs, and calculates Pb-210 atmospheric flux to aid CF model validation.

Required libraries include jsonlite, readxl, writexl, ggplot2, scales, segmented, and rstudioapi.

[Example data](https://app.geosamples.org/sample/igsn/10.58052/IEJDH0002) included here were kindly provided by [Josh Himmelstein](https://github.com/joshimmel).
