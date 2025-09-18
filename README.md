# Pb210Modeler
An R script for performing classical Pb-210 age-depth modeling (CF/CRS, CFCS, and CA/CIC).

This script is an R implementation of the popular spreadsheet-based Pb-210 modeling solutions standardized in Sanchez-Cabeza & Ruiz-Fernández, 2012. Calculations are performed automatically, guided by a TRUE/FALSE user input dialogue.

Pb210Modeler accepts both Alpha and Gamma activity data (dmp/g), and uses dry bulk density (g/cm^3) for mass-related calculations. Activity and mass data should be provided as separate Excel (.xlsx) files. 

Any missing dry bulk density and supported Pb-210 activities are calculated via linear interpolation, while missing total Pb-210 activity is calculated via exponential interpolation. Missing uncertainties are calculated via error propagation. 

Pb210Modeler performs automatic background determination for Alpha data with the changepoints package, optional user-guided background determination, user-guided surface active zone (SAZ) determination, automatic CFCS model fitting with the segmented package, and optional manual CFCS model fitting. Initial activity C(0) and DBD(0) at the surface are extrapolated by the best-fit linear models (as determined by adjusted R-square) calculated from points 1-3 through 1-10. Data tables (.csv) and accompanying plots (.pdf) are saved to model-specific folders.

Pb210Modeler also checks  for model assumption violations due to a SAZ, missing inventory, or age inversions. It also calculates Pb-210 atmospheric flux to aid model validation.

Required libraries include readxl, writexl, ggplot2, scales, segmented, dplyr, zoo, and changepoint.

Example data included here were kindly provided by Josh Himmelstein (https://app.geosamples.org/sample/igsn/10.58052/IEJDH0001).
