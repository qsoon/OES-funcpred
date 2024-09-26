# OES-funcpred

This repo contains the source code for a prediction method based on a functional approach
using multivariate functional partial least squares regression, coupled with dimension reduction and a novel outlier detection technique via functional independent component analysis.


## Description

- Code
  - `1.function_define.R` is a code for functions used for the implementation.
  - `data_generation.R` is a code for the generation of simulation dataset.
  - `outlier_detection.R` is a code for outlier detection using FICA.
  - `prediction.R` is a code for prediction using MFPLS regression.


- Data
  - `wvcoefs.RData` contains wavelet coefficients to generate the underlying data X<sub>0</sub> for scenario (iv) in the simulation study. 
  - `stationary_true.RData`, `ns_stationary_true.RData`, `sc_true.RData`, and `wavelet_true.RData` contain the underlying data X<sub>0</sub>'s for scenarios (i), (ii), (iii), and (iv) in the simulation study. 

## Code overview
We provide the code used for the simulation study in this paper.