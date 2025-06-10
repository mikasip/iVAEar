## Identifiable Autoregressive Variational Autoencoders for Nonlinear and Nonstationary Spatio-Temporal Blind Source Separation
This repository contains the code and supplementary materials for reproducing the results presented in paper "Identifiable Autoregressive Variational Autoencoders for Nonlinear and Nonstationary Spatio-Temporal Blind Source Separation" published in ECMLPKDD 2025.

## Repository Structure
```
├── data/                                   # Data files
├── ar_mismatch_sim/                        # Files for ar order mismatch simulations
├── hyper_param_search/                     # Files for case study hyper param search
├── simulations/                            # Files for the main simulations
├── athens_air_pollution_case_study.r       # A file to reproduce the results of the case study
├── Supplementary_Material_iVAEar.pdf       # Supplementary material
└── README.md                               # This file
```

## Dependencies

The main methods are in R package ```NonlinearBSS```. To ensure that the right versions of the methods are used for reproducing the results, the package can be downloaded in R as follows:
```
devtools::install_github("mikasip/NonlinearBSS@4e3e8e3c94b9858a6b467c1e965f9f6f7da470bd")
```

## Data

For the original data used in the case study, see [1].

## References

[1] Angelis, G.F., Emvoliadis, A., Theodorou, T.I., Zamichos, A., Drosou, A., Tzo-
varas, D.: Regional datasets for air quality monitoring in European cities. In:
IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Sym-
posium. pp. 6875–6880 (2024)
