# DKOsimR (Double-CRISPR Knockout Simulation in R)
## Description
DKOsimR is an R package for running Double-CRISPR Knockout Simulation (DKOsim). DKOsim is a simulation framework designed to simulate growth-based dual knockout CRISPR screens. It allows users and investigators to efficiently reproduce synthetic data where both the single gene fitness effect and the interaction of gene pairs can be pre-specified.

## Installation
DKOsimR is an open-source R package. To install DKOsimR and all required dependencies via Github, you may use the following commands in R:
```
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("yuegu-phd/DKOsimR")
devtools::install(dependencies = TRUE)
```

## Getting Started
Load the package and explore its tutorial documents in vignette
```
library(DKOsimR)
devtools::install(build_vignettes = TRUE)
vignette("DKOsimR") # see tutorial on how to generate synthetic CRISPR data using DKOsimR
browseVignettes("DKOsimR") # see the source code in tutorial
```

It might take a few minutes to build vignette. Alternatively, you may view the prebuilt vignette file after installing the package.
```
DKOsimR::open_dkosim_vignette_pdf()
```

## Overview of Study Design
![DKOsim workflow](assets/images/overview.jpg)

## Example Workflow


## References
- Gu, Y., Hart, T., Leon-Novelo, L., Shen, J.P.. Double-CRISPR Knockout Simulation (DKOsim): A Monte-Carlo Randomization System to Model Cell Growth Behavior and Infer the Optimal Library Design for Growth-Based Double Knockout Screens. bioRxiv 2025.09.11.675497. DOI: 10.1101/2025.09.11.675497.
- Shen, J., Zhao, D., Sasik, R. et al. Combinatorial CRISPR–Cas9 screens for de novo mapping of genetic interactions. Nat Methods 14, 573–576 (2017). DOI: 10.1038/nmeth.4225.
