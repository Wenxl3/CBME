# CBME (Simulation)

This project incorporates the simulation codes for "Covariate balancing with measurement error".

## Packages

To successfully reproduce the results, please check whether the following R packages are available:

```R
library(MASS)
library(ebal)
library(WeightIt)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(gtools)
library(stats)
library(nleqslv)
library(rootSolve)
library(openxlsx)
```

* ## Simulations in the the main paper

  * To create FIGURE 1 in "Section 5.1 Impact of measurement error" of the main paper, run the R scripts located in the folder `01BiasAnalysis_ASMD`.
 
  * To create TABLE 1 in "5.2 ATT estimation with measurement error correction" of the main paper, run the R scripts located in the folder `02Simulation`.
 
  * To create FIGURE 2 in "5.3 ASMD comparison" of the main paper, run the R scripts located in the folder `03ASMDComparison4Simulation`.

* ## Other simulations in the Supplementary Material

  * To create Table 1 in Section "Impact of measurement error: naive analysis for IPW, CBPS and SBW" of the Supplementary Material, run the R scripts located in the folder `06SuppSim_NaiveAnalysis`.
 
  * To create Table 2 in Section "Simulation results for corrected CBPS" of the Supplementary Material, run the R scripts located in the folder `07SuppSim_CBPSME`.
 
  * To create Table 3 and Figure 1 in Section "Correction for EB on estimating the ATE" of the Supplementary Material, run the R scripts located in the folder `10SuppSim_ATEv2`.
 
  * To create Table 4 in Section "Impact of misspecification of the measurement error model on the different methods" of the Supplementary Material, run the R scripts located in the folder `08SuppSim_DifferentMEModels`.
 
  * To Figure 2 and Figure 3 in Section "Potential multicollinearity among covariates" of the Supplementary Material, run the R scripts located in the folder `09SuppSim_Multicollinearity`.
