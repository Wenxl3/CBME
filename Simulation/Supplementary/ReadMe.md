# R scripts for the supplementary material
The R scripts to conduct simulation experiments and generate all figures and tables in the supplementary material.

Structure each folder with two subdirectories: `Code`, containing R scripts for generating figures and tables; and `Result`, containing the corresponding output figures and tables.

## Files:

- **TAB S1:**

  The Monte Carlo simulation to estimate the average treatment effect (ATE) using naive and corrected covariate balancing propensity score (CBPS) methods.

- **TAB S2 & FIG S1:**
  
  The Monte Carlo simulation to evaluate the performance of corrected entropy balancing methods for ATE estimation
  - **Table S2:** Numerical ATE estimates from naive and corrected entropy balancing methods. 
  - **Figure S1:** Trends in bias, standard deviation (SD), and mean squared error (MSE) of the ATE estimates as measurement error variance increases.

- **TAB S3:**

  The Monte Carlo simulation to assess the impact of measurement error on three naive weighting methods: inverse probability weighting (IPW), over-identified covariate balancing propensity score (CBPS; Imai & Ratkovic, 2014), and stable balancing weights (SBW; Zubizarreta, 2015).
  
- **TAB S4:**

  The Monte Carlo simulation to assess the robustness of the proposed corrected estimators under alternative (misspecified) measurement error models, reporting ATT estimates for naive EB, CEB, BCEB, CEB-HL, and CEB-HW methods.
  
- **FIG S2:**

  The Monte Carlo simulation to assess the robustness of the proposed corrected methods to multicollinearity among covariates in the presence of measurement error.
