# Bayesian-Analysis-of-Varicella-Vaccination-Coverage
A Bayesian study of childhood varicella (chickenpox) vaccination coverage in three US regions (Georgia, Wisconsin and Mississippi) across birth cohorts 2011–2020. The project moves from closed-form conjugate analysis to MCMC-based logistic regression, and uses the posterior to answer practical questions: Is coverage changing over time? Do regions differ? Was 2020 unusual?

[Read the documentation →](https://maryamshojaeiee-hub.github.io/DNA-Sequence-Assembly-with-de-Bruijn-Graphs/documentation.pdf) · [View the code →](https://maryamshojaeiee-hub.github.io/DNA-Sequence-Assembly-with-de-Bruijn-Graphs/main_codes.py) 

Estimated vaccination coverage by region and birth year

# Data
Survey counts of children who received at least one dose of the varicella vaccine, by region and birth year (2011–2020): number vaccinated (`Vaccinated`) out of the number surveyed (`Sample_Size`). Each region–year cell has roughly 180–330 children. The data are included directly in the R script.

# Methods

| Step | Approach |
|---|---|
| Coverage per region and year | Beta–Binomial conjugate model, posterior derived analytically |
| Prior sensitivity | Non-informative `Beta(1, 1)` vs. informative `Beta(50, 2)` (prior mean ≈ 0.96) |
| Overall time trend | Bayesian logistic regression, `logit(π) = β₀ + β₁·(Year − 2015)`, vague Normal priors |
| Region-specific trends | Separate intercept and slope per region, `logit(πᵢⱼ) = β₀ᵢ + β₁ᵢ·(Yearⱼ − 2015)` |
| Sampling | JAGS via `rjags`, 3 chains, 7,000 burn-in + 10,000 iterations per chain |
| Convergence | Trace and density plots, autocorrelation plots, Gelman–Rubin diagnostic (PSRF) |
| Model checking | Posterior predictive intervals for 2020 counts (model fitted on 2011–2019) |
| Regional comparison | Posterior distribution of coverage ratios between regions in 2019 |

Year is centred at 2015 to make the intercept interpretable and to improve MCMC mixing.

# Key findings
1. Prior choice barely matters. With 180–330 children per cell, posterior means under the flat and informative priors differ by at most about two percentage points. Coverage sits around 87–94% throughout.
2. No overall trend. In the pooled model the slope is essentially zero (β₁ = 0.001, 95% CrI [−0.029, 0.032]), with coverage steady at about 90%.
3. Regions tell different stories. Allowing region-specific trends reveals diverging directions:

| Region | Slope β₁ (mean) | 95% CrI | P(increasing) |
|---|---|---|---|
| Georgia | −0.018 | [−0.073, 0.036] | 0.26 |
| Wisconsin | −0.034 | [−0.088, 0.020] | 0.11 |
| Mississippi | 0.050 | [−0.002, 0.103] | **0.97** |
   
4. 2020 was in line with expectations. Observed 2020 counts in all three regions fall inside the 95% posterior predictive intervals from the 2011–2019 model.
5. 2019 comparison. Mississippi most likely had higher coverage than Wisconsin in 2019 (P(Wisconsin > Mississippi) ≈ 0.02), while the other differences are small and uncertain.

# How to run
Requires R and JAGS installed on your system.

```r
install.packages(c("rjags", "runjags", "coda", "ggplot2", "dplyr", "mcmcplots"))
source("codes.R")
```

# Tools
R, JAGS / BUGS, rjags, runjags, coda, ggplot2, dplyr
