# Bayesian-Analysis-of-Varicella-Vaccination-Coverage
A Bayesian analysis of childhood varicella (chickenpox) vaccination coverage in three US regions (Georgia, Wisconsin and Mississippi) for birth cohorts 2011–2020. The project goes from closed-form conjugate analysis to Bayesian logistic regression fitted with MCMC in JAGS, and uses the posterior to answer practical questions: Is coverage changing over time? Do the regions follow different trends? Was the 2020 cohort unusual?

This was a group project for the Bayesian Inference course at Hasselt University (2024–2025). The report follows the structure of the assignment and answers ten questions in order.

[Read the full report →](https://github.com/maryamshojaeiee-hub/Bayesian-Analysis-of-Varicella-Vaccination-Coverage/blob/main/Report_Project%20Assignment%20Bayesian%20Inference.pdf) · [View the code →](https://github.com/maryamshojaeiee-hub/Bayesian-Analysis-of-Varicella-Vaccination-Coverage/blob/main/Final%20code.R)

Estimated vaccination coverage by region and birth year

# Data
Survey counts of children who received at least one dose of the varicella vaccine, by region and birth year (2011–2020): number vaccinated (`Vaccinated`) out of the number surveyed (`Sample_Size`). There are 30 region–year observations, each with 177–334 children. The data are included directly in the R script.

# What the report covers

| Question | Topic | Method |
|---|---|---|
| Q1 | Coverage per region and birth year, and the effect of the prior | Beta–Binomial conjugate model; flat `Beta(1, 1)` vs. informative `Beta(50, 2)` prior |
| Q2 | Change in coverage over birth years 2011–2019 | Bayesian logistic regression, `logit(π) = β₀ + β₁·(Year − 2015)`, written in BUGS |
| Q3 | MCMC sampling and convergence | JAGS via `rjags`, 3 chains, 7,000 burn-in + 10,000 iterations; trace, autocorrelation and Gelman–Rubin diagnostics |
| Q4 | Posterior densities and summaries of β₀ and β₁ | Posterior means, SDs and 95% credible intervals |
| Q5 | Coverage per birth year from the regression | Comparison with the analytical estimates from Q1 |
| Q6 | Region-specific trends | Separate intercept and slope per region, `logit(πᵢⱼ) = β₀ᵢ + β₁ᵢ·(Yearⱼ − 2015)` |
| Q7 | Probability of increasing coverage per region | Posterior probability P(β₁ᵢ > 0) |
| Q8 | Estimated coverage with uncertainty vs. observed data | Posterior means and 95% credible bands per region and year |
| Q9 | Is the 2020 cohort in line with earlier years? | Posterior predictive intervals for 2020 counts |
| Q10 | Pairwise comparison of regions in 2019 | Posterior distribution of coverage ratios and P(ratio > 1) |

Year is centred at 2015 to make the intercept interpretable and to improve MCMC mixing.

# Key findings
1. Prior choice barely matters. With 180–330 children per cell, posterior means under the flat and informative priors differ by at most about two percentage points. Coverage sits around 87–94% throughout.
2. No overall trend. In the pooled model the slope is essentially zero (β₁ = 0.001, 95% CrI [−0.029, 0.032]), with coverage steady at about 90%.
3. The regions follow different trends. With a separate trend per region:

   | Region | Slope β₁ (mean) | 95% CrI | P(increasing) |
   |---|---|---|---|
   | Georgia | −0.018 | [−0.073, 0.036] | 0.26 |
   | Wisconsin | −0.034 | [−0.088, 0.020] | 0.11 |
   | Mississippi | 0.050 | [−0.002, 0.103] | **0.97** |

   Opposite trends cancel out in the pooled model, which is why it shows no change.
   
4. 2020 was in line with expectations. Observed 2020 counts in all three regions fall inside the 95% posterior predictive intervals from the 2011–2019 model.
5. Regional comparison in 2019. We compared 2019 coverage between each pair of regions using the ratio of their coverage (a ratio above 1 means the first region had higher coverage):

- **Regional comparison in 2019.** We compared 2019 coverage between each pair of regions using the ratio of their coverage (a ratio above 1 means the first region had higher coverage):

  | Comparison | Ratio (mean) | 95% CrI | P(ratio > 1) |
  |---|---|---|---|
  | Georgia / Wisconsin | 1.026 | [0.988, 1.067] | 0.91 |
  | Georgia / Mississippi | 0.989 | [0.958, 1.021] | 0.25 |
  | Wisconsin / Mississippi | 0.964 | [0.929, 0.999] | 0.02 |

  All ratios are close to 1, so coverage differences between regions were small (a few percent at most). The clearest difference is Wisconsin vs. Mississippi: Wisconsin's coverage was about 3.6% lower, its credible interval lies just below 1, and there is only a 2% probability that Wisconsin had higher coverage. Georgia likely had slightly higher coverage than Wisconsin (probability 0.91), while Georgia and Mississippi cannot be clearly separated.

# How to run
Requires R and JAGS installed on your system.

```r
install.packages(c("rjags", "runjags", "coda", "ggplot2", "dplyr", "mcmcplots"))
source("codes.R")
```

# Tools
R, JAGS / BUGS, rjags, runjags, coda, ggplot2, dplyr

Group project, Bayesian Inference course, Hasselt University (2024–2025)
