# Bayesian-Analysis-of-Varicella-Vaccination-Coverage
A Bayesian study of childhood varicella (chickenpox) vaccination coverage in three US regions (Georgia, Wisconsin and Mississippi) across birth cohorts 2011–2020. The project moves from closed-form conjugate analysis to MCMC-based logistic regression, and uses the posterior to answer practical questions: Is coverage changing over time? Do regions differ? Was 2020 unusual?

[Read the documentation →](https://maryamshojaeiee-hub.github.io/DNA-Sequence-Assembly-with-de-Bruijn-Graphs/documentation.pdf) · [View the code →](https://maryamshojaeiee-hub.github.io/DNA-Sequence-Assembly-with-de-Bruijn-Graphs/main_codes.py) · [View the tests →](https://maryamshojaeiee-hub.github.io/DNA-Sequence-Assembly-with-de-Bruijn-Graphs/test.py) 

Estimated vaccination coverage by region and birth year

# Data
Survey counts of children who received at least one dose of the varicella vaccine, by region and birth year (2011–2020): number vaccinated (`Vaccinated`) out of the number surveyed (`Sample_Size`). Each region–year cell has roughly 180–330 children. The data are included directly in the R script.

# Methods

Year is centred at 2015 to make the intercept interpretable and to improve MCMC mixing.

# Key findings

# How to run
Requires R and JAGS installed on your system.

```r
install.packages(c("rjags", "runjags", "coda", "ggplot2", "dplyr", "mcmcplots"))
source("codes.R")
```

# Tools
R, JAGS / BUGS, rjags, runjags, coda, ggplot2, dplyr
