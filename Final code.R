
library(rjags)
library(runjags)
library(coda)
library(ggplot2)
library(readxl)
library(coda)
library(dplyr)
library(mcmcplots)

#     Question 1 

data <- data.frame(
  Region = rep(c("Georgia", "Wisconsin", "Mississippi"), each = 10),
  Year = rep(2011:2020, times = 3),
  Vaccinated = c(
    196,248,261,252,276,311,265,246,251,165,     # Georgia
    207,205,212,195,231,246,215,214,197,156,     # Wisconsin
    171,208,190,215,243,276,290,242,304,161      # Mississippi
  ),
  Sample_Size = c(
    219,270,276,284,306,334,292,282,273,188,     # Georgia
    225,226,235,224,262,275,238,241,224,177,     # Wisconsin
    198,230,217,239,272,307,321,276,324,181      # Mississippi
  )
)
data






# this for flat prior with alpha =1 and beta=1 

flat_results <- data %>%
  mutate(
    post_alpha   = 1 + Vaccinated,
    post_beta    = 1 + (Sample_Size - Vaccinated),
    mean_flat    = (post_alpha / (post_alpha + post_beta)),
    mean_flat2 = round(mean_flat, 2),
    var_flat = (post_alpha*post_beta)/ (((post_alpha+post_beta)^2)*(post_alpha+post_beta+1)),
    ci_flat_low  = qbeta(0.025, post_alpha, post_beta),
    ci_flat_high = qbeta(0.975, post_alpha, post_beta)
  )
# informative prior when alpha =90 betat =10

inf_results <- data %>%
  mutate(
    post_alpha   = 50 + Vaccinated,
    post_beta    = 2 + (Sample_Size - Vaccinated),
    mean_inf     = post_alpha / (post_alpha + post_beta),
    mean_inf_2 = round(mean_inf, 2),
    var_inform = (post_alpha*post_beta)/ (((post_alpha+post_beta)^2)*(post_alpha+post_beta+1)),
    ci_inf_low   = qbeta(0.025, post_alpha, post_beta),
    ci_inf_high  = qbeta(0.975, post_alpha, post_beta)
  )
# here we are making table to see or compare mean and CI for flat and informative prior
all_results <- flat_results %>%
  select(Region,Year,Vaccinated,Sample_Size, mean_flat2,var_flat, ci_flat_low, ci_flat_high) %>%
  left_join(
    inf_results %>%
      select(Region,Year, mean_inf_2, var_inform, ci_inf_low, ci_inf_high),
    by = c("Region","Year")
  )
View(all_results)
print(all_results)
#---------------------------------------------------------------------
# question 2

# filtering years making it from 2011 to 2019
data_q2 <- data %>% 
  filter(Year >= 2011, Year <= 2019)


dataList <- list(
  Y    = data_q2$Vaccinated,
  N    = data_q2$Sample_Size,
  Year = data_q2$Year,
  Nobs = nrow(data_q2)
)

#model.inits <- list(
# list(beta0 = -1, beta1 = 0.5),
#list(beta0 = 0, beta1 = 0),
#list(beta0=1 , beta1 = -0.5)

model.inits <- list(
  list(beta0 = -1, beta1 = 0, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1001),
  list(beta0 = 0,   beta1 = 0,   .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1002),
  list(beta0 = 1,   beta1 = 0,.RNG.name = "base::Wichmann-Hill", .RNG.seed = 1003)
)
#model.inits <- list(beta0=0,beta1=0)
#file.show("C:\\Users\\basma lusta\\Desktop\\Bayesian\\the project\\BUGS.txt")
#In practice, we’ll often run 3 or 4 chains in parallel. we discard the first 
#bit of each chain (burn-in) and then use the remaining samples from all chains
#together to estimate means, intervals, and check diagnostics.

modelString <- "
model {
  for (n in 1:Nobs) {
    # 1) Likelihood
    Y[n] ~ dbin(pi[n], N[n])
    # 2) Logit link
    logit(pi[n]) <- beta0 + beta1 * (Year[n] - 2015)
  }

  # Priors
  beta0 ~ dnorm(0, 1.0E-6)
  beta1 ~ dnorm(0, 1.0E-6)
}
"

#---- How we ran the chains
#txt file: contains the model in bugs language and we run it via jags
jags <- jags.model(textConnection(modelString),
                   data = dataList,
                   inits = model.inits,
                   n.chains = 3)

#model {
  # Loop over all region–year cells (n = 1,…,Nobs)
 # for (n in 1:Nobs) {
  #  # 1) Likelihood: number vaccinated ~ Binomial(total, coverage)
   # Y[n]    ~ dbin(pi[n], N[n])
    ## 2) Linear predictor on the logit scale
    #logit(pi[n]) <- beta0 + beta1 *  (Year[n] - 2015)
    
#  }
  
  # Flat (nearly non‐informative) priors on the regression coefficients:
 # beta0 ~ dnorm(0, 1.0E-6)   # Normal(0, variance=1e6)
#  beta1 ~ dnorm(0, 1.0E-6)   # likewise
#}


#-------------------------------------------------------------
# (question 3 )
# run  or generate MCMC samples
# the burn-in

update(jags, n.iter=7000)
# Draw samples
samples <- coda.samples(
  jags,
  variable.names=c("beta0","beta1"),
  n.iter=10000, thin=1
)
#----------------------------------------------------

# History plot & posterior distributions & autocorrelation plot
plot(samples, trace=TRUE, density = TRUE)   
  

sample.mcmc <- as.mcmc.list(samples)


# Checking correlations
autocorr.plot(sample.mcmc)
autocorr.plot(samples)

# Convergence tests
gelman.diag(sample.mcmc)
gelman.plot(sample.mcmc,ask=FALSE)
#----------------------------------------------------------
# Question 4
# Combine multiple chains to get overall HPD interval
# WE CAN HAVE POSTERIOR SUMMARY AND DENSITY PLOT FROM HERE

library(runjags)
model.mcmc.combined<-combine.mcmc(samples)
View(model.mcmc.combined)

HPDinterval(model.mcmc.combined)
summary(model.mcmc.combined)
plot(model.mcmc.combined, trace=TRUE, density = TRUE)  
plot(samples, trace=TRUE, density = TRUE)  


#Question5
library(dplyr)

# 1) Compute posterior mean coefficients
#    We take the mean of the MCMC draws for beta0 and beta1
#  it is averaging over all 30 000 posterior draws of each coefficient to 
#give you the posterior mean of 𝛽0 𝛽1 ,
# this line is like we are taking 30000 B0 and 30000 B1. then we are taking
# the average of 30000 B0 and 30000 B1

# put all results for 300000 b0 and 30000 b1  from 3 chains in one matrix to esase handling 
# then take the average value for b0 column, and b1 column
post_coefs <- colMeans(as.matrix(samples)[, c("beta0","beta1")])
# post_coefs["beta0"], post_coefs["beta1"]
print(post_coefs)

#  Compute the fitted coverage (mean only) for each year 2011–2019
years <- 2011:2019

# plogis(x) computes inverse logit = e^x \1+e^x to make the transformation 
# to have it in term pf probablity 
# This uses invlogit(mean beta0 + mean beta1 * centered year)
# we will have 9 means , 1 for each year

logistic_means <- plogis(
  # ex : 2.204 + 0.0012(2017-2015)= 2.2064
  
  # transformation inverse log to:  e^2.2064 \ 1+ e^2.2064= 0.900822
  post_coefs["beta0"] + post_coefs["beta1"] * (years - 2015)
)
# 5his will give 9 estimate one for each year 
print(logistic_means)
# This uses invlogit(mean beta0 + mean beta1 * centered year)

# Compute analytical Beta‐posterior means per year
#    We average the ‐specific means_flat, mean_info from Q1 for each birth year

analytical_means_flat <- flat_results %>%
  filter(Year %in% years) %>%          # only 2011–2019
  # we are grouping by year, to take the mean of the 3 cities in one year
  # ( this like averging the 3 cities that belong to same year)
  # because the question said per year only and not per  year and region
  group_by(Year) %>%
  summarize(analytic_mean = mean(mean_flat2)) %>%
  # this mean extract only analytic mean from the flat_result table and group it by year
  pull(analytic_mean)

analytical_means_info <- inf_results %>%
  filter(Year %in% years) %>%          # only 2011–2019
  # we are grouping by year, to take the mean of the 3 cities in one year
  # ( this like averging the 3 cities that belong to same year)
  # because the question said per year only and not per  year and region
  group_by(Year) %>%
  summarize(analytic_mean = mean(mean_inf_2)) %>%
  # this mean extract only analytic mean from the informative_result table and group it by year
  pull(analytic_mean)

#  Combine into a single table
# we need to combine the results now 
compare_df <- data.frame(
  Year              = years,
  Logistic_Mean     = round(logistic_means,    4),
  Analytical_Mean_Flat   = round(analytical_means_flat,  4),
  Analytical_Mean_info   = round(analytical_means_info,  4)
)

# --- 5) Print
print(compare_df)


#Question 6
data_q2 <- data %>% 
  filter(Year >= 2011, Year <= 2019)


dataList2 <- list(
  Y    = data_q2$Vaccinated,
  N    = data_q2$Sample_Size,
  Year = data_q2$Year,
  Region = as.numeric(factor(data_q2$Region, levels = c("Georgia", "Wisconsin", "Mississippi"))),
  Nobs = nrow(data_q2),
  Nregion = length(unique(data_q2$Region))
)


model.inits2 <- function() {
  list(
    beta0 = rnorm(3, 0, 0.1), 
    beta1 = rnorm(3, 0, 0.1)    
  )
}

modelString2 <- "
model {
  for (n in 1:Nobs) {
    # 1) Likelihood
    Y[n] ~ dbin(p[n], N[n])
    # 2) Logit link
    logit(p[n]) <- beta0[Region[n]] + beta1[Region[n]] * (Year[n] - 2015)
  }

  # Priors
  for (r in 1:Nregion) {
    beta0[r] ~ dnorm(0, 1.0E-6)
    beta1[r] ~ dnorm(0, 1.0E-6)
  }
}
"
set.seed(12345)
jags2 <- jags.model(textConnection(modelString2),
                   data = dataList2,
                   inits = model.inits2,
                   n.chains = 3)

update(jags2, n.iter=7000)

samples2 <- coda.samples(
  jags2,
  variable.names=c("beta0","beta1"),
  n.iter=10000, thin=1
)

#plot(samples2, trace = TRUE, xlim = c(0, 7000))
plot(samples2, trace=TRUE, density = TRUE)   

sample2.mcmc <- as.mcmc.list(samples2)
par(mfrow = c(1, 2))

for (i in 1:3) {
  # Construct parameter names like beta0[1], beta1[1], etc.
  param1 <- paste0("beta0[", i, "]")
  param2 <- paste0("beta1[", i, "]")
  
  sub <-sample2.mcmc[, c(param1, param2), drop=FALSE]
  
  gelman.plot(sub)
  
  #readline("Press [Enter] to see next pair...")
}

gelman.diag(sample2.mcmc)
gelman.plot(sample2.mcmc,ask=FALSE)
autocorr.plot(sample2.mcmc)

model.mcmc.combined <- combine.mcmc(samples2)
HPDinterval(model.mcmc.combined)
summary(sample2.mcmc)


#Question 7
posterior_samples <- as.matrix(samples2)
# For each location, compute P(beta1[i] > 0)
prob_increase <- c(
  Georgia = mean(posterior_samples[, "beta1[1]"] > 0),
  Wisconsin = mean(posterior_samples[, "beta1[2]"] > 0),
  Mississippi = mean(posterior_samples[, "beta1[3]"] > 0)
)
round(prob_increase, 4)

# Question 8
observed_coverage <- data %>%
  filter(Year >= 2011 & Year <= 2019) %>%
  mutate(
    regions = as.numeric(factor(Region, levels = c("Georgia", "Wisconsin", "Mississippi"))),
    Observed = Vaccinated / Sample_Size
  )

years <- unique(observed_coverage$Year)

coverage_preds <- do.call(rbind, lapply(1:3, function(i) {
  sapply(years, function(y) {
    logit <- posterior_samples[, paste0("beta0[", i, "]")] + posterior_samples[, paste0("beta1[", i, "]")] * (y - 2015)
    p <- 1 / (1 + exp(-logit))
    c(mean = mean(p), lower = quantile(p, 0.025), upper = quantile(p, 0.975))
  }) %>%
    t() %>%
    as.data.frame() %>%
    mutate(Year = years, Region = c("Georgia", "Wisconsin", "Mississippi")[i])
}))

ggplot(coverage_preds, aes(x = Year, y = mean, color = Region, fill = Region)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = `lower.2.5%`, ymax = `upper.97.5%`), alpha = 0.2, color = NA) +
  geom_point(data = observed_coverage, aes(x = Year, y = Observed, color = Region), size = 2, shape = 21) +
  facet_wrap(~Region) +
  scale_x_continuous(breaks = 2011:2019) + 
  labs(y = "Vaccination Coverage Proportion", x = "Year") +
  theme_minimal()+ theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12), 
    strip.text.x = element_text(size = 11, face = "bold")
  )

# Question 9
# Define 2020 sample sizes for each region
sample_sizes_2020 <- c(Georgia = 188, Wisconsin = 177, Mississippi = 181)
observed_2020 <- c(Georgia = 165, Wisconsin = 156, Mississippi = 161)

# Region names
regions <- c("Georgia", "Wisconsin", "Mississippi")

# Run posterior predictive simulation for each region
predictive_results <- lapply(1:3, function(i) {
  region_name <- regions[i]
  size_2020 <- sample_sizes_2020[region_name]
  obs_2020 <- observed_2020[region_name]
  
  # Get posterior predictive probabilities
  logit_p <- posterior_samples[, paste0("beta0[", i, "]")] + 
    posterior_samples[, paste0("beta1[", i, "]")] * (2020 - 2015)
  p <- 1 / (1 + exp(-logit_p))
  
  # Simulate vaccinated counts
  y_pred <- rbinom(length(p), size = size_2020, prob = p)
  
  # Summary statistics
  data.frame(
    Region = region_name,
    Observed = obs_2020,
    Mean = median(y_pred),
    Lower = quantile(y_pred, 0.025),
    Upper = quantile(y_pred, 0.975),
    stringsAsFactors = FALSE
  )
})

# Combine into a data frame
pred_summary <- do.call(rbind, predictive_results)
print(pred_summary)



# Questin 10
invlogit <- function(x) 1 / (1 + exp(-x))
year_centered <- 2019 - 2015

# Compute posterior probabilities for 2019
p_G <- invlogit(posterior_samples[, "beta0[1]"] + (year_centered * posterior_samples[, "beta1[1]"]))  # Georgia
p_W <- invlogit(posterior_samples[, "beta0[2]"] + (year_centered * posterior_samples[, "beta1[2]"]))  # Mississippi
p_M <- invlogit(posterior_samples[, "beta0[3]"] +(year_centered * posterior_samples[, "beta1[3]"]))  # Wisconsin

# Ratios
Georgia_Wisconsin <- p_G / p_W
Georgia_Mississippi <- p_G / p_M
Wisconsin_Mississippi <- p_W / p_M

# Summarise
summary_stats <- function(x) {
  c(mean = mean(x), 
    lower = quantile(x, 0.025), 
    upper = quantile(x, 0.975))
}

results <- rbind(
  GA_vs_WI = summary_stats(Georgia_Wisconsin),
  GA_vs_MS = summary_stats(Georgia_Mississippi),
  WI_vs_MS = summary_stats(Wisconsin_Mississippi)
)

round(results, 3)


# Compute posterior probabilities that the ratio > 1
Pr_ratio_gt_1 <- c(
  GA_vs_WI = mean(Georgia_Wisconsin > 1),
  GA_vs_MS = mean(Georgia_Mississippi > 1),
  WI_vs_MS = mean(Wisconsin_Mississippi > 1)
)

# Display rounded results
round(Pr_ratio_gt_1, 2)

