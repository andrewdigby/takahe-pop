
# Produced by ChatGPT
# ---------------------------
# Install and load packages
# ---------------------------
install.packages("nimble")
library(nimble)

# ---------------------------
# Step 1: Simulate example data
# ---------------------------
set.seed(123)

# Parameters
n_years <- 5
n_adults <- 20
n_nests <- 8
true_S <- rep(0.85, n_years-1)   # adult survival
true_J <- rep(0.5, n_years-1)    # juvenile survival
true_lambda <- rep(2, n_years-1) # chicks per nest
true_N <- numeric(n_years)
true_N[1] <- 50

# Population process
for(t in 1:(n_years-1)) {
  R <- true_lambda[t] * n_nests * true_J[t]
  true_N[t+1] <- rpois(1, true_N[t]*true_S[t] + R)
}

# Telemetry survival data
adult_survival <- matrix(NA, nrow=n_adults, ncol=n_years-1)
for(t in 1:(n_years-1)) {
  adult_survival[,t] <- rbinom(n_adults, 1, true_S[t])
}

# Nest productivity and yearling data
chicks <- numeric(n_nests*(n_years-1))
yearlings <- numeric(n_nests*(n_years-1))
nest_year_index <- numeric(n_nests*(n_years-1))
idx <- 1
for(t in 1:(n_years-1)) {
  for(i in 1:n_nests) {
    nest_year_index[idx] <- t
    chicks[idx] <- rpois(1, true_lambda[t])
    yearlings[idx] <- rbinom(1, chicks[idx], true_J[t])
    idx <- idx + 1
  }
}

# Counts
count <- true_N

# ---------------------------
# Step 2: Define IPM in NIMBLE
# ---------------------------
ipm_code <- nimbleCode({
  for(t in 1:(n_years-1)) {
    S[t] ~ dbeta(2, 2)          # Adult survival
    J[t] ~ dbeta(2, 2)          # Juvenile survival
    lambda[t] ~ dnorm(2, 1)     # Mean chicks per nest
  }
  
  N[1] ~ dpois(init_N)
  for(t in 1:(n_years-1)) {
    R[t] <- lambda[t] * nests[t] * J[t]
    N[t+1] ~ dpois(N[t]*S[t] + R[t])
  }
  
  # Adult survival likelihood
  for(i in 1:n_adults) {
    for(t in 1:(n_years-1)) {
      y[i,t] ~ dbern(S[t])
    }
  }
  
  # Nest productivity and juvenile survival
  for(i in 1:n_nests_total) {
    chicks[i] ~ dpois(lambda[year_nest[i]])
    yearlings[i] ~ dbinom(J[year_nest[i]], chicks[i])
  }
  
  # Count data
  for(t in 1:n_years) {
    count[t] ~ dpois(N[t])
  }
})

# ---------------------------
# Step 3: Prepare data and constants
# ---------------------------
constants <- list(
  n_years = n_years,
  n_adults = n_adults,
  n_nests_total = length(chicks),
  nests = rep(n_nests, n_years-1),
  year_nest = nest_year_index
)

data <- list(
  y = adult_survival,
  chicks = chicks,
  yearlings = yearlings,
  count = count,
  init_N = count[1]
)

inits <- list(
  S = rep(0.8, n_years-1),
  J = rep(0.5, n_years-1),
  lambda = rep(2, n_years-1),
  N = count
)

# ---------------------------
# Step 4: Build and run the model
# ---------------------------
ipm_model <- nimbleModel(code = ipm_code, data = data, constants = constants, inits = inits)
ipm_conf <- configureMCMC(ipm_model, monitors = c("S", "J", "lambda", "N"))
ipm_mcmc <- buildMCMC(ipm_conf)

Cmodel <- compileNimble(ipm_model)
Cmcmc <- compileNimble(ipm_mcmc, project = ipm_model)
samples <- runMCMC(Cmcmc, niter = 10000, nburnin = 2000, thin = 10)

# ---------------------------
# Step 5: Summarize results
# ---------------------------
print(summary(samples))

# plot

# Load required libraries
library(coda)
library(ggplot2)

# Convert samples to mcmc object
mcmc_samples <- as.mcmc(samples)

# Quick summary
summary(mcmc_samples)

# Trace plots for key parameters
par(mfrow=c(2,2))
traceplot(mcmc_samples[,grep('S', colnames(mcmc_samples))], main='Traceplot: Adult Survival (S)')
traceplot(mcmc_samples[,grep('J', colnames(mcmc_samples))], main='Traceplot: Juvenile Survival (J)')
traceplot(mcmc_samples[,grep('lambda', colnames(mcmc_samples))], main='Traceplot: Chicks per Nest (lambda)')
traceplot(mcmc_samples[,grep('N', colnames(mcmc_samples))], main='Traceplot: Population Size (N)')

# Density plots using ggplot2
posterior_df <- as.data.frame(mcmc_samples)

# Example: Plot posterior for S[1]
ggplot(posterior_df, aes(x=`S[1]`)) +
  geom_density(fill='skyblue', alpha=0.5) +
  labs(title='Posterior Distribution: Adult Survival Year 1', x='S[1]', y='Density')

# Example: Plot posterior for J[1]
ggplot(posterior_df, aes(x=`J[1]`)) +
  geom_density(fill='orange', alpha=0.5) +
  labs(title='Posterior Distribution: Juvenile Survival Year 1', x='J[1]', y='Density')

# Example: Plot posterior for lambda[1]
ggplot(posterior_df, aes(x=`lambda[1]`)) +
  geom_density(fill='green', alpha=0.5) +
  labs(title='Posterior Distribution: Chicks per Nest Year 1', x='lambda[1]', y='Density')
