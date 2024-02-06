library(rstan)
load("benchmarks/FactAnal.rda")

# Change the value of k to change the model
stan_FA <<- stan(file = "benchmarks/FactAnal.stan",
                     data = list(k = 3, n = n, d = d, y = y),
                     warmup = 1000, iter = 10000, chains = 1, 
                     verbose = FALSE) 

# Default plot
plot(stan_FA)

# Basic plots of some of the challenging marginal(s) and
# bivariate distributions

par(mfrow=c(2,2))
samples <- stan_FA@sim$samples[[1]]
plot(density(samples$`beta_offdiag[9]`))
plot(samples$`beta_offdiag[9]`,
     samples$`beta_offdiag[12]`)
plot(samples$`beta_offdiag[9]`,
     samples$`beta_offdiag[11]`)
plot(samples$`beta_offdiag[6]`,
     samples$`beta_offdiag[11]`)
