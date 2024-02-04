source("benchmarks/banana.R")
source("sampler/HMC.R")
source("evaluation/ess_metrics.R")

# List containing exact samples, plus functions for log_p and grad_log_p
# (Don't worry about the 'iter' error)
banana_info <- get_banana(mu = 1, a = 1, b = matrix(5,2,2))

## Evaluate sampling using HMC
set.seed(1)
x0 <- banana_info$samples(1)
samples <- HMC_sampler(1000, x0, banana_info$log_p, banana_info$grad_log_p, 15, 0.05)
exact_samples <- banana_info$samples(1000)

## Plotting pairs of the distribution
par(mfrow = c(length(x0),length(x0)), mar = c(1,1,1,1))
for( i in 1:length(x0)){
  for(j in 1:length(x0)){
    if( i < j){
      plot(exact_samples[,i], exact_samples[,j], pch=20,col = 1)
      points(samples[,i], samples[,j], pch=20,col = 2)
    } else{
      plot.new()
    }
  }
}
legend("bottomleft", legend = c("Exact", "HMC"), col = c(1,2), lwd=2)
par(mfrow = c(1,1), mar =rep(4,4))

## Evaluate using the Bouncy Particle Sampler
library(ccpdmp)
## 15000 = num grad eval from HMC
b <- bps(15000, dnlogpi = function(x, i){return(-banana_info$grad_log_p(x)[i])}, x0 = x0, poly_order = 3)
bps_pos <- t(b$positions)

par(mfrow = c(length(x0),length(x0)), mar = c(1,1,1,1))
for( i in 1:length(x0)){
  for(j in 1:length(x0)){
    if( i < j){
      plot(exact_samples[,i], exact_samples[,j], pch=20,col = 1)
      points(bps_pos[,i], bps_pos[,j], pch=20,col = 3)
      points(samples[,i], samples[,j], pch=20,col = 2)
    } else{
      plot.new()
    }
  }
}
legend("bottomleft", legend = c("Exact", "HMC", "BPS"), col = c(1,2,3), lwd=2)
par(mfrow = c(1,1), mar = rep(4,4))

bps_samples <- gen_samples(nsample = 1000,
                           b$positions, b$times)

## ESS metrics

# Multivariate ESS
multiESS(samples) # HMC
multiESS(t(bps_samples$samples)) # BPS

# ESS of log posterior
get_logpESS(t(bps_samples$samples),banana_info$log_p)
get_logpESS(samples,banana_info$log_p)


