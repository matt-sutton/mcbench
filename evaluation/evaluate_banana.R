source("benchmarks/banana.R")
source("sampler/HMC.R")
source("sampler/RW.R")
source("evaluation/ess_metrics.R")

## Evaluate using Random walk (NOT using stan):
banana_info <- get_banana(mu = 1, a = 1, b = matrix(5,2,2), use_stan = F)

set.seed(2)
x0 <- banana_info$samples(nsamp = 1)
exact_samples <- banana_info$samples(1000)

## Evaluate sampling using RW
rw_samples <- RW_sampler(1500, x0, banana_info$log_p,eps = 0.09)
par(mfrow = c(length(x0),length(x0)), mar = c(1,1,1,1))
for( i in 1:length(x0)){
  for(j in 1:length(x0)){
    if( i < j){
      plot(exact_samples[,i], exact_samples[,j], pch=20,col = 1)
      points(rw_samples[,i], rw_samples[,j], pch=20,col = 2)
    } else{
      plot.new()
    }
  }
}
legend("bottomleft", legend = c("Exact", "RW"), col = c(1,2), lwd=2)
par(mfrow = c(1,1), mar =rep(4,4))


# List containing exact samples, plus functions for log_p and grad_log_p
# (Don't worry about the 'iter' error)
banana_info <- get_banana(mu = 1, a = 1, b = matrix(5,2,2))

## Evaluate sampling using HMC
hmc_samples <- HMC_sampler(1000, x0, banana_info$log_p, banana_info$grad_log_p, 15, 0.05)

## Plotting pairs of the distribution
par(mfrow = c(length(x0),length(x0)), mar = c(1,1,1,1))
for( i in 1:length(x0)){
  for(j in 1:length(x0)){
    if( i < j){
      plot(exact_samples[,i], exact_samples[,j], pch=20,col = 1)
      points(hmc_samples[,i], hmc_samples[,j], pch=20,col = 2)
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
      points(hmc_samples[,i], hmc_samples[,j], pch=20,col = 2)
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
get_logpESS(t(bps_samples$samples),banana_info$log_p)
get_logpESS(hmc_samples,banana_info$log_p)
get_logpESS(rw_samples,banana_info$log_p)


