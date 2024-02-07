source("benchmarks/banana.R")
source("sampler/HMC.R")
source("sampler/RW.R")
source("evaluation/ess_metrics.R")

banana_info <- get_banana(mu = 1, a = 1, b = matrix(5,2,2))

set.seed(2)
x0 <- banana_info$samples(nsamp = 1)
exact_samples <- banana_info$samples(1000)

banana_info$log_p(x0)
banana_info$grad_log_p(x0)
banana_info$samples(nsamples = 100)

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

stan_fit <<- stan("benchmarks/banana.stan", 
                  data = banana_info$stan_dat, warmup = 100, iter = 1000, chains = 1, 
                  verbose = FALSE) 
stan_samples <- extract(stan_fit)$X

## Evaluate sampling using HMC
hmc_samples <- HMC_sampler(1000, 1+0*x0, banana_info$log_p, banana_info$grad_log_p, 15, 0.05)

## Evaluate using the Bouncy Particle Sampler
library(ccpdmp)
## 15000 = num grad eval from HMC
b <- bps(15000, dnlogpi = function(x, i){return(-banana_info$grad_log_p(x)[i])}, x0 = x0, poly_order = 3)
bps_pos <- t(b$positions)

bps_samples <- t(gen_samples(nsample = 1000,
                           b$positions, b$times)$samples)

pchv = 20
par(mfrow = c(length(x0),length(x0)), mar = c(1,1,1,1))
for( i in 1:length(x0)){
  for(j in 1:length(x0)){
    if( i < j){
      plot(exact_samples[,i], exact_samples[,j], pch=pchv,col = 1)
      points(bps_samples[,i], bps_samples[,j], pch=pchv,col = 3)
      points(stan_samples[,i], stan_samples[,j], pch=pchv,col = 4)
      points(hmc_samples[,i], hmc_samples[,j], pch=pchv,col = 2)
      points(rw_samples[,i], rw_samples[,j], pch=pchv,col = 5)
    } else{
      plot.new()
    }
  }
}
legend("bottomleft", legend = c("Exact", "NUTS", "HMC", "BPS","RW"), col = c(1,4,2,3,5), lwd=2)
par(mfrow = c(1,1), mar = rep(4,4))


## ESS metrics
get_logpESS(t(bps_samples$samples),banana_info$log_p)
get_logpESS(hmc_samples,banana_info$log_p)
get_logpESS(rw_samples,banana_info$log_p)


