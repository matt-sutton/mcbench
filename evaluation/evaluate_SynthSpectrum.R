source("benchmarks/synthSpectrum.R")
source("sampler/HMC.R")
source("evaluation/ess_metrics.R")

# List containing samples from the prior, plus functions for log_p and grad_log_p
# (Don't worry about the 'iter' error)
spectrum_info <- get_spectrum()

## Neal's HMC kernel ----
## (approx. 5 min. for 50k iterations)
set.seed(1)
x0 <- spectrum_info$samples(1)
system.time(samples <- HMC_sampler(5e4, x0, 
                                   spectrum_info$log_p, spectrum_info$grad_log_p, 15, 0.05))

## ESS metrics

# Multivariate ESS
get_multiESS(samples) # HMC

# ESS of log posterior
get_logpESS(samples,spectrum_info$log_p)

# RStan ----
library(rstan)
mod <- stan_model("benchmarks/spectraModel.stan", model_name="RamanSpectrum")
dat = list(nWL=length(wvNum), nPK=3, y=y, x=wvNum)
system.time(samp <- sampling(mod, data=dat, iter=5e4))
summary(samp, pars=c("loc","amp","sca","sigma"))
traceplot(samp, pars=c("loc")) # multi-modality
