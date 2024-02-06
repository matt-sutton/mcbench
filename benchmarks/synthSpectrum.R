# simulate some spectroscopic data
set.seed(123)
wvNum <- seq(100,800,by=5)
peak1 <- 15*pi*dcauchy(wvNum, loc=250, scale=15)*100
peak2 <- 25*sqrt(2*pi)*dnorm(wvNum, mean=350, sd=25)*50
peak3 <- 20*pi*dcauchy(wvNum, loc=450, scale=20)*125

trueSD = 15
y <- peak1 + peak2 + peak3 + rnorm(length(wvNum), sd=trueSD)
nPK <- 3
nWL <- length(wvNum)

sample_priors <- function(nsamples = 10000) {
  n <- nPK*3 + 1
  X <- matrix(0, nsamples, n)
  for (i in 1:nsamples) {
    X[i,1:nPK] <- sort(runif(nPK, min(wvNum), max(wvNum))) # peak locations
    X[i,nPK + 1:nPK] <- runif(nPK, 0, max(y))              # amplitudes
    X[i,2*nPK + 1:nPK] <- rlnorm(nPK, meanlog=log(11.6) - (0.4^2)/2, sdlog=0.4)
    X[i,n] <- 1/sqrt(rgamma(1, shape=5, rate=1125))
  }
  return(X)
}

get_spectrum <- function() {
  dat = list(nWL=length(wvNum), nPK=3, y=y, x=wvNum)
  require("rstan")
  stan_fit <<- stan("benchmarks/spectraModel.stan", model_name="RamanSpectrum",
                    data = dat, warmup = 0, iter = 0, chains = 1, 
                    verbose = FALSE) 
  f <- function(x){log_prob(stan_fit, x)}
  g <- function(x){grad_log_prob(stan_fit, x)}
  sample_fun <- function(nsamp=10000){
    return(sample_priors(nsamples = nsamp))
  }
  return(list(
    log_p = f,
    grad_log_p = g,
    samples = sample_fun
  ))
}
