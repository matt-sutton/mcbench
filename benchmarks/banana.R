#########
## Hybrid Rosenbrock distribution (multidimensional banana)
## Pagani F, Wiegand M, Nadarajah S. An n-dimensional Rosenbrock distribution for 
## Markov chain Monte Carlo testing. Scand J Statist. 2022; 49: 657–680. https://doi.org/10.1111/sjos.12532
## 
## 
## Description:
## A challenging problem where curvature can be made arbitrarily hard to sample
## particularly challenging for sampling the tails. Admits exact sampling with 
## tractable normalising constant. 
##
## 
## This file returns get_banana a function that gives a list with:
##  -- log_p: the log probability of the distribution
##  -- grad_log_p: the grad of the log_p
##  -- samples: a function returning a set of sampled generated exactly
## 
## NOTE: A warning is given from stan by default as no sampling is occurring
##       See also evaluations for an example use of this problem.
#########

get_banana <- function(mu = 1, a = 1, b = matrix(5,2,2), use_stan = T){
  n2 = nrow(b); n1 = ncol(b)
  if(n1 < 2) stop("need more than 2 columns for b")
  dat = list(mu=mu, a=a,n1=n1,n2=n2,b = b)
  
  if(use_stan){
    require("rstan")
    ## Stan code for the Hybrid Rosenbrock
    stanmodelcode <- "
    data {
      real mu;
      real a;
      int <lower=0> n1;
      int <lower=0> n2;
      matrix[n2,n1] b;
    }
    parameters {
      vector[n1*n2 + 1] X;
    }
    model {
      X[n1*n2 + 1] ~ normal(mu, 1/(2*a));
      for ( j in 1:n2 ) {
        X[(j-1)*n1+1] ~ normal(X[n1*n2 + 1]^2, 1/(2*b[j,1]));
        for ( i in 2:n1 ) {
          X[(j-1)*n1+i] ~ normal(X[(j-1)*n1+i-1]^2, 1/(2*b[j,i]));
        }
      }
    }
    "
    stan_banana <<- stan(model_code = stanmodelcode,
                      data = dat, warmup = 0, iter = 0, chains = 1, 
                      verbose = FALSE) 
    f <- function(x){log_prob(stan_banana, x)}
    g <- function(x){grad_log_prob(stan_banana, x)}
  } else {
    
    ## Log likelihood for Hybrid Rosenbrock
    f <- function(X){
      n = n1*n2 + 1
      lp = dnorm(X[n], mu, 1/(2*a), log = T)
      for ( j in 1:n2 ) {
        lp = lp + dnorm(X[(j-1)*n1+1], X[n1*n2 + 1]^2, 1/(2*b[j,1]), log = T);
        for ( i in 2:n1 ) {
          lp = lp + dnorm(X[(j-1)*n1+i], X[(j-1)*n1+i-1]^2, 1/(2*b[j,i]), log = T);
        }
      }
      return(lp)
    }
    g <- NULL
  }

  sample_fun <- function(nsamples=10000){
    if(n1 < 2) stop("need more than 2 columns for b")
    n <- n2*n1 + 1
    X <- matrix(0, nsamples, n)
    X[,n] <- rnorm(nsamples, mean = mu, sd = 1/(2*a))
    for ( j in 1:n2 ) {
      X[,(j-1)*n1+1] <- rnorm(nsamples,X[,n]^2, 1/(2*b[j,1])); # pos [j,1] = (j-1)*n1+1
      for ( i in 2:n1 ) {
        X[,(j-1)*n1+i] <- rnorm(nsamples,X[,(j-1)*n1+i-1]^2, 1/(2*b[j,i])); # pos [j,i] = (j-1)*n1+i
      }
    }
    return(X)
  }
  
  return(list(
    log_p = f,
    grad_log_p = g,
    samples = sample_fun
  ))
}

