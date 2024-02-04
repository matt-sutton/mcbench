# Multivariate ESS 
get_multiESS <- function(samples){
  return(mcmcse::multiESS(samples))
}

# Min ESS 
get_minESS <- function(samples){
  return(min(mcmcse::ess(samples)))
}

# log_p ESS 
get_logpESS <- function(samples, log_p){
  lp_samp <- apply(samples, MARGIN = 1, log_p)
  return(mcmcse::ess(lp_samp))
}

