## Rosenbluth-Teller-Metropolis-Hastings Kernel with a Random Walk
## 
RW_kernel = function (log_p, epsilon, current_q) {
  q = current_q
  q = q + epsilon * rnorm(length(q),0,1)
  
  # Evaluate log prob
  current_lp = log_p(current_q)
  proposed_lp = log_p(q)
  # Accept or reject the state 
  if (runif(1) < exp(proposed_lp-current_lp))
  {
    return (q) # accept
  }
  else
  {
    return (current_q) # reject
  }
}

RW_sampler <- function(nsample, x, l_prob, eps){
  samples <- matrix(0, nrow = nsample, ncol = length(x))
  samples[1,] <- x
  
  for( i in 2:nsample){
    x <- RW_kernel(l_prob, eps, x)
    samples[i, ] <- x
  }
  return(samples)
}

