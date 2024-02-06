//////
// Factor analysis example for exchange data
// 
// Original model from
// Lopes, H. F. and West, M. (2004).
// “Bayesian model assessment in factor analysis.”
// Statistica Sinica, 14(1):41–67.
// Data from
// West, M. and Harrison, J. (1997). "Bayesian forecasting
// and dynamic models."" New York: Springer-Verlag New York.
// 
// Code by Leah South
// 
// Description:
// This is typically considered as a model choice problem with the
// number of factors being k=1, k=2 or k=3.
// The 1-factor model is 12-dimensional and appears roughly Gaussian
// so it is easy to sample from.
// The 2-factor model is multimodal with well-
// separated modes (note the additional one discovered in 
// South, L. F., Pettitt, A. N., & Drovandi, C. C. (2019).
// "Sequential Monte Carlo samplers with independent Markov chain
// Monte carlo proposals." Bayesian Analysis", 14(3):753-776)
// The 3-factor model exhibits complex dependencies and
// multimodality
//////

data {
  int <lower=0> k;
  int <lower=0> n;
  int <lower=0> d;
  matrix[n,d] y;
}

parameters {
  vector<lower=0>[d] lambda;
  vector<lower=0>[k] beta_diag;
  vector[d*k - k*(k+1)%/%2] beta_offdiag;
}

model {
  // Priors
  for (j in 1:k){
    beta_diag[j] ~ normal(0,1) T[0,];
  }
  for (j in 1:d){
    lambda[j] ~ inv_gamma(1.1,0.05);
  }
  for (j in 1:(d*k - k*(k+1)%/%2)){
    beta_offdiag[j] ~ normal(0,1);
  }
  
  // Getting lower-diagonal beta matrix
  matrix[d,k] beta = rep_matrix(0,d,k);
  int counter = 0;
  for (j in 1:k) {
    beta[j,j] = beta_diag[j];
    for (i in (j+1):d) {
      counter = counter + 1;
      beta[i,j] = beta_offdiag[counter];
    }
  }
  
  // Getting covariance and writing form for likelihood
  matrix[d,d] Omega = diag_matrix(lambda) + beta*beta';
  for ( i in 1:n ) {
    y[i,:] ~ multi_normal(rep_vector(0,d), Omega);
  }
  
}
