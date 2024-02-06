functions {
  vector spectralDensity(real[] loc, real[] sca, real[] amp,
                         vector wav, int nWV, int nPK) {
    vector[nWV] mu;
    for (i in 1:nWV) {
      mu[i] = 0;
      for (j in 1:nPK) {
        mu[i] += sca[j]*pi()*exp(cauchy_lpdf(wav[i] | loc[j], sca[j]))*amp[j];
      }
    }
    return mu;
  }
}

// The input data is a vector 'y' of length 'nWL'.
data {
  int<lower=1> nWL; # number of wavenumbers
  int<lower=0> nPK; # number of peaks
  vector[nWL] y;    # spectral density
  vector[nWL] x;    # wavenumbers
}
transformed data {
  real minWL = min(x);
  real maxWL = max(x);
  real maxY = max(y);
}

// The parameters of the spectral density functions
parameters {
  real<lower=minWL, upper=maxWL> loc[nPK]; # location
  real<lower=0, upper=maxY> amp[nPK];      # amplitude
  real<lower=0> sca[nPK];                  # scale
  real<lower=0> sigma;
}
transformed parameters {
  vector[nWL] mu = spectralDensity(loc, sca, amp, x, nWL, nPK);
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (j in 1:nPK) {
    loc[j] ~ uniform(minWL, maxWL);
    amp[j] ~ uniform(0, maxY);
    sca[j] ~ lognormal(log(11.6) - (0.4^2)/2, 0.4);
  }
  sigma ~ inv_gamma(5, 1125);
  y ~ normal(mu, sigma);
}
