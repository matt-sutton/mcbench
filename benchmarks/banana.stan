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
