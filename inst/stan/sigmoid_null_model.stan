//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] Abundance;
  vector[N] Temperature;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma;

  real<lower=0,upper=1> p;
  real<lower=0> k;
  real<lower=0> m;

}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  sigma ~ normal(0,0.2);

  p ~ normal(0.1,0.1);
  k ~ normal(500,200);
  m ~ normal(60,20);


  for(i in 1:N){
      Abundance[i] ~ normal((1-p)/(1 + exp(-k*(1/Temperature[i] - 1/m))) + p, sigma);
    }

}
