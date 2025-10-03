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
  int<lower=0,upper=1> Treatment[N];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> sigma;
  
  real<lower=0,upper=1> p_Tmt;
  real<lower=0> k_Tmt;
  real<lower=0> m_Tmt;
  
  real<lower=0,upper=1> p_Ctrl;
  real<lower=0> k_Ctrl;
  real<lower=0> m_Ctrl;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  sigma ~ normal(0,0.2);
  
  p_Tmt ~ normal(0.1,0.1);
  k_Tmt ~ normal(500,200);
  m_Tmt ~ normal(60,20);
  
  p_Ctrl ~ normal(0.1,0.1);
  k_Ctrl ~ normal(500,200);
  m_Ctrl ~ normal(60,20);
  
  for(i in 1:N){
    if(Treatment[i]==1){
      Abundance[i] ~ normal((1-p_Tmt)/(1 + exp(-k_Tmt*(1/Temperature[i] - 1/m_Tmt))) + p_Tmt, sigma);
    }else if(Treatment[i]==0){
      Abundance[i] ~ normal((1-p_Ctrl)/(1 + exp(-k_Ctrl*(1/Temperature[i] - 1/m_Ctrl))) + p_Ctrl, sigma);
    }
  }
}
generated quantities {
  real ATE;
  
  real Abundance_doTmt;
  real Abundance_doCtrl; 
  
  Abundance_doTmt = 0;
  Abundance_doCtrl = 0;
  
  for(i in 1:N){
    
    Abundance_doTmt += normal_rng((1-p_Tmt)/(1 + exp(-k_Tmt*(1/Temperature[i] - 1/m_Tmt))) + p_Tmt, sigma);
    Abundance_doCtrl += normal_rng((1-p_Ctrl)/(1 + exp(-k_Ctrl*(1/Temperature[i] - 1/m_Ctrl))) + p_Ctrl, sigma);
  }
  
  Abundance_doTmt /=N;
  Abundance_doCtrl /=N; 
  
  ATE = Abundance_doTmt - Abundance_doCtrl;
}
