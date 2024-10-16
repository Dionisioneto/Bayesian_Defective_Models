

data {
  int<lower=0> N;                        
  array[N] real time;
  array[N] int<lower=0, upper=1> delta; 
}
  
parameters {
  real alpha;
  real<lower=0> beta;
  real<lower=0> lambda;
}
  
model {
  // Prioris
  alpha ~ normal(0,10);
  beta ~ gamma(0.001,0.001);
  lambda ~ gamma(0.001,0.001);

  // Definição da verossimilhança manualmente
  real exp_term;  // Declare the variable here

  for (i in 1:N) {
    // Calculate exp_term for each i
    exp_term = exp((-beta / alpha) * (exp(alpha * time[i]) - 1));

    target += delta[i] * log(lambda * beta * exp(alpha * time[i]) * exp_term) -
              delta[i] * log((1 - (1 - lambda) * exp_term)^2) +
              (1 - delta[i]) * log(lambda * exp_term) -
              (1 - delta[i]) * log(1 - (1 - lambda) * exp_term);
  }
}

