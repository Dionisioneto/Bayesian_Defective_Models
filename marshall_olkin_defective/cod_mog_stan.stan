
// Código em Stan para a estimação Bayesiana na distribuição Marshall-Olkin Gompertz

data {
  int<lower=0> N;                        
  array[N] real time;
  array[N] int<lower=0, upper=1> delta; 
  real lambda;
}
  
parameters {
  real alpha;
  real<lower=0> beta;
  //real<lower=0> lambda;
}
  
model {
  // Prioris
    alpha ~ normal(0,10);
    beta ~ cauchy(0,5);
    //lambda ~ cauchy(0,5);
  
    
  // Definição da verossimilhança manualmente
  for (i in 1:N) {
    // Verossimilhança para cada observação i
    target += delta[i] * log(beta * lambda * exp((beta - beta * exp(alpha * time[i])) / alpha) + alpha * time[i]) -
              delta[i] * log((lambda - (lambda - 1) * exp((beta - beta * exp(alpha * time[i])) / alpha))^2) +
              (1 - delta[i]) * (log(lambda * exp(-(beta / alpha) * (exp(alpha * time[i]) - 1)))) -
              (1 - delta[i]) * (log(1 - (1 - lambda) * exp((-beta / alpha) * (exp(alpha * time[i]) - 1))));
  }
}

