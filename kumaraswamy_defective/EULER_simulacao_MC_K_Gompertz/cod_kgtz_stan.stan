
  data{
  int<lower=0> N;
  vector[N] time;
  vector[N] delta;
  }
  
  parameters{
  real alpha;
  real<lower=0> beta;
  real<lower=0> psi;
  real<lower=0> kappa;
  }
  
  model{
  
  // prioris
  
  alpha ~ normal(-1,10);
  beta ~ gamma(0.1,0.1);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.1,0.1);
  
  // Definição manual da função de verossimilhança
  for(i in 1:N){
    target += delta[i]*(log(kappa*psi*beta*exp(alpha*time[i])*exp((beta-beta*exp(alpha*time[i]))/alpha)*(1-exp((beta-beta*exp(alpha*time[i]))/alpha))^(psi-1) * (1-(1-(exp((beta-beta*exp(alpha*time[i]))/alpha)))^psi)^(kappa-1))) +
              (1 - delta[i])*(log((1 - (1 - exp((beta-beta*exp(alpha*time[i]))/alpha))^psi)^kappa));
  }
}

