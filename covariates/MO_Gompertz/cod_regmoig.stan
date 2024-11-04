
  data {
  int<lower=0> N; // Número de oaservações
  int<lower=0> p; // Número de covariáveis para alpha
  int<lower=0> s; // Número de covariáveis para beta

  vector[N] time;  // Tempo de falha
  vector[N] delta; // Indicador de censura
  
  matrix[N, p] xa; // Matriz de covariáveis para alpha (N x p)
  matrix[N, s] xb; // Matriz de covariáveis para beta (N x s)
  }
  
  parameters {
  vector[p] as; // Coeficientes para alpha
  vector[s] bs; // Coeficientes para beta
  real<lower=0> lambda; 
  }
  
  transformed parameters {
  vector[N] alpha = xa * as; // Cálculo de alpha
  vector[N] beta = exp(xb * bs); // Cálculo de beta
  }
  
  model {
    // Priori para os coeficientes em alpha
    for (j in 1:p) {
      as[j] ~ normal(0, 10);
    }
    
    // Priori para os coeficientes em beta
    for (j in 1:s) {
      bs[j] ~ normal(0, 10);
    }
    
    // Priori para lambda
    lambda ~ gamma(0.01, 0.01);
    
    for (i in 1:N) {
    // funções
      real z1 = (-1 + alpha[i] * time[i]) / sqrt(beta[i] * time[i]);
      real z2 = (-1 - alpha[i] * time[i]) / sqrt(beta[i] * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
    
    // log-verossimilhança
    target += delta[i] * log(lambda*exp((-1/(2*beta[i]*time[i]))*(1-alpha[i]*time[i])^2)) -
              delta[i] * log(sqrt(2*pi()*beta[i]*time[i]^3) * (lambda + (1-lambda)*(phi1 + exp(2*alpha[i]/beta[i])*phi2))^2) +
              (1-delta[i])*log(lambda-lambda*(phi1 + exp(2*alpha[i]/beta[i])*phi2)) -
              (1-delta[i])*log(lambda+(1-lambda)*(phi1 + exp(2*alpha[i]/beta[i])*phi2));
    }
  }

