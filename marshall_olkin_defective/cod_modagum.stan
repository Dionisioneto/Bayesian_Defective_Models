
  data {
    int<lower=0> N;                         // tamanho amostral
    array[N] real time;                    //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1> theta;
  real<lower=0> lambda;
  }

  model {
  // Prioris
  alpha ~ gamma(0.1,0.1);
  beta ~ gamma(0.1,0.1);
  theta ~ beta(100,100);
  lambda ~ gamma(0.1,0.1);

  // Definição manual da log-verossimilhança
  for (i in 1:N){
  real denom = lambda + (((1-lambda)*theta*beta)/(beta + theta*time[i]^(-alpha)));

  target += delta[i]*log(lambda*alpha*beta*theta^2*time[i]^(-alpha-1)) -
             delta[i]*log((beta+theta*time[i]^(-alpha))^2 * denom^2) +
             (1-delta[i])*log(lambda*(1 - (theta*beta)/(beta+theta*time[i]^(-alpha)))) -
             (1-delta[i])*log(denom);
  }
}


