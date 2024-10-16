
  data {
  int<lower=0> N;                         // tamanho amostral
    array[N] real time;                    //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1> theta;
  real<lower=0> u;
  real<lower=0> r;
  }

  model {
  // Prioris
  alpha ~ gamma(0.25,0.25);
  beta ~ gamma(0.25,0.25);
  theta ~ beta(1,1);
  u ~ gamma(0.25,0.25);
  r ~ gamma(0.25,0.25);

  // Definição manual da log-verossimilhança
  for (i in 1:N){
  real Ft = (theta*beta)/(beta + theta*time[i]^(-alpha));

  target += delta[i]*log((u*r*alpha*beta*theta^2*time[i]^(-alpha-1))/(beta + theta*time[i]^(-alpha))^2) +
            delta[i]*log(Ft^(r-1)) + delta[i]*log((1 - Ft^r)^(u-1)) +
            (1-delta[i])*log((1-Ft^r)^u);
  }
}

