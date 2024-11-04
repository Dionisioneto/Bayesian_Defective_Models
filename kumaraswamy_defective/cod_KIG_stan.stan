
  data {
    int<lower=0> N;                        // tamanho amostral
    vector[N] time;                        // tempo de falha observado
     array[N] int<lower=0, upper=1> delta;        // indicador do evento
  }

  parameters {
    real alpha;
    real<lower=0> beta;
    real<lower=0> psi;
    real<lower=0> kappa;
  }

  model {
    // prioris

  alpha ~ normal(-1,10);
  beta ~ gamma(0.01,0.01);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.01,0.1);

    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real exp_term = exp(2 * alpha / beta);

      // Verossimilhança do evento observado
      target += delta[i] * (log(kappa) + log(psi) - 0.5 * log(2 * beta * pi() * time[i]^3) - (1 - alpha * time[i])^2 / (2 * beta * time[i])
              + (psi - 1) * log(phi1 + exp_term * phi2));

      // Verossimilhança para censura
      target += delta[i] * (kappa - 1) * log(1 - (phi1 + exp_term * phi2)^psi);
      target += (1 - delta[i]) * kappa * log(1 - (phi1 + exp_term * phi2)^psi);
    }
  }

