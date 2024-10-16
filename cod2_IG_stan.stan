
  data {
    int<lower=0> N;                        // tamanho amostral
    array[N] real time;                   //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
    real alpha;
    real<lower=0> beta;
  }

  model {
    // Prioris
    alpha ~ normal(0,50);
    beta ~ gamma(0.001,0.001);

    // Definição da log-verossimilhança manualmente
    for (i in 1:N) {
      target += delta[i] * log((1 / sqrt(2 * pi() * beta * time[i]^3)) *
                exp(-(1 / (2 * beta * time[i])) * pow(1 - alpha * time[i], 2))) +
                (1 - delta[i]) * log(1 - (Phi((alpha * time[i] - 1) / sqrt(beta * time[i])) +
                exp(2 * alpha / beta) * Phi((-alpha * time[i] - 1) / sqrt(beta * time[i]))));
    }
  }

