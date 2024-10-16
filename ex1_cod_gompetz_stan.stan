

// Estimação Bayesiana na distribuição Gompertz

  data {
    int<lower=0> N;                    // tamanho amostral
    array[N] real time;                // tempo de falha observado
    array[N] int<lower=0, upper=1> delta;  // indicador da falha
  }

  parameters {
    real alpha;                   
    real<lower=0> beta;
  }

  model {
    // Priori para beta: p(beta) = 1/beta
  target += -log(beta);  // Log-priori para beta (Jeffreys)

    // Definição da verossimilhança manualmente
    for (i in 1:N) {
      // Verossimilhança para cada observação i
      target +=  log(beta)*delta[i] + alpha*delta[i]*time[i] - (beta / alpha) * (exp(alpha * time[i]) - 1);
    }
  }


