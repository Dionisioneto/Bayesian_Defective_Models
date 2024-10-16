

data{
  int<lower=0> N;  // Número de observações
  array[N] real y;      // Dados observados
}
parameters{
  real mu;              // Parâmetro da média
  real<lower=0> sigma;  // Parâmetro do desvio-padrão
}
model{
y ~ normal(mu, sigma); // Especificação dos dados assumirem uma distribuição normal
}

