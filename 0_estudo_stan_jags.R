### ---
### 0. Estimação Bayesiana dos sofwares Stan e Jags 
### ---

## ---
## Passo 1: Instalação das bilbiotecas e pacotes na máquina
## ---

# stan
#install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))

# JAGS
# baixar o JAGS do site: https://sourceforge.net/projects/mcmc-jags/
# instalar o pacote R2JAGS

#install.packages("R2jags")

library(rstan)
library(R2jags)

## ---
## Passo 2: Exemplo de estimação de um modelo gaussiano
## ---

# 2.1: Gerando uma amostra de uma distribuição normal

n = 1000
mu = 12.98
sigma2 = 5.35

s.normal = rnorm(n=n,mean=mu,sd=sqrt(sigma2))


# 2.2: Realizando a estimação via stan

# 2.2.1: Escrita do código para a especificação da amostra aleatória, parâmetros e modelo dos dados.

codigo_normal <- "

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
"

## Passar o código escrito para um file stan 
writeLines(codigo_normal, con = "codigo_normal.stan")

# 2.2.2: Criação da lista dos dados
# Deve especificar o tamanho amostral e os dados amostrais

data_list = list(
    N = length(s.normal),
    y = s.normal
)


# 2.2.3: Utilização da função stan para o Hamitolniano Monte Carlo, 
# sem especificação da distribuição a priori

normbayes = stan(file = "codigo_normal.stan",        ## Escrita do código (2.2.1) 
                 data = data_list,                 ## Lista de dados (2.2.2)
                 iter = 5000,                      ## Número de iterações
                 chains = 1)                       ## Quantidade de cadeias 


# 2.2.4: Resultados da estimação

print(normbayes) # sumário da estimação Bayesiana

# Valores estimados para a média e variancia

posterior_samples = extract(normbayes)

mu_samples = posterior_samples$mu
sigma2_samples = posterior_samples$sigma^2


# traceplots

plot(1:length(mu_samples),mu_samples,type="l")
plot(1:length(sigma2_samples),sigma2_samples,type="l")

# já realiza o burn-in!


# autocorrelation plot

acf(mu_samples)
acf(sigma2_samples)

# As autocorrelações estão perfeitas

# visualização da distribuição a posteriori

hist(mu_samples)
abline(v=mu,lwd=2,col="red")


hist(sigma2_samples)
abline(v=sigma2,lwd=2,col="red")


# 2.2.3: Utilização da função stan para o Hamitolniano Monte Carlo, 
# com especificação da distribuição a priori.

# No problema de Behens-Fisher, é comum de se utilizar a priori
# em dois estágios Normal-Gama Inversa.

# Dessa forma, para este problema de Inferência, devemos modificar o 
# texto script no stan.


cod_normal_BF <- "

data{
  int<lower=0> N;  // Número de observações
  array[N] real y;      // Dados observados
}

parameters{
  real mu;              // Parâmetro da média dos dados
  real<lower=0> sigma;  // Parâmetro do desvio-padrão dos dados
  real theta;            // Parâmetro da média de mu
  real<lower=0> tau;    // Desvio-padrão de mu
}

model{

// [Prioris]
tau ~ inv_gamma(0.001, 0.001); // Priori não informativa para tau
mu ~ normal(theta,sqrt(tau));  // Priori para mu, com média theta e d.p tau

// [Verossimilhança]
y ~ normal(mu, sigma); // Especificação dos dados assumirem uma distribuição normal
}
"

## Passar o código escrito para um file stan 
writeLines(cod_normal_BF, con = "cod_normal_BF.stan")

# 2.2.2: Criação da lista dos dados
# Deve especificar o tamanho amostral e os dados amostrais

data_list = list(
  N = length(s.normal),
  y = s.normal
)


# 2.2.3: Utilização da função stan para o Hamitolniano Monte Carlo, 
# sem especificação da distribuição a priori

normbayes_BF = stan(file = "cod_normal_BF.stan",        ## Escrita do código (2.2.1) 
                 data = data_list,                 ## Lista de dados (2.2.2)
                 iter = 5000,                      ## Número de iterações
                 chains = 2)                       ## Quantidade de cadeias 


# 2.2.4: Resultados da estimação

print(normbayes_BF) # sumário da estimação Bayesiana



# Valores estimados para a média e variancia

posterior_samples = extract(normbayes_BF)

mu_samples = posterior_samples$mu
sigma2_samples = posterior_samples$sigma^2


# traceplots

plot(1:length(mu_samples),mu_samples,type="l")
plot(1:length(sigma2_samples),sigma2_samples,type="l")

# já realiza o burn-in!


# autocorrelation plot

acf(mu_samples)
acf(sigma2_samples)

# As autocorrelações estão perfeitas

# visualização da distribuição a posteriori

hist(mu_samples)
abline(v=mu,lwd=2,col="red")


hist(sigma2_samples)
abline(v=sigma2,lwd=2,col="red")



# 2.3: Realizando a estimação via JAGS
# Parece o mesmo processo do stan, porém é preciso realizar
# algumas mudanças na sintaxe

# no código, precisa ser sigma^-2, pois o jags modela a 
# precisão ao invés da variância

cod_jags_normal = "
model
  { # Verossimilhança
  for(i in 1:n){y[i] ~ dnorm(mu,sigma^-2)}
  
    # Prioris
    mu ~ dunif(-1000,1000)
    sigma ~ dunif(0,1000)
  }

"

model_data = list(n=n,y=s.normal)
model_params = c("mu", "sigma")


normjags = jags(data = model_data,
                parameters.to.save = model_params,
                model.file = textConnection(cod_jags_normal),
                n.chains = 4,
                n.iter = 5000,
                n.burnin = 500,
                n.thin=2,
                DIC = FALSE)

## analise

print(normjags)
traceplot(normjags)


## amostras
mcmc_samples = as.mcmc(normjags)

# Extraindo amostras específicas
mu_samples_jags = unlist(mcmc_samples[, "mu"])
sigma_samples_jags = unlist(mcmc_samples[, "sigma"])

# Calculando as amostras de variância
sigma2_samples_jags = sigma_samples_jags^2

# Exibir as primeiras amostras da média e variância
hist(mu_samples_jags)
abline(v=mu,lwd=2,col="red")

hist(sigma2_samples_jags)
abline(v=sigma2,lwd=2,col="red")



















