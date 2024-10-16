### ---
### 1. Estimação Bayesiana em modelos de Regressão
### ---


## Bibliotecas
library(rstan)
library(R2jags)

## ---
## Passo 1: Modelo de Regressão Linear Gaussiano
## ---

## 1.1: Geração de dados
## yi = (beta0 + beta1 xi1 + beta2 xi2 + beta3) + ei
## ei ~ Normal(0,sigma^2)

set.seed(123)

n = 2346
sigma2 = 1
betas = c(1.4,0.5,-0.67,0.87)
x1 = rnorm(n=n,mean=10,sd=2.9)
x2 = rnorm(n=n,mean=1,sd=0.87)
x3 = rbinom(n=n,size=1,prob=0.56)

X = cbind(1,x1,x2,x3)
mus = X %*% betas 


dados = rnorm(n=n,mean = mus, sd = sqrt(sigma2))

## 1.2: Estimação dentro do stan
## data list, código e transcrever o código.

data_list = list(n = n, 
                 y = dados,
                 x1 = x1, x2 = x2, x3 = x3)


cod_regnorm_stan = "

// Modelo de Regressão Linear no Stan

data {
  int<lower=0> n;
  vector[n] x1;         // Listando as covariáveis
  vector[n] x2;
  vector[n] x3;
  vector[n] y;         // Listando a variável resposta

}

parameters {
  real beta_0;         // Intercepto
  real beta_1;        // Efeitos
  real beta_2;
  real beta_3;
  real<lower=0> sigma;  // Desvio-padrão do erro
}

model {
  // Priors
  beta_0 ~ normal(0,100);
  beta_1 ~ normal(0,100);
  beta_2 ~ normal(0,100);
  beta_3 ~ normal(0,100);
  sigma ~ gamma(0.01,0.01);
  
  
  y ~ normal(beta_0 + beta_1 * x1 + beta_2 * x2 + beta_3 * x3, sigma);
}

"

## Transcrever o código escrito para um file stan 
writeLines(cod_regnorm_stan, con = "cod_regnorm_stan.stan")


# 1.3: Utilização da função stan para o Hamitolniano Monte Carlo.

regnormstan = stan(file = "cod_regnorm_stan.stan",        ## Escrita do código (2.2.1) 
                 data = data_list,                 ## Lista de dados (2.2.2)
                 iter = 5000,                      ## Número de iterações
                 chains = 1)                       ## Quantidade de cadeias 



print(regnormstan) # sumário da estimação Bayesiana

matrix(c(betas,sqrt(sigma2)))



## ---
## Passo 2: Processo de estimação do JAGS
## ---

## O JAGS utiliza do mesmo data list que escrevemos para o stan;
## dessa forma, iremos apenas escrever o código, transcrevelo e
## utilizar a função jags model.


## 2.1: O código JAGS

cod_jags_regnorm = "
  model{
  # Prioris
  beta_0 ~ dnorm(0, 0.001)
  beta_1 ~ dnorm(0, 0.001)
  beta_2 ~ dnorm(0, 0.001)
  beta_3 ~ dnorm(0, 0.001)
  sigma ~ dgamma(0.001,0.001)
  
  # Verossimilhança
  for(i in 1:n){
  mu[i] = beta_0 + beta_1*x1[i] + beta_2*x2[i] + beta_3*x3[i]
  y[i] ~ dnorm(mu[i], sigma^-2)
  }
  
  }
"

params_jags_regnorm = c("beta_0", "beta_1", "beta_2", "beta_3", "sigma")


regnorm_jags = jags(data = data_list,
                parameters.to.save = params_jags_regnorm,
                model.file = textConnection(cod_jags_regnorm),
                n.chains = 2,
                n.iter = 5000,
                n.burnin = 500,
                n.thin=2,
                DIC = FALSE)

print(regnorm_jags)
traceplot(regnorm_jags)

matrix(c(betas,sqrt(sigma2)))


## amostras
mcmc_samples = as.mcmc(regnorm_jags)

# Extraindo amostras específicas
b0_samples_jags = unlist(mcmc_samples[, "beta_0"])
b1_samples_jags = unlist(mcmc_samples[, "beta_1"])
b2_samples_jags = unlist(mcmc_samples[, "beta_2"])
b3_samples_jags = unlist(mcmc_samples[, "beta_3"])
sigma_samples_jags = unlist(mcmc_samples[, "sigma"])

# Calculando as amostras de variância
sigma2_samples_jags = sigma_samples_jags^2

# Exibir as primeiras amostras da média e variância
hist(b0_samples_jags)
abline(v=betas[1],lwd=2,col="red")

hist(b1_samples_jags)
abline(v=betas[2],lwd=2,col="red")

hist(b2_samples_jags)
abline(v=betas[3],lwd=2,col="red")

hist(b3_samples_jags)
abline(v=betas[4],lwd=2,col="red")

hist(sigma2_samples_jags)
abline(v=sigma2,lwd=2,col="red")



## ---
## Passo 2: Modelo de Regressão Logística
## ---

1/(1+exp(-mus))
























