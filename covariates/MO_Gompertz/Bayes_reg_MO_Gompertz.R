### -----------
### Estudo inicial do modelo 
### de regressão MO Gompertz
### Geração de dados e estimação Bayesiana
### -----------

library(pacman)
p_load(survival,rstan)

# Step 1. Base functions

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}


Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*t)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}

# MO Gompertz

ftmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*ft_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = (1-(1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))^2
  return(num/dem)
}

Ftmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*St_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = 1 - ((1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))
  return(1-(num/dem))
}

Stmo_gompertz = function(t,alpha,beta,lambda){
  num = lambda*St_Gompertz(t=t,alpha=alpha,beta=beta)
  dem = 1 - ((1-lambda)*St_Gompertz(t=t,alpha=alpha,beta=beta))
  return(num/dem)
}


## inserir xa e xb sem o vetor de 1's do intercepto
gen.cure.reg.mog = function(n,as,bs,xa,xb,l){
  finv.mog = function(t,alpha,beta,lambda,unif){Ftmo_gompertz(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.mog = function(alpha,beta,lambda,unif){t=uniroot(finv.mog,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  p_base = exp(b/a)
  p_cure = (l*p_base)/(l*p_base+1-p_base)
  
  t = rep(NA,n)
  
  for(i in 1:n){
    rmi = rbinom(n=1,size=1,prob=1-p_cure[i,])
    uni = runif(n=1,0,max=1-p_cure[i,])
    t[i] = ifelse(rmi==0, Inf, ger.mog(alpha=a[i,],beta=b[i,],lambda=l,unif=uni))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta,xa,xb))
}

# Step 2. Generate data

## Parameters information
n=5000

# first (15%): 
# a0 = c(-1, -0.5, -0.8)
# b0 = c(-1, 2, 4)
# l0 = 0.5

# second (45%)
a0 = c(-1.2, 0.5, 0.2)
b0 = c(-1.1, 1.5, 0.9)
l0 = 2


## covariates data
xa1 = rbinom(n=n,size=1,prob=0.7); xa2 = runif(n=n,0,1)
xb1 = rbinom(n=n,size=1,prob=0.5); xb2 = runif(n=n,0,1)

xa0 = cbind(xa1, xa2)
xa0 = cbind(xb1, xb2)

## generate data given coeficients, covariates data and lambda given.
dados.reg.mog = gen.cure.reg.mog(n=n,as=a0,bs=b0,xa=xa0,xb=xa0,l=l0) ## well done!
dados.reg.mog = as.data.frame(dados.reg.mog)


## Avaliando o dataset gerado
mean(dados.reg.mog$delta)
a = cbind(1,xa0) %*% a0
b = exp(cbind(1,xa0) %*% b0)
p = (l0*exp(b/a))/(l0*exp(b/a)+1-exp(b/a))

summary(a)
summary(b)
summary(p)

# Step 3. bayesian estimation via stan (Hamiltonian Monte Carlo)

# 3.1: Stan code

cod_regmog = "
data {
  int<lower=0> N; // Número de observações
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
  
  // Definição da função log-verossimilhança
  for (i in 1:N) {
    target += delta[i] * log(beta[i] * lambda * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i] + alpha[i] * time[i])) -
              delta[i] * log((1 - ((1-lambda) * exp((beta[i] - beta[i] * exp(alpha[i] * time[i])) / alpha[i])))^2) +
              (1 - delta[i]) * (log(lambda * exp(-(beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))) -
              (1 - delta[i]) * (log(1 - ((1 - lambda) * exp((-beta[i] / alpha[i]) * (exp(alpha[i] * time[i]) - 1)))));
  }
}
"


## Transcrever o código escrito para um file stan 
writeLines(cod_regmog, con = "cod_regmog.stan")


# 3.2: Data and stan function

data_regmog = list(N = dim(dados.reg.mog)[1],
                   p = dim(cbind(1,dados.reg.mog[,3:4]))[2],
                   s = dim(cbind(1,dados.reg.mog[,5:6]))[2],
                   time = dados.reg.mog[,1],
                   delta = dados.reg.mog[,2],
                   xa = cbind(1,dados.reg.mog[,3:4]),
                   xb = cbind(1,dados.reg.mog[,5:6]))

## 3.3: Run the Bayesian model

regmogfit = stan(file = 'cod_regmog.stan', data = data_regmog,
              chains = 1, iter=2000, warmup=300)

## 3.4: Verificação de todas as cadeias

regmogfit_post_samples = extract(regmogfit)

par(mfrow=c(3,3))

# as
plot(regmogfit_post_samples$as[,1], type='l')
abline(h=a0[1], col='red', lwd=2)

plot(regmogfit_post_samples$as[,2], type='l')
abline(h=a0[2], col='red', lwd=2)

plot(regmogfit_post_samples$as[,3], type='l')
abline(h=a0[3], col='red', lwd=2)

#bs 
plot(regmogfit_post_samples$bs[,1], type='l')
abline(h=b0[1], col='red', lwd=2)

plot(regmogfit_post_samples$bs[,2], type='l')
abline(h=b0[2], col='red', lwd=2)

plot(regmogfit_post_samples$bs[,3], type='l')
abline(h=b0[3], col='red', lwd=2)

#lambda
plot(regmogfit_post_samples$lambda, type='l')
abline(h=l0, col='red', lwd=2)


par(mfrow=c(1,1))
plot(regmogfit_post_samples$lambda, type='l')
abline(h=l0, col='red', lwd=2)

acf(regmogfit_post_samples$as[,1])
acf(regmogfit_post_samples$as[,2])
acf(regmogfit_post_samples$as[,3])

acf(regmogfit_post_samples$bs[,1])
acf(regmogfit_post_samples$bs[,2])
acf(regmogfit_post_samples$bs[,3])

acf(regmogfit_post_samples$lambda)
#dev.off()

## posterior mean
ma0=mean(regmogfit_post_samples$as[,1])
ma1=mean(regmogfit_post_samples$as[,2])
ma2=mean(regmogfit_post_samples$as[,3])
mb0=mean(regmogfit_post_samples$bs[,1])
mb1=mean(regmogfit_post_samples$bs[,2])
mb2=mean(regmogfit_post_samples$bs[,3])
mlamb=mean(regmogfit_post_samples$lambda)
post_mean=c(ma0,ma1,ma2,mb0,mb1,mb2,mlamb)

## sd mean
sda0=sd(regmogfit_post_samples$as[,1])
sda1=sd(regmogfit_post_samples$as[,2])
sda2=sd(regmogfit_post_samples$as[,3])
sdb0=sd(regmogfit_post_samples$bs[,1])
sdb1=sd(regmogfit_post_samples$bs[,2])
sdb2=sd(regmogfit_post_samples$bs[,3])
sdlamb=sd(regmogfit_post_samples$lambda)
post_sd=c(sda0,sda1,sda2,sdb0,sdb1,sdb2,sdlamb)

## credible interval
q025a0 = quantile(regmogfit_post_samples$as[,1],probs=0.025);q975a0 = quantile(regmogfit_post_samples$as[,1],probs=0.975)
q025a1 = quantile(regmogfit_post_samples$as[,2],probs=0.025);q975a1 = quantile(regmogfit_post_samples$as[,2],probs=0.975)
q025a2 = quantile(regmogfit_post_samples$as[,3],probs=0.025);q975a2 = quantile(regmogfit_post_samples$as[,3],probs=0.975)
q025b0 = quantile(regmogfit_post_samples$bs[,1],probs=0.025);q975b0 = quantile(regmogfit_post_samples$bs[,1],probs=0.975)
q025b1 = quantile(regmogfit_post_samples$bs[,2],probs=0.025);q975b1 = quantile(regmogfit_post_samples$bs[,2],probs=0.975)
q025b2 = quantile(regmogfit_post_samples$bs[,3],probs=0.025);q975b2 = quantile(regmogfit_post_samples$bs[,3],probs=0.975)
q025lamb = quantile(regmogfit_post_samples$lambda,probs=0.025);q975lamb = quantile(regmogfit_post_samples$lambda,probs=0.975)

qinf025 = c(q025a0,q025a1,q025a2,q025b0,q025b1,q025b2,q025lamb)
qsup975 = c(q975a0,q975a1,q975a2,q975b0,q975b1,q975b2,q975lamb)

# Verificação de ajuste
params=c(a0,b0,l0)

## relative bias
rbias = (post_mean-params)/params*100
rbias

## coverage probability

cpost95 = as.integer((params >= qinf025) & (params <= qsup975))
cpost95


















