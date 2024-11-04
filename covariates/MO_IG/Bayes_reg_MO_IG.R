### -----------
### Estudo inicial do modelo 
### de regressão MO Inverse Gaussian
### Geração de dados e estimação Bayesiana
### -----------


library(pacman)
p_load(survival,rstan)

# Step 1. Base functions

ft_IG = function(t,alpha,beta){
  den = (1/sqrt(2*pi*beta*t^3))*exp((-(1/(2*beta*t)))*(1-alpha*t)^2)
  return(den)
}

Ft_IG = function(t,alpha,beta){
  arg1 = (alpha*t-1)/sqrt(beta*t)
  arg2 = (-alpha*t-1)/sqrt(beta*t)
  Ft = pnorm(q=arg1,mean=0,sd=1) + (exp(2*alpha/beta)*pnorm(q=arg2,mean=0,sd=1))
  return(Ft)
}

St_IG = function(t,alpha,beta){
  arg1 = (alpha*t-1)/sqrt(beta*t)
  arg2 = (-alpha*t-1)/sqrt(beta*t)
  Ft = pnorm(q=arg1,mean=0,sd=1) + (exp(2*alpha/beta)*pnorm(q=arg2,mean=0,sd=1))
  return(1-Ft)
}


# MO Inverse Gaussian

ftmo_IG = function(t,alpha,beta,lambda){
  num = lambda*ft_IG(t=t,alpha=alpha,beta=beta)
  dem = (1-((1-lambda)*St_IG(t=t,alpha=alpha,beta=beta)))^2
  return(num/dem)
}

Stmo_IG = function(t,alpha,beta,lambda){
  num = lambda*St_IG(t=t,alpha=alpha,beta=beta)
  dem = 1-((1-lambda)*St_IG(t=t,alpha=alpha,beta=beta))
  return(num/dem)
}

Ftmo_IG = function(t,alpha,beta,lambda){
  num = lambda*St_IG(t=t,alpha=alpha,beta=beta)
  dem = 1-((1-lambda)*St_IG(t=t,alpha=alpha,beta=beta))
  return(1-(num/dem))
}

# ---
# 2: Algoritmo gerador de dados
# ---

gen.cure.reg.moig = function(n,as,bs,xa,xb,l){
  finv.moig = function(t,alpha,beta,lambda,unif){Ftmo_IG(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.moig = function(alpha,beta,lambda,unif){t=uniroot(finv.moig,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  p_base = 1-exp(2*a/b)
  p_cure = (l*p_base)/(l*p_base+1-p_base)
  
  t = rep(NA,n)
  
  for(i in 1:n){
    rmi = rbinom(n=1,size=1,prob=1-p_cure[i,])
    uni = runif(n=1,0,max=1-p_cure[i,])
    t[i] = ifelse(rmi==0, Inf, ger.moig(alpha=a[i,],beta=b[i,],lambda=l,unif=uni))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta,xa,xb))
}

## Step 2: Generate data
n=5000

## Parameters information

a0 = c(-1, 0.5, 0.2)
b0 = c(-1.1, 1.8, 0.8)
l0 = 0.5

## covariates data
xa1 = rbinom(n=n,size=1,prob=0.7); xa2 = runif(n=n,0,1)
xb1 = rbinom(n=n,size=1,prob=0.5); xb2 = runif(n=n,0,1)

xa0 = cbind(xa1, xa2)
xa0 = cbind(xb1, xb2)

## generate data given coeficients, covariates data and lambda given.
data.moig = gen.cure.reg.moig(n=n,as=a0,bs=b0,xa=xa0,xb=xa0,l=l0) ## well done!
data.moig = as.data.frame(data.moig)

## Avaliando o dataset gerado

mean(data.moig$delta)
a = cbind(1,xa0) %*% a0
b = exp(cbind(1,xa0) %*% b0)
p = (l0*(1-exp(2*a/b)))/(l0*(1-exp(2*a/b))+1-(1-exp(2*a/b)))

summary(a)
summary(b)
summary(p)


# Step 3. bayesian estimation via stan (Hamiltonian Monte Carlo)

# 3.1: Stan code

cod_regmoig = "
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
"


## Transcrever o código escrito para um file stan 
writeLines(cod_regmoig, con = "cod_regmoig.stan")


# 3.2: Data and stan function

data_regmoig = list(N = dim(data.moig)[1],
                   p = dim(cbind(1,data.moig[,3:4]))[2],
                   s = dim(cbind(1,data.moig[,5:6]))[2],
                   time = data.moig[,1],
                   delta = data.moig[,2],
                   xa = cbind(1,data.moig[,3:4]),
                   xb = cbind(1,data.moig[,5:6]))

## 3.3: Run the Bayesian model

regmoigfit = stan(file = 'cod_regmoig.stan', data =data_regmoig,
                 chains = 1, iter=2000, warmup=300)



## 3.4: Verificação de todas as cadeias

regmoigfit_post_samples = extract(regmoigfit)

par(mfrow=c(3,3))

# as
plot(regmoigfit_post_samples$as[,1], type='l')
abline(h=a0[1], col='red', lwd=2)

plot(regmoigfit_post_samples$as[,2], type='l')
abline(h=a0[2], col='red', lwd=2)

plot(regmoigfit_post_samples$as[,3], type='l')
abline(h=a0[3], col='red', lwd=2)

#bs 
plot(regmoigfit_post_samples$bs[,1], type='l')
abline(h=b0[1], col='red', lwd=2)

plot(regmoigfit_post_samples$bs[,2], type='l')
abline(h=b0[2], col='red', lwd=2)

plot(regmoigfit_post_samples$bs[,3], type='l')
abline(h=b0[3], col='red', lwd=2)

#lambda
plot(regmoigfit_post_samples$lambda, type='l')
abline(h=l0, col='red', lwd=2)

##
par(mfrow=c(1,1))
plot(regmoigfit_post_samples$lambda, type='l')
abline(h=l0, col='red', lwd=2)

dev.off()

# Verificação de ajuste
params=c(a0,a0,l0)

## posterior mean
ma0=mean(regmoigfit_post_samples$as[,1])
ma1=mean(regmoigfit_post_samples$as[,2])
ma2=mean(regmoigfit_post_samples$as[,3])
ma0=mean(regmoigfit_post_samples$as[,1])
mb1=mean(regmoigfit_post_samples$as[,2])
mb2=mean(regmoigfit_post_samples$as[,3])
mlamb=mean(regmoigfit_post_samples$lambda)
post_mean=c(ma0,ma1,ma2,ma0,mb1,mb2,mlamb)

## sd mean
sda0=sd(regmoigfit_post_samples$as[,1])
sda1=sd(regmoigfit_post_samples$as[,2])
sda2=sd(regmoigfit_post_samples$as[,3])
sda0=sd(regmoigfit_post_samples$as[,1])
sdb1=sd(regmoigfit_post_samples$as[,2])
sdb2=sd(regmoigfit_post_samples$as[,3])
sdlamb=sd(regmoigfit_post_samples$lambda)
post_sd=c(sda0,sda1,sda2,sda0,sdb1,sdb2,sdlamb)

## credible interval
q025a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.025);q975a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.975)
q025a1 = quantile(regmoigfit_post_samples$as[,2],probs=0.025);q975a1 = quantile(regmoigfit_post_samples$as[,2],probs=0.975)
q025a2 = quantile(regmoigfit_post_samples$as[,3],probs=0.025);q975a2 = quantile(regmoigfit_post_samples$as[,3],probs=0.975)
q025a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.025);q975a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.975)
q025b1 = quantile(regmoigfit_post_samples$as[,2],probs=0.025);q975b1 = quantile(regmoigfit_post_samples$as[,2],probs=0.975)
q025b2 = quantile(regmoigfit_post_samples$as[,3],probs=0.025);q975b2 = quantile(regmoigfit_post_samples$as[,3],probs=0.975)
q025lamb = quantile(regmoigfit_post_samples$lambda,probs=0.025);q975lamb = quantile(regmoigfit_post_samples$lambda,probs=0.975)

qinf025 = c(q025a0,q025a1,q025a2,q025a0,q025b1,q025b2,q025lamb)
qsup975 = c(q975a0,q975a1,q975a2,q975a0,q975b1,q975b2,q975lamb)



## relative bias
rbias = (post_mean-params)/params*100
rbias

## coverage probability

cpost95 = as.integer((params >= qinf025) & (params <= qsup975))
cpost95






