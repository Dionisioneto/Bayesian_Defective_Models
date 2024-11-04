### -----------
### Estudo inicial do modelo 
### de regressão MO inversa Gaussiana
### Geração de dados e estimação via Máxima Verossimilhança
### -----------

library(pacman)
p_load(survival)

# 1. Base functions

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


# t: Failure time. 
# alpha, beta: Shape parameters.
# lambda: MO parameter.


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

## Parameters information
n=5000

a0 = c(-1, 0.5, 0.2)
b0 = c(-1.1, 1.8, 0.8)
l0 = 0.5

## covariates data

xa1 = rbinom(n=n,size=1,prob=0.7); xa2 = runif(n=n,0,1)
xb1 = rbinom(n=n,size=1,prob=0.5); xb2 = runif(n=n,0,1)

xa0 = cbind(xa1, xa2)
xb0 = cbind(xb1, xb2)

## generate data given coeficients, covariates data and lambda given.
data.moig = gen.cure.reg.moig(n=n,as=a0,bs=b0,xa=xa0,xb=xb0,l=l0) ## well done!
data.moig = as.data.frame(data.moig)

## Avaliando o dataset gerado

mean(data.moig$delta)
a = cbind(1,xa0) %*% a0
b = exp(cbind(1,xa0) %*% b0)
p = (l0*(1-exp(2*a/b)))/(l0*(1-exp(2*a/b))+1-(1-exp(2*a/b)))

summary(a)
summary(b)
summary(p)


## ---
## Maximum likelihood estimation
## BFGS algorithm in optim() function
## ---

## Step 1: Build the log-likelihood function

loglik_regmoig = function(par,time,delta,xa,xb){
  #parameters
  as = par[1:(dim(xa)[2]+1)] # b0, b1, b2
  bs = par[(dim(xa)[2]+2):(dim(xa)[2]+1+dim(xb)[2]+1)] # d0, d1, d2
  l = par[(dim(xa)[2]+1+dim(xb)[2]+1)+1] # l0
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  
  den = as.vector(ftmo_IG(t=time,alpha=a,beta=b,lambda=l))
  surv = as.vector(Stmo_IG(t=time,alpha=a,beta=b,lambda=l))
  
  log_lik = sum(delta*log(den) + (1-delta)*log(surv))
  return(-1*log_lik)
}

## Step 2: Choose the initial steps and maximize using optim

params = c(a0,b0,l0)
#chutes = rep(1,4)
start = params

#start = c(-1,0.1,0.1,-1,2,1,1)

maxlike_regmoig = optim(par = start,
                       fn = loglik_regmoig,
                       gr = NULL,
                       hessian = T,
                       method = "BFGS",
                       time = data.moig$t2,
                       delta = data.moig$delta,
                       xa=cbind(data.moig$xa1,data.moig$xa2),
                       xb=cbind(data.moig$xb1,data.moig$xb2))


params
## point estimation (MLE estimate)
maxlike_regmoig$par

## relative bias 
(maxlike_regmoig$par-params)/params

## standard error
se = sqrt(diag(solve(maxlike_regmoig$hessian)))
coverage = se*1.96

## coverage probability

inf95 = maxlike_regmoig$par - coverage
sup95 = maxlike_regmoig$par + coverage
cp95 = as.integer((params >= inf95) & (params <= sup95))
cp95




































