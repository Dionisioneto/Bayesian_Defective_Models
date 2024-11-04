### -----------
### Geração de dados e MLE 
### de regressão MO Gompertz
### Geração de dados e estimação via Máxima Verossimilhança
### -----------

library(pacman)
p_load(survival)

# 1. Base functions

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

## Parameters information
n=5000

# first (15%): 
# a0 = c(-1, -0.8, 0.1)
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
xb0 = cbind(xb1, xb2)

## generate data given coeficients, covariates data and lambda given.
data = gen.cure.reg.mog(n=n,as=a0,bs=b0,xa=xa0,xb=xb0,l=l0) ## well done!
data = as.data.frame(data)

## Avaliando o dataset gerado
mean(data$delta)

a = cbind(1,xa0) %*% a0
b = exp(cbind(1,xa0) %*% b0)
p = (l0*exp(b/a))/(l0*exp(b/a)+1-exp(b/a))

summary(a)
summary(b)
summary(p)

## ---
## Maximum likelihood estimation
## BFGS algorithm in optim() function
## ---

## Step 1: Build the log-likelihood function

loglik_regmog = function(par,time,delta,xa,xb){
  #parameters
  as = par[1:(dim(xa)[2]+1)] # b0, b1, b2
  bs = par[(dim(xa)[2]+2):(dim(xa)[2]+1+dim(xb)[2]+1)] # d0, d1, d2
  l = par[(dim(xa)[2]+1+dim(xb)[2]+1)+1] # l0
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  
  den = as.vector(ftmo_gompertz(t=time,alpha=a,beta=b,lambda=l))
  surv = as.vector(Stmo_gompertz(t=time,alpha=a,beta=b,lambda=l))
  
  log_lik = sum(delta*log(den) + (1-delta)*log(surv))
  return(-1*log_lik)
}


## Step 2: Choose the initial steps and maximize using optim

params = c(a0,b0,l0)

start = params-0.25
#start = params
#start = c(-1,0.1,0.1,-1,2,1,1)

maxlike_regmog = optim(par = start,
                        fn = loglik_regmog,
                        gr = NULL,
                        hessian = T,
                        method = "BFGS",
                        time = data$t2,
                        delta = data$delta,
                        xa=cbind(data$xa1,data$xa2),
                        xb=cbind(data$xb1,data$xb2))


params
## point estimation (MLE estimate)
maxlike_regmog$par

## relative bias 
(maxlike_regmog$par-params)/params

## standard error
se = sqrt(diag(solve(maxlike_regmog$hessian)))
coverage = se*1.96

## coverage probability

inf95 = maxlike_regmog$par - coverage
sup95 = maxlike_regmog$par + coverage
cp95 = as.integer((params >= inf95) & (params <= sup95))
cp95








