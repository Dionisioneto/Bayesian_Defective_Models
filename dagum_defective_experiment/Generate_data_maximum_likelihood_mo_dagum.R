### ------
### Fase 1: Estudos das funções do modelo 
###  Distribuição defectiva dagum
### ------

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(survival, survminer)

### ------ 
### 1. Geração de dados do modelo Dagum defeituoso ====
### ------

# 1.1: Funções Importantes (Dagum).

# t: Failure time. 
# alpha, beta, theta
# theta in (0,1)


ft_dagum = function(t,alpha,beta,theta){
  den = (alpha*beta*theta^2)*t^(-alpha-1)/(beta+theta*t^(-alpha))^2
  return(den)
}

Ft_dagum = function(t,alpha,beta,theta){
  dist = theta*beta/(beta + theta*t^(-alpha))
  return(dist)
}


st_dagum = function(t,alpha,beta,theta){
  dist = theta*beta/(beta + theta*t^(-alpha))
  return(1-dist)
}

## Funções
## 1.3 Marshall-Olkin Dagum

ft_modagum = function(t,alpha,beta,theta,lambda){
  den = (lambda*ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta))/(1-(1-lambda)*st_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^2
  return(den)
}


Ft_modagum = function(t,alpha,beta,theta,lambda){
  surv = lambda*st_dagum(t=t,alpha=alpha,beta=beta,theta=theta)/(1-(1-lambda)*st_dagum(t=t,alpha=alpha,beta=beta,theta=theta))
  return(1-surv)
}

st_modagum = function(t,alpha,beta,theta,lambda){
  surv = lambda*st_dagum(t=t,alpha=alpha,beta=beta,theta=theta)/(1-(1-lambda)*st_dagum(t=t,alpha=alpha,beta=beta,theta=theta))
  return(surv)
}


### ------ 
### 2. Algoritmo Geração de dados com fração de cura ====
### ------

# 1.2: Algoritmo de geração de dados sob a F(t).

gen.cure.modagum = function(n,a,b,th,lb,p){
  
  finv.modagum = function(t,alpha,beta,theta,lambda,unif){Ft_modagum(t,alpha,beta,theta,lambda)-unif}
  ger.modagum = function(alpha,beta,theta,lambda,unif){t=uniroot(finv.modagum, c(0,1000),tol=0.0000001,alpha=alpha,beta=beta,theta=theta,lambda=lambda,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0, Inf, ger.modagum(alpha=a,beta=b,theta=th,lambda=lb,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  
  t2 = pmin(t,u2) ; delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}



n=10000
a0=2
b0=2
th0=0.8
l0=1

pbase = 1-th0; p0 = (l0*pbase)/(l0*pbase + 1 - pbase); p0

dados.modagum = gen.cure.modagum(n=n,a=a0,b=b0,th=th0,lb=l0,p=p0)

## ---
## Estimação loglik via BFGS ===
## ---

## função log-verossimilhança

loglik_modagum = function(par, time,delta){
  
  # parametros
  alpha = par[1]
  beta = par[2]
  theta = par[3]
  lambda = par[4]
  
  densidade = ft_modagum(t=time,alpha=alpha,beta=beta,theta=theta,lambda=lambda)
  sobrevivencia = st_modagum(t=time,alpha=alpha,beta=beta,theta=theta,lambda=lambda)
  
  log_vero = sum(delta*log(densidade) + (1 - delta)*log(sobrevivencia))
  return(-1*log_vero)
}



params = c(a0,b0,th0,l0)
chutes = params
#chutes = rep(1,4)

maxlike = optim(par = chutes,
                fn = loglik_modagum,
                gr = NULL,
                hessian = TRUE,
                method = "BFGS",
                time = dados.modagum[,1],
                delta = dados.modagum[,2])

params
maxlike$par

sqrt(diag(solve(maxlike$hessian)))

maxlike$value


## reparametrização

## Reparametrização da verossimilhança

## alpha1 = log(alpha); beta1 = log(beta)
## theta1 = exp(theta)/(1 + exp(theta))
## u1 = log(u); r1 = log(r) 

ft1_dagum = function(t,alpha,beta,theta){
  alpha1 = exp(alpha); beta1 = exp(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  den = (alpha1*beta1*theta1^2)*t^(-alpha1-1)/(beta1+theta1*t^(-alpha1))^2
  return(den)
}


Ft1_dagum = function(t,alpha,beta,theta){
  alpha1 = exp(alpha); beta1 = exp(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  dist = theta1*beta1/(beta1 + theta1*t^(-alpha1))
  return(dist)
}


st1_dagum = function(t,alpha,beta,theta){
  alpha1 = exp(alpha); beta1 = exp(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  dist = theta1*beta1/(beta1 + theta1*t^(-alpha1))
  return(1-dist)
}


## funções na marshall olkin com reparametrização
ft1_modagum = function(t,alpha,beta,theta,lambda){
  lambda1=log(lambda)
  den = (lambda1*ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta))/(1-(1-lambda1)*st1_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^2
  return(den)
}


Ft1_modagum = function(t,alpha,beta,theta,lambda){
  lambda1=log(lambda)
  surv = lambda1*st1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)/(1-(1-lambda1)*st1_dagum(t=t,alpha=alpha,beta=beta,theta=theta))
  return(1-surv)
}

st1_modagum = function(t,alpha,beta,theta,lambda){
  lambda1=log(lambda)
  surv = lambda1*st1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)/(1-(1-lambda1)*st1_dagum(t=t,alpha=alpha,beta=beta,theta=theta))
  return(surv)
}



## ---
## Estimação loglik via BFGS ===
## ---

## função log-verossimilhança

loglik1_modagum = function(par, time,delta){
  
  # parametros
  alpha = par[1]
  beta = par[2]
  theta = par[3]
  lambda = par[4]
  
  densidade = ft1_modagum(t=time,alpha=alpha,beta=beta,theta=theta,lambda=lambda)
  sobrevivencia = st1_modagum(t=time,alpha=alpha,beta=beta,theta=theta,lambda=lambda)
  
  log_vero = sum(delta*log(densidade) + (1 - delta)*log(sobrevivencia))
  return(-1*log_vero)
}



params = c(a0,b0,th0,l0)
chutes = params
#chutes = rep(0.8,4)

maxlike = optim(par = chutes,
                fn = loglik1_modagum,
                gr = NULL,
                hessian = TRUE,
                method = "BFGS",
                time = dados.modagum[,1],
                delta = dados.modagum[,2])

params
maxlike$par

exp(maxlike$par[1]);exp(maxlike$par[2])
maxlike$par

sqrt(diag(solve(maxlike$hessian)))

maxlike$value








