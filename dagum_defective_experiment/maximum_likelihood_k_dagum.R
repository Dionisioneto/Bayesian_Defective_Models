### ------
### Fase 3: Estimação Frequentista. 
###  Distribuição defectiva Kumaraswamy dagum
### ------

library(pacman)
p_load(survival)

## Funções de base do modelo

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

## 1.2 Kumaraswamy Dagum

ft_kdagum = function(t,alpha,beta,theta, u, r){
  den = u*r*ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)*(Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^(r-1)*(1-(Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^r)^(u-1)
  return(den)
}

Ft_kdagum = function(t,alpha,beta,theta, u, r){
  surv = (1 - (Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^r)^u
  return(1-surv)
}

st_kdagum = function(t,alpha,beta,theta, u, r){
  surv = (1 - (Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta))^r)^u
  return(surv)
}

## função log-verossimilhança

loglik_kdagum = function(par, time,delta){
  
  # parametros
  alpha = par[1]
  beta = par[2]
  theta = par[3]
  u = par[4]
  r = par[5]
  
  densidade = ft_kdagum(t=time, alpha=alpha,beta=beta,theta=theta,u=u,r=r)
  sobrevivencia = st_kdagum(t=time, alpha=alpha,beta=beta,theta=theta,u=u,r=r)
  
  log_vero = sum(delta*log(densidade) + (1 - delta)*log(sobrevivencia))
  return(-1*log_vero)
}




# Algoritmo de geração de dados sob a F(t) Kdagum ===

gen.cure.kdagum = function(n,a,b,th,u,r,p){
  
  finv.kdagum = function(t,alpha,beta,theta,u,r,unif){Ft_kdagum(t,alpha,beta,theta,u,r)-unif}
  ger.kdagum = function(alpha,beta,theta,u,r,unif){t=uniroot(finv.kdagum, c(0,1000),tol=0.0000001,alpha=alpha,beta=beta,theta=theta,u=r,r=r,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0, Inf, ger.kdagum(alpha=a,beta=b,theta=th,r=r,u=u,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  
  t2 = pmin(t,u2) ; delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}


n=1000
a0=2
b0=1
th0=0.45
u0=0.5# ou 0.5, 2
r0=0.5# ou 0.5, 2

pbase = 1-th0; p0 = (1-(1-pbase)^r0)^u0; p0

dados.kdagum = gen.cure.kdagum(n=n,a=a0,b=b0,th=th0,u=u0,r=r0,p=p0)


## Verificando na curva de Kaplan-Meier 

survival_object = Surv(dados.kdagum[,1], dados.kdagum[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)
st_t_grid =st_kdagum(t=t_grid,alpha=a0,beta=b0,theta=th0, u=u0, r=r0)

lines(t_grid,st_t_grid, lwd=2, col = "deeppink")
abline(h=p0, lwd=2, col='steelblue')







## Estimação loglik via BFGS ===

params = c(a0,b0,th0,u0,r0)
chutes = params
#chutes = c(1,1.2,1,1,1)

maxlike = optim(par = chutes,
                fn = loglik_kdagum,
                hessian = TRUE,
                method = "BFGS",
                time = dados.kdagum[,1],
                delta = dados.kdagum[,2])

params
maxlike$par

sqrt(diag(solve(maxlike$hessian)))

maxlike$value














## Reparametrização da verossimilhança

## alpha1 = log(alpha); beta1 = log(beta)
## theta1 = exp(theta)/(1 + exp(theta))
## u1 = log(u); r1 = log(r) 

ft1_dagum = function(t,alpha,beta,theta){
  alpha1 = log(alpha); beta1 = log(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  den = (alpha1*beta1*theta1^2)*t^(-alpha1-1)/(beta1+theta1*t^(-alpha1))^2
  return(den)
}


Ft1_dagum = function(t,alpha,beta,theta){
  alpha1 = log(alpha); beta1 = log(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  dist = theta1*beta1/(beta1 + theta1*t^(-alpha1))
  return(dist)
}


st1_dagum = function(t,alpha,beta,theta){
  alpha1 = log(alpha); beta1 = log(beta)
  theta1 = exp(theta)/(1 + exp(theta))
  
  dist = theta1*beta1/(beta1 + theta1*t^(-alpha1))
  return(1-dist)
}

## 1.2 Kumaraswamy Dagum

ft1_kdagum = function(t,alpha,beta,theta, u, r){
  u1=exp(u);r1=exp(r)
  
  den = u1*r1*ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)*Ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^(r1-1)*(1-Ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r1)^(u1-1)
  return(den)
}

Ft1_kdagum = function(t,alpha,beta,theta, u, r){
  u1=exp(u);r1=exp(r)
  
  surv = (1 - Ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r1)^u1
  return(1-surv)
}

st1_kdagum = function(t,alpha,beta,theta, u, r){
  u1=exp(u);r1=exp(r)
  
  surv = (1 - Ft1_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r1)^u1
  return(surv)
}



## função log-verossimilhança

loglik1_kdagum = function(par, time,delta){
  
  # parametros
  alpha = par[1]
  beta = par[2]
  theta = par[3]
  u = par[4]
  r = par[5]
  
  densidade = ft1_kdagum(t=time, alpha=alpha,beta=beta,theta=theta,u=u,r=r)
  sobrevivencia = st1_kdagum(t=time, alpha=alpha,beta=beta,theta=theta,u=u,r=r)
  
  log_vero = sum(delta*log(densidade) + (1 - delta)*log(sobrevivencia))
  return(-log_vero)
}



# Algoritmo de geração de dados sob a F(t) Kdagum ===

n=10000
a0=2
b0=0.8
th0=0.55
u0=2# ou 0.5, 2
r0=2# ou 0.5, 2

pbase = 1-th0; p0 = (1-(1-pbase)^r0)^u0; p0

dados.kdagum = gen.cure.kdagum(n=n,a=a0,b=b0,th=th0,u=u0,r=r0,p=p0)


## Estimação loglik via BFGS ===
parametros = c(a0,b0,th0,u0,r0)
#chutes = rep(0.1,5)
chutes = parametros

maxlike1 = optim(par = chutes,
                fn = loglik1_kdagum,
                gr = NULL,
                hessian = TRUE,
                method = "BFGS",
                time = dados.kdagum[,1],
                delta = dados.kdagum[,2])


parametros
exp(maxlike1$par)




