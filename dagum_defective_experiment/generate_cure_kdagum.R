### ------
### Fase 2: Geração de dados com fração de cura dos modelos. 
###  Distribuição defectiva dagum
### ------

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(survival, survminer)

### ------ 
### 1. Funções base do modelo Dagum defeituoso ====
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
  den = u*r*ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)*Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^(r-1)*(1-Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r)^(u-1)
  return(den)
}

Ft_kdagum = function(t,alpha,beta,theta, u, r){
  surv = (1 - Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r)^u
  return(1-surv)
}

st_kdagum = function(t,alpha,beta,theta, u, r){
  surv = (1 - Ft_dagum(t=t,alpha=alpha,beta=beta,theta=theta)^r)^u
  return(surv)
}


### ------ 
### 2. Algoritmo Geração de dados com fração de cura ====
### ------

# 1.2: Algoritmo de geração de dados sob a F(t).

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
b0=4
th0=0.35
u0=0.5# ou 0.5, 2
r0=0.5# ou 0.5, 2

pbase = 1-th0; p0 = (1-(1-pbase)^r0)^u0; p0

dados.kdagum = gen.cure.kdagum(n=n,a=a0,b=b0,th=th0,u=u0,r=r0,p=p0)


## Verificando a curva de Kaplan-Meier

survival_object = Surv(dados.kdagum[,1], dados.kdagum[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.kdagum[,1])+10,length=200)
lines(t_grid,st_kdagum(t=t_grid,alpha=a0,beta=b0,theta=th0, u=u0, r=r0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="red")










