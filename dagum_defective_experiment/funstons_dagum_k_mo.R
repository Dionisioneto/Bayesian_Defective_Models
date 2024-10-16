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

## ---
## 2. Estudo do comportamento das distribuições
## ---

a0 = 0.5
b0 = 2
th0 = 0.4


time = seq(0,2,length=500)

# densidade e sobrevivência [Dagum]
plot(time,ft_dagum(t=time,alpha=a0,beta=b0,theta=th0), type='l')
plot(time,st_dagum(t=time,alpha=a0,beta=b0,theta=th0), type='l',
     ylim=c(0,1))


# densidade e sobrevivência [Dagum - Kumaraswamy]

a0 = 0.5
b0 = 2
th0 = 0.4
u0 = 1.2
r0 = 0.3

time = seq(0,2,length=500)


plot(time,ft_kdagum(t=time,alpha=a0,beta=b0,theta=th0,u=u0,r=r0), type='l')
plot(time,st_kdagum(t=time,alpha=a0,beta=b0,theta=th0,u=u0,r=r0), type='l',
     ylim=c(0,1))

# densidade e sobrevivência [Dagum - Marshall Olkin]
a0 = 0.5
b0 = 2
th0 = 0.4
l0 = 0.2

time = seq(0,2,length=500)

plot(time,ft_modagum(t=time,alpha=a0,beta=b0,theta=th0,lambda=l0), type='l')
plot(time,st_modagum(t=time,alpha=a0,beta=b0,theta=th0,lambda=l0), type='l',
     ylim=c(0,1))
