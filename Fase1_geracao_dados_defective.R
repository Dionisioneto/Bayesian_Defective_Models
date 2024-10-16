### ------
### Fase 1: Geração de dados com fração de cura dos modelos. 
###  1. () 2. () 3. () 4. () 5. () 6. ()
### ------

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(survival, survminer)

### ------ 
### 1. Geração de dados do modelo Gompertz defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------


# 1.1: Funções Importantes.

# t: Failure time. 
# alpha, beta: Shape parameters.

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*beta)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}


Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}

# 1.2: Algoritmo de geração de dados sob a F(t).

gen.cure.gompertz = function(n,a,b,p){
  rm = rbinom(n=n,size=1,prob=1-p)
  
  t=rep(NA,n)
  
  for(i in 1:n){
    t[i]=ifelse(rm[i]==0, Inf, 
                log((-(alpha/beta)*log(1-runif(n=1,min=0,max=1-p))) + 1)*(1/alpha))
    
  }
  
  t_finite = ifelse(t==Inf,0,t)
  
  u2 = runif(n=n,0,max(t_finite))
  
  t2 = pmin(t,u2) ; delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
} 

## Determine some values
n = 10000

alpha=-1
beta=1
p = exp(beta/alpha)

data.gompertz = gen.cure.gompertz(n=n,a=alpha,b=beta,p=p)
colnames(data.gompertz) = c("tempo", "delta")

## Verificando na curva de Kaplan-Meier 

survival_object = Surv(data.gompertz[,1], data.gompertz[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,10,by=0.01)
st_t_grid = St_Gompertz(t=t_grid,alpha=alpha,beta=beta)

lines(t_grid,st_t_grid, lwd=2, col = "deeppink")



### ------ 
### 2. Geração de dados do modelo Inversa Gaussiana defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------


# 2.1: Funções Importantes.

# t: Failure time. 
# alpha, beta: Shape parameters.

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

# estudo da distribuição
# alpha = -1
# beta = 4
# p0 = 1 - exp(2*alpha/beta);p0
# t_grid = seq(0,6,length=100)

#plot(t_grid,ft_IG(t=t_grid,alpha=alpha,beta=beta),type='l',lwd=2)
#plot(t_grid,Ft_IG(t=t_grid,alpha=alpha,beta=beta),type='l',lwd=2)
#plot(t_grid,St_IG(t=t_grid,alpha=alpha,beta=beta),type='l',lwd=2,ylim=c(0,1))
#abline(h=p0,lwd=2,col="deeppink")


# 2.2: Algoritmo de geração de dados sob a F(t).
# Não foi possível obter a forma fechada para a t = F^-1(u).
# Temos que utilizar o uniroort para encontrar a raiz unitária!


gen.cure.IG = function(n,a,b,p){
  
  finv.IG = function(t,alpha,beta,unif){Ft_IG(t=t,alpha=alpha,beta=beta)-unif}
  ger.IG = function(alpha,beta,unif){t=uniroot(finv.IG, c(0,1000),tol=0.01,alpha=alpha,beta=beta,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0, Inf, ger.IG(alpha=a,beta=b,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  
  t2 = pmin(t,u2) ; delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}

n=500
a0 = -1
b0 = 4
p0 = 1 - exp(2*a0/b0);p0

dados.IG = gen.cure.IG(n=n,a=a0,b=b0,p=p0)
head(dados.IG)

## Verificando a curva de Kaplan-Meier

survival_object = Surv(dados.IG[,1], dados.IG[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.IG[,1])+10,length=100)
lines(t_grid,St_IG(t=t_grid,alpha=a0,beta=b0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="red")

### ------ 
### 3. Geração de dados do modelo Marshall-Olkin Gompertz defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------

# 3.1: Funções Importantes.
# Obs.: We need to run the functions St_Gompertz and ft_Gompertz before.

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*beta)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}
Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}
St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}

# t: Failure time. 
# alpha, beta: Shape parameters.
# lambda: MO parameter.

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

## Estudo gráfico

a0= -2;b0 = 2;l0 = 2

t_grid = seq(0,3,length=1000)

## densidade de probabilidade
plot(t_grid,ftmo_gompertz(t=t_grid,alpha=a0,beta=b0,lambda=l0), 
     type = 'l', lwd=2)

pg=exp(b0/a0)
p0=(l0*pg)/(l0*pg+1-pg)


plot(t_grid,Stmo_gompertz(t=t_grid,alpha=a0,beta=b0,lambda=l0), 
     type = 'l', lwd=2, ylim=c(0,1))

abline(h=p0,lwd=2,col="deeppink")

# 3.2: Algoritmo de geração de dados sob a F(t).
# Não foi possível obter a forma fechada para a t = F^-1(u).
# Temos que utilizar o uniroort para encontrar a raiz unitária!


gen.cure.mog = function(n,a,b,l,p){
  finv.mog = function(t,alpha,beta,lambda,unif){Ftmo_gompertz(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.mog = function(alpha,beta,lambda,unif){t=uniroot(finv.mog,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0, Inf,  ger.mog(alpha=a,beta=b,lambda=l,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}

n = 300
a0= -2;b0 = 2;l0 = 2

pg=exp(b0/a0)
p0=(l0*pg)/(l0*pg+1-pg); p0

dados.mog = gen.cure.mog(n=n,a=a0,b=b0,l=l0,p=p0)
head(dados.mog)

## Verificando a curva de Kaplan-Meier

survival_object = Surv(dados.mog[,1], dados.mog[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.mog[,1])+10,length=100)
lines(t_grid,Stmo_gompertz(t=t_grid,alpha=a0,beta=b0,lambda=l0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="red")


### ------ 
### 4. Geração de dados do modelo Marshall-Olkin Gaussiano Inverso defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------

# 4.1: Funções Importantes.
# Obs.: We need to run the functions ft_IG and St_IG before.

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


# 4.2: Algoritmo de geração de dados sob a F(t).
# Não foi possível obter a forma fechada para a t = F^-1(u).
# Temos que utilizar o uniroot para encontrar a raiz unitária!

gen.cure.moig = function(n,a,b,l,p){
  finv.moig = function(t,alpha,beta,lambda,unif){Ftmo_IG(t=t,alpha=alpha,beta=beta,lambda=lambda)-unif}
  ger.moig = function(alpha,beta,lambda,unif){t=uniroot(finv.moig,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i] = ifelse(rm[i]==0, Inf,ger.moig(alpha=a,beta=b,lambda=l,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}


n = 300
a0= -2;b0 =10;l0 = 2

pg= 1 - exp(2*a0/b0)
p0=(l0*pg)/(l0*pg+1-pg); p0

dados.moig = gen.cure.moig(n=n,a=a0,b=b0,l=l0,p=p0)
head(dados.mog)


## Verificando a curva de Kaplan-Meier
survival_object = Surv(dados.moig[,1], dados.moig[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.moig[,1])+10,length=100)
lines(t_grid,Stmo_IG(t=t_grid,alpha=a0,beta=b0,lambda=l0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="deeppink")



### ------ 
### 5. Geração de dados do modelo Kumaraswamy Gompertz defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------

# 5.1: Funções Importantes.
# Obs.: We need to run the functions ft_Gompertz and St_Gompertz.

# t: Failure time. 
# alpha, beta: Shape parameters.

ft_Kgompertz = function(t,alpha,beta,kappa,psi){
  den1 = kappa*psi*ft_Gompertz(t=t,alpha=alpha,beta=beta)
  den2 = Ft_Gompertz(t=t,alpha=alpha,beta=beta)^(psi-1)*(1-Ft_Gompertz(t=t,alpha=alpha,beta=beta)^(psi))^(kappa-1)
  return(den1*den2)
}

Ft_Kgompertz = function(t,alpha,beta,kappa,psi){
  surv = (1 - Ft_Gompertz(t=t,alpha=alpha,beta=beta)^(psi))^(kappa)
  return(1-surv)
}

st_Kgompertz = function(t,alpha,beta,kappa,psi){
  surv = (1 - Ft_Gompertz(t=t,alpha=alpha,beta=beta)^(psi))^(kappa)
  return(surv)
}

a0=-2
b0=2
psi0=2
k0=2

t_grid = seq(0,2.5,length=1000)
st_Kg = st_Kgompertz(t=t_grid,alpha=a0,beta=b0,kappa=k0,psi=psi0)

plot(t_grid,st_Kg, lwd=2,type='l', ylim=c(0,1))

pb = exp(b0/a0)
p0 = (1-(1-pb)^psi0)^k0; p0

abline(h=p0, lwd=2, col = "purple")


# 5.2: Algoritmo de geração de dados sob a F(t).
# Temos a forma fechada para a t = F^-1(u)!


gen.cure.kgz = function(n,a,b,k,ps,p){
  finv.kgz = function(t,alpha,beta,kappa,psi,unif){Ft_Kgompertz(t=t,alpha=alpha,beta=beta,kappa=kappa,psi=psi)-unif}
  ger.kgz = function(alpha,beta,kappa,psi,unif){t=uniroot(finv.kgz,c(0,1000),tol=0.01,alpha=alpha,beta=beta,kappa=kappa,psi=psi,unif=unif);return(t$root)}
  
  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i]=ifelse(rm[i]==0, Inf, ger.kgz(alpha=a,beta=b,kappa=k,psi=ps,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}

n=300
a0=-1
b0=1
psi0=2
k0=2

pb = exp(b0/a0)
p0 = (1-(1-pb)^psi0)^k0; p0


dados.kgz=gen.cure.kgz(n=n,a=a0,b=b0,k=k0,ps=psi0,p=p0)
head(dados.kgz)

## Verificando a curva de Kaplan-Meier
survival_object = Surv(dados.kgz[,1], dados.kgz[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.kgz[,1])+10,length=100)
lines(t_grid,st_Kgompertz(t=t_grid,alpha=a0,beta=b0,psi=psi0,kappa=k0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="deeppink")


### ------ 
### 6. Geração de dados do modelo Kumaraswamy Gaussiano Inverso defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------

# 6.1: Funções Importantes.
# Obs.: We need to run the functions ft_IG and St_IG before.

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

st_IG = function(t,alpha,beta){
  arg1 = (alpha*t-1)/sqrt(beta*t)
  arg2 = (-alpha*t-1)/sqrt(beta*t)
  Ft = pnorm(q=arg1,mean=0,sd=1) + (exp(2*alpha/beta)*pnorm(q=arg2,mean=0,sd=1))
  return(1-Ft)
}

# t: Failure time. 
# alpha, beta: Shape parameters.

ft_KIG = function(t,alpha,beta,kappa,psi){
  den1 = kappa*psi*ft_IG(t=t,alpha=alpha,beta=beta)
  den2 = (Ft_IG(t=t,alpha=alpha,beta=beta))^(psi-1)*(1-Ft_IG(t=t,alpha=alpha,beta=beta)^(psi))^(kappa-1)
  
  return(den1*den2)
}

Ft_KIG = function(t,alpha,beta,kappa,psi){
  surv = (1 - Ft_IG(t=t,alpha=alpha,beta=beta)^(psi))^(kappa)
  return(1-surv)
}

st_KIG = function(t,alpha,beta,kappa,psi){
  surv = (1 - Ft_IG(t=t,alpha=alpha,beta=beta)^(psi))^(kappa)
  return(surv)
}


a0=-0.5
b0=3
psi0=0.8
k0=0.8

t_grid = seq(0,20,length=1000)
st_kig = st_KIG(t=t_grid,alpha=a0,beta=b0,kappa=k0,psi=psi0)

plot(t_grid,st_kig, lwd=2,type='l', ylim=c(0,1))

pb = 1 - exp(2*a0/b0)
p0 = (1-(1-pb)^psi0)^k0; p0

abline(h=p0, lwd=2, col = "purple")


# 6.2: Algoritmo de geração de dados sob a F(t).
# Temos a forma fechada para a t = F^-1(u)!


gen.cure.kig = function(n,a,b,k,ps,p){
  finv.kig = function(t,alpha,beta,kappa,psi,unif){Ft_KIG(t=t,alpha=alpha,beta=beta,kappa=kappa,psi=psi)-unif}
  ger.kig = function(alpha,beta,kappa,psi,unif){t=uniroot(finv.kig,c(0,1000),tol=0.001,alpha=alpha,beta=beta,kappa=kappa,psi=psi,unif=unif);return(t$root)}

  rm = rbinom(n=n,size=1,prob=1-p)
  un = runif(n=n,0,max=1-p)
  t = rep(NA,n)
  
  for(i in 1:n){
    t[i]=ifelse(rm[i]==0, Inf, ger.kig(alpha=a,beta=b,kappa=k,psi=ps,unif=un[i]))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
}


n=300
a0=-0.5
b0=3
psi0=0.8
k0=0.8

pb = 1 - exp(2*a0/b0)
p0 = (1-(1-pb)^psi0)^k0; p0

dados.kig=gen.cure.kig(n=n,a=a0,b=b0,k=k0,ps=psi0,p=p0)
head(dados.kig)

## Verificando a curva de Kaplan-Meier
survival_object = Surv(dados.kig[,1], dados.kig[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.kig[,1])+10,length=100)
lines(t_grid,st_KIG(t=t_grid,alpha=a0,beta=b0,kappa=k0,psi=psi0),lwd=2,col="steelblue")
abline(h=p0,lwd=2,col="deeppink")






