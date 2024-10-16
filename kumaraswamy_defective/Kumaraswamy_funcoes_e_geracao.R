## ------
## Geração de dados dos
## Modelos na família Kumaraswamy
## ------

### ------ 
### 1. Funções base e Geração do Modelo Kumaraswamy Gompertz====
### ------


## 1.1 Distribuição Gompertz 

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}


Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*beta)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}


## 1.2 Distribuição gussiana Inversa

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

## 1.3 Distribuição Kumaraswamy Gompertz

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


## 1.4 Distribuição Kumaraswamy Gaussiana Inversa

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


## Algoritmo de Geração de dados ====


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




