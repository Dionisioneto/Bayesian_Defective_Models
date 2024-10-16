### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 5. Modelo Kumaraswamy Gompertz;
### ------

# 1. Modelo Gompertz;
# 2. Modelo Gaussiano-inverso;
# 3. Modelo Marshall-Olkin Gompertz;
# 4. Modelo Marshall-Olkin Gaussiano-inverso;
# [5. Modelo Kumaraswamy Gompertz];
# 6. Modelo Kumaraswamy Gaussiano inverso.

# -------------------------------------------------------------------------------
library(pacman)
p_load(survival,rstan, R2jags,rjags)

# ---

### ------ 
### 5. Geração de dados do modelo Kumaraswamy Gompertz defeituoso ====
### para fração de cura da tese do Ricardo Ferreira
### ------

# 5.1: Funções Importantes.
# Obs.: We need to run the functions ft_Gompertz and St_Gompertz.

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

n=800
a0kg=-1.25
b0kg=1
psi0kg=2
k0kg=2

pgkg = exp(b0kg/a0kg)
p0kg = (1-(1-pgkg)^psi0kg)^k0kg; p0kg

dados.kgz=gen.cure.kgz(n=n,a=a0kg,b=b0kg,k=k0kg,ps=psi0kg,p=p0kg)
head(dados.kgz)

## Verificando na curva de Kaplan-Meier 

survival_object = Surv(dados.kgz[,1], dados.kgz[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)
st_t_grid = st_Kgompertz(t=t_grid,alpha=a0kg,
                         beta=b0kg,kappa=k0kg,psi=psi0kg)

lines(t_grid,st_t_grid, lwd=2, col = "deeppink")
abline(h=p0kg, lwd=2, col='steelblue')


# ---
# 5.3 Programação no stan  ====
# ---


## Portanto, no objeto dados.mog, 
## temos os dados do tempo de falha e o indicador de censura.

## Vamos atribuir uma priori para alpha, beta e lambda não informativas.
## Revisar o código, pois o amaostrador não está rodando nos valores iniciais

cod_kgtz_stan = "
  data{
  int<lower=0> N;
  array[N] real time;
  array[N] int<lower=0, upper=1> delta;
  }
  
  parameters{
  real alpha;
  real<lower=0> beta;
  real<lower=0> psi;
  real<lower=0> kappa;
  }
  
  model{
  
  // prioris
  
  alpha ~ normal(-1,10);
  beta ~ gamma(0.001,0.001);
  kappa ~ gamma(0.25,0.25);
  psi ~ gamma(0.25,0.25);
  
  // Definição manual da função de verossimilhança
  for(i in 1:N){
    target += delta[i]*(log(kappa*psi*beta*exp(alpha*time[i])*exp((beta-beta*exp(alpha*time[i]))/alpha)*(1-exp((beta-beta*exp(alpha*time[i]))/alpha))^(psi-1) * (1-(1-(exp((beta-beta*exp(alpha*time[i]))/alpha)))^psi)^(kappa-1))) +
              (1 - delta[i])*(log((1 - (1 - exp((beta-beta*exp(alpha*time[i]))/alpha))^psi)^kappa));
  }
}
"


## Transcrever o código escrito para um file stan 
writeLines(cod_kgtz_stan, con = "cod_kgtz_stan.stan")



## Organizando os dados [data list]
#dados.kgz

data_kgz = list(N = dim(dados.kgz)[1], 
                 time = dados.kgz[,1],
                 delta = dados.kgz[,2])


## Compilar e rodar o modelo
kgzfit = stan(file = 'cod_kgtz_stan.stan', data = data_kgz, 
               chains = 1, iter = 2000, warmup = 200)

a0kg;b0kg;psi0kg;k0kg
summary(kgzfit)$summary



kgzfit_post_samples = extract(kgzfit)


plot(kgzfit_post_samples$alpha, type='l')
abline(h=a0kg,col="red", lwd=2)

plot(kgzfit_post_samples$beta, type='l')
abline(h=b0kg,col="red", lwd=2)

plot(kgzfit_post_samples$psi, type='l')
abline(h=psi0kg,col="red", lwd=2)

plot(kgzfit_post_samples$kappa, type='l')
abline(h=psi0kg,col="red", lwd=2)

## estimativas pontuais
# mean_alpha = mean(mogfit_post_samples$alpha)
# mean_beta = mean(mogfit_post_samples$beta)
# mean_lambda= mean(mogfit_post_samples$lambda)








