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

# primeiro: Alpha = -1, Beta = 5, psi = 5, kappa = 0.5;
# primeiro: Alpha = -2, Beta = 10, psi = 1.2, kappa = 0.5;
# segundo: Alpha = -1, Beta = 1, psi = 2, kappa = 2.



n=800
a0kg=-1
b0kg=5
psi0kg=5
k0kg=0.5

pgkg = exp(b0kg/a0kg)
p0kg = (1-(1-pgkg)^psi0kg)^k0kg; p0kg

dados.kgz=gen.cure.kgz(n=n,a=a0kg,b=b0kg,k=k0kg,ps=psi0kg,p=p0kg)
#head(dados.kgz)

1-mean(dados.kgz[,2])

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
  vector[N] time;
  vector[N] delta;
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
  beta ~ gamma(0.1,0.1);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.1,0.1);

  // Definição manual da função de verossimilhança
  for(i in 1:N){
    target += delta[i]*(log(kappa*psi*beta*exp(alpha*time[i])*exp((beta-beta*exp(alpha*time[i]))/alpha)*(1-exp((beta-beta*exp(alpha*time[i]))/alpha))^(psi-1) * (1-(1-(exp((beta-beta*exp(alpha*time[i]))/alpha)))^psi)^(kappa-1))) +
              (1 - delta[i])*(log((1 - (1 - exp((beta-beta*exp(alpha*time[i]))/alpha))^psi)^kappa));
  }
}
"


# cod_kgtz_stan = "
#   data{
#   int<lower=0> N;
#   array[N] real time;
#   array[N] int<lower=0, upper=1> delta;
#   }
# 
#   parameters {
#       real alpha;
#       real log_beta;
#       real log_psi;
#       real log_kappa;
#     }
# 
#     transformed parameters{
#     real<lower=0> beta = exp(log_beta);
#     real<lower=0> psi = exp(log_psi);
#     real<lower=0> kappa = exp(log_kappa);
#     }
# 
#     model {
#     alpha ~ normal(-1,10);
#     log_beta ~ normal(0,10);
#     log_kappa ~ cauchy(0,10);
#     log_psi ~ cauchy(0,10);
#     
#      // Definição manual da função de verossimilhança
#   for(i in 1:N){
#     target += delta[i]*(log(kappa*psi*beta*exp(alpha*time[i])*exp((beta-beta*exp(alpha*time[i]))/alpha)*(1-exp((beta-beta*exp(alpha*time[i]))/alpha))^(psi-1) * (1-(1-(exp((beta-beta*exp(alpha*time[i]))/alpha)))^psi)^(kappa-1))) +
#               (1 - delta[i])*(log((1 - (1 - exp((beta-beta*exp(alpha*time[i]))/alpha))^psi)^kappa));
#     }
#   }
# 
# "


## Transcrever o código escrito para um file stan 
writeLines(cod_kgtz_stan, con = "cod_kgtz_stan.stan")


## Organizando os dados [data list]
#dados.kgz

data_kgz = list(N = dim(dados.kgz)[1], 
                 time = dados.kgz[,1],
                 delta = dados.kgz[,2])


## Compilar e rodar o modelo

# init_fun <- function() {
#   list(alpha = runif(1, -3, -2),       # No constraints on alpha
#        beta = runif(1, 5, 6),       # Positive constraint
#        kappa = runif(1, 0.4, 0.6),
#        psi = runif(1, 0.4, 0.6))     # Positive constraint
# }

# a0kg=-2
# b0kg=6
# psi0kg=0.5
# k0kg=0.5

kgzfit = stan(file = 'cod_kgtz_stan.stan', data = data_kgz, 
               chains = 1, iter = 2000, warmup = 200)

 
kgzfit_post_samples = extract(kgzfit)

par(mfrow=c(2,2))
plot(kgzfit_post_samples$alpha, type='l', ylab = "alpha")
abline(h=a0kg,col="red", lwd=2)

plot(kgzfit_post_samples$beta, type='l', ylab = "beta")
abline(h=b0kg,col="red", lwd=2)

plot(kgzfit_post_samples$psi, type='l', ylab = "psi")
abline(h=psi0kg,col="red", lwd=2)

plot(kgzfit_post_samples$kappa, type='l', ylab = "kappa")
abline(h=k0kg,col="red", lwd=2)

par(mfrow=c(1,1))

# acf(kgzfit_post_samples$alpha)
# acf(kgzfit_post_samples$beta)
# acf(kgzfit_post_samples$psi)
# acf(kgzfit_post_samples$kappa)

## estimativas pontuais
## para a configuração de prioris não informativas
## alpha ~ normal(-1,10);
## beta ~ gamma(0.25,0.25);
## kappa ~ gamma(0.25,0.25);
## psi ~ gamma(0.25,0.25);
## A estimação pontual pela mediana esyá bem melhor que a média.

mean_alpha = median(kgzfit_post_samples$alpha)
mean_beta = median(kgzfit_post_samples$beta)
mean_psi = median(kgzfit_post_samples$psi)
mean_kappa = median(kgzfit_post_samples$kappa)

a0kg;b0kg;psi0kg;k0kg
summary(kgzfit)$summary


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

median(kgzfit_post_samples$beta)
a0kg;b0kg;psi0kg;k0kg
mean_alpha;mean_beta;mean_kappa;mean_psi

t_grid = seq(0,100,by=0.01)
st_t_est = st_Kgompertz(t=t_grid,alpha=mean_alpha,
                         beta=mean_beta,kappa=mean_kappa,psi=mean_psi)

# st_t_est  = st_Kgompertz(t=t_grid,alpha=a0kg,
#                          beta=b0kg,kappa=k0kg,psi=psi0kg)

lines(t_grid,st_t_est, lwd=2, col = "deeppink")




fc_base_est = exp(mean_beta/mean_alpha)
fc_est = (1-(1-fc_base_est)^mean_psi)^mean_kappa

abline(h=fc_est, lwd=2, col='steelblue')
text(x = 8, y = fc_est- 0.05, 
     labels = bquote(hat(p) == .(round(fc_est, 4))))



# pgkg = exp(b0kg/a0kg)
# p0kg = (1-(1-pgkg)^psi0kg)^k0kg; p0kg






