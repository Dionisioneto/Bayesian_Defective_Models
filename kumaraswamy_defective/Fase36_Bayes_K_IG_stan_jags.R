### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 6. Modelo Kumaraswamy Gaussiano inverso.
### ------

# 1. Modelo Gompertz;
# 2. Modelo Gaussiano-inverso;
# 3. Modelo Marshall-Olkin Gompertz;
# 4. Modelo Marshall-Olkin Gaussiano-inverso;
# 5. Modelo Kumaraswamy Gompertz;
# [6. Modelo Kumaraswamy Gaussiano inverso].

# -------------------------------------------------------------------------------
library(pacman)
p_load(survival,rstan, R2jags,rjags)

# ---
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


# primeiro: Alpha = -0.2, Beta = 2, psi = 2.2, kappa = 2;
# segundo: Alpha = -1, Beta = 10, psi = 1.4, kappa = 0.6.

n=800
a0kig=-0.2
b0kig=2
psi0kig=2.2
k0kig=2


pbkig = 1 - exp(2*a0kig/b0kig)
p0kig = (1-(1-pbkig)^psi0kig)^k0kig; p0kig

dados.kig = gen.cure.kig(n=n,a=a0kig,b=b0kig,k=k0kig,ps=psi0kig,p=p0kig)

1-mean(dados.kig[,2])

## Verificando na curva de Kaplan-Meier 

survival_object = Surv(dados.kig[,1], dados.kig[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,400,by=0.01)
st_t_grid = st_KIG(t=t_grid,alpha=a0kig,beta=b0kig,
                   kappa=k0kig,psi=psi0kig)

lines(t_grid,st_t_grid, lwd=2, col = "deeppink")
abline(h=p0kig, lwd=2, col='steelblue')




# ---
# 6.2 Programação no stan  ====
# ---

## Portanto, no objeto dados.IG, temos os dados do tempo de falha e o indicador de censura.

# Com o código feito, vamos estudar o estudo de prioris
# alpha ~ normal(-1,10);
# beta ~ gamma(0.25,0.25);
# kappa ~ gamma(0.25,0.25);
# psi ~ gamma(0.25,0.25);

cod_KIG_stan = "
  data {
    int<lower=0> N;                        // tamanho amostral
    vector[N] time;                        // tempo de falha observado
     array[N] int<lower=0, upper=1> delta;        // indicador do evento
  }

  parameters {
    real alpha;
    real<lower=0> beta;
    real<lower=0> psi;
    real<lower=0> kappa;
  }

  model {
    // prioris

  alpha ~ normal(-1,10);
  beta ~ gamma(0.1,0.1);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.01,0.1);

    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real exp_term = exp(2 * alpha / beta);

      // Verossimilhança do evento observado
      target += delta[i] * (log(kappa) + log(psi) - 0.5 * log(2 * beta * pi() * time[i]^3) - (1 - alpha * time[i])^2 / (2 * beta * time[i])
              + (psi - 1) * log(phi1 + exp_term * phi2));

      // Verossimilhança para censura
      target += delta[i] * (kappa - 1) * log(1 - (phi1 + exp_term * phi2)^psi);
      target += (1 - delta[i]) * kappa * log(1 - (phi1 + exp_term * phi2)^psi);
    }
  }
"

## versão considerando os parametros na reta
# cod_KIG_stan = "
#   data {
#     int<lower=0> N;                        // tamanho amostral
#     vector[N] time;                        // tempo de falha observado
#      array[N] int<lower=0, upper=1> delta;        // indicador do evento
#   }
# 
#   parameters {
#     real alpha;
#     real log_beta;
#     real log_psi;
#     real log_kappa;
#   }
#   
#   transformed parameters{
#   real<lower=0> beta = exp(log_beta);  
#   real<lower=0> psi = exp(log_psi);
#   real<lower=0> kappa = exp(log_kappa);
#   }
# 
#   model {
#   alpha ~ normal(-1,4);
#   log_beta ~ normal(0,4);
#   log_kappa ~ normal(0,4);
#   log_psi ~ normal(0,4);
#     
#     for (i in 1:N) {
#       real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
#       real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
#       real phi1 = Phi(z1);
#       real phi2 = Phi(z2);
#       real exp_term = exp(2 * alpha / beta);
# 
#       // Verossimilhança do evento observado
#       target += delta[i] * (log(kappa) + log(psi) - 0.5 * log(2 * beta * pi() * time[i]^3) - (1 - alpha * time[i])^2 / (2 * beta * time[i])
#               + (psi - 1) * log(phi1 + exp_term * phi2));
# 
#       // Verossimilhança para censura
#       target += delta[i] * (kappa - 1) * log(1 - (phi1 + exp_term * phi2)^psi);
#       target += (1 - delta[i]) * kappa * log(1 - (phi1 + exp_term * phi2)^psi);
#     }
#   }
# "



## Transcrever o código escrito para um file stan 
writeLines(cod_KIG_stan, con = "cod_KIG_stan.stan")

#dados.kig
## Organizando os dados [data list]

data_kig = list(N = dim(dados.kig)[1], 
               time = dados.kig[,1],
               delta = dados.kig[,2])

## chutes iniciais
#init_vals = list(list(alpha = -0.5, beta = 0.5, kappa = 0.5, psi = 0.5))

## Compilar e rodar o modelo
kigfit = stan(file = 'cod_KIG_stan.stan', data = data_kig, 
              chains = 1, iter = 2000, warmup = 200)

a0kig;b0kig;psi0kig;k0kig
summary(kigfit)$summary


kigfit_post_samples = extract(kigfit)

par(mfrow=c(2,2))
plot(kigfit_post_samples$alpha, type='l', ylab="Alpha")
abline(h=a0kig,col="red", lwd=2)

plot(kigfit_post_samples$beta, type='l', ylab="Beta")
abline(h=b0kig,col="red", lwd=2)

plot(kigfit_post_samples$psi, type='l', ylab="Psi")
abline(h=psi0kig,col="red", lwd=2)

plot(kigfit_post_samples$kappa, type='l', ylab="Kappa")
abline(h=k0kig,col="red", lwd=2)

par(mfrow=c(1,1))


## verificar as correlações
# acf(kigfit_post_samples$alpha)
# acf(kigfit_post_samples$beta)
# acf(kigfit_post_samples$psi)
# acf(kigfit_post_samples$kappa)


survival_object = Surv(dados.kig[,1], dados.kig[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,400,by=0.01)
st_t_est = st_KIG(t=t_grid,
                   alpha=mean(kigfit_post_samples$alpha),
                   beta=mean(kigfit_post_samples$beta),
                   kappa=mean(kigfit_post_samples$kappa),
                   psi=mean(kigfit_post_samples$psi))

pbkig_est = 1 - exp(2*mean(kigfit_post_samples$alpha)/mean(kigfit_post_samples$beta))
p0kig_est = (1-(1-pbkig_est)^mean(kigfit_post_samples$psi))^mean(kigfit_post_samples$kappa);


lines(t_grid,st_t_est, lwd=2, col = "deeppink")

abline(h=p0kig_est, lwd=2, col='steelblue')
text(x = 40, y = p0kig_est- 0.05, 
     labels = bquote(hat(p) == .(round(p0kig_est, 4))))



plot(dgamma(x=seq(0,100,0.1), shape = 0.25, rate = 0.25),type='l')












