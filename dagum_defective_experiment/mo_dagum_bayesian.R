### ------
### Fase 4: Estimação bayesiana via MCMC
###  Distribuição defectiva dagum
### ------

if(!require(pacman)) install.packages("pacman"); library(pacman)
p_load(survival, survminer, rstan, R2jags)

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



n=2000
a0=2
b0=5
th0=0.7 ## Entre 0 e 1
l0=1.2

pbase = 1-th0; p0 = (l0*pbase)/(l0*pbase + 1 - pbase); p0

dados.modagum = gen.cure.modagum(n=n,a=a0,b=b0,th=th0,lb=l0,p=p0)

# Verificando na curva de Kaplan-Meier 

surv_obj = Surv(dados.modagum[,1], dados.modagum[,2])
km_fit = survfit(surv_obj  ~ 1)

par(mfrow=c(1,1))
plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)
st_modagum_t_grid = st_modagum(t=t_grid,alpha=a0,beta=b0,theta=th0,lambda=l0)


lines(t_grid,st_modagum_t_grid, lwd=2, col = "deeppink")
abline(h=p0, lwd=2, col='steelblue')



## Estimação via stan


cod_modagum = "
  data {
    int<lower=0> N;                         // tamanho amostral
    array[N] real time;                    //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1> theta;
  real<lower=0> lambda;
  }

  model {
  // Prioris
  alpha ~ gamma(0.1,0.1);
  beta ~ gamma(0.1,0.1);
  theta ~ beta(1,1);
  lambda ~ gamma(0.1,0.1);

  // Definição manual da log-verossimilhança
  for (i in 1:N){
  real denom = lambda + (((1-lambda)*theta*beta)/(beta + theta*time[i]^(-alpha)));

  target += delta[i]*log(lambda*alpha*beta*theta^2*time[i]^(-alpha-1)) -
             delta[i]*log((beta+theta*time[i]^(-alpha))^2 * denom^2) +
             (1-delta[i])*log(lambda*(1 - (theta*beta)/(beta+theta*time[i]^(-alpha)))) -
             (1-delta[i])*log(denom);
  }
}

"

## Transcrever o código escrito para um file stan 
writeLines(cod_modagum, con = "cod_modagum.stan")

#dados.kig
## Organizando os dados [data list]

data_cod_modagum = list(N = dim(dados.modagum)[1], 
                        time = dados.modagum[,1],
                        delta = dados.modagum[,2])


## Compilar e rodar o modelo
modagum = stan(file = 'cod_modagum.stan', data = data_cod_modagum, 
              chains = 1, iter = 2000, warmup = 200)



# a0;b0;th0;l0
# summary(modagum)$summary

modagum_post_samples = extract(modagum)

plot(modagum_post_samples$alpha, type='l')
abline(h=a0,col="red", lwd=2)

plot(modagum_post_samples$beta, type='l')
abline(h=b0,col="red", lwd=2)

plot(modagum_post_samples$theta, type='l')
abline(h=th0,col="red", lwd=2)

plot(modagum_post_samples$lambda, type='l')
abline(h=th0,col="red", lwd=2)

a0;b0;th0;l0
summary(modagum)$summary




fbaseest_modagum = 1-mean(modagum_post_samples$theta); 
fc_est_modgaum = (mean(modagum_post_samples$lambda)*fbaseest_modagum)/(mean(modagum_post_samples$lambda)*fbaseest_modagum + 1 - fbaseest_modagum)

fc_est_modgaum

t_grid = seq(0,100,by=0.01)

surv_est_modagum = st_modagum(t=t_grid,
                   alpha = mean(modagum_post_samples$alpha),
                   beta=mean(modagum_post_samples$beta),
                   theta = mean(modagum_post_samples$theta),
                   lambda = mean(modagum_post_samples$lambda))


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)

lines(t_grid,surv_est_modagum, lwd=2, col = "deeppink")
abline(h=fc_est_modgaum, lwd=2, col='steelblue')




