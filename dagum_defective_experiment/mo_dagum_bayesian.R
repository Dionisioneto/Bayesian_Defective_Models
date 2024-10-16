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



n=1000
a0=0.8
b0=1.5
th0=0.45
l0=0.5

pbase = 1-th0; p0 = (l0*pbase)/(l0*pbase + 1 - pbase); p0

dados.modagum = gen.cure.modagum(n=n,a=a0,b=b0,th=th0,lb=l0,p=p0)
#View(dados.modagum)


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
  alpha ~ gamma(0.01,0.01);
  beta ~ gamma(0.25,0.25);
  theta ~ beta(1,1);
  lambda ~ gamma(0.25,0.25);

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

## log-verossimilhanca - Marshall-Olkin Dagum

# denom = lambda + (((1-lambda)*theta*beta)/(beta + theta*time[i]^(-alpha)))
# 
# delta[i]*log(lambda*alpha*beta*theta^2*time^(-alpha-1)) -
#   delta[i]*log((beta+theta*time[i]^(-alpha))^2 * denom^2) +
#   (1-delta[i])*log(lambda*(1 - (theta*beta)/(beta+theta*t^(-alpha)))) -
#   (1-delta[i])*log(denom)




# cog_modagum2 = "
# 
# data {
#   int<lower=0> N;                         // Tamanho amostral
#   array[N] real time;                      // Tempo de falha observado
#   array[N] int<lower=0, upper=1> delta;    // Indicador do evento
# }
# 
# parameters {
#   real alpha0;              // alpha0 = log(alpha)
#   real beta0;               // beta0 = log(beta)
#   real theta0;              // theta0 = log(theta/(1-theta))
#   real<lower=0> lambda0;     // lambda0 = log(lambda)
# }
# 
# transformed parameters {
#   real<lower=0> alpha = exp(alpha0);                // Inversa da transformação logarítmica
#   real<lower=0> beta = exp(beta0);                  // Inversa da transformação logarítmica
#   real<lower=0, upper=1> theta = inv_logit(theta0); // Inversa da transformação logístico (para restringir a [0,1])
#   real<lower=0> lambda = exp(lambda0)
# }
# 
# model {
#   // Prioris no espaço transformado
#   alpha0 ~ normal(0, 10);  // Log(alpha)
#   beta0 ~ normal(0, 10);   // Log(beta)
#   theta0 ~ normal(0, 10);  // Logit(theta) (transformação de Beta(1,1) para logit normal)
#   lambda0 ~ normal(0, 10);
# 
#   // Definição manual da log-verossimilhança
#   for (i in 1:N) {
#     real F1 = 1 - (theta * beta / (beta + theta * time[i]^(-alpha)));
#     
#     target += delta[i] * log((lambda * alpha * beta * theta^2 * time[i]^(-alpha-1)) / (beta + theta * time[i]^(-alpha))^2 * (1 - (1 - lambda) * F1)^(-1)) 
#               + (1 - delta[i]) * log(lambda * F1 / (1 - (1 - lambda) * F1));
#   }
# }
# 
# "

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



a0;b0;th0;l0
summary(modagum)$summary


modagum_post_samples = extract(modagum)


plot(modagum_post_samples$alpha, type='l')
abline(h=a0,col="red", lwd=2)

plot(modagum_post_samples$beta, type='l')
abline(h=b0,col="red", lwd=2)

plot(modagum_post_samples$theta, type='l')
abline(h=th0,col="red", lwd=2)

plot(modagum_post_samples$lambda, type='l')
abline(h=th0,col="red", lwd=2)


