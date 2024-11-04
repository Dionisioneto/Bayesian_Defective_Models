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
a0=1.2
b0=0.5
th0=0.75
u0=0.5# ou 0.5, 2
r0=0.5# ou 0.5, 2

pbase = 1-th0; p0 = (1-(1-pbase)^r0)^u0; p0

dados.kdagum = gen.cure.kdagum(n=n,a=a0,b=b0,th=th0,u=u0,r=r0,p=p0)
#View(dados.kdagum)




## ---
## Estimação via stan
## ---

cod_kdagum = "
  data {
  int<lower=0> N;                         // tamanho amostral
    array[N] real time;                    //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0,upper=1> theta;
  real<lower=0> u;
  real<lower=0> r;
  }

  model {
  // Prioris
  alpha ~ gamma(0.25,0.25);
  beta ~ gamma(0.25,0.25);
  theta ~ beta(1,1);
  u ~ gamma(0.25,0.25);
  r ~ gamma(0.25,0.25);

  // Definição manual da log-verossimilhança
  for (i in 1:N){
  real Ft = (theta*beta)/(beta + theta*time[i]^(-alpha));

  target += delta[i]*log((u*r*alpha*beta*theta^2*time[i]^(-alpha-1))/(beta + theta*time[i]^(-alpha))^2) +
            delta[i]*log(Ft^(r-1)) + delta[i]*log((1 - Ft^r)^(u-1)) +
            (1-delta[i])*log((1-Ft^r)^u);
  }
}
"

#3 log-verossimilhança da K-Dagum

# Ft = (theta*beta)/(beta + theta*time[i]^(-alpha))
# 
# delta[i]*log((u*r*alpha*beta*theta^2*time[i]^(-alpha-1))/(beta + theta*time[i]^(-alpha))^2) +
#   delta[i]*log(Ft^(r-1)) + delta[i]*log((1 - Ft^r)^(u-1)) +
#   (1-delta[i])*log((1-Ft^r)^u)





## Transcrever o código escrito para um file stan
writeLines(cod_kdagum, con = "cod_kdagum.stan")

#dados.kig
## Organizando os dados [data list]

data_cod_kdagum = list(N = dim(dados.kdagum)[1],
                        time = dados.kdagum[,1],
                        delta = dados.kdagum[,2])

# init_values = list(
#   list(theta = 1, beta = 1, alpha = 1, r = 1, u = 1)
# )

## Compilar e rodar o modelo
mod_kdagum = stan(file = 'cod_kdagum.stan', data = data_cod_kdagum,
               chains = 1, iter = 2000, warmup = 300)
               
#,init = init_values)

a0;b0;th0;u0;r0
summary(mod_kdagum)$summary


kdagum_post_samples = extract(mod_kdagum)


plot(kdagum_post_samples$alpha, type='l')
abline(h=a0,col="red", lwd=2)

plot(kdagum_post_samples$beta, type='l')
abline(h=b0,col="red", lwd=2)

plot(kdagum_post_samples$theta, type='l')
abline(h=th0,col="red", lwd=2)

plot(kdagum_post_samples$u, type='l')
abline(h=u0,col="red", lwd=2)

plot(kdagum_post_samples$r, type='l')
abline(h=r0,col="red", lwd=2)



# 
# ## versão 2
# 
# cod2_kdagum = "
# data {
#   int<lower=0> N;                         // tamanho amostral
#   array[N] real time;                    //  tempo de falha observado
#   array[N] int<lower=0, upper=1> delta;  // indicador do evento
# }
# 
# parameters {
#   real alpha1;                            // log(alpha)
#   real beta1;                             // log(beta)
#   real theta1;                            // logit(theta), ou seja, log(theta / (1 - theta))
#   real u1;
#   real r1;
# }
# 
# transformed parameters {
#   real<lower=0> alpha = exp(alpha1);      // alpha = exp(alpha1)
#   real<lower=0> beta = exp(beta1);        // beta = exp(beta1)
#   real<lower=0,upper=1> theta = exp(theta1) / (1 + exp(theta1));  // theta = logistic(theta1)
#   real<lower=0> u = exp(u1);        // beta = exp(beta1)
#   real<lower=0> r = exp(r1);        // beta = exp(beta1)
# }
# 
# model {
#   // Prioris não-informativas (normais de grande variância)
#   alpha1 ~ normal(0, 10);                 // Prior não-informativa para alpha1
#   beta1 ~ normal(0, 10);                  // Prior não-informativa para beta1
#   theta1 ~ normal(0, 10);                 // Prior não-informativa para theta1
#   u1 ~ normal(0, 10);
#   r1 ~ normal(0, 10);
# 
#   // Definição manual da log-verossimilhança
#   for (i in 1:N) {
#     real F1 = (theta * beta / (beta + theta * pow(time[i], -alpha)));
# 
#     target += delta[i] * log(u * r * alpha * beta * theta^2 * (pow(time[i], -alpha-1) / pow(beta + theta * pow(time[i], -alpha), 2)) * pow(F1, r-1) * pow(1 - pow(F1, r), u-1))
#             + (1 - delta[i]) * log(pow(1 - F1, u));
#   }
# }
# "
# ## Transcrever o código escrito para um file stan 
# writeLines(cod2_kdagum, con = "cod2_kdagum.stan")
# 
# 
# # Exemplo de inicialização em R:
# init_values <- list(
#   list(alpha1 = 0.0, beta1 = 0.0, theta1 = 0.0, u = 1.0, r = 1.0)
# )
# 
# fit <- stan(file = 'cod2_kdagum.stan', data =data_cod_kdagum, init = init_values, iter = 2000, chains = 1)
# 
# 
# 
# 
# kdagum_post_samples = extract(fit)
# 
# 
# plot(kdagum_post_samples$alpha, type='l')
# abline(h=a0,col="red", lwd=2)
# 
# plot(kdagum_post_samples$beta, type='l')
# abline(h=b0,col="red", lwd=2)
# 
# plot(kdagum_post_samples$theta, type='l')
# abline(h=th0,col="red", lwd=2)
# 
# plot(kdagum_post_samples$u, type='l')
# abline(h=u0,col="red", lwd=2)
# 
# plot(kdagum_post_samples$r, type='l')
# abline(h=r0,col="red", lwd=2)
