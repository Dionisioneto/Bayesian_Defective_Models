### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 2. Modelo Gaussiano-inverso
### ------

# 1. Modelo Gompertz;
# [2. Modelo Gaussiano-inverso];
# 3. Modelo Marshall-Olkin Gompertz;
# 4. Modelo Marshall-Olkin Gaussiano-inverso;
# 5. Modelo Kumaraswamy Gompertz;
# 6. Modelo Kumaraswamy Gaussiano inverso.

# -------------------------------------------------------------------------------
library(pacman)
p_load(survival,rstan, R2jags)

# ---

# 3.1: Geração de dados


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

n = 800
a0ig = -1
b0ig = 4
p0ig = 1 - exp(2*a0ig/b0ig);p0ig

dados.IG = gen.cure.IG(n=n,a=a0ig,b=b0ig,p=p0ig)
#head(dados.IG)


# ---
# 3.2 Programação no stan  ====
# ---

## Portanto, no objeto dados.IG, temos os dados do tempo de falha e o indicador de censura.

## Vamos atribuir uma priori para alpha e beta não informativas.


# cod_IG_stan = "
#   data {
#     int<lower=0> N;                        // tamanho amostral
#     array[N] real time;                   //  tempo de falha observado
#     array[N] int<lower=0, upper=1> delta; //   indicador do evento
#   }
# 
#   parameters {
#     real alpha;
#     real<lower=0> beta;
#   }
# 
#   model {
#     target += -log(beta); // Priori para beta: p(beta) = 1/beta
# 
#     // Definição da log-verossimilhança manualmente
#     for (i in 1:N) {
#       target += delta[i] * log((1 / sqrt(2 * pi() * beta * time[i]^3)) *
#                 exp(-(1 / (2 * beta * time[i])) * pow(1 - alpha * time[i], 2))) +
#                 (1 - delta[i]) * log(1 - (Phi((alpha * time[i] - 1) / sqrt(beta * time[i])) +
#                 exp(2 * alpha / beta) * Phi((-alpha * time[i] - 1) / sqrt(beta * time[i]))));
#     }
#   }
# "

## Com o código feito, vamos estudar o estudo de prioris 
## alpha ~ Normal(mu0,sigma20) e beta~gamma(a2,b2), diferente
## da abordagem desenvolvida pelo Ricardo.

cod2_IG_stan = "
  data {
    int<lower=0> N;                        // tamanho amostral
    array[N] real time;                   //  tempo de falha observado
    array[N] int<lower=0, upper=1> delta; //   indicador do evento
  }

  parameters {
    real alpha;
    real<lower=0> beta;
  }

  model {
    // Prioris
    alpha ~ normal(0,100);
    beta ~ gamma(0.001,0.001);

    // Definição da log-verossimilhança manualmente
    for (i in 1:N) {
      target += delta[i] * log((1 / sqrt(2 * pi() * beta * time[i]^3)) *
                exp(-(1 / (2 * beta * time[i])) * pow(1 - alpha * time[i], 2))) +
                (1 - delta[i]) * log(1 - (Phi((alpha * time[i] - 1) / sqrt(beta * time[i])) +
                exp(2 * alpha / beta) * Phi((-alpha * time[i] - 1) / sqrt(beta * time[i]))));
    }
  }
"



## Transcrever o código escrito para um file stan 
writeLines(cod2_IG_stan, con = "cod2_IG_stan.stan")


## Organizando os dados [data list]

data_IG = list(N = dim(dados.IG)[1], 
              time = dados.IG[,1],
              delta = dados.IG[,2])


## Compilar e rodar o modelo
migfit = stan(file = 'cod2_IG_stan.stan', data = data_IG, 
            chains = 1, iter = 5000, warmup = 500)

a0ig;b0ig
migfit



migfit_post_samples = extract(migfit)


plot(migfit_post_samples$alpha, type='l')
abline(h=a0ig,col="red", lwd=2)

plot(migfit_post_samples$beta, type='l')
abline(h=b0ig,col="red", lwd=2)



## estimativas pontuais
mean_alphaig = mean(migfit_post_samples$alpha)
mean_betaig =mean(migfit_post_samples$beta)

## Verificando a curva de Kaplan-Meier
survival_object = Surv(dados.IG[,1], dados.IG[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


# Plotanto o modelo
t_grid = seq(0,30,by=0.01)
st_t_grid = St_IG(t=t_grid,alpha=mean_alphaig,beta=mean_betaig)
lines(t_grid,st_t_grid, lwd=2, col = "deeppink")

## Fração de cura calculada
fc_IG_stan = 1 - exp(2*mean_alphaig/mean_betaig)
abline(h=fc_IG_stan,lwd=2,col="purple")

## Para os modelo defeituoso inversa gaussiana as prioris 
## alpha ~ Normal(mu0,sigma20) e beta~gamma(a2,b2) estão funcionando.

# ---
# 3.3 Programação no Jags  ====
# ---

library(rjags)

## Código para o modelo

cod_IG_jags = "
model {
  pi = 3.141593
  for (i in 1:n) {
    logftig[i]=log(1/sqrt(2*pi*beta*pow(time[i], 3))*exp(-1/(2*beta*time[i])*pow(1 - alpha*time[i], 2)))
    logstig[i]=log(1 - (pnorm((alpha*time[i] - 1)/sqrt(beta*time[i]),0,1) + exp(2*alpha/beta) * pnorm((-alpha*time[i]-1)/sqrt(beta*time[i]),0,1)))
    
    # Definição da log-verossimilhança, utilizando o truque de zeros
    p[i] = 1000 - delta[i]*logftig[i] - (1-delta[i]) * logstig[i]
    zeros[i] ~ dpois(p[i])
  }
  # Priori não informativas para alpha e beta (assumindo independência)
  alpha ~ dunif(-10, 10)
  beta ~ dunif(0, 10)
}
"


## Especificação dos dados a serem utilizados

data.jags.ig = list(n=nrow(dados.IG), delta = dados.IG[,2], 
                    time = dados.IG[,1], zeros = rep(0,nrow(dados.IG)))


## nome dos parametros a serem salvos
p.ig.jags = c("alpha", "beta")
i.ig.jags = function(){list(alpha=rnorm(1),beta=rgamma(1,1,1))}

migtz.jags = jags.model(data=data.jags.ig, file = textConnection(cod_IG_jags),
                       n.chains = 2, inits = i.ig.jags)


## run the model for 1000 burn-in simulations
update(migtz.jags, 1000)

## run the model for 5000 aditional simulations to keep one in 10.
res.migtz.jags = coda.samples(migtz.jags, variable.names = p.ig.jags,
                            n.iter = 5000, n.thin=10)

summary(res.migtz.jags)
a0;b0

samples.migtz.jags = as.mcmc(do.call(rbind, res.migtz.jags))

## Análise gráfica de alpha
plot(samples.migtz.jags[,1])
a0

## Análise gráfica de beta
plot(samples.migtz.jags[,2])
b0






























