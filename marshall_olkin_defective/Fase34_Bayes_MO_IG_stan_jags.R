### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 4. Modelo Marshall-Olkin Gaussiano-inverso;
### ------

# 1. Modelo Gompertz;
# 2. Modelo Gaussiano-inverso;
# 3. Modelo Marshall-Olkin Gompertz;
# [4. Modelo Marshall-Olkin Gaussiano-inverso];
# 5. Modelo Kumaraswamy Gompertz;
# 6. Modelo Kumaraswamy Gaussiano inverso.

## O stan não está conseguindo amostrar!!!
# -------------------------------------------------------------------------------
library(pacman)
p_load(survival,rstan, R2jags,rjags)

# ---


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
  ger.moig = function(alpha,beta,lambda,unif){t=uniroot(finv.moig,c(0,1000),tol=0.0001,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
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


n = 1000
a0moig=-2;b0moig=10;l0moig = 2

pgmoig= 1 - exp(2*a0moig/b0moig)
p0moig=(l0moig*pgmoig)/(l0moig*pgmoig+1-pgmoig); p0moig

dados.moig = gen.cure.moig(n=n,a=a0moig,b=b0moig,l=l0moig,p=p0moig)
head(dados.moig)


## Verificando na curva de Kaplan-Meier 

survival_object = Surv(dados.moig[,1], dados.moig[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)
st_t_grid = Stmo_IG(t=t_grid,alpha=a0moig,beta=b0moig,
                   lambda=l0moig)


lines(t_grid,st_t_grid, lwd=2, col = "deeppink")
abline(h=p0moig, lwd=2, col='steelblue')




# ---
# 4.3 Programação no stan  ====
# ---


## Portanto, no objeto dados.mog, 
## temos os dados do tempo de falha e o indicador de censura.

## Vamos atribuir uma priori para alpha, beta e lambda não informativas.
## Revisar o código, pois o amaostrador não está rodando nos valores iniciais

cod_moig_stan = "
  data {
   int<lower=0> N;                        // tamanho amostral
   vector[N] time;                        // tempo de falha observado
   array[N] int<lower=0, upper=1> delta;        // indicador do evento
  }
  
  parameters {
    real alpha;
    real<lower=0> beta;
    real<lower=0> lambda;
  }
  
  model {
    
    alpha~normal(-1,10);
    beta~gamma(0.25,0.25);
    lambda~gamma(0.25,0.25);
    
    // Definição manual da função de verossimilhança
    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real exp_term = exp(2 * alpha / beta);

    // Verossimilhança do evento observado e censurado
	
	  target += delta[i] * log(lambda*exp((-1/(2*beta*time[i]))*(1-alpha*time[i])^2)) -
              delta[i] * log(sqrt(2*pi()*beta*time[i]^3) * (lambda + (1-lambda)*(phi1 + exp(2*alpha/beta)*phi2))^2) +
              (1-delta[i])*log(lambda-lambda*(phi1 + exp(2*alpha/beta)*phi2)) -
              (1-delta[i])*log(lambda+(1-lambda)*(phi1 + exp(2*alpha/beta)*phi2));
    }
  }
"

# # real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
# # real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
# # real phi1 = Phi(z1);
# # real phi2 = Phi(z2);
# 
# delta[i] * log(lambda*exp((-1/(2*beta*time[i]))*(1-alpha*time[i])^2)) -
#   delta[i] * log(sqrt(2*pi()*beta*t^3) * (lambda + (1-lambda)*(phi1 + exp(2*alpha/beta)*phi2))^2) +
#   (1-delta[i])*log(lambda-lambda*(phi1 + exp(2*alpha/beta)*phi2)) -
#   (1-delta[i])*log(lambda+(1-lambda)*(phi1 + exp(2*alpha/beta)*phi2))

## Transcrever o código escrito para um file stan 
writeLines(cod_moig_stan, con = "cod_moig_stan.stan")


## Organizando os dados [data list]

data_moig = list(N = dim(dados.moig)[1], 
               time = dados.moig[,1],
               delta = dados.moig[,2])

## Definindo os chutes como uma lista
#init_values = list(alpha = a0moig, beta = b0moig, lambda = l0moig)

## Compilar e rodar o modelo
moigfit = stan(file = 'cod_moig_stan.stan', data = data_moig, 
              chains = 1, iter = 10000, warmup = 1000)
              
#, init = list(init_values))

a0moig;b0moig;l0moig
summary(moigfit)$summary



moigfit_post_samples = extract(moigfit)


plot(moigfit_post_samples$alpha, type='l')
abline(h=a0moig,col="red", lwd=2)

plot(moigfit_post_samples$beta, type='l')
abline(h=b0moig,col="red", lwd=2)

plot(moigfit_post_samples$lambda, type='l')
abline(h=l0moig,col="red", lwd=2)



# ---
# 4.4 Programação no jags  ====
# ---


## Código para o modelo

cod_moig_jags = "
model {
  for (i in 1:n) {
    # Definição da função de densidade de falha logaritmada
    logftmoig1[i] = log(lambda * exp(-1*pow(alpha * time[i] - 1, 2) / (2 * beta * time[i])))     
    logftmoig2[i] = log(sqrt(beta * pow(time[i], 3)) *
                        pow((lambda - 1) * pnorm((alpha * time[i] - 1) / sqrt(beta * time[i]), 0, 1) +
                            (lambda - 1) * exp(2 * alpha / beta) *
                            pnorm(-1*(alpha * time[i] + 1) / sqrt(beta * time[i]), 0, 1) - lambda, 2))

    # Definição da função de sobrevivência logaritmada
    logstmoig1[i] = log(lambda *
                        (pnorm((alpha * time[i] - 1) / sqrt(beta * time[i]), 0, 1) +
                         exp(2 * alpha / beta) *
                         pnorm(-1*(alpha * time[i] + 1) / sqrt(beta * time[i]), 0, 1) - 1)) 
    logstmoig2[i] = log(-1*(lambda - 1) * pnorm((alpha * time[i] - 1) / sqrt(beta * time[i]), 0, 1) +
                         (lambda - 1) * exp(2 * alpha / beta) *
                         (pnorm((alpha * time[i] + 1) / sqrt(beta * time[i]), 0, 1) + 1) - 1)

    # Definição da log-verossimilhança utilizando o truque de zeros
    p[i] = 10000 - delta[i] * (logftmoig1[i] - logftmoig2[i]) -
                     (1 - delta[i]) * (logstmoig1[i] - logstmoig2[i])
    zeros[i] ~ dpois(p[i])
  }

  # Priori não informativas para alpha, beta e lambda (assumindo independência)
  alpha ~ dunif(-10, 10)
  beta ~ dgamma(10, 1)
  lambda ~ dgamma(2, 1)
}

"

## Especificação dos dados a serem utilizados

data.jags.moig = list(n=nrow(dados.moig), delta = dados.moig[,2],
                      time = dados.moig[,1], zeros = rep(0,nrow(dados.moig)))


## nome dos parametros a serem salvos
p.moig.jags = c("alpha", "beta", "lambda")
i.moig.jags = function(){list(alpha=rnorm(1),beta=rgamma(1,1,1),lambda=rgamma(1,1,1))}

moig.jags = jags.model(data=data.jags.moig, file = textConnection(cod_moig_jags),
                        n.chains = 2, inits = i.moig.jags)




lambda=l0moig;alpha=a0moig;beta=b0moig
timecens = dados.moig[which(dados.moig[,2]==0),1]


## Revisar o código, as funções log estão retornando valores negativos!!

log(-1*(lambda - 1) * pnorm((alpha * timecens - 1) / sqrt(beta * timecens), 0, 1) +
      (lambda - 1) * exp(2 * alpha / beta) *
      (pnorm((alpha * timecens + 1) / sqrt(beta * timecens), 0, 1) + 1) - 1)



log(lambda *
      (pnorm((alpha * log(lambda *
                        (pnorm((alpha * timecens - 1) / sqrt(beta * timecens), 0, 1) +
                         exp(2 * alpha / beta) *
                         pnorm(-1*(alpha * timecens + 1) / sqrt(beta * timecens), 0, 1) - 1)) - 1) / sqrt(beta * timecens), 0, 1) +
         exp(2 * alpha / beta) *
         pnorm(-1*(alpha * timecens + 1) / sqrt(beta * timecens), 0, 1) - 1))




