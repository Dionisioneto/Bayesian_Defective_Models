### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 3. Modelo Marshall-Olkin Gompertz;
### ------

# 1. Modelo Gompertz;
# 2. Modelo Gaussiano-inverso;
# [3. Modelo Marshall-Olkin Gompertz];
# 4. Modelo Marshall-Olkin Gaussiano-inverso;
# 5. Modelo Kumaraswamy Gompertz;
# 6. Modelo Kumaraswamy Gaussiano inverso.

# -------------------------------------------------------------------------------
library(pacman)
p_load(survival,rstan, R2jags)

# ---
# 3.1 Geração de Dados  ====
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

# 3.1: Geração de dados

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

n=1000
a0mog=-1.2;b0mog=2;l0mog=0.8
pbmog=exp(b0mog/a0mog); p0mog=(l0mog*pbmog)/(l0mog*pbmog+1-pbmog)
p0mog

dados.mog = gen.cure.mog(n=n,a=a0mog,b=b0mog,l=l0mog,p=p0mog)

prop.table(table(dados.mog[,2]))

## Verificando na curva de Kaplan-Meier 

survival_object = Surv(dados.mog[,1], dados.mog[,2])
km_fit = survfit(survival_object ~ 1)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)


t_grid = seq(0,100,by=0.01)
st_t_grid = Stmo_gompertz(t=t_grid,alpha=a0mog,beta=b0mog,
                    lambda=l0mog)


lines(t_grid,st_t_grid, lwd=2, col = "deeppink")
abline(h=p0mog, lwd=2, col='steelblue')





# ---
# 3.2 Programação no stan  ====
# ---


## Portanto, no objeto data.gompertz, 
## temos os dados do tempo de falha e o indicador de censura.

## Vamos atribuir uma priori para alpha, beta e lambda não informativas.
## Observacao: As prioris estão afetando a verossimilhanca!

## Com o código feito, vamos estudar o estudo de prioris 
## alpha ~ Normal(mu0,sigma20); beta~gamma(a2,b2); lambda~HalfCauchy(0,25)

cod_mog_stan = "

data {
  int<lower=0> N;                        
  array[N] real time;
  array[N] int<lower=0, upper=1> delta; 
}
  
parameters {
  real alpha;
  real<lower=0> beta;
  real<lower=0> lambda;
}
  
model {
  // Prioris
  alpha ~ normal(-1,10);
  beta ~ gamma(0.25,0.25);
  lambda ~ gamma(0.25,0.25);


  for (i in 1:N) {
    // Calculate exp_term for each i
    real exp_term = exp((-beta / alpha) * (exp(alpha * time[i]) - 1));

    target += delta[i] * log(lambda * beta * exp(alpha * time[i]) * exp_term) -
              delta[i] * log((1 - (1 - lambda) * exp_term)^2) +
              (1 - delta[i]) * log(lambda * exp_term) -
              (1 - delta[i]) * log(1 - (1 - lambda) * exp_term);
  }
}
"


## Transcrever o código escrito para um file stan 
writeLines(cod_mog_stan, con = "cod_mog_stan.stan")

## Organizando os dados [data list]


data_mog = list(N = dim(dados.mog)[1],
                time = dados.mog[,1],
                delta = dados.mog[,2])


## Compilar e rodar o modelo

## Compilar e rodar o modelo
mogfit = stan(file = 'cod_mog_stan.stan', data = data_mog, 
            chains = 1, iter = 10000, warmup = 1000)

a0mog;b0mog;l0mog
mogfit


mogfit_post_samples = extract(mogfit)


plot(mogfit_post_samples$alpha, type='l')
abline(h=a0mog,col="red", lwd=2)

plot(mogfit_post_samples$beta, type='l')
abline(h=b0mog,col="red", lwd=2)

plot(mogfit_post_samples$lambda, type='l')
abline(h=l0mog,col="red", lwd=2)


## estimativas pontuais
mean_alphamog = mean(mogfit_post_samples$alpha)
mean_betamog = mean(mogfit_post_samples$beta)
mean_lambdamog= mean(mogfit_post_samples$lambda)

## Verificando a curva de Kaplan-Meier
survival_object = Surv(dados.mog[,1], dados.mog[,2])
km_fit = survfit(survival_object ~ 1)

plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

t_grid = seq(0,max(dados.mog[,1])+10,length=100)
lines(t_grid,
      Stmo_gompertz(t=t_grid,alpha=mean_alphamog,beta=mean_betamog,lambda=mean_lambdamog),lwd=2,col="steelblue")

fc_base = exp(mean_betamog/mean_alphamog)
fc_mog = (mean_lambdamog*fc_base)/(mean_lambdamog*fc_base+1-fc_base);fc_mog

abline(h=fc_mog,lwd=2,col="red")



# ---
# 3.3 Programação no jags  ====
# ---










