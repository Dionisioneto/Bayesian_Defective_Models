### ------
### Fase 3: Estimação Bayesiana dos modelos MCMC via STAN e JAGS.
### 1. Modelo Gompertz
### ------

# [1. Modelo Gompertz];
# 2. Modelo Gaussiano-inverso;
# 3. Modelo Marshall-Olkin Gompertz;
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

gen.cure.gompertz = function(n,a,b,p){
  rm = rbinom(n=n,size=1,prob=1-p)

  t=rep(NA,n)
  
  for(i in 1:n){
    t[i]=ifelse(rm[i]==0, Inf, 
                log((-(a/b)*log(1-runif(n=1,min=0,max=1-p))) + 1)*(1/a))
    
  }
  
  t_finite = ifelse(t==Inf,0,t)
  
  u2 = runif(n=n,0,max(t_finite))
  
  t2 = pmin(t,u2) ; delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta))
} 

## Determine some values
n = 800

a0g=-0.8
b0g=0.5
p0g = exp(beta/alpha) ## escolher a proporção da fração de cura de acordo exp(beta/alpha)


data.gompertz = gen.cure.gompertz(n=n,a=a0g,b=b0g,p=p0g)
colnames(data.gompertz) = c("tempo", "delta")

# ## Verificando na curva de Kaplan-Meier 
# survival_object = Surv(data.gompertz[,1], data.gompertz[,2])
# km_fit = survfit(survival_object ~ 1)
# 
# 
# plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
#      main = "Curva de Kaplan-Meier", conf.int = F)
# 
# # Plotanto o modelo
# t_grid = seq(0,10,by=0.01)
# st_t_grid = St_Gompertz(t=t_grid,alpha=alpha,beta=beta)
# lines(t_grid,st_t_grid, lwd=2, col = "deeppink")

# ---
# 3.2 Programação no stan  ====
# ---

## Portanto, no objeto data.gompertz, temos os dados do tempo de falha e o indicador de censura.

## Vamos atribuir uma priori para alpha e beta não informativas.


# ex1_cod_gompertz_stan = "
# 
# // Estimação Bayesiana na distribuição Gompertz
# 
#   data {
#     int<lower=0> N;                    // tamanho amostral
#     array[N] real time;                // tempo de falha observado
#     array[N] int<lower=0, upper=1> delta;  // indicador da falha
#   }
# 
#   parameters {
#     real alpha;                   
#     real<lower=0> beta;
#   }
# 
#   model {
#     // Priori para beta: p(beta) = 1/beta
#   target += -log(beta);  // Log-priori para beta (Jeffreys)
# 
#     // Definição da verossimilhança manualmente
#     for (i in 1:N) {
#       // Verossimilhança para cada observação i
#       target +=  log(beta)*delta[i] + alpha*delta[i]*time[i] - (beta / alpha) * (exp(alpha * time[i]) - 1);
#     }
#   }
# 
# "

## Com o código feito, vamos estudar o estudo de prioris 
## alpha ~ Normal(mu0,sigma20) e beta~gamma(a2,b2), diferente
## da abordagem desenvolvida pelo Ricardo.

cod2_gompertz_stan = "

// Estimação Bayesiana na distribuição Gompertz

  data {
    int<lower=0> N;                    // tamanho amostral
    array[N] real time;                // tempo de falha observado
    array[N] int<lower=0, upper=1> delta;  // indicador da falha
  }

  parameters {
    real alpha;                   
    real<lower=0> beta;
  }

  model {
  alpha ~ normal(0,100);
  beta ~ gamma(0.001,0.001);

    // Definição da verossimilhança manualmente
    for (i in 1:N) {
      // Verossimilhança para cada observação i
      target +=  log(beta)*delta[i] + alpha*delta[i]*time[i] - (beta / alpha) * (exp(alpha * time[i]) - 1);
    }
  }

"

## Transcrever o código escrito para um file stan 
writeLines(cod2_gompertz_stan, con = "cod2_gompertz_stan.stan")


## Organizando os dados [data list]

data_ex1_gompertz = list(N = dim(data.gompertz)[1], 
                         time = data.gompertz[,1],
                         delta = data.gompertz[,2])


## Compilar e rodar o modelo
gfit = stan(file = 'cod2_gompertz_stan.stan', data = data_ex1_gompertz, 
            chains = 1, iter = 5000, warmup = 500)

a0g;b0g
gfit

gfit_post_samples = extract(gfit)


plot(gfit_post_samples$alpha, type='l')
abline(h=a0g,col="red", lwd=2)

plot(gfit_post_samples$beta, type='l')
abline(h=b0g,col="red", lwd=2)


## estimativas pontuais
mean_alphag = mean(gfit_post_samples$alpha)
mean_betag =mean(gfit_post_samples$beta)


plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

# Plotanto o modelo
t_grid = seq(0,10,by=0.01)
st_t_grid = St_Gompertz(t=t_grid,alpha=mean_alphag,beta=mean_betag)
lines(t_grid,st_t_grid, lwd=2, col = "deeppink")

## Fração de cura calculada
fc_stan = exp(mean_betag/mean_alphag)
abline(h=fc_stan,lwd=2,col="purple")


## com o uso de prioris não informativas, deu tudo certo!

# ---
# 3.3 Programação no JAGS  ====
# ---

library(rjags)

## Código para o modelo

cod2_jags_gompertz = "
model{
for(i in 1:n){
	logftgtz[i] = log(beta) + alpha*time[i] - (beta/alpha)*(exp(alpha*time[i]) - 1)
	logstgtz[i] = -(beta/alpha)*(exp(alpha*time[i]) - 1)
	
	## Definição da log-verossimilhança, utilizando o truque de zeros.
	phi[i] = 100000 - delta[i]*logftgtz[i] - (1-delta[i])*logstgtz[i]
	zeros[i] ~ dpois(phi[i])
}
## Prioris para alpha e beta (assumindo independência)
	alpha ~ dunif(-20, 20)
  beta ~ dunif(0, 20)
}
"

## Especificação dos dados a serem utilizados

data.jags.gtz = list(n=nrow(data.gompertz), delta = data.gompertz[,2],
                     time = data.gompertz[,1], zeros = rep(0,nrow(data.gompertz)))


## nome dos parametros a serem salvos
p.gtz.jags = c("alpha", "beta")
i.gtz.jags = function(){list(alpha=rnorm(1),beta=rgamma(1,0.01, 0.01))}

mgtz.jags = jags.model(data = data.jags.gtz, file = textConnection(cod2_jags_gompertz),
                       n.chains = 3, inits = i.gtz.jags)

## run the model for 1000 burn-in simulations
update(mgtz.jags, 1000)

## run the model for 5000 aditional simulations to keep one in 10.
res.gtz.jags = coda.samples(mgtz.jags, variable.names = p.gtz.jags,
                            n.iter = 5000, n.thin=10)


result = as.mcmc(do.call(rbind, res.gtz.jags))


## estimativas das médias a posteriori
colMeans(result)

## Estudo do percurso da cadeia
plot(result[,1])
plot(result[,2])


## estimativas pontuais

mean_jags_alpha = mean(result[,1])
mean_jags_beta = mean(result[,2])


## Visualização gráfica
plot(km_fit, xlab = "Tempo", ylab = "Probabilidade de Sobrevivência",
     main = "Curva de Kaplan-Meier", conf.int = F)

# Plotanto o modelo
t_grid = seq(0,10,by=0.01)
st_t_grid = St_Gompertz(t=t_grid,alpha=mean_jags_alpha,beta=mean_jags_beta)
lines(t_grid,st_t_grid, lwd=2, col = "deeppink")

## Fração de cura calculada
fc_stan = exp(mean_jags_beta/mean_jags_alpha)
abline(h=fc_stan,lwd=2,col="purple")


















