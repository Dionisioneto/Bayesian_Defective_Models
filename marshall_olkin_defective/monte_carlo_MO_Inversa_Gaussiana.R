## Simulação Monte Carlo
## Modelo Marshall-Olkin Gaussiana Inversa

## N = 50, 100, 1.000 e 10.000.
## Censoring(%) = 15%, 45%

# ----
# Functions and pacakages
# ----

library(pacman)
p_load(survival,rstan, R2jags)
source("https://raw.githubusercontent.com/Dionisioneto/Bayesian_Defective_Models/refs/heads/master/marshall_olkin_defective/MO_funcoes_e_geracao.R")

## ---
## Funções importantes ===
## ---

bias = function(true,est){(est-true)/true*100}


## ---
## Scenarios: 
## ---
# primeiro: Alpha = -1.2, Beta = 2, Lambda = 0.8
# segundo: Alpha = -2, Beta = 5, Lambda = 2


n.replicas = 1000

## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4
n_amostral = c(50, 100, 1000, 10000)
#n_amostral = 1000

moig_mc_alpha = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))
moig_mc_beta = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))
moig_mc_lambda = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## colunas: 3 valores do effective sample size (n_eff)
moig_mc_eff = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## colunas: 3 valores do Rhat
moig_mc_Rhat = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## Censura

moig_mc_cens = array(data=0, dim = c(n.replicas, 1, length(n_amostral)))

## ajuste dos parametros
a0moig=-2;b0moig=10;l0moig = 2

pgmoig= 1 - exp(2*a0moig/b0moig)
p0moig=(l0moig*pgmoig)/(l0moig*pgmoig+1-pgmoig)
p0moig


for(j in 1:length(n_amostral)){
  for(i in 1:n.replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    dados.moig = gen.cure.moig(n=n_amostral[j],a=a0moig,b=b0moig,l=l0moig,p=p0moig)
    
    data_moig = list(N = dim(dados.moig)[1], 
                     time = dados.moig[,1],
                     delta = dados.moig[,2])
    
    
    ## Compilar e rodar o modelo
    moigfit = stan(file = 'cod_moig_stan.stan', data = data_moig, 
                   chains = 1, iter = 1000, warmup = 100)
    
    moigfit_post_samples = extract(moigfit)
    mean_alpha = mean(moigfit_post_samples$alpha); sd_alpha = sd(moigfit_post_samples$alpha)
    mean_beta = mean(moigfit_post_samples$beta); sd_beta = sd(moigfit_post_samples$beta)
    mean_lambda = mean(moigfit_post_samples$lambda); sd_lambda = sd(moigfit_post_samples$lambda)
    
    ## Intervalos de credibilidade
    
    sumario_moigfit = summary(moigfit)
    
    check_alpha = sum(a0moig >= sumario_moigfit$summary[1,"2.5%"] && a0moig <= sumario_moigfit$summary[1,"97.5%"]) 
    check_beta = sum(b0moig >= sumario_moigfit$summary[2,"2.5%"] && b0moig <= sumario_moigfit$summary[2,"97.5%"]) 
    check_lambda = sum(l0moig >= sumario_moigfit$summary[3,"2.5%"] && l0moig <= sumario_moigfit$summary[3,"97.5%"]) 
    
    ## salvando valores
    
    moig_mc_alpha[i,1,j] = mean_alpha; moig_mc_alpha[i,2,j] = sd_alpha; moig_mc_alpha[i,3,j] = bias(est=mean_alpha,true=a0moig)
    moig_mc_beta[i,1,j] = mean_beta; moig_mc_beta[i,2,j] = sd_beta; moig_mc_beta[i,3,j] = bias(est=mean_beta,true=b0moig)
    moig_mc_lambda[i,1,j] = mean_lambda; moig_mc_lambda[i,2,j] = sd_lambda; moig_mc_lambda[i,3,j] = bias(est=mean_lambda,true=l0moig)
    
    ## salvando intervalso de credibilidade
    moig_mc_alpha[i,4,j] = check_alpha; moig_mc_beta[i,4,j] = check_beta; moig_mc_lambda[i,4,j] = check_lambda
    
    ## salvando diagnóstico das cadeias
    moig_mc_eff[i,1:4,j] = sumario_moigfit$summary[,"n_eff"]
    moig_mc_Rhat[i,1:4,j] = sumario_moigfit$summary[,"Rhat"]
    
    ## salvando o percentual censura dos dados
    moig_mc_cens[i,1,j] = 1-(sum(dados.moig[,2])/n_amostral[j])
  }
}

moig_mc_alpha
moig_mc_beta
moig_mc_lambda

moig_mc_eff
moig_mc_Rhat


colMeans(moig_mc_alpha)
colMeans(moig_mc_beta)
colMeans(moig_mc_lambda)



