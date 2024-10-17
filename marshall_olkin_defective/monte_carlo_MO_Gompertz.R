## Simulação Monte Carlo
## Modelo Marshall-Olkin Gompertz

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
## 1. Scenario: 
## ---
# primeiro: Alpha = -1.2, Beta = 2, Lambda = 0.8 
# segundo: Alpha = -3, Beta = 4, Lambda = 2

n.replicas = 1000

## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4
n_amostral = c(50, 100, 1000, 10000)
#n_amostral = 1000

mo_mc_alpha = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))
mo_mc_beta = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))
mo_mc_lambda = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## colunas: 3 valores do effective sample size (n_eff)
mo_mc_eff = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## colunas: 3 valores do Rhat
mo_mc_Rhat = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)))

## Censura

mog_mc_cens = array(data=0, dim = c(n.replicas, 1, length(n_amostral)))


## ajuste dos parametros
a0mog=-2;b0mog=5;l0mog=2
pbmog=exp(b0mog/a0mog); p0mog=(l0mog*pbmog)/(l0mog*pbmog+1-pbmog)
p0mog

dados.mog = gen.cure.mog(n=n,a=a0mog,b=b0mog,l=l0mog,p=p0mog)

for(j in 1:length(n_amostral)){
  for(i in 1:n.replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    dados.mog = gen.cure.mog(n=n_amostral[j],a=a0mog,b=b0mog,l=l0mog,p=p0mog)
    
    data_mog = list(N = dim(dados.mog)[1],
                    time = dados.mog[,1],
                    delta = dados.mog[,2])
    
    ## Compilar e rodar o modelo
    mogfit = stan(file = 'cod_mog_stan.stan', data = data_mog, 
                  chains = 1, iter = 10000, warmup = 1000)
    
    mogfit_post_samples = extract(mogfit)
    mean_alpha = mean(mogfit_post_samples$alpha); sd_alpha = sd(mogfit_post_samples$alpha)
    mean_beta = mean(mogfit_post_samples$beta); sd_beta = sd(mogfit_post_samples$beta)
    mean_lambda = mean(mogfit_post_samples$lambda); sd_lambda = sd(mogfit_post_samples$lambda)
  
    ## Intervalos de credibilidade
    
    sumario_mogfit = summary(mogfit)
    
    check_alpha = sum(a0mog >= sumario_mogfit$summary[1,"2.5%"] && a0mog <= sumario_mogfit$summary[1,"97.5%"]) 
    check_beta = sum(b0mog >= sumario_mogfit$summary[2,"2.5%"] && b0mog <= sumario_mogfit$summary[2,"97.5%"]) 
    check_lambda = sum(l0mog >= sumario_mogfit$summary[3,"2.5%"] && l0mog <= sumario_mogfit$summary[3,"97.5%"]) 
    
    ## salvando valores
    
    mo_mc_alpha[i,1,j] = mean_alpha; mo_mc_alpha[i,2,j] = sd_alpha; mo_mc_alpha[i,3,j] = bias(est=mean_alpha,true=a0mog)
    mo_mc_beta[i,1,j] = mean_beta; mo_mc_beta[i,2,j] = sd_beta; mo_mc_beta[i,3,j] = bias(est=mean_beta,true=b0mog)
    mo_mc_lambda[i,1,j] = mean_lambda; mo_mc_lambda[i,2,j] = sd_lambda; mo_mc_lambda[i,3,j] = bias(est=mean_lambda,true=l0mog)
    
    ## salvando intervalso de credibilidade
    mo_mc_alpha[i,4,j] = check_alpha; mo_mc_beta[i,4,j] = check_beta; mo_mc_lambda[i,4,j] = check_lambda
  
    ## salvando diagnóstico das cadeias
    mo_mc_eff[i,1:4,j] = sumario_mogfit$summary[,"n_eff"]
    mo_mc_Rhat[i,1:4,j] = sumario_mogfit$summary[,"Rhat"]
    
    ## salvando os percentuais de censua nos dados
    mog_mc_cens[i,1,j] = 1-(sum(dados.mog[,2])/n_amostral[j])
  }
}

mo_mc_alpha
mo_mc_beta
mo_mc_lambda

mo_mc_eff
mo_mc_Rhat

colMeans(mo_mc_alpha)
colMeans(mo_mc_beta)
colMeans(mo_mc_lambda)



# 
# ## 2. Scenario: 
# 
# # sumario_mogfit = summary(mogfit)
# # 
# # sum(a0mog >= sumario_mogfit$summary[1,"2.5%"] && a0mog <= sumario_mogfit$summary[1,"97.5%"]) 
# # sum(b0mog >= sumario_mogfit$summary[2,"2.5%"] && a0mog <= sumario_mogfit$summary[2,"97.5%"]) 
# # sum(l0mog >= sumario_mogfit$summary[3,"2.5%"] && a0mog <= sumario_mogfit$summary[3,"97.5%"]) 
# 
# # fc_base_2_5 = exp(sumario_mogfit$summary[2,"2.5%"]/sumario_mogfit$summary[1,"2.5%"])
# # fc_mog_2_5 = (sumario_mogfit$summary[3,"2.5%"]*fc_base_2_5)/(sumario_mogfit$summary[3,"2.5%"]*fc_base_2_5+1-fc_base_2_5)
# # 
# # fc_base_97_5 = exp(sumario_mogfit$summary[2,"97.5%"]/sumario_mogfit$summary[1,"97.5%"])
# # fc_mog_97_5 = (sumario_mogfit$summary[3,"97.5%"]*fc_base_97_5)/(sumario_mogfit$summary[3,"97.5%"]*fc_base_97_5+1-fc_base_97_5)
# # 
# # sum(p0mog >= fc_mog_2_5 && p0mog <= fc_mog_97_5)
# # 
# # sum(p0mog >= sumario_mogfit$summary[1,"2.5%"] && a0mog <= sumario_mogfit$summary[1,"97.5%"]) 
# 
# 
# 
# sumari$summary[1,"2.5%"]
# sumari$summary[1,"97.5%"]  
# 
# sumari$summary[,"n_eff"]
# sumari$summary[,"Rhat"]
#   
# # Alpha = -2, Beta = 5, Lambda = 2, p = 0.4172 
# 
# 

















