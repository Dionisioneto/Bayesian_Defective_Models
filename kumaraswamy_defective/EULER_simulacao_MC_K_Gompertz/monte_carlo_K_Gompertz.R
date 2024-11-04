## Simulação Monte Carlo
## Modelo Kumaraswamy Gompertz

## N = 50, 100, 1.000 e 10.000.
## Censoring(%) = 15%, 45%

# ----
# Functions and packages
# ----


library(rstan)
source("https://raw.githubusercontent.com/Dionisioneto/Bayesian_Defective_Models/refs/heads/master/kumaraswamy_defective/Kumaraswamy_funcoes_e_geracao.R")

cod_kgtz_stan = "
  data{
  int<lower=0> N;
  vector[N] time;
  vector[N] delta;
  }
  
  parameters{
  real alpha;
  real<lower=0> beta;
  real<lower=0> psi;
  real<lower=0> kappa;
  }
  
  model{
  
  // prioris
  
  alpha ~ normal(-1,10);
  beta ~ gamma(0.1,0.1);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.1,0.1);
  
  // Definição manual da função de verossimilhança
  for(i in 1:N){
    target += delta[i]*(log(kappa*psi*beta*exp(alpha*time[i])*exp((beta-beta*exp(alpha*time[i]))/alpha)*(1-exp((beta-beta*exp(alpha*time[i]))/alpha))^(psi-1) * (1-(1-(exp((beta-beta*exp(alpha*time[i]))/alpha)))^psi)^(kappa-1))) +
              (1 - delta[i])*(log((1 - (1 - exp((beta-beta*exp(alpha*time[i]))/alpha))^psi)^kappa));
  }
}
"


## Transcrever o código escrito para um file stan 
writeLines(cod_kgtz_stan, con = "cod_kgtz_stan.stan")



## ---
## Funções importantes ===
## ---

bias = function(true,est){(est-true)/true*100}

## ---
## Scenarios: 
## ---
# primeiro: Alpha = -2, Beta = 10, psi = 1.2, kappa = 0.5;
# segundo: Alpha = -1, Beta = 1, psi = 2, kappa = 2.


n.replicas = 500
## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4
n_amostral = c(50, 100, 1000, 10000)

## Tabelas

kg_mc_alpha = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                    dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kg_mc_beta = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                   dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kg_mc_psi = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                  dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kg_mc_kappa = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                    dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))


## colunas: 4 valores do effective sample size (n_eff)
kg_mc_eff = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                  dimnames = list(1:n.replicas,c("alpha", "beta", "psi","kappa"),n_amostral))

## colunas: 4 valores do Rhat
kg_mc_Rhat = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                   dimnames = list(1:n.replicas,c("alpha", "beta", "psi","kappa"),n_amostral))

## Censura

kg_mc_cens = array(data=0, dim = c(n.replicas, 1, length(n_amostral)))

# primeiro: Alpha = -2, Beta = 10, psi = 1.2, kappa = 0.5;
# segundo: Alpha = -1, Beta = 1, psi = 2, kappa = 2.


a0kg=-2; b0kg=10; psi0kg=1.2; k0kg=0.5

pgkg = exp(b0kg/a0kg)
p0kg = (1-(1-pgkg)^psi0kg)^k0kg; p0kg


# ## cutes de valores inicias aleatórios
# # init_fun = function() {
# #   list(alpha = runif(1, -2.55,-2), 
# #        beta = runif(1, 9,11), 
# #        psi = runif(1, 1,2),
# #       kappa = runif(1, 0,1))
# # }
# 
init_val = list(list(alpha = -1,
               beta = 1,
               psi = 1,
               kappa = 1))



for(j in 1:length(n_amostral)){
  for(i in 1:n.replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    dados.kg = gen.cure.kgz(n=n_amostral[j],a=a0kg,b=b0kg,k=k0kg,ps=psi0kg,p=p0kg)
    
    data_kg = list(N = dim(dados.kg)[1], 
                    time = dados.kg[,1],
                    delta = dados.kg[,2])
    
    ## Compilar e rodar o modelo
    kgfit = stan(file = 'cod_kgtz_stan.stan', data = data_kg, 
                  chains = 1, iter = 8000, warmup = 1000, init=init_val)
    
    kgfit_post_samples = extract(kgfit)
    median_alpha = median(kgfit_post_samples$alpha); sd_alpha = sd(kgfit_post_samples$alpha)
    median_beta = median(kgfit_post_samples$beta); sd_beta = sd(kgfit_post_samples$beta)
    median_psi = median(kgfit_post_samples$psi); sd_psi = sd(kgfit_post_samples$psi)
    median_kappa = median(kgfit_post_samples$kappa); sd_kappa = sd(kgfit_post_samples$kappa)
    
    ## Intervalos de credibilidade
    
    sumario_kgfit = summary(kgfit)
    
    check_alpha = sum(a0kg >= sumario_kgfit$summary[1,"2.5%"] && a0kg <= sumario_kgfit$summary[1,"97.5%"]) 
    check_beta = sum(b0kg >= sumario_kgfit$summary[2,"2.5%"] && b0kg <= sumario_kgfit$summary[2,"97.5%"]) 
    check_psi = sum(psi0kg >= sumario_kgfit$summary[3,"2.5%"] && psi0kg <= sumario_kgfit$summary[3,"97.5%"]) 
    check_kappa = sum(k0kg >= sumario_kgfit$summary[4,"2.5%"] && k0kg <= sumario_kgfit$summary[4,"97.5%"]) 
    
    ## salvando valores
    
    kg_mc_alpha[i,1,j] = median_alpha; kg_mc_alpha[i,2,j] = sd_alpha; kg_mc_alpha[i,3,j] = bias(est=median_alpha,true=a0kg)
    kg_mc_beta[i,1,j] = median_beta; kg_mc_beta[i,2,j] = sd_beta; kg_mc_beta[i,3,j] = bias(est=median_beta,true=b0kg)
    kg_mc_psi[i,1,j] = median_psi; kg_mc_psi[i,2,j] = sd_psi; kg_mc_psi[i,3,j] = bias(est=median_psi,true=psi0kg)
    kg_mc_kappa[i,1,j] = median_kappa; kg_mc_kappa[i,2,j] = sd_kappa; kg_mc_kappa[i,3,j] = bias(est=median_kappa,true= k0kg)
    
    ## salvando intervalso de credibilidade
    kg_mc_alpha[i,4,j] = check_alpha; kg_mc_beta[i,4,j] = check_beta; kg_mc_psi[i,4,j] = check_psi; kg_mc_kappa[i,4,j] = check_kappa
    
    ## salvando diagnóstico das cadeias
    kg_mc_eff[i,1:4,j] = sumario_kgfit$summary[1:4,"n_eff"]
    kg_mc_Rhat[i,1:4,j] = sumario_kgfit$summary[1:4,"Rhat"]
    
    ## salvando o percentual censura dos dados
    kg_mc_cens[i,1,j] = 1-(sum(dados.kg[,2])/n_amostral[j])
    
    
    ## Salvando os resultados em um arquivo .csv
    if(i==n.replicas){
      write.csv2(kg_mc_alpha[,,j], file = paste("resumos/","kg_mc_alpha_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_beta[,,j], file = paste("resumos/","kg_mc_beta_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_psi[,,j], file = paste("resumos/","kg_mc_psi_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_kappa[,,j], file = paste("resumos/","kg_mc_kappa_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_eff[,,j], file = paste("resumos/","kg_mc_eff_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_Rhat[,,j], file = paste("resumos/","kg_mc_Rhat_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kg_mc_cens[,,j], file = paste("resumos/","kg_mc_cens_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      
    }
  }
}

# kg_mc_alpha
# kg_mc_beta
# kg_mc_psi
# kg_mc_kappa
# 
# kg_mc_eff
# kg_mc_Rhat
# 
# a0kg; b0kg; psi0kg; k0kg
# colMeans(kg_mc_alpha)
# colMeans(kg_mc_beta)
# colMeans(kg_mc_psi)
# colMeans(kg_mc_kappa)

