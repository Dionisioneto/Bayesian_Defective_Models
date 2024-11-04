## Simulação Monte Carlo
## Modelo Kumaraswamy Gaussiana Inversa

## N = 50, 100, 1.000 e 10.000.
## Censoring(%) = 15%, 45%


# ----
# Functions and packages
# ----

library(rstan)
source("https://raw.githubusercontent.com/Dionisioneto/Bayesian_Defective_Models/refs/heads/master/kumaraswamy_defective/Kumaraswamy_funcoes_e_geracao.R")

cod_KIG_stan = "
  data {
    int<lower=0> N;                         // tamanho amostral
    vector[N] time;                        // tempo de falha observado
    vector[N] delta;                      // indicador do evento
  }

  parameters {
    real alpha;
    real<lower=0> beta;
    real<lower=0> psi;
    real<lower=0> kappa;
  }

  model {
    // prioris
  
  alpha ~ normal(-1,10);
  beta ~ gamma(0.1,0.1);
  kappa ~ gamma(0.1,0.1);
  psi ~ gamma(0.1,0.1);
    
    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real exp_term = exp(2 * alpha / beta);

      // Verossimilhança do evento observado
      target += delta[i] * (log(kappa) + log(psi) - 0.5 * log(2 * beta * pi() * time[i]^3) - (1 - alpha * time[i])^2 / (2 * beta * time[i])
              + (psi - 1) * log(phi1 + exp_term * phi2));

      // Verossimilhança para censura
      target += delta[i] * (kappa - 1) * log(1 - (phi1 + exp_term * phi2)^psi);
      target += (1 - delta[i]) * kappa * log(1 - (phi1 + exp_term * phi2)^psi);
    }
  }
"



## Transcrever o código escrito para um file stan 
writeLines(cod_KIG_stan, con = "cod_KIG_stan.stan")



## ---
## Funções importantes ===
## ---

bias = function(true,est){(est-true)/true*100}


## Tabelas e informações

n.replicas = 500
## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4

n_amostral = c(50, 100, 1000, 10000)

kig_mc_alpha = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                     dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kig_mc_beta = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                    dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kig_mc_psi = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                   dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

kig_mc_kappa = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                     dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))


## colunas: 4 valores do effective sample size (n_eff)
kig_mc_eff = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                   dimnames = list(1:n.replicas,c("alpha", "beta", "psi","kappa"),n_amostral))

## colunas: 4 valores do Rhat
kig_mc_Rhat = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                    dimnames = list(1:n.replicas,c("alpha", "beta", "psi","kappa"),n_amostral))

## Censura

kig_mc_cens = array(data=0, dim = c(n.replicas, 1, length(n_amostral)))

## ---
## Scenarios: 
## ---
# primeiro: Alpha = -0.20, Beta = 2, psi = 2.2, kappa = 2;
# segundo: Alpha = -1, Beta = 10, psi = 1.4, kappa = 0.6.


a0kig=-0.20; b0kig=2; psi0kig=2.2; k0kig=2

pbkig = 1 - exp(2*a0kig/b0kig)
p0kig = (1-(1-pbkig)^psi0kig)^k0kig; p0kig



for(j in 1:length(n_amostral)){
  for(i in 1:n.replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    dados.kig = gen.cure.kig(n=n_amostral[j],a=a0kig,b=b0kig,k=k0kig,ps=psi0kig,p=p0kig)
    
    data_kig = list(N = dim(dados.kig)[1], 
                   time = dados.kig[,1],
                   delta = dados.kig[,2])
    
    ## Compilar e rodar o modelo
    kigfit = stan(file = 'cod_KIG_stan.stan', data = data_kig, 
                 chains = 1, iter = 8000, warmup = 1000)
    
    kigfit_post_samples = extract(kigfit)
    mean_alpha = mean(kigfit_post_samples$alpha); sd_alpha = sd(kigfit_post_samples$alpha)
    mean_beta = mean(kigfit_post_samples$beta); sd_beta = sd(kigfit_post_samples$beta)
    mean_psi = mean(kigfit_post_samples$psi); sd_psi = sd(kigfit_post_samples$psi)
    mean_kappa = mean(kigfit_post_samples$kappa); sd_kappa = sd(kigfit_post_samples$kappa)
    
    ## Intervalos de credibilidade
    
    sumario_kigfit = summary(kigfit)
    
    check_alpha = sum(a0kig >= sumario_kigfit$summary[1,"2.5%"] && a0kig <= sumario_kigfit$summary[1,"97.5%"]) 
    check_beta = sum(b0kig >= sumario_kigfit$summary[2,"2.5%"] && b0kig <= sumario_kigfit$summary[2,"97.5%"]) 
    check_psi = sum(psi0kig >= sumario_kigfit$summary[3,"2.5%"] && psi0kig <= sumario_kigfit$summary[3,"97.5%"]) 
    check_kappa = sum(k0kig >= sumario_kigfit$summary[4,"2.5%"] && k0kig <= sumario_kigfit$summary[4,"97.5%"]) 
    
    ## salvando valores
    
    kig_mc_alpha[i,1,j] = mean_alpha; kig_mc_alpha[i,2,j] = sd_alpha; kig_mc_alpha[i,3,j] = bias(est=mean_alpha,true=a0kig)
    kig_mc_beta[i,1,j] = mean_beta; kig_mc_beta[i,2,j] = sd_beta; kig_mc_beta[i,3,j] = bias(est=mean_beta,true=b0kig)
    kig_mc_psi[i,1,j] = mean_psi; kig_mc_psi[i,2,j] = sd_psi; kig_mc_psi[i,3,j] = bias(est=mean_psi,true=psi0kig)
    kig_mc_kappa[i,1,j] = mean_kappa; kig_mc_kappa[i,2,j] = sd_kappa; kig_mc_kappa[i,3,j] = bias(est=mean_kappa,true= k0kig)
    
    ## salvando intervalso de credibilidade
    kig_mc_alpha[i,4,j] = check_alpha; kig_mc_beta[i,4,j] = check_beta; kig_mc_psi[i,4,j] = check_psi; kig_mc_kappa[i,4,j] = check_kappa
    
    ## salvando diagnóstico das cadeias
    kig_mc_eff[i,1:4,j] = sumario_kigfit$summary[1:4,"n_eff"]
    kig_mc_Rhat[i,1:4,j] = sumario_kigfit$summary[1:4,"Rhat"]
    
    ## salvando o percentual censura dos dados
    kig_mc_cens[i,1,j] = 1-(sum(dados.kig[,2])/n_amostral[j])
    
    ## Salvando os resultados em um arquivo .csv
    if(i==n.replicas){
      write.csv2(kig_mc_alpha[,,j], file = paste("resumos/","kig_mc_alpha_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_beta[,,j], file = paste("resumos/","kig_mc_beta_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_psi[,,j], file = paste("resumos/","kig_mc_psi_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_kappa[,,j], file = paste("resumos/","kig_mc_kappa_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_eff[,,j], file = paste("resumos/","kig_mc_eff_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_Rhat[,,j], file = paste("resumos/","kig_mc_Rhat_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(kig_mc_cens[,,j], file = paste("resumos/","kig_mc_cens_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
      
    }
  }
}

# kig_mc_alpha
# kig_mc_beta
# kig_mc_psi
# kig_mc_kappa
# 
# kig_mc_eff
# kig_mc_Rhat
# 
# a0kig; b0kig; psi0kig; k0kig
# colMeans(kig_mc_alpha)
# colMeans(kig_mc_beta)
# colMeans(kig_mc_psi)
# colMeans(kig_mc_kappa)




