## Simulação Monte Carlo
## Modelo Marshall-Olkin Gaussiana Inversa

## N = 50, 100, 1.000 e 10.000.
## Censoring(%) = 15%, 45%

# ----
# Functions and pacakages
# ----

library(rstan)
source("https://raw.githubusercontent.com/Dionisioneto/Bayesian_Defective_Models/refs/heads/master/marshall_olkin_defective/MO_funcoes_e_geracao.R")

## ---
## Funções importantes ===
## ---

bias = function(true,est){(est-true)/true*100}

## Codigo Stan

cod_moig_stan = "
  data {
   int<lower=0> N;                        // tamanho amostral
   vector[N] time;                       // tempo de falha observado
   vector[N] delta;                      // indicador do evento
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
## Transcrever o código escrito para um file stan 
writeLines(cod_moig_stan, con = "cod_moig_stan.stan")

## ---
## Scenarios: 
## ---
# primeiro: Alpha = -1.2, Beta = 2, Lambda = 0.8
# segundo: Alpha = -2, Beta = 5, Lambda = 2


n.replicas = 500

## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4
n_amostral = c(50, 100, 1000, 10000)


moig_mc_alpha = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                      dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))
moig_mc_beta = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                     dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))
moig_mc_lambda = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                       dimnames = list(1:n.replicas,c("mean", "sd", "bias", "cp"),n_amostral))

## colunas: 3 valores do effective sample size (n_eff)
moig_mc_eff = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                    dimnames = list(1:n.replicas,c("alpha", "beta", "lambda", "lp"),n_amostral))

## colunas: 3 valores do Rhat
moig_mc_Rhat = array(data=0, dim = c(n.replicas, ncols_mc, length(n_amostral)),
                     dimnames = list(1:n.replicas,c("alpha", "beta", "lambda", "lp"),n_amostral))

## Censura

moig_mc_cens = array(data=0, dim = c(n.replicas, 1, length(n_amostral)))

## ajuste dos parametros
a0moig=-2;b0moig=10;l0moig = 2

pgmoig= 1 - exp(2*a0moig/b0moig)
p0moig=(l0moig*pgmoig)/(l0moig*pgmoig+1-pgmoig)
p0moig


for(j in 1:length(n_amostral)){
  for(i in 1:n.replicas){
    sucesso = FALSE  # Flag para indicar se a iteração foi bem-sucedida
    
    while(!sucesso){
      tryCatch({
        print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
        
        dados.moig = gen.cure.moig(n=n_amostral[j],a=a0moig,b=b0moig,l=l0moig,p=p0moig)
        
        data_moig = list(N = dim(dados.moig)[1], 
                         time = dados.moig[,1],
                         delta = dados.moig[,2])
        
        
        ## Compilar e rodar o modelo
        moigfit = stan(file = 'cod_moig_stan.stan', data = data_moig, 
                       chains = 1, iter = 8000, warmup = 1000)
        
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
        
        ## Salvando os resultados em um arquivo .csv
        if(i==n.replicas){
          write.csv2(moig_mc_alpha[,,j], file = paste("resumos/","moig_mc_alpha_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          write.csv2(moig_mc_beta[,,j], file = paste("resumos/","moig_mc_beta_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          write.csv2(moig_mc_lambda[,,j], file = paste("resumos/","moig_mc_lambda_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          write.csv2(moig_mc_eff[,,j], file = paste("resumos/","moig_mc_eff_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          write.csv2(moig_mc_Rhat[,,j], file = paste("resumos/","moig_mc_Rhat_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          write.csv2(moig_mc_cens[,,j], file = paste("resumos/","moig_mc_cens_", "rep_",n.replicas, "_n_",n_amostral[j], ".csv", sep=""))
          
        }
        
        sucesso = TRUE # Marca sucesso se não houver erro
      }, error = function(e) {
        message(paste("Erro na iteração", i, "com n_amostral =", n_amostral[j], ":", e$message))
        # A iteração será repetida devido ao `while (!sucesso)`
      })
    }
  }
}




# Dentro da pasta, deve haver o codigo stan e o a pasta de resumos para salvar. 


# moig_mc_alpha
# moig_mc_beta
# moig_mc_lambda

# moig_mc_eff
# moig_mc_Rhat

# colMeans(moig_mc_alpha)
# colMeans(moig_mc_beta)
# colMeans(moig_mc_lambda)



