### -----------
### Monte Carlo 
### de regressão MO Gompertz
### Geração de dados e estimação Bayesiana
### -----------

library(pacman)
p_load(survival,rstan)

#setwd('C:/Users/dioni/OneDrive - University of São Paulo/Doutorado em Estatística/2024.2/2_Topicos_Avancados_de_Pesquisa_I/R_code/covariates/MO_Gompertz')
## 1. Base functions ====

St_Gompertz = function(t,alpha,beta){
  st = exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(st)
}


Ft_Gompertz = function(t,alpha,beta){
  Ft = 1 - exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(Ft)
}

ft_Gompertz = function(t,alpha,beta){
  ft = beta*exp(alpha*t)*exp(-(beta/alpha)*(exp(alpha*t)-1))
  return(ft)
}

# MO Gompertz

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


## 2. Generate data ====

## inserir xb e xd sem o vetor de 1's do intercepto
gen.cure.reg.mog = function(n,as,bs,xa,xb,l){
  finv.mog = function(t,alpha,beta,lambda,unif){Ftmo_gompertz(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.mog = function(alpha,beta,lambda,unif){t=uniroot(finv.mog,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  p_base = exp(b/a)
  p_cure = (l*p_base)/(l*p_base+1-p_base)
  
  t = rep(NA,n)
  
  for(i in 1:n){
    rmi = rbinom(n=1,size=1,prob=1-p_cure[i,])
    uni = runif(n=1,0,max=1-p_cure[i,])
    t[i] = ifelse(rmi==0, Inf, ger.mog(alpha=a[i,],beta=b[i,],lambda=l,unif=uni))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta,xa,xb))
}


## 4. Tables to save data ====

n_replicas = 500                        ## number of Monte Carlo replicates  
ncols_mc = 4                             ## numero of information in each replicate (mean,sd,bias,cp)
n_amostral = c(50, 100, 500, 1000, 5000) ## sample sizes to study

## Tables

## as = (a0, a1, a2)
mog_mc_a0a1a2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "sd", "bias", "cp"),3),n_amostral))


## bs = (b0, b1, b2)
mog_mc_b0b1b2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "sd", "bias", "cp"),3),n_amostral))

## lambda

mog_mc_lambda = array(data=0, dim = c(n_replicas, ncols_mc, length(n_amostral)),
                      dimnames = list(1:n_replicas,c("mean", "sd", "bias", "cp"),n_amostral))

## Censura

mog_mc_cens = array(data=0, dim = c(n_replicas, 1, length(n_amostral)))


## ---
## Scenarios: 
## ---

# first (15%): 
# a0 = c(-1.2, 0.5, 0.2)
# b0 = c(-0.5, 1.5, 1.2)
# l0 = 0.5

# second (45%)
a0 = c(-1.2, 0.5, 0.2)
b0 = c(-1.1, 1.5, 0.9)
l0 = 2

params=c(a0,b0,l0)

for(j in 1:length(n_amostral)){
  for(i in 1:n_replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    ## covariates data
    xa1 = rbinom(n=n_amostral[j],size=1,prob=0.7); xa2 = runif(n=n_amostral[j],0,1)
    xb1 = rbinom(n=n_amostral[j],size=1,prob=0.5); xb2 = runif(n=n_amostral[j],0,1)
    xa0 = cbind(xa1, xa2)
    xb0 = cbind(xb1, xb2)
    
    dados.reg.mog = gen.cure.reg.mog(n=n_amostral[j],as=a0,bs=b0,xa=xa0,xb=xb0,l=l0)
    dados.reg.mog = as.data.frame(dados.reg.mog)
    
    # Data and stan function
    
    data_regmog = list(N = dim(dados.reg.mog)[1],
                       p = dim(cbind(1,dados.reg.mog[,3:4]))[2],
                       s = dim(cbind(1,dados.reg.mog[,5:6]))[2],
                       time = dados.reg.mog[,1],
                       delta = dados.reg.mog[,2],
                       xa = cbind(1,dados.reg.mog[,3:4]),
                       xb = cbind(1,dados.reg.mog[,5:6]))
    
    ## Run the Bayesian model
    regmogfit = stan(file = 'cod_regmog.stan', data = data_regmog,
                     chains = 1, iter=1000, warmup=5000)
    
    ## Extract posterior samples
    regmogfit_post_samples = extract(regmogfit)
    
    ## posterior mean 
    ma0=mean(regmogfit_post_samples$as[,1])
    ma1=mean(regmogfit_post_samples$as[,2])
    ma2=mean(regmogfit_post_samples$as[,3])
    mb0=mean(regmogfit_post_samples$bs[,1])
    mb1=mean(regmogfit_post_samples$bs[,2])
    mb2=mean(regmogfit_post_samples$bs[,3])
    mlamb=mean(regmogfit_post_samples$lambda)
    post_mean=c(ma0,ma1,ma2,mb0,mb1,mb2,mlamb)
    
    ## posterior standard deviation
    sda0=sd(regmogfit_post_samples$as[,1])
    sda1=sd(regmogfit_post_samples$as[,2])
    sda2=sd(regmogfit_post_samples$as[,3])
    sdb0=sd(regmogfit_post_samples$bs[,1])
    sdb1=sd(regmogfit_post_samples$bs[,2])
    sdb2=sd(regmogfit_post_samples$bs[,3])
    sdlamb=sd(regmogfit_post_samples$lambda)
    post_sd=c(sda0,sda1,sda2,sdb0,sdb1,sdb2,sdlamb)
    
    ## relative bias
    rbias = (post_mean-params)/params*100
    
    ## coverage probability (95%)
    q025a0 = quantile(regmogfit_post_samples$as[,1],probs=0.025);q975a0 = quantile(regmogfit_post_samples$as[,1],probs=0.975)
    q025a1 = quantile(regmogfit_post_samples$as[,2],probs=0.025);q975a1 = quantile(regmogfit_post_samples$as[,2],probs=0.975)
    q025a2 = quantile(regmogfit_post_samples$as[,3],probs=0.025);q975a2 = quantile(regmogfit_post_samples$as[,3],probs=0.975)
    q025b0 = quantile(regmogfit_post_samples$bs[,1],probs=0.025);q975b0 = quantile(regmogfit_post_samples$bs[,1],probs=0.975)
    q025b1 = quantile(regmogfit_post_samples$bs[,2],probs=0.025);q975b1 = quantile(regmogfit_post_samples$bs[,2],probs=0.975)
    q025b2 = quantile(regmogfit_post_samples$bs[,3],probs=0.025);q975b2 = quantile(regmogfit_post_samples$bs[,3],probs=0.975)
    q025lamb = quantile(regmogfit_post_samples$lambda,probs=0.025);q975lamb = quantile(regmogfit_post_samples$lambda,probs=0.975)
    
    qinf025 = c(q025a0,q025a1,q025a2,q025b0,q025b1,q025b2,q025lamb)
    qsup975 = c(q975a0,q975a1,q975a2,q975b0,q975b1,q975b2,q975lamb)
    
    cpost95 = as.integer((params >= qinf025) & (params <= qsup975))
    cpost95
    
    ## salvando valores
    
    mog_mc_a0a1a2[i,1,j] = post_mean[1]; mog_mc_a0a1a2[i,2,j] = post_sd[1]; mog_mc_a0a1a2[i,3,j] = rbias[1]; mog_mc_a0a1a2[i,4,j] = cpost95[1]
    mog_mc_a0a1a2[i,5,j] = post_mean[2]; mog_mc_a0a1a2[i,6,j] = post_sd[2]; mog_mc_a0a1a2[i,7,j] = rbias[2]; mog_mc_a0a1a2[i,8,j] = cpost95[2]
    mog_mc_a0a1a2[i,9,j] = post_mean[3]; mog_mc_a0a1a2[i,10,j] = post_sd[3]; mog_mc_a0a1a2[i,11,j] = rbias[3]; mog_mc_a0a1a2[i,12,j] = cpost95[3]
    
    mog_mc_b0b1b2[i,1,j] = post_mean[4]; mog_mc_b0b1b2[i,2,j] = post_sd[4]; mog_mc_b0b1b2[i,3,j] = rbias[4]; mog_mc_b0b1b2[i,4,j] = cpost95[4]
    mog_mc_b0b1b2[i,5,j] = post_mean[5]; mog_mc_b0b1b2[i,6,j] = post_sd[5]; mog_mc_b0b1b2[i,7,j] = rbias[5]; mog_mc_b0b1b2[i,8,j] = cpost95[5]
    mog_mc_b0b1b2[i,9,j] = post_mean[6]; mog_mc_b0b1b2[i,10,j] = post_sd[6]; mog_mc_b0b1b2[i,11,j] = rbias[6]; mog_mc_b0b1b2[i,12,j] = cpost95[6]
    
    
    mog_mc_lambda[i,1,j] = post_mean[7]; mog_mc_lambda[i,2,j] = post_sd[7]; mog_mc_lambda[i,3,j] = rbias[7]; mog_mc_lambda[i,4,j] = cpost95[7]
    
    
    ## salvando o percentual censura dos dados
    mog_mc_cens[i,,j] = 1-(sum(dados.reg.mog[,2])/n_amostral[j])
    
    ## Salvando os resultados em um arquivo .csv
    if(i==n_replicas){
      write.csv2(mog_mc_a0a1a2[,,j], file = paste("resumos_Bayes_MC_reg_MO_Gompertz/","mog_mc_a0a1a2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(mog_mc_b0b1b2[,,j], file = paste("resumos_Bayes_MC_reg_MO_Gompertz/","mog_mc_b0b1b2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(mog_mc_lambda[,,j], file = paste("resumos_Bayes_MC_reg_MO_Gompertz/","mog_mc_lambda_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(mog_mc_cens[,,j], file = paste("resumos_Bayes_MC_reg_MO_Gompertz/","mog_mc_cens_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
    }
    
  }
}









