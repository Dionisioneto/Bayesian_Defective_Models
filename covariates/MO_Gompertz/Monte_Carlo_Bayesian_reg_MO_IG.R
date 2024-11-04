### -----------
### Monte Carlo 
### de regressão MO inverse Gaussian
### Geração de dados e estimação Bayesiana
### -----------


#setwd('C:/Users/dioni/OneDrive - University of São Paulo/Doutorado em Estatística/2024.2/2_Topicos_Avancados_de_Pesquisa_I/R_code/covariates/MO_IG')

## ---
## 1. Base functions ====
## ---

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


# MO Inverse Gaussian

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

# ---
# 2: Algoritmo gerador de dados
# ---

gen.cure.reg.moig = function(n,as,bs,xa,xb,l){
  finv.moig = function(t,alpha,beta,lambda,unif){Ftmo_IG(t=t,alpha=alpha,beta=beta,lambda=lambda) - unif}
  ger.moig = function(alpha,beta,lambda,unif){t=uniroot(finv.moig,c(0,1000),tol=0.01,alpha=alpha,beta=beta,lambda=lambda,unif=unif);return(t$root)}
  
  ## individual parameters
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  p_base = 1-exp(2*a/b)
  p_cure = (l*p_base)/(l*p_base+1-p_base)
  
  t = rep(NA,n)
  
  for(i in 1:n){
    rmi = rbinom(n=1,size=1,prob=1-p_cure[i,])
    uni = runif(n=1,0,max=1-p_cure[i,])
    t[i] = ifelse(rmi==0, Inf, ger.moig(alpha=a[i,],beta=b[i,],lambda=l,unif=uni))
  }
  
  t_finite = ifelse(t==Inf,0,t)
  u2 = runif(n=n,0,max(t_finite))
  t2 = pmin(t,u2); delta = ifelse(t<u2,1,0)
  
  return(cbind(t2,delta,xa,xb))
}

## ---
## 3. Tables to save data ====
## ---

n_replicas = 500                         ## number of Monte Carlo replicates  
ncols_mc = 4                             ## numero of information in each replicate (mean,sd,bias,cp)
n_amostral = c(50, 100, 500, 1000, 5000) ## sample sizes to study

## Tables

## as = (a0, a1, a2)
moig_mc_a0a1a2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "sd", "bias", "cp"),3),n_amostral))


## bs = (b0, b1, b2)
moig_mc_b0b1b2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "sd", "bias", "cp"),3),n_amostral))

## lambda

moig_mc_lambda = array(data=0, dim = c(n_replicas, ncols_mc, length(n_amostral)),
                      dimnames = list(1:n_replicas,c("mean", "sd", "bias", "cp"),n_amostral))

## Censura

moig_mc_cens = array(data=0, dim = c(n_replicas, 1, length(n_amostral)))


## Scenarios
## Parameters information

a0 = c(-1, 0.5, 0.2)
b0 = c(-1.1, 1.8, 0.8)
l0 = 0.5

params=c(a0,b0,l0)


for(j in 1:length(n_amostral)){
  for(i in 1:n_replicas){
    print(paste("Realização Monte Carlo: ", i, ", tamanho amostral:", n_amostral[j]))
    
    ## covariates data
    xa1 = rbinom(n=n_amostral[j],size=1,prob=0.7); xa2 = runif(n=n_amostral[j],0,1)
    xb1 = rbinom(n=n_amostral[j],size=1,prob=0.5); xb2 = runif(n=n_amostral[j],0,1)
    xa0 = cbind(xa1, xa2)
    xb0 = cbind(xb1, xb2)
    
    dados.reg.moig = gen.cure.reg.moig(n=n_amostral[j],as=a0,bs=b0,xa=xa0,xb=xb0,l=l0)
    dados.reg.moig = as.data.frame(dados.reg.moig)
    
    # Data and stan function
    
    data_regmoig = list(N = dim(dados.reg.moig)[1],
                       p = dim(cbind(1,dados.reg.moig[,3:4]))[2],
                       s = dim(cbind(1,dados.reg.moig[,5:6]))[2],
                       time = dados.reg.moig[,1],
                       delta = dados.reg.moig[,2],
                       xa = cbind(1,dados.reg.moig[,3:4]),
                       xb = cbind(1,dados.reg.moig[,5:6]))
    
    ## Run the Bayesian model
    regmoigfit = stan(file = 'cod_regmoig.stan', data = data_regmoig,
                     chains = 1, iter=5000, warmup=1000)
    
    ## Extract posterior samples
    regmoigfit_post_samples = extract(regmoigfit)
    
    ## posterior mean 
    ma0=mean(regmoigfit_post_samples$as[,1])
    ma1=mean(regmoigfit_post_samples$as[,2])
    ma2=mean(regmoigfit_post_samples$as[,3])
    mb0=mean(regmoigfit_post_samples$bs[,1])
    mb1=mean(regmoigfit_post_samples$bs[,2])
    mb2=mean(regmoigfit_post_samples$bs[,3])
    mlamb=mean(regmoigfit_post_samples$lambda)
    post_mean=c(ma0,ma1,ma2,mb0,mb1,mb2,mlamb)
    
    ## posterior standard deviation
    sda0=sd(regmoigfit_post_samples$as[,1])
    sda1=sd(regmoigfit_post_samples$as[,2])
    sda2=sd(regmoigfit_post_samples$as[,3])
    sdb0=sd(regmoigfit_post_samples$bs[,1])
    sdb1=sd(regmoigfit_post_samples$bs[,2])
    sdb2=sd(regmoigfit_post_samples$bs[,3])
    sdlamb=sd(regmoigfit_post_samples$lambda)
    post_sd=c(sda0,sda1,sda2,sdb0,sdb1,sdb2,sdlamb)
    
    ## relative bias
    rbias = (post_mean-params)/params*100
    
    ## coverage probability (95%)
    q025a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.025);q975a0 = quantile(regmoigfit_post_samples$as[,1],probs=0.975)
    q025a1 = quantile(regmoigfit_post_samples$as[,2],probs=0.025);q975a1 = quantile(regmoigfit_post_samples$as[,2],probs=0.975)
    q025a2 = quantile(regmoigfit_post_samples$as[,3],probs=0.025);q975a2 = quantile(regmoigfit_post_samples$as[,3],probs=0.975)
    q025b0 = quantile(regmoigfit_post_samples$bs[,1],probs=0.025);q975b0 = quantile(regmoigfit_post_samples$bs[,1],probs=0.975)
    q025b1 = quantile(regmoigfit_post_samples$bs[,2],probs=0.025);q975b1 = quantile(regmoigfit_post_samples$bs[,2],probs=0.975)
    q025b2 = quantile(regmoigfit_post_samples$bs[,3],probs=0.025);q975b2 = quantile(regmoigfit_post_samples$bs[,3],probs=0.975)
    q025lamb = quantile(regmoigfit_post_samples$lambda,probs=0.025);q975lamb = quantile(regmoigfit_post_samples$lambda,probs=0.975)
    
    qinf025 = c(q025a0,q025a1,q025a2,q025b0,q025b1,q025b2,q025lamb)
    qsup975 = c(q975a0,q975a1,q975a2,q975b0,q975b1,q975b2,q975lamb)
    
    cpost95 = as.integer((params >= qinf025) & (params <= qsup975))
    cpost95
    
    ## salvando valores
    
    moig_mc_a0a1a2[i,1,j] = post_mean[1]; moig_mc_a0a1a2[i,2,j] = post_sd[1]; moig_mc_a0a1a2[i,3,j] = rbias[1]; moig_mc_a0a1a2[i,4,j] = cpost95[1]
    moig_mc_a0a1a2[i,5,j] = post_mean[2]; moig_mc_a0a1a2[i,6,j] = post_sd[2]; moig_mc_a0a1a2[i,7,j] = rbias[2]; moig_mc_a0a1a2[i,8,j] = cpost95[2]
    moig_mc_a0a1a2[i,9,j] = post_mean[3]; moig_mc_a0a1a2[i,10,j] = post_sd[3]; moig_mc_a0a1a2[i,11,j] = rbias[3]; moig_mc_a0a1a2[i,12,j] = cpost95[3]
    
    moig_mc_b0b1b2[i,1,j] = post_mean[4]; moig_mc_b0b1b2[i,2,j] = post_sd[4]; moig_mc_b0b1b2[i,3,j] = rbias[4]; moig_mc_b0b1b2[i,4,j] = cpost95[4]
    moig_mc_b0b1b2[i,5,j] = post_mean[5]; moig_mc_b0b1b2[i,6,j] = post_sd[5]; moig_mc_b0b1b2[i,7,j] = rbias[5]; moig_mc_b0b1b2[i,8,j] = cpost95[5]
    moig_mc_b0b1b2[i,9,j] = post_mean[6]; moig_mc_b0b1b2[i,10,j] = post_sd[6]; moig_mc_b0b1b2[i,11,j] = rbias[6]; moig_mc_b0b1b2[i,12,j] = cpost95[6]
    
    
    moig_mc_lambda[i,1,j] = post_mean[7]; moig_mc_lambda[i,2,j] = post_sd[7]; moig_mc_lambda[i,3,j] = rbias[7]; moig_mc_lambda[i,4,j] = cpost95[7]
    
    
    ## salvando o percentual censura dos dados
    moig_mc_cens[i,,j] = 1-(sum(dados.reg.moig[,2])/n_amostral[j])
    
    ## Salvando os resultados em um arquivo .csv
    if(i==n_replicas){
      write.csv2(moig_mc_a0a1a2[,,j], file = paste("resumos_Bayes_MC_reg_MO_IG/","moig_mc_a0a1a2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_b0b1b2[,,j], file = paste("resumos_Bayes_MC_reg_MO_IG/","moig_mc_b0b1b2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_lambda[,,j], file = paste("resumos_Bayes_MC_reg_MO_IG/","moig_mc_lambda_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_cens[,,j], file = paste("resumos_Bayes_MC_reg_MO_IG/","moig_mc_cens_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
    }
    
  }
}







