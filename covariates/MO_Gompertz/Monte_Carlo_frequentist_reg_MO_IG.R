### -----------
### Monte Carlo MLE 
### de regressão MO inverse Gaussian
### Geração de dados e estimação via Máxima Verossimilhança
### -----------

# 1. Base functions

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

# ---
# 2: Algoritmo gerador de dados ====

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

## Parameters information
n=5000

a0 = c(-1, 0.5, 0.2)
b0 = c(-1.1, 1.8, 0.8)
l0 = 0.5


## Step 1: Build the log-likelihood function

loglik_regmoig = function(par,time,delta,xa,xb){
  #parameters
  as = par[1:(dim(xa)[2]+1)] # b0, b1, b2
  bs = par[(dim(xa)[2]+2):(dim(xa)[2]+1+dim(xb)[2]+1)] # d0, d1, d2
  l = par[(dim(xa)[2]+1+dim(xb)[2]+1)+1] # l0
  
  a = cbind(1,xa) %*% as
  b = exp(cbind(1,xb) %*% bs)
  
  den = as.vector(ftmo_IG(t=time,alpha=a,beta=b,lambda=l))
  surv = as.vector(Stmo_IG(t=time,alpha=a,beta=b,lambda=l))
  
  log_lik = sum(delta*log(den) + (1-delta)*log(surv))
  return(-1*log_lik)
}


## 4. Tables to save data ====

n_replicas = 500                        ## number of Monte Carlo replicates  
ncols_mc = 4                             ## numero of information in each replicate (mean,sd,bias,cp)
n_amostral = c(50, 100, 500, 1000, 5000) ## sample sizes to study

## Tables

## as = (a0, a1, a2)
moig_mc_a0a1a2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "se", "bias", "cp"),3),n_amostral))


## bs = (b0, b1, b2)
moig_mc_b0b1b2 = array(data=0, dim = c(n_replicas, ncols_mc*3, length(n_amostral)),
                      dimnames = list(1:n_replicas,rep(c("mean", "se", "bias", "cp"),3),n_amostral))

## lambda

moig_mc_lambda = array(data=0, dim = c(n_replicas, ncols_mc, length(n_amostral)),
                      dimnames = list(1:n_replicas,c("mean", "se", "bias", "cp"),n_amostral))

## Censura

moig_mc_cens = array(data=0, dim = c(n_replicas, 1, length(n_amostral)))


## ---
## Scenarios: 
## ---

## (Parameters information (45%)

a0 = c(-1, 0.5, 0.2)
b0 = c(-1.1, 1.8, 0.8)
l0 = 0.5

params = c(a0,b0,l0)
start = params-0.15


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
    
    # MLE estimator via optim
    maxlike_regmoig = optim(par = start,
                           fn = loglik_regmoig,
                           gr = NULL,
                           hessian = T,
                           method = "BFGS",
                           time = dados.reg.moig$t2,
                           delta = dados.reg.moig$delta,
                           xa=cbind(dados.reg.moig$xa1,dados.reg.moig$xa2),
                           xb=cbind(dados.reg.moig$xb1,dados.reg.moig$xb2))
    
    ## point estimation (MLE estimate)
    mles = maxlike_regmoig$par
    
    ## relative bias 
    rbias = (mles-params)/params
    
    ## standard error
    se = sqrt(diag(solve(maxlike_regmoig$hessian)))
    coverage = se*1.96
    
    inf95 = maxlike_regmoig$par - coverage
    sup95 = maxlike_regmoig$par + coverage
    cp95 = as.integer((params >= inf95) & (params <= sup95))
    
    
    ## salvando valores
    
    moig_mc_a0a1a2[i,1,j] = mles[1]; moig_mc_a0a1a2[i,2,j] = se[1]; moig_mc_a0a1a2[i,3,j] = rbias[1]; moig_mc_a0a1a2[i,4,j] = cp95[1]
    moig_mc_a0a1a2[i,5,j] = mles[2]; moig_mc_a0a1a2[i,6,j] = se[2]; moig_mc_a0a1a2[i,7,j] = rbias[2]; moig_mc_a0a1a2[i,8,j] = cp95[2]
    moig_mc_a0a1a2[i,9,j] = mles[3]; moig_mc_a0a1a2[i,10,j] = se[3]; moig_mc_a0a1a2[i,11,j] = rbias[3]; moig_mc_a0a1a2[i,12,j] = cp95[3]
    
    moig_mc_b0b1b2[i,1,j] = mles[4]; moig_mc_b0b1b2[i,2,j] = se[4]; moig_mc_b0b1b2[i,3,j] = rbias[4]; moig_mc_b0b1b2[i,4,j] = cp95[4]
    moig_mc_b0b1b2[i,5,j] = mles[5]; moig_mc_b0b1b2[i,6,j] = se[5]; moig_mc_b0b1b2[i,7,j] = rbias[5]; moig_mc_b0b1b2[i,8,j] = cp95[5]
    moig_mc_b0b1b2[i,9,j] = mles[6]; moig_mc_b0b1b2[i,10,j] = se[6]; moig_mc_b0b1b2[i,11,j] = rbias[6]; moig_mc_b0b1b2[i,12,j] = cp95[6]
    
    
    moig_mc_lambda[i,1,j] = mles[7]; moig_mc_lambda[i,2,j] = se[7]; moig_mc_lambda[i,3,j] = rbias[7]; moig_mc_lambda[i,4,j] = cp95[7]
    
    
    ## salvando o percentual censura dos dados
    moig_mc_cens[i,,j] = 1-(sum(dados.reg.moig[,2])/n_amostral[j])
    
    ## Salvando os resultados em um arquivo .csv
    if(i==n_replicas){
      write.csv2(moig_mc_a0a1a2[,,j], file = paste("resumos_freq_MC_reg_MO_IG/","moig_mc_a0a1a2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_b0b1b2[,,j], file = paste("resumos_freq_MC_reg_MO_IG/","moig_mc_b0b1b2_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_lambda[,,j], file = paste("resumos_freq_MC_reg_MO_IG/","moig_mc_lambda_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
      write.csv2(moig_mc_cens[,,j], file = paste("resumos_freq_MC_reg_MO_IG/","moig_mc_cens_", "rep_",n_replicas, "_n_",n_amostral[j], ".csv", sep=""))
    }
  }
}



