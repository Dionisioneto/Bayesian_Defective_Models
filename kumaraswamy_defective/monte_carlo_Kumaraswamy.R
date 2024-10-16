## Simulação Monte Carlo
## Modelo Marshall-Olkin Gompertz

## N = 50, 100, 1.000 e 10.000.
## Censoring(%) = 15%, 45%

# functions
#source("")

## ---
## Funções importantes ===
## ---

bias = function(true,est){(est-true)/true*100}

## 1. Scenario: 

# Alpha = -1.2, Beta = 2, Lambda = 0.8, p = 0.1570 

n.replicas = 500

## colunas: a média a posteriori, desvio padrão a posteriori, bias, coverage  
ncols_mc = 4
ncols_mc_medidas = 4


m_mc = matrix(data=0, nrow = n.replicas, ncol = ncols_mc)


# for(i in 1:n.replicas){
#   
# }

## 2. Scenario: 

sumari = summary(mogfit)
  
sumari$summary[,"n_eff"]
sumari$summary[,"Rhat"]
  
# Alpha = -2, Beta = 5, Lambda = 2, p = 0.4172 






















