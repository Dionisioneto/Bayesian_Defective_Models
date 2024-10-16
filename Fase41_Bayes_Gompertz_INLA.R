### ------
### Fase 4: Estimação Bayesiana dos modelos dentro do INLA.
### 1. Modelo Gompertz
### ------

# [1. Modelo Gompertz];
# 2. Modelo Gaussiano-inverso;
# 3. Modelo Marshall-Olkin Gompertz;
# 4. Modelo Marshall-Olkin Gaussiano-inverso;
# 5. Modelo Kumaraswamy Gompertz;
# 6. Modelo Kumaraswamy Gaussiano inverso.

options(timeout=600)

install.packages("INLA", repos=c(getOption("repos"), INLA="inla.r-inla-download.org/R/stable"), dep=TRUE) 
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("graph", "Rgraphviz"), dep=TRUE)


library(INLA)
















