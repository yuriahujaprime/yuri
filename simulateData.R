# Data simulator for LASSO Gauss Copula
# Simulates datasets (X,Y,beta) assuming a naive correlation structure in X and user-specified correlation structure in beta
# Yuri Ahuja
# Last modified 8/21/2019


library(glmnet)
library(ggplot2)
library(plotly)
library(tidyr)
library(doParallel)
library(extraDistr)
library(mvtnorm)
library(bindata)

registerDoParallel(4)

source("utils.r")
source("load_kmer.r")

setwd("~/Documents/Prime Discoveries")


expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# Simulates data
# Inputs:
# N = number of individuals; k = number of features (i.e. kmers)
# nnz = number of nonzero betas to impose sparsity
# rho = correlation between pairwise betas ([beta1,beta2], [beta3,beta4], etc.)
# Outputs:
# X = Nxk matrix of features (i.e. kmer counts)
# Y = N-dimensional vector of outcomes (i.e. disease status)
# beta = k-dimensional dimension of regression coefficients (note that this does NOT include the intercept coefficient)
simulate <- function(N,k,nnz,rho=0.9){
  nnz <- 2*round(nnz/2)
  means = rmultinom(1,100*k,rep(1/k,k))
  
  sigma = matrix(0,k,k)
  for (i in 1:(k/4)){
    l = 4*i-3
    u =4*i
    sigma[l:u,l:u] = 0.25*(100*k*.25^6)^2
  }
  diag(sigma) = (100*k*.25^6)^2
  
  X <- rmvnorm(N,mean=means,sigma=sigma)
  
  betas <- rep(0,k)
  beta_sigma <- matrix(0,nnz,nnz)
  for (i in 1:(nnz/2)){
    beta_sigma[(2*i-1),2*i] <- beta_sigma[2*i,(2*i-1)] <- rho * (4/nnz)^2
  }
  diag(beta_sigma) <- (4/nnz)^2
  betas[1:nnz] <- rmvnorm(1,rep(0,nnz),beta_sigma)
  
  Yprob <- X %*% betas
  Yprob = expit(Yprob - mean(Yprob))
  Y <- rbinom(N,1,Yprob)
  
  return(list("X"=X,"Y"=Y,"beta"=betas))
}


# Example of use
data = simulate(100,64,8,0.9)
X = data$X
Y = data$Y
betas = data$beta
write.csv(X,"X.csv",row.names=FALSE)
write.csv(Y,"Y.csv",row.names=FALSE)
write.csv(betas,"betas.csv",row.names=FALSE)
