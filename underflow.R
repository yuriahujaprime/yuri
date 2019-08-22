# Workbook for experiments regarding dealing with underflow in mbx data
# Eric should have everything in here already but please contact me if any questions arise
# Yuri Ahuja
# Last modified 8/21/19


install.packages('VGAM')
install.packages('dirmult')
install.packages('MGLM')

library(ggplot2)
library(plotly)
library(dplyr)
library(ggpubr)
library(plyr)
library(tidyr)
library(VGAM)
library(MASS)
library(dirmult)
library(MGLM)
library(MASS)

setwd("~/Documents/Prime Discoveries")

probiotic_taxa = c('g__Propionibacterium', 'g__Lactobacillus', 'g__Lactococcus',
                   'g__Bifidobacterium', 'g__Roseburia', 'g__Faecalibacterium',
                   'g__Ruminococcus', 'g__Akkermansia', 'g__Aldercreutzia')
tax_level = 'genus'

valdata = read.table(paste0("~/Documents/Prime Discoveries/Validation/", tax_level, "_level_table.tsv"), sep='\t', header=TRUE, row.names=1, check.names=FALSE, comment.char='', skip=1)


# Data preprocessing

taxa_names = c()
for (rowname in rownames(valdata)) {
  taxa = strsplit(rowname, ';')[[1]]
  
  tax = taxa[length(taxa)]
  for (i in length(taxa):1) {
    if (nchar(taxa[i]) <= 3)
      tax = taxa[i-1]
  }
  taxa_names = c(taxa_names, tax)
}

mbx04_replicates = data.frame(tax=taxa_names)
mbx04_replicates = cbind(mbx04_replicates, valdata[,substr(names(valdata), 1, 5) == 'MBX04'])
mbx04_replicates = cbind(mbx04_replicates,valdata[,61:70])
rownames(mbx04_replicates) = rownames(valdata)
apply(mbx04_replicates[,2:ncol(mbx04_replicates)], 1, prod)
a = mbx04_replicates[(apply(mbx04_replicates[,2:ncol(mbx04_replicates)], 1, prod) == 0) & (rowSums(mbx04_replicates[,2:ncol(mbx04_replicates)]) > 0),]

mbx04_replicates_probiotic = mbx04_replicates[mbx04_replicates$tax %in% probiotic_taxa,]

z = apply(mbx04_replicates[,2:ncol(mbx04_replicates)], 1, sd)
z = z[z != 0]

x = mbx04_replicates[mbx04_replicates$tax == 'f__Ruminococcaceae',]
x = as.vector(as.matrix(x[2:ncol(x)]))

y = mbx04_replicates[mbx04_replicates$tax == 'o__Clostridiales',]
y = as.vector(as.matrix(y[2:ncol(y)]))


# BETA BINOMIAL MOM INFERENCE

n <- 8e4 # Number of raw FASTA reads per sample
m1 <- n * as.numeric(rowMeans(mbx04_replicates[,-1]))
m2 <- m1^2 + n^2 * as.numeric(apply(mbx04_replicates[,-1],1,var))
alpha <- (n*m1 - m2) / (n*(m2/m1 - m1 - 1) + m1)
beta <- (n - m1) * (n - m2/m1) / (n*(m2/m1 - m1 - 1) + m1)
BBMOM_params <- rbind(alpha,beta) # 2 x K matrix of (alpha,beta) estimators for each bacterium


# Quick and dirty means of interpolating raw read count matrix (X) from proportions (mbx04_replicates_num)

mbx04_replicates_num <- as.matrix(mbx04_replicates[,-1])
multipliers <- apply(mbx04_replicates_num,2,function(vec){2/min(vec[vec>0])})
multipliers[multipliers>120000] <- multipliers[multipliers>120000]/2
X <- multipliers * t(mbx04_replicates_num)
X <- X[,colSums(X)>0]


# DIRICHLET-MULTINOMIAL VS MULTINOMIAL VS NEGATIVE BINOMIAL MLE INFERENCE

set.seed(123)
mnFit <- MGLMfit(X,dist="MN")
dmFit <- MGLMfit(X,dist="DM")
dmFit2 <- dirmult(X)
gdmFit <- MGLMfit(X,dist="GDM")
#gdmFit <- MGLMfit(X[,-underdispersed],dist="GDM")
negmnFit <- MGLMfit(X,dist="NegMN")

nEffective <- 1/dmFit2$theta
dmFitCoefs <- dmFit2$gamma


# Plot of fitted alpha value versus proportion of reads observed to be 0 for bacteria with a mix of zeros and non-zeros
atLeastOneZero <- which(apply(X,2,function(vec){any(vec==0)}))
propNonZeros <- apply(X[,atLeastOneZero],2,function(vec){mean(vec>0)})

png('DM_alpha_vs_propNonZero.png')
plot(propNonZeros, dmFit2$gamma[atLeastOneZero], xlab="Proportion of Reads =/= 0", ylab="Dirichlet Multinomial Fitted alpha")
dev.off()


# BETA-BINOMIAL MLE

bbFit <- sapply(1:ncol(X),function(i){
  print(i)
  fit <- vglm(cbind(X[,i],rowSums(X[,-i])) ~ 1, betabinomial)
  coefs <- Coef(fit)
  alpha <- coefs[1] * (1-coefs[2])/coefs[2]
  beta <- (1-coefs[1]) * (1-coefs[2])/coefs[2]
  c(alpha,beta)
})


# How much does actual variance deviate from that expected under MN/DM model assumptions?

n <- mean(rowSums(X))
alphas <- dmFit2$gamma
pis <- X/rowSums(X)
pi_hats <- colMeans(pis)
observed_vars <- as.vector(apply(pis,2,var))

mn_expected_vars <- as.vector(pi_hats * (1-pi_hats) / n)
underdispersed <- which(observed_vars-mn_expected_vars <= 0)

dm_expected_vars <- as.vector(sapply(1:ncol(X),function(i){
  1/n*alphas[i]/sum(alphas) * (1 - alphas[i]/sum(alphas)) * ((n+sum(alphas))/(1+sum(alphas)))
}))

alphas_sc10000 <- alphas*10000/sum(alphas)
dm_expected_vars_sc10000 <- as.vector(sapply(1:ncol(X),function(i){
  1/n*alphas_sc10000[i]/sum(alphas_sc10000) * (1 - alphas_sc10000[i]/sum(alphas_sc10000)) * ((n+sum(alphas_sc10000))/(1+sum(alphas_sc10000)))
}))


strepto <- grep('g__Streptococcus',colnames(X))

# Plot of expected (under multinomial) versus observed variance
summary(lm(log(observed_vars)~log(mn_expected_vars)))
png("MN_expected_vs_observed_variance.png")
plot(log(mn_expected_vars),log(observed_vars),xlab="Expected Variance Under Multinomial",ylab="Observed Variance")
lines(seq(-30,0),seq(-30,0),col='blue')
lines(log(mn_expected_vars), fitted(lm(log(observed_vars)~log(mn_expected_vars))), col='red',add=T)
#lines(rep(log(mn_expected_vars)[strepto],31),seq(-30,0), col='black')    # Uncomment to identify point corresponding to 'g__Streptococcus'
dev.off()


# Plot of expected (under dirichlet-multinomial) versus observed variance
summary(lm(log(observed_vars)~log(dm_expected_vars)))
png("DM_expected_vs_observed_variance.png")
plot(log(dm_expected_vars),log(observed_vars),xlab="Expected Variance Under Dirichlet-Multinomial",ylab="Observed Variance")
lines(seq(-30,0),seq(-30,0),col='blue')
lines(log(dm_expected_vars), fitted(lm(log(observed_vars)~log(dm_expected_vars))), col='red',add=T)
#lines(rep(log(dm_expected_vars)[strepto],31),seq(-30,0), col='black')    # Uncomment to identify point corresponding to 'g__Streptococcus'
dev.off()

