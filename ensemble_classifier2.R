# LASSO/Random Forest Ensemble Classifier
# Implements classifier ensembling LASSO-regularized logistic regression and random forest
# Yuri Ahuja
# Last modified 8/21/2019


library(glmnet)
library(ggplot2)
library(plotly)
library(tidyr)
library(doParallel)
library(randomForest)

registerDoParallel(4)

source("utils.r")
source("load_kmer.r")
source("load_cohort.r")


# Specify studies to read into memory
studies = c('amerigut','pascal','hmp','ibd','crc258','crc280','crc445','crc607','crc678')
load_kmer(studies=studies)


# Choose datasets/classes to load for subsequent experimentation
datasets = list(list(crc258.kmer.healthy, crc258.kmer.cancer),
                list(crc280.kmer.healthy, crc280.kmer.cancer, crc280.kmer.adenoma),
                list(crc445.kmer.healthy, crc445.kmer.cancer),
                list(crc607.kmer.healthy, crc607.kmer.cancer),
                list(crc678.kmer.healthy, crc678.kmer.cancer))
classes = c("healthy", "cancer", "adenoma")
plot.title = NULL
legend.title = "pheno"
class.order = 1:2
mycolors = c(primeyellow, primegreen)
fam = 'multinomial'

features = list()
labels = list()
for (i in 1:length(datasets)) {
  study = studies[i]
  features[[study]] = matrix(NA, nrow=0, ncol=ncol(datasets[[i]][[1]]))
  labels[[study]] = c()
  for (j in 1:length(datasets[[i]])){
    features[[study]] = rbind(features[[i]], datasets[[i]][[j]])
    labels[[study]] = c(labels[[study]], rep(j, nrow(datasets[[i]][[j]])))
  }
  labels[[study]] = as.factor(labels[[study]])
}


# Cross-validates and runs LASSO-regularized logistic regression and random forest (see example below)
# Automatically determines from number of unique labels whether problem is binomial or multinomial
# Inputs:
# features = nPatients x nKmers matrix of patient-level kmer frequencies
# labels = nPatients-dimensional vector of patient disease labels (can be factors or integers)
# Outputs:
# lr = (binary) nPatients-dimensional vector or (multiclass) nPatients x nClasses matrix of LASSO-predicted disease probabilities
# rf = (binary) nPatients-dimensional vector or (multiclass) nPatients x nClasses matrix of random forest-predicted disease probabilities
# results = (binary) 3x2 matrix of mean accuracies and AUCs or (multiclass) 3x1 matrix of mean accuracies for LASSO, random forest, and ensemble classifiers
runExperimentEnsNew <- function(features, labels){
  cv.logistic = cv.glmnet(features, labels, family=fam, type.measure='class', nfolds=length(labels), alpha=1, grouped=FALSE, keep=FALSE, parallel=TRUE)
  lambda = cv.logistic$lambda.min
  rf = randomForest(features, labels, ntree=1000)
  
  # Binary outcome case
  if (length(unique(labels)) == 2){
    probs_rf = rf$votes[,2]
    
    probs_lr = rep(0,length(labels))
    for (i in 1:length(labels)) {
      print(paste0("Progress ",i,"/",length(labels)))
      
      features.train = features[-i,]
      labels.train = labels[-i]
      features.test = features[i,]
      labels.test = labels[i]
      
      logistic = glmnet(features.train, labels.train, family=fam, lambda=lambda, alpha=1)
      probs_lr[i] = predict(logistic,t(as.matrix(features.test)), type="response")
    }
    
    probs_ens = 0.5*(probs_lr+probs_rf)
    
    results = matrix(0,3,2)
    rownames(results) = c("Lasso","RandFor","Ensemble")
    colnames(results) = c("Accuracy", "AUC")
    labels_numeric = as.numeric(labels)-1
    results[1,1] = mean(labels_numeric == round(probs_lr))
    results[1,2] = auc(labels_numeric, probs_lr)
    results[2,1] = mean(labels_numeric == round(probs_rf))
    results[2,2] = auc(labels_numeric, probs_rf)
    results[3,1] = mean(labels_numeric == round(probs_ens))
    results[3,2] = auc(labels_numeric, probs_ens)
  }
  
  # Multiclass outcome case
  else{
    probs_rf = rf$votes[,]
    
    probs_lr = matrix(nrow=length(labels), ncol=length(classes))
    for (i in 1:length(labels)) {
      print(paste0(i,"/",length(labels)))
      
      features.train = features[-i,]
      labels.train = labels[-i]
      features.test = features[i,]
      labels.test = labels[i]
      
      logistic = glmnet(features.train, labels.train, family=fam, lambda=lambda, alpha=1)
      probs_lr[i,] = as.vector(predict(logistic, t(as.matrix(features.test)), type='response'))
    }
    
    results = matrix(0,3,1)
    rownames(results) = c("Lasso","RandFor","Ensemble")
    colnames(results) = c("Accuracy")
    labels_numeric = as.numeric(labels)
    results[1,1] = mean(labels_numeric == apply(probs_lr, 1, which.max))
    results[2,1] = mean(labels_numeric == apply(probs_rf, 1, which.max))
    results[3,1] = mean(labels_numeric == apply(0.5 * (probs_lr + probs_rf), 1, which.max))
  }

  return(list("lr"=probs_lr, "rf"=probs_rf, "results"=results))
}


# Example use
experiment <- runExperimentEnsNew(features,labels)
results <- ibd_experiment$results
print(results)


# Experiments with different relative weightings of LASSO predicted probabilities and random forest predicted probabilities (see example below)
# in generating ensemble predicted probabilities
# alpha = weight of logistic regression prediction (Ensemble = alpha*LASSO + (1-alpha)*RF)
# Inputs:
# experiment: unaltered output of runExperimentEnsNew function
# labels: nPatients-dimensional vector of disease labels (again can be factors or integers)
# Output:
# (binary) 2 x 11 matrix with accuracy and AUC on rows, alphas = {0.1,0.2,...,1.0} on columns
# (multiclass) (nClasses+1) x 11 matrix with accuracy and AUC for each outcome on rows, alphas = {0.1,0.2,...,1.0} on columns
relative_weighting <- function(experiment,labels){
  alphas <- seq(0,1,.1)
  
  # Multiclass outcome case
  if (length(unique(labels)) > 2){
    result <- matrix(nrow=(length(unique(labels))+1), ncol=length(alphas))
    rownames(result) <- c("Accuracy",paste(as.character(levels(labels)),"AUC"))
    colnames(result) <- alphas
    labels_numeric <- as.numeric(labels)

    for (i in 1:length(alphas)){
      alpha <- alphas[i]
      prob <- alpha*experiment$lr + (1-alpha)*experiment$rf
      result[1,i] <- mean(labels_numeric == apply(prob, 1, which.max))
      for (j in 1:length(unique(labels))){
        result[(j+1),i] <- auc(labels_numeric==j,prob[,j])
      }
    }
  }
  
  # Binary outcome case
  else{
    result <- matrix(nrow=2, ncol=length(alphas))
    colnames(result) <- alphas
    rownames(result) <- c("Accuracy","AUC")
    labels_numeric <- as.numeric(labels)
    
    for (i in 1:length(alphas)){
      alpha <- alphas[i]
      prob <- alpha*experiment$lr + (1-alpha)*experiment$rf
      result[1,i] <- mean(labels_numeric-1 == round(prob))
      result[2,i] <- auc(labels_numeric==2,prob)
    }
  }
  
  return(result)
}


# Example use
relative_weighting(experiment,labels)


