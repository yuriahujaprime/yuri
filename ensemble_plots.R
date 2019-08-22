# Code for plotting results from classifiers (ensemble and originals) on various diseases (CD, UC, CC etc.)
# Yuri Ahuja
# Last modified sometime in July 2019
# Not really intended for general use - contact me if you want to use this

library(ggplot2)
library(reshape2)
library(ggpubr)


# Plot ensemble results

ibd_results <- relative_weighting(ibd_experiment,labels_ibd)
hmp_results <- rbind(relative_weighting(hmp_experiment,labels_hmp),0)
pascal_results <- relative_weighting(pascal_experiment,labels_pascal)

ibdacc <- data.frame(alpha=rep(colnames(ibd_results),2), Accuracy=c(ibd_results[1,],pascal_results[1,]), Source=c(rep('IBD',11),rep('Pascal',11)))
crohnauc <- data.frame(alpha=rep(colnames(ibd_results),3), AUC=c(ibd_results[2,],pascal_results[2,],hmp_results[2,]), Source=c(rep('IBD',11),rep('Pascal',11),rep('HMP',11)))
ucauc <- data.frame(alpha=rep(colnames(ibd_results),2), AUC=c(ibd_results[3,],pascal_results[3,]), Source=c(rep('IBD',11),rep('Pascal',11)))
canceracc <- data.frame(alpha=rep(colnames(ibd_results),3), Accuracy=c(cancer_weighting_experiment$crc280[1,],cancer_weighting_experiment$crc607[1,],cancer_weighting_experiment$crc678[1,]), Source=c(rep('CRC 280',11),rep('CRC 607',11), rep('CRC 678',11)))
cancerauc <- data.frame(alpha=rep(colnames(ibd_results),3), AUC=c(cancer_weighting_experiment$crc280[2,],cancer_weighting_experiment$crc607[2,],cancer_weighting_experiment$crc678[2,]), Source=c(rep('CRC 280',11),rep('CRC 607',11), rep('CRC 678',11)))
adenomaauc <- data.frame(alpha=colnames(ibd_results), AUC=c(cancer_weighting_experiment$crc280[3,]))

altogether <- data.frame(alpha=rep(colnames(ibd_results),9),
                         out=c(ibd_results[2,],pascal_results[2,],hmp_results[2,], ibd_results[3,],pascal_results[3,],
                               cancer_weighting_experiment$crc280[2,],cancer_weighting_experiment$crc607[2,],cancer_weighting_experiment$crc678[2,],
                               cancer_weighting_experiment$crc280[3,]),
                         line=as.vector(sapply(1:9,function(i){rep(i,11)})),
                         disease=c(rep('Crohns AUC',33), rep('UC AUC',22), rep('Carcinoma AUC',33), rep('Adenoma AUC',11)))


ibd_dataframe <- data.frame(alpha=colnames(ibd_results), Accuracy=ibd_results[1,], AUCCD=ibd_results[2,], AUCUC=ibd_results[3,])
hmp_dataframe <- data.frame(alpha=colnames(ibd_results), Accuracy=hmp_results[1,], AUCCD=hmp_results[2,])
pascal_dataframe <- data.frame(alpha=colnames(pascal_results), Accuracy=pascal_results[1,], AUCCD=pascal_results[2,], AUCUC=pascal_results[3,])
crc258_dataframe <- data.frame(alpha=colnames(cancer_weighting_experiment$crc258), Accuracy=cancer_weighting_experiment$crc258[1,], AUCCarcinoma=cancer_weighting_experiment$crc258[2,])
crc280_dataframe <- data.frame(alpha=colnames(cancer_weighting_experiment$crc280), Accuracy=cancer_weighting_experiment$crc280[1,], AUCCarcinoma=cancer_weighting_experiment$crc280[2,], AUCAdenoma=cancer_weighting_experiment$crc280[3,])
crc445_dataframe <- data.frame(alpha=colnames(cancer_weighting_experiment$crc445), Accuracy=cancer_weighting_experiment$crc445[1,], AUCCarcinoma=cancer_weighting_experiment$crc445[2,])
crc445_dataframe <- data.frame(alpha=colnames(cancer_weighting_experiment$crc607), Accuracy=cancer_weighting_experiment$crc607[1,], AUCCarcinoma=cancer_weighting_experiment$crc607[2,])


# IBD accuracy

p_ibdacc <- ggplot(ibdacc, aes(x=alpha,y=Accuracy,color=Source,group=Source)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='Accuracy', labs='Data Source')
p_ibdacc

p_crohnauc <- ggplot(crohnauc, aes(x=alpha,y=AUC,color=Source,group=Source)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='AUROC', labs='Data Source')
p_crohnauc
  
p_ucauc <- ggplot(ucauc, aes(x=alpha,y=AUC,color=Source,group=Source)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='AUROC', labs='Data Source')
p_ucauc  

p_canceracc <- ggplot(canceracc, aes(x=alpha,y=Accuracy,color=Source,group=Source)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='Accuracy', labs='Data Source')
p_canceracc

p_cancerauc <- ggplot(cancerauc, aes(x=alpha,y=AUC,color=Source,group=Source)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='AUROC', labs='Data Source')
p_cancerauc

p_adenomaauc <- ggplot(adenomaauc, aes(x=alpha,y=AUC,group=1)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='AUROC')
p_adenomaauc  

p_altogether <- ggplot(altogether, aes(x=alpha,y=out,color=disease,group=line)) +
  geom_line() +
  labs(x='Proportion Logistic Regression', y='AUROC', labs='Disease/Measure') +
  scale_fill_brewer(palette="Set1")
p_altogether
  
  