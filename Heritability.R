#!/usr/bin/env Rscript
# File: Heritability.R
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

library(data.table)
library(ggplot2)

load("Resistomics_normalized_removal_outliers.Rdata")
##Resampling 10 times for the 90% of total dataset
for (j in 1:10){
  select_train <- sample(246, 246*0.9)
  tmp <- mydata1[select_train, ]
  mydata_tmp<-tmp[,-c(1:13)]
  mydata_tmp$FID<-tmp$ID
  mydata_tmp$IID<-c(rep(1,221))
  mydata_tmp<-mydata_tmp[,c(146,145,1:144)]
  mydata_tmp$FID <- sub("^", "BE", mydata_tmp$FID )
  pheno_j<-paste(j,"pheno_",sep = "_")
  write.table(mydata_tmp,paste(pheno_j,"txt",sep = "."), quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  GCTA_qcov<-mydata_tmp[,c(1:2)]
  GCTA_qcov$age<-tmp[,13]
  GCTA_qcov1<-paste(j,"GCTA_qcov",sep = "_")
  write.table(GCTA_qcov,paste(GCTA_qcov1,"txt",sep = "."), quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  GCTA_cov<-mydata_tmp[,c(1:2)]
  GCTA_cov$age<-tmp[,8]
  GCTA_cov1<-paste(j,"GCTA_cov",sep = "_")
  write.table(GCTA_cov,paste(GCTA_cov1,"txt",sep = "."), quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  GCTA_ID<-mydata_tmp[,c(1:2)]
  GCTA_ID1<-paste(j,"GCTA_ID",sep = "_")
  write.table(GCTA_ID,paste(GCTA_ID1,"txt",sep = "."), quote = FALSE,row.names = FALSE,col.names = FALSE)
}

##permutation 1000 times for each immunophenotypes
data0<-mydata1[,c(8,13)]
mydata1<-mydata1[,c(25:168)]

for (j in c(1:144)){
  data0$V500<-mydata1[,j]
  nsim <- 1000
  res <- numeric(nsim) ## set aside space for results
  tmp_table <- matrix(nrow =246 ,
                           ncol = 1000)
  tmp_table<-data.frame(tmp_table)
  for (i in 1:nsim) {
    ## standard approach: scramble response value
    perm <- sample(nrow(data0))
    bdat <- transform(data0,colonies=V500[perm])
    tmp_table[, i] <- bdat[,4]
  }
  tmp_table$FID<-tbl3$V2
  tmp_table$IID<-c(rep(1,246))
  tmp_table<-tmp_table[,c(1002,1001,1:1000)]
  tmp_table$FID <- sub("^", "BE", tmp_table$FID )
  write.table(tmp_table,paste(j,"txt",sep = "."), quote = FALSE,row.names = FALSE,col.names = FALSE)
}


##compare real estimate and permuted one
path<-setwd("Resistomics\\H2\\perm")
fileNames = list.files(path=path,pattern="*.txt", full.names = TRUE)
tmp_table <- matrix(nrow =144,
                         ncol = 3)
tmp_table<-data.frame(tmp_table)
colnames(tmp_table)<-c("ID","mean_perm","min_top50")
for(i in 1:length(fileNames)){
  data<-read.table(fileNames[i],header=F)
  tmp_table[i,2]<-mean(data$V2)
  n <- 5
  tmp<-data[data$V2 > quantile(data$V2,prob=1-n/100),]
  tmp_table[i,3]<-min(tmp$V2)
  tmp_table[i,1]<-fileNames[i]
}
write.csv(tmp_table, 'h2_perum.csv', quote = FALSE)

