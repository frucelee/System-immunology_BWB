#!/usr/bin/env Rscript
# File: machine_learning.R
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

library(data.table)
library(purrrlyr)
library(dplyr)
library(caret)
library(glmnet)

load("Resistomics_data.Rdata")
path<-setwd("S1/")
fileNames = list.files(path=path,pattern=".txt", full.names = TRUE)
data_al<-t(data.frame(fileNames))
colnames(data_al)<-data_al[1,]
colnames(data_al) <- sub('S1/mc_', '', colnames(data_al))
colnames(data_al) <- gsub('_loco_S1_genotype.txt','',colnames(data_al))
names(fileNames)<-colnames(data_al)
length(colnames(data_al))
result <- matrix(nrow =104,
                     ncol = 3)
result<-data.frame(result)
colnames(result )<-c("r","RMSE","p")
for (j in c(names(fileNames))){
  file= paste('S1/mc_', j, sep='')
  file<-paste(file,"_loco_S1_genotype.txt",sep="")
  gumlm<-fread(file,header=T)
  write.csv(gumlm, 'tmp.csv', quote = FALSE)
  data_snp <- read.csv('tmp.csv',header = TRUE)
  jj<-as.numeric(j)
  mydata<-mydata[data_snp$IID,]
  tmp<-mydata[,c(3:24)]
  data_snp$V500<-tmp$V500
  data_snp$ID<-mydata$ID
  data_snp <- data_snp[,-1]
  data_snp<-na.omit(data_snp)
  tmp_data <- matrix(nrow =1 ,ncol = 3)
  mc_test <- read.csv('S1_sample_test_left_.csv',header = TRUE)
  mc_test$ID <- sub("^", "BE", mc_test$ID)
  for (i in 1){
    mc_train <- data_snp[!c(data_snp$ID%in%mc_test$ID), ]
    mc_test <- data_snp[data_snp$ID%in%mc_test$ID, ]
    mc_train<-subset(mc_train, select = - c(ID,IID))
    mc_test<-subset(mc_test, select = - c(ID,IID))
    lambda.grid <- seq(0, 100)
    alpha.grid <- seq(0, 0.5, length = 6)
    trnCtrl = trainControl(
      method = "repeatedcv",
      number = 5,
      repeats = 5)
    srchGrd = expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)
    my.train <- train(x = as.matrix(mc_train[,-which( colnames(mc_train)=="V500" )]),
                      y = mc_train$V500,
                      method = "glmnet",
                      tuneGrid = srchGrd,
                      trControl = trnCtrl,
                      standardize = FALSE,
                      maxit = 1000000)
    ##my.train$finalModel
    ##my.train$results
    ##my.train$bestTune
    prediction<-predict(my.train,mc_test[,-which( colnames(mc_test)=="V500" )])
    cortest<-cor.test(method="spearman",mc_test$V500,prediction)
    r<-cortest$estimate[[1]]
    tmp_data[i, 1] <- r
    RMSE<-sqrt(mean((mc_test$V500 - prediction)^2))
    tmp_data[i, 2] <- RMSE
    tmp_data[i, 3] <- cortest$p.value
  }
  tmp_data[tmp_data == 0] <- NA
  result[jj,1]<-tmp_data[1, 1]
  result[jj,2]<-tmp_data[1, 2]
  result[jj,3]<-tmp_data[1, 3]
}
write.csv(result,"T1_MC_SNP_only_results.csv",quote = FALSE)