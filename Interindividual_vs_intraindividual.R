#!/usr/bin/env Rscript
# File: Interdividual vs intraindividual.R
# Author: Shifang Li
# E-mail: fruceslee@gmail.com
rm(list=ls())
set.seed(1234)
load("subset_T1T2.Rdata")
# Create a new matrix with 28 rows and 4 columns:
IntravsInter_results <- matrix(nrow=144, ncol =4)
# Rename column and row names of the new matrix:
colnames(IntravsInter_results) <- c("p_of_model","r2_of_model",
                       "prop_attributed_to_sampling", "prop_attributed_to_ID")
rownames(IntravsInter_results) <- colnames(data[,c(15:158)])
data_sub<-data[,c(3,6,15:158)]
# Estimate the relative r2 contributions of inter-individual variation and intra-individual variation over time(T1 and T2)
# for each trait:
library(relaimpo)
library(RColorBrewer)
library(ellipse)
library(gplots)
library(vegan)
for(i in c(3:146)){
  model_vac02 <- lm(data_sub[,i] ~ data_sub[,"Sampling"]
                    + data_sub[,"ID"],na.action="na.exclude")
  IntravsInter_results[c(i - 3 + 1), 2] <- summary(model_vac02)$r.squared
  IntravsInter_results[c(i - 3 + 1), 1] <- pf(summary(model_vac02)$fstatistic[1],
                                 summary(model_vac02)$fstatistic[2],
                                 summary(model_vac02)$fstatistic[3], lower.tail = F)
  {
    ca01 <- calc.relimp(model_vac02, diff = T, rela = T)
    IntravsInter_results[c(i - 3 + 1), 3] <- ca01$lmg[grepl(pattern = "Sampling",
                                               x = names(ca01$lmg))]
    IntravsInter_results[c(i - 3 + 1), 4] <- ca01$lmg[grepl(pattern = "ID",
                                               x = names(ca01$lmg))]
    rm(ca01)
  }
  rm(model_vac02)
  cat(i, sep = "\t")
}
# Plotting all traits:
library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
IntravsInter_results<-data.frame(IntravsInter_results)
IntravsInter_results$p_of_model_adj<- p.adjust(IntravsInter_results$p_of_model, method = "bonferroni",
                                  n = length(IntravsInter_results$p_of_model))
write.csv(IntravsInter_results, 'IntravsInter_results.csv', quote = FALSE)
summary(IntravsInter_results$prop_attributed_to_ID)
summary(IntravsInter_results$prop_attributed_to_sampling)
summary(IntravsInter_results$r2_of_model)
IntravsInter_results$X<-row.names(IntravsInter_results)
names(IntravsInter_results)[names(IntravsInter_results)=="prop_attributed_to_ID"]="Interindividual"
names(IntravsInter_results)[names(IntravsInter_results)=="prop_attributed_to_sampling"]="Intraindividual"
nba.m <- melt(IntravsInter_results,id = c('X',"p_of_model","r2_of_model","p_of_model_adj"))
ggplot(nba.m, aes(x=X, y=value, fill=factor(variable, levels=c("Interindividual", "Intraindividual")))) +coord_flip()+
  geom_bar(stat="identity")+scale_fill_manual(values=c('lightgray','black'))+theme_classic()+theme(axis.text.x = element_text(angle = 90, size =0.6,hjust = 0.5, vjust = 0.5))+labs(x = '', y = 'Proportion of R2',fill="")+theme(legend.position="top")

IntravsInter_results$p_of_model_adj<-log10(IntravsInter_results$p_of_model_adj)
IntravsInter_results$p_of_model_adj<-IntravsInter_results$p_of_model_adj*(-1)
P1<-ggplot(data = IntravsInter_results,mapping = aes(x=X,y=r2_of_model)) + geom_point(shape=7,size=3)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+labs(x = '', y = 'R2')+ylim(0,1)
P2<-ggplot(data = IntravsInter_results,mapping = aes(x=X,y=p_of_model_adj)) + geom_point(shape=7,size=3)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+labs(x = '', y = '-log10 adjusted P')+ geom_hline(aes(yintercept= 1.30))
library(gridExtra)          
grid.arrange(P1,P2)
P2<-ggplot(data = IntravsInter_results,mapping = aes(x=X,y=p_of_model_adj)) + geom_point(shape=7,size=1.5)+theme_classic()+ theme(axis.text.x = element_blank())+labs(x = '', y = '-log10 adjusted P')+ geom_hline(aes(yintercept= 1.30))
P1<-ggplot(data = IntravsInter_results,mapping = aes(x=X,y=r2_of_model)) + geom_point(shape=7,size=1.5)+theme_classic()+ theme(axis.text.x = element_blank())+labs(x = '', y = 'R2')+ylim(0,1)

# Correlation matrix
rm(list=ls())
load("subset_T1T2.Rdata")
data[,15:285] <- log(data[15:285],2)
data$Sampling<-factor(data$Sampling,levels = c("1","2"),labels = c("T1","T2"))
data_sub<-data[,c(3,6,15:158)]
data_sub<-na.omit(data_sub)
data_sub1<-subset(data_sub,Sampling=="T1")
data_sub2<-subset(data_sub,Sampling=="T2")
ID<-intersect(data_sub1$ID,data_sub2$ID)
data_sub1<-data_sub1[data_sub1$ID%in%ID,]
data_sub2<-data_sub2[data_sub2$ID%in%ID,]
data_sel<-rbind(data_sub1,data_sub2)
cell.cor <- cor(t(data_sel[,c(3:146)]), use = "pairwise.complete.obs", method = "spearman") ##for all the immunophenotypes
rownames(cell.cor) <- paste0(data_sel$ID, "V", data_sel$Sampling)
colnames(cell.cor) <- paste0(data_sel$ID, "V", data_sel$Sampling)
# Distance matrix:
d <- dist(cell.cor)
# MDS scaling:
ff1 <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# ff$points stores the x y coordinates, so we pull them out for ease:
ff <- ff1$points
ff <- as.data.frame(ff)
# Measure distances: make a 'pair' factor:
ff$Sampling<-data_sel$Sampling
ff$ID1<-data_sel$ID
library(ggplot2)
p<-ggplot(data=ff,aes(x=V1,y=V2,group=Sampling)) + geom_line(aes(group=ID1),size=0.1) + geom_point(aes(shape=Sampling,colour=Sampling),size=2)+theme_bw() + theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  scale_color_manual(values = pal <- c("#B916A5","#0096DB","#9F59D6","#FF2C11","#0158A3",
                                               "#EDE100","#01B7FF","#FF9DA0","#581EA8","#009292",
                                               "#3C5CD4","#FF8907","#44C84E","#FE2A77","#C73E3A","#D9AB42"))+   ## 
  theme(axis.line = element_line(colour = "black"))+labs(x="PRINCIPAL COORDINATE 1",y="PRINCIPAL COORDINATE 2",title ="")
p
library(ggExtra)
ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
wilcox.test(ff[ff$Sampling %in% c("T1"),"V1"],ff[ff$Sampling %in% c("T2"),"V1"],exact=F, paired = T)
wilcox.test(ff[ff$Sampling %in% c("T1"),"V2"],ff[ff$Sampling %in% c("T2"),"V2"],exact=F, paired = T)