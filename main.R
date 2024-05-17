#!/usr/bin/env Rscript
# File: main.R
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

library(psych)
library(ggplot2)
library(reshape2)
data<-read.csv("Resistomics\\data_corrected_age_season.csv",header=T)
Igs<-data[,c(1:5)]
Cyto_serum<-data[,c(6:8)]
FACS<-data[,c(9:29)]
Stim_cyto<-data[,c(30:144)]
##correlation among stimulated cytokines, similar used to estimate other correlations.
phylum_corr <-corr.test(Stim_cyto, use = "pairwise",
                        method = "spearman",adjust="none")
Corr_cyto_R<-phylum_corr$r
Corr_cyto_P<-phylum_corr$p
for (i in c(1:115)){Corr_cyto_P[,i]<-p.adjust(Corr_cyto_P[,i],method = "fdr",n=length(Corr_cyto_P[,i]))}
#write.csv(Corr_cyto_P, 'p_corr_stim_cells.csv', quote = FALSE)
#write.csv(Corr_cyto_R, 'r_corr_stim_cells.csv', quote = FALSE)
Corr_cyto_R<-data.frame(Corr_cyto_R)
Corr_cyto_P<-data.frame(Corr_cyto_P)
Corr_cyto_R$X<-rownames(Corr_cyto_R)
Corr_cyto_P$X<-rownames(Corr_cyto_P)
Corr_cyto_R <- melt(Corr_cyto_R,id = c('X'))
Corr_cyto_P <- melt(Corr_cyto_P,id = c('X'))
Corr_cyto_R$value1<-Corr_cyto_P$value
names(Corr_cyto_R)[names(Corr_cyto_R)=="value"]="R"
names(Corr_cyto_R)[names(Corr_cyto_R)=="value1"]="P"

anno<-read.csv("Resistomics\\anno.csv",header=T)
Corr_cyto_R<-phylum_corr$r
anno<-data.frame(anno[30:144,])
anno<-anno[,-c(1:2)]
rownames(anno) = colnames(Corr_cyto_R)
library(pheatmap)
newCols <- colorRampPalette(grDevices::rainbow(length(unique(anno$Type))))
annoCol <- newCols(length(unique(anno$Type)))

annoCol <-c("#B916A5","#0096DB","#9F59D6","#FF2C11","#0158A3")
names(annoCol) <- unique(anno$Type)

annoCol1 <-c("blue","green")
names(annoCol1) <- unique(anno$Type1)

annoCol2 <-c("#FFFC9B","#FF8907","#924900","#01B7FF","#3C5CD4","#44C84E")
names(annoCol2) <- unique(anno$Type2)

annoCol3 <-c("#FE2A77","#009393")
names(annoCol3) <- unique(anno$Type3)
annoCol4 <- list(Type = annoCol,Type1 = annoCol1,Type2 = annoCol2,Type3 = annoCol3)
range <- max(abs(Corr_cyto_R))
hm_1<-pheatmap(Corr_cyto_R,cluster_row = T,annotation_colors = annoCol4,annotation_row = anno,symm=T,breaks = seq(-range, range, length.out = 100),fontsize_row=4, fontsize_col=4,angle_col =90,fontsize_number = 6, number_color = "white", cellwidth =4, cellheight =4,cluster_col = T,color = colorRampPalette(c("blue", "white","red"))(100))

##NMDS
cell.cor <- cor(data[,c(30:144)], use = "pairwise.complete.obs", method = "spearman")
d <- dist(cell.cor)
ff1 <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
ff2 <- as.data.frame(ff1$points)
library(ggplot2)
ff2$Type<-anno[,2]
ff2$Type0<-anno[,4]
ff2$Type1<-anno[,3]
p<-ggplot(data=ff2,aes(x=V1,y=V2,group=Type))+ geom_point(aes(colour=Type1,shape=Type),size=2.5)+theme_bw() + theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +scale_color_manual(values = pal <- c("#01B7FF","#3C5CD4","#44C84E","#FE2A77","#924900","#FF8907","red","red","blue"))+   ## 
  theme(axis.line = element_line(colour = "black"))+labs(x="NMDS1",y="NMDS2",title ="",size=13)+ stat_ellipse(aes(color = Type1), level = 0.95, linetype = 2, show.legend = FALSE)+ scale_shape_manual(values = c(17,19))  #scale_shape_manual(values=c(18, 19))
p
library(ggExtra)
ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
library(PMCMR)
library(PMCMRplus)
All3$Type <- as.factor(All3$Type)
ff2$Type1 <- as.factor(ff2$Type1)
PMCMRplus::kwAllPairsDunnTest(ff2$V2,ff2$Type1,p.adjust.method ="fdr")
result1<-c()
for (i in c(1:2)){
  list<-data.frame(colnames(ff2))
  fit1<- wilcox.test(ff2[ff2$Type %in% c("WWBC"),c(list[i,])],ff2[ff2$Type %in% c("PBMC"),c(list[i,])],exact=F, paired = F)
  result1<-rbind(result1,fit1$p.value) #fit1$p.value
}
#write.csv(result1, 'NMDS_PBMC_WWBC.csv', quote = FALSE)

### variance explained using lm
### functions, eg., select features associated with PBMC_IL1B_LPS
immpheno.assoc<-function(dat){
  res<-c()
  for(i in colnames(dat)){
    cor.tmp<-cor.test(input.pheno$PBMC_IL1B_LPS,method = "spearman")
    tmp.res<-data.frame(mt=i,cor=cor.tmp$estimate,p=cor.tmp$p.value)
    res<-rbind(res,tmp.res) 
  }
  return(res)
}

### select features that are highly correlated with each other to remove colinearity
high.assoc<-function(dat){
  id<-colnames(dat)
  res.cor<-c()
  for(i in 1:(length(id)-1)){
    for(j in (i+1):length(id)){
      cor.tmp<-cor.test(dat[,i],dat[,j],method = "spearman")
      tmp.res<-data.frame(f1=colnames(dat)[i],f2=colnames(dat)[j],cor=cor.tmp$estimate,p=cor.tmp$p.value)
      res.cor<-rbind(res.cor,tmp.res)
    }
  }
  return(res.cor)
}

### remove highly correlated features and get independent features
indep_cov<-function(cor.pheno,cor.feature){
  cor.feature$p1<-cor.pheno$p[match(cor.feature$f1,cor.pheno$mt)]
  cor.feature$p2<-cor.pheno$p[match(cor.feature$f2,cor.pheno$mt)]
  
  check<-cor.feature[which(abs(cor.feature$cor)>0.4),]
  if(dim(check)[1]==0){
    mt.remove<-c()
  }else{
    mt.remove<-c()
    for(i in 1:dim(check)[1]){
      ifelse(check$p1[i]>check$p1[i],id.remove<-check$f1[i],id.remove<-check$f2[i])
      mt.remove<-c(mt.remove,id.remove)
    }
  }
  return(mt.remove)
}

### Calculate the variance explained
variance_lm <- function(y, x, return.type="r.squared") {
  if (!is.null(x)) {
    if (class(x) == "numeric") {
      model <- lm(y ~ x)
    } else {
      model <- lm(y ~ ., data=x)
    }
    if (return.type=="r.squared") {
      return (summary(model)$adj.r.squared)
    } else if (return.type == "model") {
      return (model)
    }
  }
}

###distance matrix‒based variance estimation
tmp<-na.omit(data_sel)
inBC <- vegdist(scale(tmp[,-c(1:21)]),method = "euclidean")
adonisResults <- NULL
for (i in c(1:21)) {
  ad <- adonis(inBC ~ tmp[,i],permutations=1000)
  aov_table <- ad$aov.tab
  oneRow <- data.frame(Var=i,
                       DF=aov_table[1,1],
                       SumsOfSqs=aov_table[1,2],
                       MeanSqs=aov_table[1,3],
                       FModel=aov_table[1,4],
                       R2=aov_table[1,5],
                       pval=aov_table[1,6],
                       FDR.BH=NA,
                       Significant=NA)
  adonisResults <- rbind.data.frame(adonisResults,oneRow)
}
rownames(adonisResults) = colnames(tmp[,c(1:21)])
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant="No"
adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"
adonisResults <- subset(adonisResults,Significant=="Yes")