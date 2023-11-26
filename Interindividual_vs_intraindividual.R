#Interdividual vs intraindividual
rm(list=ls())
set.seed(1234)
#mydata <- read.csv("D:\\meeting\\datanew\\all_T1T2.csv",header=T)
#mydata2<- mydata[mydata$Sampling %in% c("2"),]
#namesID<-data.frame(mydata2[,3])
#write.table(namesID, 'namesID.txt', quote = FALSE,row.names = FALSE)
#tbl2<-read.csv("D:\\meeting\\datanew\\all_T1T2.csv",header=T)
#tbl2<- tbl2[tbl2$Sampling %in% c("1"),]
#tbl3<-read.table("namesID.txt",header=T)
#data<-tbl2[match(tbl3$ID, tbl2$ID),]
mydata <- read.csv("D:\\meeting\\datanew\\all_T1T2.csv",header=T)
mydata1<-mydata[,]
mydata1$ID<-factor(mydata1$ID)
mydata1$Sampling<-factor(mydata1$Sampling)
inormal <- function(x)#
{
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

mydata1[,15:285] <- log(mydata1[15:285],2)
for (i in 27:55){
  mydata1[,i]<-inormal(mydata1[,i])
}
mydata1<- mydata1[mydata1$Sampling %in% c("1","2"),]
mydata2<-data.frame(mydata1)
# Create a new matrix with 28 rows and 4 columns:
mydata3 <- matrix(nrow=144, ncol =4)
# Rename column and row names of the new matrix:
colnames(mydata3) <- c("p_of_model","r2_of_model",
                       "prop_attributed_to_sampling", "prop_attributed_to_ID")
#dim(mydata2[,c(14:15,20:22,24:42,46:50,56:66,68,70,72,73:83,85,87,89,106,107:108,110:117,119,121,123,124:134,136,138,140,141:142,144,146:151,153,155,157,158,159,161:168,170,172,179,180,185,191,192,193,195:202,206,208,209:219,221,223,225,227,232,234,236,242)])
test<-mydata2[,c(14:15,20:22,24:42,46:50,56:66,68,70,72,73:83,85,87,89,106,107:108,110:117,119,121,123,124:134,136,138,140,141:142,144,146:151,153,155,157,158,159,161:168,170,172,179,180,185,191,192,193,195:202,206,208,209:219,221,223,225,227,232,234,236,242)]

rownames(mydata3) <- colnames(mydata2[,c(14:15,20:22,24:42,46:50,56:66,68,70,72,73:83,85,87,89,106,107:108,110:117,119,121,123,124:134,136,138,140,141:142,144,146:151,153,155,157,158,159,161:168,170,172,179,180,185,191,192,193,195:202,206,208,209:219,221,223,225,227,232,234,236,242)])
mydata4<-mydata2[,c(3,6,14:15,20:22,24:42,46:50,56:66,68,70,72,73:83,85,87,89,106,107:108,110:117,119,121,123,124:134,136,138,140,141:142,144,146:151,153,155,157,158,159,161:168,170,172,179,180,185,191,192,193,195:202,206,208,209:219,221,223,225,227,232,234,236,242)]
# Estimate the relative r2 contributions of inter-individual variation and
# intra-individual variation over time(day 0 and day 7)
# for each trait:
library(relaimpo)
library(RColorBrewer)
library(ellipse)
library(gplots)
library(vegan)
for(i in c(3:146)){
  model_vac02 <- lm(mydata4[,i] ~ mydata4[,"Sampling"]
                    + mydata4[,"ID"],na.action="na.exclude")
  mydata3[c(i - 3 + 1), 2] <- summary(model_vac02)$r.squared
  mydata3[c(i - 3 + 1), 1] <- pf(summary(model_vac02)$fstatistic[1],
                                 summary(model_vac02)$fstatistic[2],
                                 summary(model_vac02)$fstatistic[3], lower.tail = F)
  {
    ca01 <- calc.relimp(model_vac02, diff = T, rela = T)
    mydata3[c(i - 3 + 1), 3] <- ca01$lmg[grepl(pattern = "Sampling",
                                               x = names(ca01$lmg))]
    mydata3[c(i - 3 + 1), 4] <- ca01$lmg[grepl(pattern = "ID",
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
mydata3<-data.frame(mydata3)
mydata3$p_of_model_adj<- p.adjust(mydata3$p_of_model, method = "bonferroni",
                                  n = length(mydata3$p_of_model))
write.csv(mydata3, 'nature immunology.csv', quote = FALSE)
mydata3 <- read.csv("nature immunology.csv",header=T)
summary(mydata3$prop_attributed_to_ID)
summary(mydata3$prop_attributed_to_sampling)
summary(mydata3$r2_of_model)
mydata3$X<-row.names(mydata3)
names(mydata3)[names(mydata3)=="prop_attributed_to_ID"]="Interindividual"
names(mydata3)[names(mydata3)=="prop_attributed_to_sampling"]="Intraindividual"
nba.m <- melt(mydata3,id = c('X',"p_of_model","r2_of_model","p_of_model_adj"))
ggplot(nba.m, aes(x=X, y=value, fill=factor(variable, levels=c("Interindividual", "Intraindividual")))) +coord_flip()+
  geom_bar(stat="identity")+scale_fill_manual(values=c('lightgray','black'))+theme_classic()+theme(axis.text.x = element_text(angle = 90, size =0.6,hjust = 0.5, vjust = 0.5))+labs(x = '', y = 'Proportion of R2',fill="")+theme(legend.position="top")

mydata3$p_of_model_adj<-log10(mydata3$p_of_model_adj)
mydata3$p_of_model_adj<-mydata3$p_of_model_adj*(-1)
P1<-ggplot(data = mydata3,mapping = aes(x=X,y=r2_of_model)) + geom_point(shape=7,size=3)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+labs(x = '', y = 'R2')+ylim(0,1)
P2<-ggplot(data = mydata3,mapping = aes(x=X,y=p_of_model_adj)) + geom_point(shape=7,size=3)+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))+labs(x = '', y = '-log10 adjusted P')+ geom_hline(aes(yintercept= 1.30))
library(gridExtra)          
grid.arrange(P1,P2)
P2<-ggplot(data = mydata3,mapping = aes(x=X,y=p_of_model_adj)) + geom_point(shape=7,size=1.5)+theme_classic()+ theme(axis.text.x = element_blank())+labs(x = '', y = '-log10 adjusted P')+ geom_hline(aes(yintercept= 1.30))
P1<-ggplot(data = mydata3,mapping = aes(x=X,y=r2_of_model)) + geom_point(shape=7,size=1.5)+theme_classic()+ theme(axis.text.x = element_blank())+labs(x = '', y = 'R2')+ylim(0,1)


#Code for panel F.  PCA
# Correlation matrix:
mydata <- read.csv("D:\\meeting\\datanew\\all_T1T2.csv",header=T)
mydata[,15:285] <- log(mydata[15:285],2)
mydata$Sampling<-factor(mydata$Sampling,levels = c("1","2"),labels = c("T1","T2"))

mydata4<-mydata[,c(3,6,14:15,20:22,24:42,46:50,56:66,68,72,73:83,85,89,106,107:108,110:117,119,123,124:134,136,140,141:142,144,146:151,153,157,158,159,161:168,170,179,180,185,192,193,195:202,208,209:219,221,225,227,232,234,236,242)]
mydata4<-na.omit(mydata4)
mydata41<-subset(mydata4,Sampling=="T1")
mydata42<-subset(mydata4,Sampling=="T2")
ID<-intersect(mydata41$ID,mydata42$ID)
mydata41<-mydata41[mydata41$ID%in%ID,]
mydata42<-mydata42[mydata42$ID%in%ID,]
mydata5<-rbind(mydata41,mydata42)
mydata51<-mydata5[,c(11:31)]  ##32:137 3:7 8:10 11:31 

cell.cor <- cor(t(mydata5[,c(3:137)]), use = "pairwise.complete.obs", method = "spearman")
rownames(cell.cor) <- paste0(mydata5$ID, "V", mydata5$Sampling)
colnames(cell.cor) <- paste0(mydata5$ID, "V", mydata5$Sampling)
# Distance matrix:
d <- dist(cell.cor)
# MDS scaling:
ff1 <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
# ff$points stores the x y coordinates, so we pull them out for ease:
ff <- ff1$points
ff <- as.data.frame(ff)
# Measure distances: make a 'pair' factor:
ff$Sampling<-mydata5$Sampling
ff$ID1<-mydata5$ID
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

pc1_boxplot <- ggplot(data = ff, aes(x = Sampling, y = V1)) +
  geom_boxplot(aes(color = Sampling)) +
  scale_color_manual(values = c('red', 'blue'), limits = c('T1', 'T2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position= 'none', axis.title = element_text(color = 'transparent')) +labs(x="",y="Immunological distance",title ="p=0.836")+
  coord_flip()

pc2_boxplot <- ggplot(data = ff, aes(x = Sampling, y = V2)) +
  geom_boxplot(aes(color = Sampling)) +labs(x="",y="Immunological distance",title ="p=0.01706")+
  scale_color_manual(values = c('red', 'blue'), limits = c('T1', 'T2')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.position= 'none', axis.title = element_text(color = 'transparent'))

library(grid)

grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 3)))
print(p, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1:2))
print(pc1_boxplot, vp = viewport(layout.pos.row = 3, layout.pos.col = 1:2))
print(pc2_boxplot, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 3))
