#!/usr/bin/env Rscript
# File: MOFA.R
# Author: Shifang Li
# E-mail: fruceslee@gmail.com

library("readxl") 
library(GGally) 
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(MOFA2)
mydata <- read.csv("Resistomics/Sampling/subset_T1T2.csv",header=T)
Igs_serum = data.matrix(mydata[,8:10])
Cytokines = data.matrix(mydata[,c(3:7,32:146)])
FACS =data.matrix(mydata[,c(11:31)])
Group<-mydata[,1]
# ================================================== NORMALIZATION ==================================================
Igs_serum <- log2(Igs_serum + 1)
FACS <- log2(FACS + 1)
Cytokines <- log2(Cytokines + 1)
# PRT is already in logarithmic scale
# ================================================== FILTERING ==================================================
if(F){
  
  limIQR <- 0.3
  
  # Filter for IQR (Q3-Q1)
  iqr <- apply(Igs_serum, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) Igs_serum <- Igs_serum[, -ir]
  
  iqr <- apply(Cytokines, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) Cytokines <- Cytokines[, -ir]
  
  iqr <- apply(FACS, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) FACS <- FACS[, -ir]
}
ldata <- list(Ig= t(Igs_serum),FACS = t(FACS), Cytokines = t(Cytokines))

# ================================================== CREATE MOFA OBJECT ==================================================
MOFAobject <- create_mofa(ldata)
# Visualise data structure
plot_data_overview(MOFAobject)
# ================================================== DEFINE OPTIONS ==================================================
# Data options
# - scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance (default is FALSE)
# - scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance (default is FALSE)
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- T

# Model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
# Training options
# - maxiter: number of iterations. Default is 1000.
# - convergence_mode: "fast", "medium", "slow". For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - drop_factor_threshold: minimum variance explained criteria to drop factors while training. Default is -1 (no dropping of factors)
# - gpu_mode: use GPU mode? This needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 4000
library(MOFAdata)
library(MOFA)
library(ggplot2)
library(qusage)
TrainOptions <- get_default_model_options(MOFAobject)
TrainOptions$DropFactorThreshold <- 0.02 # Automatically drop factors that explain less than 2% of variance in all omics
# ================================================== PREPARE MOFA OBJECT ==================================================

MOFAobject <- prepare_mofa(MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts) # stochastic_options = stochastic_opts

# ================================================== TRAIN THE MODEL ==================================================
# R requires the package reticulate to communicate with Python, and this is the source of most of the installation problems. 
# https://rstudio.github.io/reticulate/

library(reticulate)
# ho creato il mio ENV chiamato myR
# source /anaconda3/etc/profile.d/conda.sh
# conda activate myR
# pip install mofapy2

use_python("C:\\Users\\u232828\\AppData\\Local\\Programs\\Python\\Python39\\python.exe", required = T)
conda_list()
py_config()

outfile <- paste0(getwd(),"/MOFA/model.hdf5")
MOFAmodel <- run_mofa(MOFAobject, outfile,use_basilisk = TRUE)

# ================================================== LOAD TRAINED MODEL AND ADD METADATA ==================================================
model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))
sample_metadata <- model@samples_metadata
sample_metadata$condition<-c(rep("T1",221), rep("T2",221))
samples_metadata(model) <- sample_metadata
# ================================================== VARIANCE DECOMPOSITION ==================================================
plot_variance_explained(model, x="view", y="factor") +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 15), axis.text.y = element_text(size = 15), legend.title = element_text(size = 13))
ggsave(("MOFA/Variance Decomposition.tiff"), width = 4, height = 7)


plot_variance_explained(model, x="view", y="factor", plot_total = T)[[2]] +
  theme(axis.text.x = element_text(angle=50, vjust=1, hjust= 1, size = 19), axis.text.y = element_text(size = 17), axis.title.y = element_text(size = 15))
#calculateVarianceExplained(model, views = "all", factors = "all")

r2 <- calculate_variance_explained(model)
r2$r2_total
# ================================================== VISUALIZATION OF SINGLE FACTORS ==================================================
# Get factor values
Z <- get_factors(model, factors = 1:7 , groups = "all", as.data.frame=TRUE)
#Z$factor <- as.factor(Z$factor)
#df <- merge(Z,  cbind(sample = rownames(Igs_serum), condition = Group), by="sample")

#df$condition <- factor(df$condition, levels = c("PRNT+", "PRNT-"))
Z$group <-c(rep("T1",221), rep("T2",221))
Z$condition <-c(rep("T1",221), rep("T2",221))
df<-Z
# Beeswarm Plots v1
ggplot(Z, aes(x=factor, y=value, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  labs(x="", y="Factor value") + geom_quasirandom(size = 3, dodge.width = 0, alpha = 0.8, method = "quasirandom", width = 0.4) + 
  scale_color_manual(values=mycols) + theme_classic() + geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  scale_shape_manual(values=c(20, 18)) +
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggsave(("MOFA/Beeswarm plots.tiff"), width = 14, height = 6)

if (F){
  # Violin Plots
  plot_factor(model, factor = 1:7, color_by = "condition", dot_size = 1.5, dodge = T, legend = T, add_violin = T, violin_alpha = 0.3) + 
    scale_fill_manual(values=mycols) + theme(text = element_text(size=18)) 
  ggsave(("MOFA/Violin plots12.tiff"), width = 14, height = 6)
}


# Violin Plots
ggplot(df[df$factor %in% c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7"),], aes(x=factor, y=value, fill=condition, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") +
  geom_quasirandom(size = 1.8, dodge.width = 3, alpha = 1, method = "quasirandom", width = 0.4, show.legend = F) + 
  geom_violin(alpha= 0.2, trim = T, scale = "width", position = position_dodge(width=3), show.legend = T, col = "black", width = 2.8) +
  scale_shape_manual(values=c(20, 18)) + scale_color_manual(values=mycols) + scale_fill_manual(values=mycols) + 
  # geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  labs(x="", y="Factor value") + theme_classic() + 
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
ggplot(df[df$factor %in% c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7"),], aes(x=factor, y=value, fill=condition, col=condition)) + facet_wrap(~factor, nrow=1, scales="free_x") + 
  geom_violin(alpha= 0.2, trim = T, scale = "width", position = position_dodge(width=3), show.legend = T, col = "black", width = 2.8) +
  scale_shape_manual(values=c(20, 18)) + scale_color_manual(values=mycols) + scale_fill_manual(values=mycols) + 
  # geom_hline(yintercept=0, linetype="dashed", size=0.4, alpha=0.5) +
  labs(x="", y="Factor value") + theme_classic() + 
  theme(text = element_text(size=18), panel.border = element_rect(color="black", size=0.2, fill=NA), 
        strip.background = element_rect(colour = "black", size=0.3), panel.spacing = unit(0,"lines"),
        axis.line = element_line(colour = 'grey', size = 0.2), axis.text.x = element_text(color="white"),
        axis.ticks.x = element_line(color = "white"))
wilcox.test(df2[df2$condition %in% c("T1"),"Factor7"],df2[df2$condition %in% c("T2"),"Factor7"],exact=F, paired = F)

# ================================================== VISUALIZATION OF COMBINATION OF FACTORS ==================================================

# plot_factors(model, factor = 1:5, color_by = "condition", shape_by = "sex") + 
#  scale_fill_manual(values=mycols) + scale_color_manual(values=mycols) + theme(text = element_text(size=18))
# spread over factors
df2 <- spread(df, key="factor", value="value")
library(ggridges)
diamonds<-df2[,c(3,8)]
diamonds<-df2[,c(3,9)]
diamonds<-df2[,c(3,10)]
library(reshape2)
nba.m <- melt(diamonds,id = c("condition"))
p4<-ggplot(nba.m, aes(x=value, y=variable, fill=condition))+
  geom_density_ridges(scale=4)+
  scale_fill_cyclical(values = c("red","blue"))+
  theme_ridges(grid = FALSE)+labs(x = '', y="",title = '',main="")
p5<-ggplot(nba.m, aes(x=value, y=variable, fill=condition))+
  geom_density_ridges(scale=4)+
  scale_fill_cyclical(values = c("red","blue"))+
  theme_ridges(grid = FALSE)+labs(x = '', y="",title = '',main="")
p6<-ggplot(nba.m, aes(x=value, y=variable, fill=condition))+
  geom_density_ridges(scale=4)+
  scale_fill_cyclical(values = c("red","blue"))+
  theme_ridges(grid = FALSE)+labs(x = '', y="",title = '',main="")
library(gridExtra)          
grid.arrange(p4,p5,p6,ncol = 1)
ggpairs(df2, columns = c(4:10),
        lower = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)), 
        diag = list(continuous=GGally::wrap("densityDiag", alpha = 0.5, size = 0.2)), 
        upper = list(continuous=GGally::wrap("points", size=2, alpha = 0.6)),
        mapping = ggplot2::aes(color=condition)) + 
  scale_fill_manual(values=mycols) + scale_color_manual(values=mycols) + 
  theme(text = element_text(size=16), panel.grid.major = element_blank(),
        panel.border = element_rect(color="black", size=0.3, fill=NA),
        panel.background = element_rect(fill = "white"), strip.background = element_rect(color = "black", size=0.3, fill = "gray94"))
ggsave(("MOFA/Combinations of factors.tiff"), width = 7, height = 7)
# ================================================== VISUALIZATION OF FEATURE WEIGHTS ==================================================
#plot_weights(model, view = "FACS", factor = 2, nfeatures = 7,  scale = T, text_size = 4)
var <- model@cache$variance_explained$r2_per_factor[[1]]
var <- cbind(Factor = 1:nrow(var), var)

for(z in 1:length(ldata)){
  
  plot_top_weights(object = model, view = names(ldata)[z], factor = var[which.max(var[,z + 1]), 1], nfeatures = 10) + 
    theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
  ggsave(paste0("MOFA/Feature weights F1", var[which.max(var[,z + 1]), 1]," (", names(ldata)[z], ").tiff"), width = 10, height = 4)
}

z <- 1
plot_top_weights(object = model, view = names(ldata)[z], factor = 1, nfeatures = 10) + 
  theme(text = element_text(size=15), strip.background = element_rect(color = "black", fill = "gray94"))
ggsave(paste0("MOFA/Feature weights F", 3," (", names(ldata)[z], ").tiff"), width = 10, height = 4)
#weights <- get_weights(model)
weights <- get_weights(model, views = 3, factors = 2)
#weights$FACS
C1<-data.frame(weights$Cytokines)
C1$names<-row.names(C1)
C1$Factor61<-abs(C1$Factor2)
C1 <-C1[order(C1[,3], decreasing = TRUE), ]
C1<-C1[c(1:10),]
C2<-data.frame(weights$FACS)
C2$names<-row.names(C2)
C2$Factor61<-abs(C2$Factor2)
C2 <-C2[order(C2[,3], decreasing = TRUE), ]
C2<-C2[c(1:10),]
C3<-weights$Igs_serum
p1<-ggplot(C1,aes(x = reorder(names,Factor61),Factor61)) + 
  geom_point(size=4,color="blue") + 
  geom_segment(aes(x=names,xend=names,y=0,yend=Factor61),
               size=1, color="blue") +
  theme_minimal()+theme_bw() +theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(colour = "black")) +
  coord_flip()

p2<-ggplot(C2,aes(x = reorder(names,Factor61),Factor61)) + 
  geom_point(size=4,color="red") + 
  geom_segment(aes(x=names,xend=names,y=0,yend=Factor61),
               size=1, color="red") +
  theme_minimal()+theme_bw() +theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(colour = "black")) +
  coord_flip()
library(gridExtra)          
grid.arrange(p1,p2,ncol = 2)