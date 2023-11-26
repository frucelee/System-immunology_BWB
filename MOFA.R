##MOFA analysis
##MOFA analysis
library("readxl") # read excel
library(GGally) # ggpairs() 
library(tidyr)
library(ggplot2)
#install.packages('devtools')
#devtools::install_github('VPetukhov/ggrastr')
library(ggbeeswarm) # geom_quasirandom()

# pip install mofapy2 (on terminal)
# devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
# https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md
# https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html
# (Optional) set up reticulate connection with Python
# reticulate::use_python("/Users/ricard/anaconda3/envs/base_new/bin/python", required = T)
## Warning: replacing previous import 'DelayedArray::pnorm' by 'stats::pnorm'
## when loading 'MOFA2'
library(MOFA2)
mydata <- read.csv("D:\\datanew\\all_T1T2.csv",header=T)
mydata <- mydata[,c(1:15,20:22,24:42,46:50,56:66,68,72,73:83,85,89,106,107:108,110:117,119,123,124:134,136,140,141:142,144,146:151,153,157,158,159,161:168,170,179,180,185,192,193,195:202,208,209:219,221,225,227,232,234,236,242)]
mycols <- c('T1' = 'red', 'T2' = 'blue')
#rm(CD40Lpos, BposUS, ALL, FACS, PRTcov)
which( colnames(mydata)=="MHCII.in.lin." )
Igs = data.matrix(mydata[,19:21])
Cytokines1 =data.matrix(mydata[,c(14:18)])
Cytokines = data.matrix(mydata[,c(43:148)])
FACS =data.matrix(mydata[,c(22:42)])

Group<-mydata[,6]

# ================================================== NORMALIZATION ==================================================

Igs <- log2(Igs + 1)
FACS <- log2(FACS + 1)
Cytokines <- log2(Cytokines + 1)
Cytokines1 <- log2(Cytokines1 + 1)
# PRT is already in logarithmic scale


# ================================================== FILTERING ==================================================
if(F){
  
  limIQR <- 0.3
  
  # Filter for IQR (Q3-Q1)
  iqr <- apply(Igs, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) Igs <- Igs[, -ir]
  
  iqr <- apply(Cytokines, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) Cytokines <- Cytokines[, -ir]
  
  iqr <- apply(FACS, 2, IQR, na.rm = T)
  hist(iqr, nclass = 100, xT2="IQR", yT2="Frequency", col="blue2")
  ir <- which(iqr < limIQR)
  if(length(ir)!=0) FACS <- FACS[, -ir]
}
# ================================================== DATA LOAD ==================================================

# Multiple formats are allowed for the input data:

## -- Option 1 -- ##
# nested list of matrices, where the first index refers to the view and the second index refers to the group.
# samples are stored in the rows and features are stored in the columns.
# Missing values must be filled with NAs, including samples missing an entire view

# (...)

## -- Option 2 -- ##
# data.frame with columns ["sample","feature","view","group","value"]
# In this case there is no need to have missing values in the data.frame,
# they will be automatically filled in when creating the corresponding matrices


#ldata <- list(CD40L = t(CD40L), NK = t(NK), CD8 = t(CD8), Bcells = t(Bcells), Proteomics = t(PRT))
# ldata <- list('FACS UM' = t(cbind(CD40L, Bcells)), 'FACS CL' = t(cbind(NK, CD8)), Proteomics = t(PRT))
# ldata <- list(FACS = t(cbind(CD40L, Bcells, NK, CD8)), Proteomics = t(PRT))

ldata <- list(FACS = t(FACS), Cytokines = t(Cytokines))

# ================================================== CREATE MOFA OBJECT ==================================================

# The aim of the multi-group framework is to identify the sources of variability *within* the groups. If your aim is to find a factor that 'separates' the groups, 
# you DO NOT want to use the multi-group framework. Please see the FAQ (https://github.com/bioFAM/MOFA2#2-faq-on-the-multi-group-functionality)

# MOFAobject <- create_mofa(ldata, groups=Group)
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
# - likelihoods: likelihood per view (options are "gaussian","poisson","bernoulli"). By default, they are guessed inteIgslly. We advise users to use “gaussian??? whenever possible!
# - num_factors: number of factors. By default K=10
# - spikesT2_factors: use spike-sT2 sparsity prior in the factors? default is FALSE.
# - spikesT2_weights: use spike-sT2 sparsity prior in the weights? default is TRUE.
# - ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups.
# - ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views.
# Only change the default model options if you are familiar with the underlying mathematical model!
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 7
# model_opts$likelihoods <- c("gaussian","gaussian","gaussian","gaussian","gaussian")
# model_opts$likelihoods <- c("gaussian","gaussian","gaussian")

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
#BiocManager::install("MOFAdata")
#BiocManager::install("MOFA")
library(MOFAdata)
library(MOFA)
library(ggplot2)
library(qusage)
TrainOptions <- get_default_model_options(MOFAobject)
TrainOptions$DropFactorThreshold <- 0.02 # Automatically drop factors that explain less than 2% of variance in all omics

###  Prepare MOFA models x10 and select the model with best fit
n_inits <- 10
MOFAlist <- lapply(seq_len(n_inits), function(it) {
  TrainOptions$seed <- 2019 + it
  MOFAobject <- prepare_mofa(
    MOFAobject, data_options = data_opts, model_options = model_opts, training_options = train_opts
  )
  runMOFA(MOFAobject)
})
# (Optional) Set stochastic inference options
# Only recommended with very large sample size (>1e6) and when having access to GPUs
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.1 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)
# stochastic_opts <- get_default_stochastic_options(MOFAobject)


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

# Using a conda enviroment called myR
#use_condaenv("myR", required = TRUE)

outfile <- paste0(getwd(),"/MOFA/model.hdf5")
MOFAmodel <- run_mofa(MOFAobject, outfile,use_basilisk = TRUE)



# ================================================== LOAD TRAINED MODEL AND ADD METADATA ==================================================

model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))

sample_metadata <- model@samples_metadata
# head(sample_metadata, n=3)
#sample_metadata$Group <- Group
# sample_metadata$age <- Age
# sample_metadata$sex <- Sex
# sample_metadata
sample_metadata$condition<-c(rep("T1",221), rep("T2",221))
samples_metadata(model) <- sample_metadata
# head(model@samples_metadata, n=3)


# ================================================== VARIANCE DECOMPOSITION ==================================================

# Total variance explained per view and group
# head(model@cache$variance_explained$r2_total[[1]]) # group 1

# Variance explained for every factor in per view and group
# head(model@cache$variance_explained$r2_per_factor[[1]]) # group 1

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
#df <- merge(Z,  cbind(sample = rownames(Igs), condition = Group), by="sample")

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


# axis.text = element_blank(), , axis.ticks = element_blank()

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
C3<-weights$Igs
#weights <- get_weights(model, views = 2, factors = 2)

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

p1<-ggplot(C1,aes(names,Factor2)) + 
  geom_point(size=4,color="blue") + 
  geom_segment(aes(x=names,xend=names,y=0,yend=Factor2),
               size=1, color="gray") +
  theme_minimal()+theme_bw() +theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(colour = "black")) +
  coord_flip()

p1<-ggplot(C1,aes(names,Factor2)) + 
  geom_point(size=4,color="blue") +
  coord_flip()+ 
  geom_segment(aes(x=names,xend=names,y=0,yend=Factor2),
               size=1, color="gray")+theme_bw() + theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(colour = "black"))+labs(x="",y="Log2 cytokines levels",title ="BMC_IL-10_HKLM p=0.062")
# ================================================== VISUALIZATION OF PATTERNS IN THE DATA ==================================================

# df3 <- merge(df2, cbind(sample = rownames(freq), freq), by = "sample")
# m <- df3[, c(1:3,10,17)]
# names(m)[5] <- "CD3_Freq_of_Parent"
# ggplot(m, aes(Factor5, CD3_Freq_of_Parent, color= condition)) + geom_point(alpha = 0.4, size = 3) + geom_smooth(formula=y~x, method=lm) + scale_color_manual(values=mycols)

library(ggpubr)

plot_data_scatter(model, view = "FACS", factor = 5,
                  features = 6,           # number of features to plot (they are selected by loading)
                  add_lm = TRUE, lm_per_group = F) + scale_fill_manual(values=mycols) 
ggsave(("MOFA/Scatter plots NK_values vs Factor1_values.tiff"), width = 10, height = 6)


plot_data_scatter(model, view = "OTU", factor = 2,
                  features = 6,           # number of features to plot (they are selected by loading)
                  add_lm = TRUE,          # add linear regression
                  color_by = "Group", lm_per_group = F) + scale_fill_manual(values=mycols) 
ggsave(("MOFA/Scatter plots Proteomics_values vs Factor1_values.tiff"), width = 10, height = 6)

library(uwot)

set.seed(42)
#model <- run_umap(model)
#model <- load_model(paste0(getwd(),"/MOFA/model.hdf5"))
#sample_metadata <- model@samples_metadata
#head(sample_metadata, n=3)
#sample_metadata$condition <- Group
#sample_metadata
#samples_metadata(model) <- sample_metadata
#head(model@samples_metadata, n=3)
model <- run_tsne(model, dims = 2, perplexity = 7, verbose=F, max_iter = 5000, theta = 0.0, pca_scale = F, normalize = T)
Z <- model@dim_red[["TSNE"]]
Z$condition <-c(rep("T1",221), rep("T2",221))
ggplot(Z, aes(TSNE1, TSNE2, col=condition)) +  geom_point(size = 6, alpha = 0.8)  + scale_color_manual(values=mycols) + 
  theme_bw(base_size = 17) + theme(legend.title = element_blank(), panel.grid.major = element_blank())  + geom_line(aes(group=sample))
ggsave(("MOFA/TSE.tiff"), width = 9, height = 8)
Z$ID<-mydata$ID
p<-ggplot(data=Z,aes(x=TSNE1,y=TSNE2,group=condition)) + geom_line(aes(group=ID)) + geom_point(aes(shape=condition,colour=condition),size=3)+theme_bw() + theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +   ## 
  theme(axis.line = element_line(colour = "black"))+labs(x="TSNE1",y="TSNE2",title ="")
p
library(ggExtra)
ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
wilcox.test(Z[Z$condition %in% c("T1"),"TSNE2"],Z[Z$condition %in% c("T2"),"TSNE2"],exact=F, paired = F)

wilcox.test(df2[df2$condition %in% c("T1"),"Factor7"],df2[df2$condition %in% c("T2"),"Factor7"],exact=F, paired = F)
