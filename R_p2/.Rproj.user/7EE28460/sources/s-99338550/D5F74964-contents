library(ape)
library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)
library("PerformanceAnalytics")
library(stringr)
library(gplots)


os = Sys.info()['sysname']

source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'

emdata = read.nexus(emdatadir)
# tree plot
baleenwhale = phylo2L(emdata,error = 1e-5)
L_ext = baleenwhale$L
extantspecieslabel = baleenwhale$ESL

extantspecieslabel <- c("B.mysticetus", "E.australis"       
, "E.glacialis", "E.japonica"        
, "B.acutorostrata" ,"B.bonaerensis"  
, "B.borealis" ,"B.brydei"       
, "B.edeni" , "B.omurai"       
, "B.musculus" , "B.physalus"     
, "M.novaeangliae" , "E.robustus"     
, "C.marginata" )
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
phylo_test$tip.label <- extantspecieslabel
plot(phylo_test,show.tip.label = TRUE)


# load data
my.data <- phylo_test$tip.label
# extract numbers only
my.data.num <- as.numeric(str_extract(my.data, "[0-9]+"))
sorted.order <- order(my.data.num)
sorted.species.names <- extantspecieslabel[sorted.order]

# 1000 simulation results
filepredict_name =  paste0(dir,'predictsimTPUS.csv')
predictZ = read.csv(filepredict_name)
filepredict_DR =  paste0(dir,'predictsimDR_fixm.csv')
predictZDR = read.csv(filepredict_DR)

corr <- c()
for(tip_label in extantspecieslabel){
  corr <- cbind(corr, predictZ[,which(sorted.species.names == tip_label)])
} 

colnames(corr) <- extantspecieslabel
corr.DR <- cor(corr)

hc <- as.hclust(phylo_test) #Compulsory step as as.dendrogram doesn't have a method for phylo objects.
dend <- as.dendrogram(hc)
plot(dend, horiz=TRUE)


heatmap.2(corr.DR, Rowv=dend, Colv=dend,col=bluered)


chart.Correlation(corr.DR)
