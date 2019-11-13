library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)

source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/Est/'

emdata = read.nexus(emdatadir)

dropextinct = T
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
phylo_test$node.label <- c(1:14)
# plot(phylo_test,show.tip.label = TRUE,show.node.label = TRUE)

# writeNexus(phylo_test,file='C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/bw.nex')
# Extract the order names of simulated species
prunedL = phylo2L(phylo_test)
pl = prunedL$L
phylo_pl = DDD::L2phylo(pl)
# plot(phylo_pl,show.tip.label = TRUE,show.node.label = TRUE)

rank.sim.species <- as.numeric(gsub("t","", phylo_pl$tip.label))
sim.species.col.names = c()
for(i in c(1:15)){
  sim.species.col.names[rank.sim.species[i]] <- extantspecieslabel[i]
}




# empirical data
fileemp_name = paste0(dir,'emp.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]


obsZ_pic = 10^obsZ_emp[,2]
names(obsZ_pic) = extantspecieslabel



# calculate sim contrast
mode.label = c('TVP','TV','TVM')


simfile_TVP = paste0(dir,'predictsimTVP.csv')
simfile_TV = paste0(dir,'predictsimTV.csv')
simfile_TVM = paste0(dir,'predictsimTVM.csv')


# predict simulations
predictZ_TVP = read.csv(simfile_TVP)
predictZ_TV = read.csv(simfile_TV)
predictZ_TVM = read.csv(simfile_TVM)




sort = 0

if(sort == 0){
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  
  dimnames(predictZ_matrix_TVP)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TV)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TVM)[[2]] = sim.species.col.names
  
  
}else{
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  TVP_Z_order = t(apply(predictZ_matrix_TVP,1,order))
  TVM_Z_order = t(apply(predictZ_matrix_TVM,1,order))
  
  predictZ_matrix_TVP=t(apply(predictZ_matrix_TVP,1,sort))
  predictZ_matrix_TV=t(apply(predictZ_matrix_TV,1,sort))
  predictZ_matrix_TVM=t(apply(predictZ_matrix_TVM,1,sort))
  
  dimnames(predictZ_matrix_TVP)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TV)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TVM)[[2]] = sorted.species.labels
  
}


samplesize = nrow(predictZ_matrix_TVP)


TVP.pic <- c()
TV.pic <- c()
TVM.pic <- c()

scaled = FALSE
# empirical contrast
empirical_pic<-pic(obsZ_pic,phylo_test, scaled = scaled)

for(size in c(1:samplesize)){
  TVP.pic <- rbind(TVP.pic,pic(predictZ_matrix_TVP[size,],phylo_test, scaled = scaled))
  TV.pic <- rbind(TV.pic,pic(predictZ_matrix_TV[size,],phylo_test, scaled = scaled))
  TVM.pic <- rbind(TVM.pic,pic(predictZ_matrix_TVM[size,],phylo_test, scaled = scaled))
  
}
restructured_col <- c(2,3,4,1,5,7,6,13,14,8,9,10,11,12)
restructured_col_plottree <- c(2,3,4,1,7,6,11,12,10,9,8,14,13,5)
TVP.pic <- TVP.pic[,restructured_col]
TV.pic <- TV.pic[,restructured_col]
TVM.pic <- TVM.pic[,restructured_col]



d_all_TVP = as.data.frame(as.table(TVP.pic))[,2:3]
colnames(d_all_TVP) = c('betspecies','contrast')
d_all_TV = as.data.frame(as.table(TV.pic))[,2:3]
colnames(d_all_TV) = c('betspecies','contrast')
d_all_TVM = as.data.frame(as.table(TVM.pic))[,2:3]
colnames(d_all_TVM) = c('betspecies','contrast')

df_emp = as.data.frame((as.table(empirical_pic)))
colnames(df_emp) = c('betspecies','contrast')
ordered.emp <- empirical_pic[restructured_col_plottree]

xe<-setNames(d_all_TVP$contrast,d_all_TVP$betspecies)
emp_dots<-setNames(df_emp$contrast,df_emp$betspecies)

par(mfrow=c(1,2))
plotTree(phylo_test,mar=c(5.1,1.1,2.1,0))
# nodelabels(text=phylo_test$node.label,node=phylo_test$node.label+Ntip(phylo_test),
#            frame="none",adj=c(1.1,-0.4))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(xe~factor(names(xe),levels=restructured_col_plottree),horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(phylo_test)),at = c(2:15)-0.5,varwidth=TRUE,
        ylim=c(-1500,1500),notch=TRUE,col='#cff0da')
stripchart(emp_dots~factor(names(emp_dots),levels=restructured_col_plottree), horizontal=TRUE, 
           method = "jitter", add = TRUE, pch = 20,at = c(2:15)-0.5, col = 'red')
axis(1)
title(xlab="log(body size)")

boxplot(xe~factor(names(xe),levels=phylo_test$node.label),horizontal=TRUE,
        axes=FALSE,xlim=c(1,Ntip(phylo_test)),at = c(2:15)-0.5,varwidth=TRUE)
axis(1)
title(xlab="log(body size)")




savefile = paste0(dir,'predictcontrastimage_TVP_TV_TVM_sorted',count,'.png')
ggsave(savefile,p_finalTVM)




# phylomorphospace(phylo_test,cbind(obsZ_pic,predictZ_matrix_TVP[1,]),label="off",node.size=c(0,1))




