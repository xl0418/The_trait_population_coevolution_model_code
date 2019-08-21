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
dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/results_0724_contrast/'
dir_emp = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/result_cluster/Est/'
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
plot(phylo_test,show.tip.label = TRUE)


# empirical data
fileemp_name = paste0(dir_emp,'emp.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]
# Extract the order names of simulated species
prunedL = phylo2L(phylo_test)
pl = prunedL$L
phylo_pl = DDD::L2phylo(pl)
plot(phylo_pl,show.tip.label = TRUE,show.node.label = TRUE)

rank.sim.species <- as.numeric(gsub("t","", phylo_pl$tip.label))
sim.species.col.names = c()
for(i in c(1:15)){
  sim.species.col.names[rank.sim.species[i]] <- extantspecieslabel[i]
}

mode.label = c('TVP','TV','TVM')


simfile_TVP = paste0(dir,'predictsimTVP.csv')
simfile_TV = paste0(dir,'predictsimTV.csv')
simfile_TVM = paste0(dir,'predictsimTVM.csv')
simfile_NTVP=paste0(dir,'predictsimNTVP.csv')
simfile_NTVM=paste0(dir,'predictsimNTVM.csv')

# predict simulations
predictZ_TVP = read.csv(simfile_TVP)
predictZ_TV = read.csv(simfile_TV)
predictZ_TVM = read.csv(simfile_TVM)
predictN_TVP = read.csv(simfile_NTVP)/1e5
predictN_TVM = read.csv(simfile_NTVM)/1e5


sort = 1

if(sort == 0){
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  predictN_matrix_TVP = as.matrix(predictN_TVP)
  predictN_matrix_TVM = as.matrix(predictN_TVM)
  dimnames(predictZ_matrix_TVP)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TV)[[2]] = sim.species.col.names
  dimnames(predictZ_matrix_TVM)[[2]] = sim.species.col.names
  dimnames(predictN_matrix_TVP_sorted)[[2]] = sim.species.col.names
  dimnames(predictN_matrix_TVM_sorted)[[2]] = sim.species.col.names
  
  
  
}else{
  predictZ_matrix_TVP = as.matrix(predictZ_TVP)
  predictZ_matrix_TV = as.matrix(predictZ_TV)
  predictZ_matrix_TVM = as.matrix(predictZ_TVM)
  predictN_matrix_TVP = as.matrix(predictN_TVP)
  predictN_matrix_TVM = as.matrix(predictN_TVM)
  TVP_Z_order = t(apply(predictZ_matrix_TVP,1,order))
  TVM_Z_order = t(apply(predictZ_matrix_TVM,1,order))
  
  predictZ_matrix_TVP=t(apply(predictZ_matrix_TVP,1,sort))
  predictZ_matrix_TV=t(apply(predictZ_matrix_TV,1,sort))
  predictZ_matrix_TVM=t(apply(predictZ_matrix_TVM,1,sort))
  predictN_matrix_TVP_sorted = c()
  predictN_matrix_TVM_sorted = c()
  
  for(Nrow in c(1:nrow(predictN_matrix_TVP))){
    predictN_matrix_TVP_sorted=rbind(predictN_matrix_TVP_sorted,
                              predictN_matrix_TVP[Nrow,TVP_Z_order[Nrow,]])
    predictN_matrix_TVM_sorted=rbind(predictN_matrix_TVM_sorted,
                                     predictN_matrix_TVM[Nrow,TVM_Z_order[Nrow,]])
  }
  dimnames(predictZ_matrix_TVP)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TV)[[2]] = sorted.species.labels
  dimnames(predictZ_matrix_TVM)[[2]] = sorted.species.labels
  dimnames(predictN_matrix_TVP_sorted)[[2]] = sorted.species.labels
  dimnames(predictN_matrix_TVM_sorted)[[2]] = sorted.species.labels
  
}

samplesize = nrow(predictZ_matrix_TVP)

d_all_TVP = as.data.frame(as.table(predictZ_matrix_TVP))[,2:3]
colnames(d_all_TVP) = c('species','traitall')
d_all_TV = as.data.frame(as.table(predictZ_matrix_TV))[,2:3]
colnames(d_all_TV) = c('species','traitall')
d_all_TVM = as.data.frame(as.table(predictZ_matrix_TVM))[,2:3]
colnames(d_all_TVM) = c('species','traitall')
d_N_TVP = cbind(as.data.frame(as.table(predictN_matrix_TVP_sorted))[,2:3],'TVP')

colnames(d_N_TVP) = c('species','Pop','model')
d_N_TVM = cbind(as.data.frame(as.table(predictN_matrix_TVM_sorted))[,2:3],'TVM')
colnames(d_N_TVM) = c('species','Pop','model')

d_N = rbind(d_N_TVP,d_N_TVM)

meanpop_TVP = cbind(as.data.frame(as.table(colMeans(predictN_matrix_TVP_sorted))),'TVP')
meanpop_TVM = cbind(as.data.frame(as.table(colMeans(predictN_matrix_TVM_sorted))),'TVM')
colnames(meanpop_TVP) = c('species','Pop','model')
colnames(meanpop_TVM) = c('species','Pop','model')

df_mean_pop = rbind(meanpop_TVP,meanpop_TVM)


obsZ_mean = 10^(obsZ_emp[,2])



d_meanemp = data.frame(species=species_label, trait=obsZ_mean)


plot_tree <- ggtree(phylo_test,size=.5) +geom_tiplab(size=3.5,fontface="bold") #+xlim(0,80)


plot_sepboxplt_TVP <- facet_plot(plot_tree, panel="TVP", data=d_all_TVP, geom_boxploth, 
                             mapping = aes(x=traitall, group=label),color='#3ac569',fill= '#cff0da',outlier.colour = NULL)  + theme_tree2()

p_finalTVP <- facet_plot(plot_sepboxplt_TVP+xlim_tree(40), panel="TVP", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'black')


plot_sepboxplt_TV <- facet_plot(p_finalTVP, panel="TV", data=d_all_TV, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ),color='#f9320c',fill='#f1bbba',outlier.colour = NULL )  + theme_tree2()

p_finalTV <- facet_plot(plot_sepboxplt_TV+xlim_tree(40), panel="TV", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'black')


plot_sepboxplt_TVM <- facet_plot(p_finalTV, panel="TVM", data=d_all_TVM, geom_boxploth, 
                                mapping = aes(x=traitall, group=label),color='#6a60a9',fill='#dedcee' ,outlier.colour = NULL )  + theme_tree2()

p_finalTVM <- facet_plot(plot_sepboxplt_TVM+xlim_tree(40), panel="TVM", data=d_meanemp, geom_point, 
                        mapping = aes(x=trait, group=label ),color = 'black')+
  theme(strip.background = element_rect(colour="white", fill="white"), 
        strip.text.x = element_text(size=12, face="bold"),
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))

# 
# p_final_pop_TVP <- facet_plot(p_finalTVM, panel = 'pop', 
#                  data = df_mean_pop, geom = geom_barh,
#                  mapping =  aes(x=Pop, fill = model),width = 0.5,
#                  stat='identity' ) 
# 
# p_final_pop_TVP


lbs <- c(Tree = "Phylogenetic tree\nof baleen whales", TVP = "Trait evolution \n+ population dynamics",
         TV = "Trait evolution", TVM ="Trait evolution \n+ metabolism dynamics") # ,pop="Abundance distribution (10)")
facet_labeller(p_finalTVM, lbs)


savefile = paste0(dir,'predictimage_TVP_TV_TVM',count,'.png')
# ggsave(savefile,p_finalTVM)









