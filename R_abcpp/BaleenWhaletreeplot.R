library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
library(DDD)
os = Sys.info()['sysname']
if(os == 'Darwin'){
  source('~/Documents/GitHub/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('~/Documents/GitHub/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
  dir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/'
  
}else{
  source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
  dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'
}
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
fileemp_name = paste0(dir,'trait_labels.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]


filepredict_name =  paste0(dir,'predictsimTPUS.csv')
predictZ = read.csv(filepredict_name)
filepredict_DR =  paste0(dir,'predictsimDR_fixm.csv')
predictZDR = read.csv(filepredict_DR)
# filepredict_NH =  paste0(dir,'predictsimNHUS5s.csv')
# predictZNH = read.csv(filepredict_NH)

sort = 0

if(sort == 0){
  predictZ_matrix = as.matrix(predictZ)
  
  predictZ_matrixDR = as.matrix(predictZDR)
  
  # predictZ_matrixNH = as.matrix(predictZNH)
}else{
  predictZ_matrix = as.matrix(predictZ)
  predictZ_matrix=t(apply(predictZ_matrix,1,sort))
  predictZ_matrixDR = as.matrix(predictZDR)
  predictZ_matrixDR=t(apply(predictZ_matrixDR,1,sort))
  
  # predictZ_matrixNH = as.matrix(predictZNH)
  # predictZ_matrixNH=t(apply(predictZ_matrixNH,1,sort))
  
}




samplesize = nrow(predictZ_matrix)

dimnames(predictZ_matrix)[[2]] = sorted.species.labels
dimnames(predictZ_matrixDR)[[2]] = sorted.species.labels
# dimnames(predictZ_matrixNH)[[2]] = sorted.species.labels

d_all = as.data.frame(as.table(predictZ_matrix))[,2:3]
colnames(d_all) = c('species','traitall')
d_allDR = as.data.frame(as.table(predictZ_matrixDR))[,2:3]
colnames(d_allDR) = c('species','traitall')
# d_allNH = as.data.frame(as.table(predictZ_matrixNH))[,2:3]
# colnames(d_allNH) = c('species','traitall')

obsZ_mean = obsZ_emp[,2]

d_meanemp = data.frame(species=species_label, trait=obsZ_mean)


plot_tree <- ggtree(phylo_test)+geom_tiplab(size=4) #+xlim(0,80)

# plot_dottips = plot_tree %<+% d_meanemp+ geom_tippoint(aes(size=trait))




# plot_sepdots = facet_plot(plot_tree, panel="dot", data=d_meanemp, geom=geom_point, aes(x=trait), color='firebrick')+ theme_tree2()

plot_sepboxplt <- facet_plot(plot_tree, panel="TP Trait", data=d_all, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#D1B6E1"))

p_finalTP <- facet_plot(plot_sepboxplt, panel="TP Trait", data=d_meanemp, geom_point, 
           mapping = aes(x=trait, group=label ),color = 'red')



plot_sepboxpltDR <- facet_plot(p_finalTP, panel="DR Trait", data=d_allDR, geom_boxploth, 
                             mapping = aes(x=traitall, group=label ))  + theme_tree2()+
  theme(strip.background = element_rect(fill="#D1B6E1"))

p_finalDR <- facet_plot(plot_sepboxpltDR+xlim_tree(30), panel="DR Trait", data=d_meanemp, geom_point, 
                      mapping = aes(x=trait, group=label ),color = 'red')

p_finalDR

# plot_sepboxpltNH <- facet_plot(p_finalDR, panel="NH Trait", data=d_allNH, geom_boxploth, 
#                              mapping = aes(x=traitall, group=label ))  + theme_tree2()+
#   theme(strip.background = element_rect(fill="#D1B6E1"))

# p_final <- facet_plot(plot_sepboxpltNH+xlim_tree(45), panel="NH Trait", data=d_meanemp, geom_point, 
#                       mapping = aes(x=trait, group=label ),color = 'red')

# p_final









