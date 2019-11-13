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
plot(phylo_test,show.tip.label = TRUE)


# empirical data
fileemp_name = paste0(dir,'emp.csv')
obsZ_emp = read.csv(fileemp_name)
obsZ_emp[,1] = extantspecieslabel
species_label = extantspecieslabel
sorted.species.labels <- obsZ_emp[order(obsZ_emp[,2]),1]

mode.label = c('TVP','TV','TVM')

for(count in mode.label){
  simfile = paste0(dir,'predictsim',count,'.csv')
  
  # predict simulations
  predictZ = read.csv(simfile)
  
  
  sort = 1
  
  if(sort == 0){
    predictZ_matrix = as.matrix(predictZ)
    
    
  }else{
    predictZ_matrix = as.matrix(predictZ)
    predictZ_matrix=t(apply(predictZ_matrix,1,sort))
    
    
  }
  
  samplesize = nrow(predictZ_matrix)
  dimnames(predictZ_matrix)[[2]] = sorted.species.labels
  
  
  d_all = as.data.frame(as.table(predictZ_matrix))[,2:3]
  colnames(d_all) = c('species','traitall')
  
  
  
  obsZ_mean = 10^(obsZ_emp[,2])
  
  
  
  d_meanemp = data.frame(species=species_label, trait=obsZ_mean)
  
  
  plot_tree <- ggtree(phylo_test)+geom_tiplab(size=4) #+xlim(0,80)
  
  
  plot_sepboxplt <- facet_plot(plot_tree, panel="The total length", data=d_all, geom_boxploth, 
                               mapping = aes(x=traitall, group=label ))  + theme_tree2()+
    theme(strip.background = element_rect(fill="#99CCFF"))
  
  p_finalTP <- facet_plot(plot_sepboxplt+xlim_tree(40), panel="The total length", data=d_meanemp, geom_point, 
                          mapping = aes(x=trait, group=label ),color = 'red')
  

  p_finalTP <- p_finalTP+ggtitle(count)+
    theme(plot.title = element_text(hjust = 0.5))
  
  savefile = paste0(dir,'predictimage',count,'.png')
  ggsave(savefile,p_finalTP)
}








