library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggstance)
library(ggimage)
os = Sys.info()['sysname']
if(os == 'Darwin'){
  source('~/Documents/GitHub/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
  source('~/Documents/GitHub/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
  emdatadir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
  dir = '~/GoogleDrive/Research/Project2/planktonic_foraminifera_macroperforate/'
  
}else{
source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
dir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/'
}
emdata = read.tree(emdatadir)
# prune tree by phytools
phy_prune = fancyTree(emdata, type="droptip",tip = getExtinct(emdata),cex = 0.7)

# trait data 
setwd('C:/Liang/Code/Pro2/data')

sort=1
if(sort==0){
  fileunsort_name = 'unsort.csv'
  obsZ_read = read.csv(fileunsort_name)
  species_label = phy_prune$tip.label
}else{
  fileunsort_name = 'sort.csv'
  obsZ_read = read.csv(fileunsort_name)
  fileemp_name = 'emp.csv'
  obsZ_emp = read.csv(fileemp_name)
  species_label = obsZ_emp[,1]
}

obsZ_matrix = as.matrix(obsZ_read)
obsZ_mean = colMeans(obsZ_matrix)

samplesize = nrow(obsZ_matrix)
dimnames(obsZ_matrix)[[2]] = species_label

# Group out the species with missing data 
missingtips = obsZ_emp[which(obsZ_emp[,2]==-1),1]
missingtree <- groupOTU(phy_prune, .node=missingtips)

groupout = 1

d_all = as.data.frame(as.table(obsZ_matrix))[,2:3]
colnames(d_all) = c('species','traitall')

d_mean1 = data.frame(species=species_label, trait=obsZ_mean)
d_mean2 = data.frame(species=species_label, traitdot=obsZ_mean)

if(groupout == 1){
  plot_tree = ggtree(missingtree, aes(color=group))+
          scale_colour_manual(values=c('#285943','#A593E0'))
}else{
  plot_tree <- ggtree(phy_prune)
}

plot_dottips = plot_tree %<+% d_mean1+ geom_tippoint(aes(size=trait))

plot_sepdots = facet_plot(plot_dottips, panel="dot", data=d_mean2, geom=geom_point, aes(x=traitdot), color='firebrick')+ theme_tree2()

plot_sepboxplt <- facet_plot(plot_dottips, panel="Trait", data=d_all, geom_boxploth, 
                 mapping = aes(x=traitall, group=label ))  + theme_tree2()+
                 theme(strip.background = element_rect(fill="#D1B6E1"))
  # geom_phylopic(image="6fe75b0b-1488-4193-8523-f240c2d59575", color="#cff0da", alpha = .1, size=Inf)


plot_sepboxplt





