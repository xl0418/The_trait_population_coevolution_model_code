library(DDD)
dir = 'C:/Liang/Googlebox/Research/Project2/treesim_newexp/example'

for(i in c(1:14)){
  filename = paste0(dir,i,'/Treedata.Rdata')
  pngname = paste0(dir,i,'/testtree',i,'.png')
  load(filename)
  phy = DDD::L2phylo(result$L,dropextinct = F)
  png(filename=pngname)
  plot(phy,show.tip.label=FALSE)
  # title(paste("Phylogenetic tree",i),line=-25)
  dev.off()
}