library(phytools)
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
plot(emdata,show.tip.label = FALSE)

# Pruned tree plot
dropextinct = T
L_ext = phylo2L(emdata)
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
plot(phylo_test,show.tip.label = TRUE)

# prune tree by phytools
phy_prune = fancyTree(emdata, type="droptip",tip = getExtinct(emdata),cex = 0.7)

L = L_ext

prune = 0
if (prune==1){
  L=pruneL(L)
}
phylo_prune = DDD::L2phylo(L,dropextinct = dropextinct)
plot(phylo_prune,show.tip.label = FALSE)


time.list = c(sort(c(L[,1],L[which(L[,4]!= -1),4]),decreasing = TRUE),0)
#total number of species
num.species = nrow(L)
trait.table = matrix(0,nrow = length(time.list)-1,ncol = nrow(L)+1)
time.branching = match(L[,1],time.list)
time.end = match(L[,4],time.list)
time.end[is.na(time.end)] = length(time.list)
survival.time = cbind(time.branching,time.end)
timelist = as.data.frame(time.list)
timebranch = as.data.frame(time.branching)
timeend = as.data.frame(time.end)


for(i in 1:num.species){
  
  trait.table[,i+1][time.branching[i]:(time.end[i]-1) ] = 1
}
trait.table = rbind(trait.table,trait.table[dim(trait.table)[1],])
trait.table[,1] = time.list
existing_species_table = trait.table[-1,-1]


write.csv(timelist, file = paste0(dir,"timelist.csv"))
write.csv(timebranch, file = paste0(dir,"timebranch.csv"))
write.csv(timeend, file = paste0(dir,"timeend.csv"))
write.csv(existing_species_table, file = paste0(dir,"traittable.csv"))
write.csv(L, file = paste0(dir,"Ltable.csv"))

para_list = list(pars=pars,prune=prune,seed=seed_fun)
# 
# output = list(L=L, timelist= timelist, timebranch = timebranch, timeend = timeend,traittable = existing_species_table)
write.csv(para_list, file = paste0(dir,"para_setting.csv"))