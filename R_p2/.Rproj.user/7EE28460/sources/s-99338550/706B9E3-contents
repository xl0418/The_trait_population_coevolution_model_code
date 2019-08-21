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
  emdatadir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/slater_mcct.txt'
  dir = 'C:/Liang/Googlebox/Research/Project2/BaleenWhales/treedata/'
}
emdata = read.nexus(emdatadir)
plot(emdata,show.tip.label = FALSE)
emdata_labelchange = emdata
newlabels = c()
for(i in c(1:78)){
  newlabels=c(newlabels,paste0('t',i ))
}
emdata_labelchange$tip.label = newlabels

plot(emdata_labelchange,show.tip.label = TRUE)


# Pruned tree plot
dropextinct = T
baleenwhale = phylo2L(emdata,error = 1e-5)
L_ext = baleenwhale$L
extantspecieslabel = baleenwhale$ESL
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
phylo_test$tip.label <- extantspecieslabel
plot(phylo_test,show.tip.label = TRUE)

# prune tree by phytools
phy_prune = fancyTree(emdata, type="droptip",tip = getExtinct(emdata),cex = 0.7)

L = L_ext

prune = 1
if (prune==1){
  L=pruneL(L)
}

# correct the crown age for baleen whales
L[1,] = L[2,]
L[1,2] = 0
L[1,3] = -1
L[2,3] = 2

positive.clade = c(-2)
do = TRUE
while(do){
  negative.row = which(match(L[,2],positive.clade)>0)
  if(length(negative.row) == 0){
    break
  }
  L[negative.row,2] = - L[negative.row,2]
  positive.clade = L[negative.row,3]
  L[negative.row,3] = - L[negative.row,3]

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
write.csv(extantspecieslabel, file = paste0(dir,"extantspecieslabels.csv"))

# para_list = list(pars=pars,prune=prune,seed=seed_fun)
# # 
# # output = list(L=L, timelist= timelist, timebranch = timebranch, timeend = timeend,traittable = existing_species_table)
# write.csv(para_list, file = paste0(dir,"para_setting.csv"))