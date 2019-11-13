source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/Nindex.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/sddsim.R', echo=TRUE)
source('C:/Liang/Googlebox/Research/Project1/R_pro1/Final/event_matrix.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/Plottree_single_Pro1.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)

library(DDD)
library(MASS)
library(rgl)
library(stringr)
# library(matrixcalc)
library("reshape2")
library('Matrix')
library(plyr) 
# library(twitteR)

pars=c(0.8,0.2,300)
seed_fun=29
prune=1
result = sddsim(n=2,parsN=c(2,0),age=70,pars=pars , seed_fun = seed_fun, lambda_allo0 = 0.2, M0=0,K_fix = 1)
# dir = 'C:/Liang/PhdIntroProject2/Example/'
dir = 'C:/Liang/Googlebox/Research/Project2/treesim_newexp/largetree/'
filename = paste0(dir, 'Treedata.Rdata')
save(result,file = filename)
print(dim(result$L))
plottree(file = filename,dropextinct =F)

load(file = filename)

L = result$L
if (prune==1){
  L=pruneL(L)
}
phylo_p=L2phylo(L,dropextinct =F)
plot(phylo_p,show.tip.label = F)

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

