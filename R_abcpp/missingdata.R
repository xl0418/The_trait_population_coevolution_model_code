source('C:/Liang/Code/Pro2/R_p2/phylo2L.R', echo=TRUE)
source('C:/Liang/Code/Pro2/R_p2/pruneL.R', echo=TRUE)
library(phytools)
emdatadir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/aLb_renamed.tre'
dir = 'C:/Liang/Googlebox/Research/Project2/planktonic_foraminifera_macroperforate/'

emdata = read.tree(emdatadir)
plot(emdata,show.tip.label = FALSE)


dropextinct = T
L_ext = phylo2L(emdata)
phylo_test = DDD::L2phylo(L_ext,dropextinct = dropextinct)
plot(phylo_test,show.tip.label = TRUE)

phy_de <- fancyTree(emdata, type = "droptip", tip = getExtinct(emdata), cex = 0.7)
phy_de$tip.label

missingdata = c("Globigerinella_adamsi","Globorotalia_theyeri","Globorotalia_ungulata",
                "Globorotalia_bermudezi","Globorotalia_cavernula")
missingindex = match(missingdata,phy_de$tip.label)
missingnames = phylo_test$tip.label[missingindex]
missingspecies_fulltree = as.numeric(gsub("t", "", missingnames))

allspecies= sort(as.numeric(gsub("t", "", phylo_test$tip.label)))
missingspecies_recontree = match(missingspecies_fulltree,allspecies)


