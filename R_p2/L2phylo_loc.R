L2phylo_loc=function(L, dropextinct = T, loc = 0) 
{
  if(loc == 0){    
    L = L[order(abs(L[, 3])), 1:4]
  }
  if(loc == 1){    
    L = L[order(abs(L[, 3])), 1:5]
  }
  
  age = L[1, 1]
  L[, 1] = age - L[, 1]
  L[1, 1] = -1
  notmin1 = which(L[, 4] != -1)
  L[notmin1, 4] = age - L[notmin1, 4]
  if (dropextinct == T) {
    sall = which(L[, 4] == -1)
    tend = age
  }
  else {
    sall = which(L[, 4] >= -1)
    tend = (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  }
  L = L[, -4]
  if(loc == 0){
    linlist = cbind(L[sall, ], paste("t", abs(L[sall, 3]), sep = ""), 
                    tend)
  }
  if(loc == 1){
    linlist = cbind(L[sall, ], abs(L[sall, 3]),tend)
    
    # linlist = cbind(L[sall, ], paste("loc", abs(L[sall, 4]), sep = ""), tend)
    linlist = linlist[,-4]
  }
  done = 0
  while (done == 0) {
    j = which.max(linlist[, 1])
    daughter = linlist[j, 3]
    parent = linlist[j, 2]
    parentj = which(parent == linlist[, 3])
    parentinlist = length(parentj)
    if (parentinlist == 1) {
      spec1 = paste(linlist[parentj, 4], ":", as.numeric(linlist[parentj, 
                                                                 5]) - as.numeric(linlist[j, 1]), sep = "")
      spec2 = paste(linlist[j, 4], ":", as.numeric(linlist[j, 
                                                           5]) - as.numeric(linlist[j, 1]), sep = "")
      linlist[parentj, 4] = paste("(", spec1, ",", spec2, 
                                  ")", sep = "")
      linlist[parentj, 5] = linlist[j, 1]
      linlist = linlist[-j, ]
    }
    else {
      linlist[j, 1:3] = L[which(L[, 3] == as.numeric(parent)), 
                          1:3]
    }
    if (is.null(nrow(linlist))) {
      done = 1
    }
  }
  linlist[4] = paste(linlist[4], ":", linlist[5], ";", sep = "")
  phy = read.tree(text = linlist[4])
  tree = as.phylo(phy)
  return(tree)
}