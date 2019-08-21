seprate_vec = function(vec,rightendpoint,leftendpoint,binwidth=0.1){
  sepvec = c()
  length = (rightendpoint-leftendpoint)/binwidth
  leftend = leftendpoint
  rightend = leftendpoint+binwidth
  for(i in c(1:length)){
    num = sum(vec > leftend& vec < rightend)
    sepvec = rbind(sepvec, c(num,(leftend+rightend)/2))
    leftend = rightend
    rightend = leftend + binwidth
  }
  return(sepvec)
}
