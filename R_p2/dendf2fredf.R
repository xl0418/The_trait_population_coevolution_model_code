dendf2fredf = function(dendf, iteration, population,leftendpoint,rightendpoint,binwidth){
  finvec=c()
  for(i in c(1:iteration)){
    startp = (i-1)*population+1
    endp = i*population
    focalvec = dendf$Samples[startp:endp]
    finvec1 = seprate_vec(vec = focalvec,leftendpoint = leftendpoint, 
                          rightendpoint = rightendpoint,
                          binwidth = binwidth)
    finvec = rbind(finvec, cbind(finvec1,i))
  }
  return(finvec)
}