BHcorrection = function(p, thresh){
  sp = sort(p, index.return=TRUE)
  pp = sp$x
  ind = sp$ix
  cutoffs = which(pp < 1:length(p)/length(p) * thresh)
  if(length(cutoffs)==0){
    print("no selection"); return(NA)
  }else{
    where = 1:max(cutoffs)
    return(ind[where])
  }
}
