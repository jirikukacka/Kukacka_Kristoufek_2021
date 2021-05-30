gh07_cd_MF<-function(seed,obs,burn,smin,smax,sstep,qmin,qmax,qstep,kkk,jjj,metrics){

  out=rep(0,metrics)

  ### generate original time series
  series=gh07_cd(seed,obs,burn,kkk,jjj)
  
  ### shuffle the series randomly
  seriesshuf=sample(series,replace=FALSE)
  
  ### estimate multifractality measures
  out[1:2]=MF_DFA_strength(series,smin,smax,sstep,qmin,qmax,qstep)
  out[3:4]=MF_DFA_strength(seriesshuf,smin,smax,sstep,qmin,qmax,qstep)
  
  return(out)
}