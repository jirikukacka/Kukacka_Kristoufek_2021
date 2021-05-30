# fluctuation function for fast linear trend fits for MF-DFA
fluctuation<-function(x){
  t<-1:length(x)
  res<-x-t*sum((t-mean(t))*x)/sum((t-mean(t))^2)
  return(mean((res-mean(res))^2))
}

# MF_DFA function returns a matrix of generalized Hurst exponents H(q) for a given range of qs based on MF-DFA
# function arguments:
# - x - time series of (log-)returns
# - smin, smax - minimum and maximum scales to consider for the scaling law
# - sstep - step length for scales between smin and smax
# - qmin, qmax - minim and maximum moments q for estimation of the multifractal spectrum
# - qstep - step length for moments between qmin and qmax
MF_DFA<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  xx<-cumsum(x-mean(x))
  qcount<-(qmax-qmin)/qstep+1
  F2_s<-matrix(rep(NA,(qcount+1)*((smax-smin)/sstep+1)),ncol=qcount+1)
  H_q<-matrix(rep(NA,2*qcount),ncol=2)
  counter_s<-1
  for(s in seq(smin,smax,by=sstep)){
    cut<-length(xx)-floor(length(xx)/s)*s
    if(cut==0){
      xmat_L<-matrix(xx,nrow=s,ncol=floor(length(xx)/s))
    } else {
      xmat_L<-matrix(xx[-c((length(xx)-cut+1):length(xx))],nrow=s,ncol=floor(length(xx)/s))
    }
    F2_L<-apply(xmat_L,2,fluctuation)
    if(cut==0){
      F2_R<-NULL
    } else {
      xmat_R<-matrix(xx[-c(1:cut)],nrow=s,ncol=floor(length(xx)/s))
      F2_R<-apply(xmat_R,2,fluctuation)
    }
    F2<-c(F2_L,F2_R)
    F2_s[counter_s,1]<-s
    counter_q<-2
    for(q in seq(qmin,qmax,by=qstep)){
      if(q==0){
        F2_s[counter_s,counter_q]<-exp(mean(log(F2))/2)
      } else {
        F2_s[counter_s,counter_q]<-mean(F2^(q/2))^(1/q)  
      }
      counter_q<-counter_q+1
    }
    counter_s<-counter_s+1
  }
  counter_q<-1
  for(q in seq(qmin,qmax,by=qstep)){
    H_q[counter_q,1]<-q
    H_q[counter_q,2]<-lm(log(F2_s[,(counter_q+1)])~log(F2_s[,1]))$coefficients[2]
    counter_q<-counter_q+1
  }
  return(H_q)
}

# MF_DFA_strength function returns $\Delta H$ and $\Delta \alpha$
# the parameters hold from the MF_DFA function above
MF_DFA_strength<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  H_q<-MF_DFA(x,smin,smax,sstep,qmin,qmax,qstep)
  delta_H<-max(H_q[,2])-min(H_q[,2])
  qcount<-(qmax-qmin)/qstep+1
  alpha_f<-H_q[-1,]
  alpha_f[,1]<-H_q[-1,][,2]+H_q[-1,][,1]*(H_q[-1,][,2]-H_q[-qcount,][,2])
  alpha_f[,2]<-1+H_q[-1,][,1]*(alpha_f[,1]-H_q[-1,][,2])
  delta_alpha<-max(alpha_f[,1])-min(alpha_f[,1])
  return(c(delta_H,delta_alpha))
}