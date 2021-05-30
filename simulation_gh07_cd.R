### load packages and functions
rm(list=ls(all=TRUE))
setwd(dirname(parent.frame(2)$ofile))   # comment on servers

# install.packages("doParallel")
# install.packages("foreach")
# install.packages("iterators")
# install.packages("iterators")

source("src/gh07_cd.R")
source("src/gh07_cd_MF.R")
source("src/MF_DFA.R")

#library(doParallel,lib.loc="/your_path/R-lib")   # possibly needed on servers
#library(foreach,lib.loc="/your_path/R-lib")   # possibly needed on servers
require(doParallel)
require(foreach)

### set up the paralellization
cl=makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()
# clusterEvalQ(cl,.libPaths("/your_path/R-lib"))   # possibly needed on servers

start=Sys.time()
print(Sys.time())

### general simulation setup
runs=10   # set 1000 to replicate results in the paper
burn=100   # set 1000 to replicate results in the paper
obs=2000   # set 20000 to replicate results in the paper
metrics=4

### parameterization and data structures
k_val=1800
j_val=2
kk=c(1/1.25^5*k_val,1/1.25^4*k_val,1/1.25^3*k_val,1/1.25^2*k_val,1/1.25^1*k_val,k_val,1.25^1*k_val,1.25^2*k_val,1.25^3*k_val,1.25^4*k_val,1.25^5*k_val)   # $\psi$ ($\alpha$ in the original code)
jj=c(1/1.25^5*j_val,1/1.25^4*j_val,1/1.25^3*j_val,1/1.25^2*j_val,1/1.25^1*j_val,j_val,1.25^1*j_val,1.25^2*j_val,1.25^3*j_val,1.25^4*j_val,1.25^5*j_val)   # $\gamma$ ($\beta$ in the original code)
est=array(0,dim=c(length(kk),length(jj),runs,metrics))

### multifractal spectrum estimation setup
smin=10
smax=obs/5
sstep=10
qmax=4
qmin=-qmax
qstep=1

### run in paralell
for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    print(c(k,j))
    estpar=foreach(i=1:runs,.combine=rbind) %dopar%
    {
      gh07_cd_MF(i,obs,burn,smin,smax,sstep,qmin,qmax,qstep,kk[k],jj[j],metrics)
    }
    est[k,j,1:runs,1:metrics]=estpar
  }
}

### export data set with results
save(est,file="gh07_cd.Rda")

print(Sys.time())
print(Sys.time()-start)

# quit(save="no")   # possibly needed on servers