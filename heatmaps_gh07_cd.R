### load packages and functions
rm(list=ls(all=TRUE))
setwd(dirname(parent.frame(2)$ofile))

source('src/heatmap.2.R')

require(gplots)
require(gtools)
plot.dendrogram = stats:::plot.dendrogram

### Gaunersdorfer and Hommes (2007) model, $\gamma$ vs. $\psi$:
k_val=1800
j_val=2
kk=c(1/1.25^5*k_val,1/1.25^4*k_val,1/1.25^3*k_val,1/1.25^2*k_val,1/1.25^1*k_val,k_val,1.25^1*k_val,1.25^2*k_val,1.25^3*k_val,1.25^4*k_val,1.25^5*k_val)   # $\psi$ ($\alpha$ in the original code)
jj=c(1/1.25^5*j_val,1/1.25^4*j_val,1/1.25^3*j_val,1/1.25^2*j_val,1/1.25^1*j_val,j_val,1.25^1*j_val,1.25^2*j_val,1.25^3*j_val,1.25^4*j_val,1.25^5*j_val)   # $\gamma$ ($\beta$ in the original code)
rows=round(as.matrix(kk),1)
cols=round(as.matrix(jj),2)
kk.name=expression(psi)
jj.name=expression(gamma)

### load data set with results
load(file="gh07_cd.Rda")
runs=length(est[1,1, ,1])
metrics=length(est[1,1,1, ])

### compute individual $\Delta \alpha$ ratios original/shuffled
ratios=array(0,dim=c(length(kk),length(jj),runs))

for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    for (i in 1:runs)
    {  
      ratios[k,j,i]=est[k,j,i,2]/est[k,j,i,4]   
    }
  }
}

### detect potentially diverging runs
count=matrix(rep(0,length(kk)*length(jj)),nrow=length(kk),ncol=length(jj))

for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    count[k,j]=(runs-sum(is.na(ratios[k,j,1:runs])))/runs
  }
}

### compute multifractality strengths
filepdf="heatmap_gh07_c_mfst.pdf"
matrix=matrix(rep(0,length(kk)*length(jj)),nrow=length(kk),ncol=length(jj))
matrix=data.frame(matrix)
rownames(matrix)=rows
colnames(matrix)=cols

for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    matrix[k,j]=(floor(mean(ratios[k,j,1:runs],na.rm=TRUE)*10))/10   # floored for graphical purposes
    if (count[k,j] < 0.95) {   # only non-diverging setups considered in graphics
      matrix[k,j] = NA
    }
  }
}

pdf(filepdf)
heatmap.2(as.matrix(matrix),Colv=FALSE,Rowv=FALSE,scale="none",
          xlab=jj.name,ylab=kk.name,
          dendrogram='none',cellnote=round(as.matrix(matrix),1),notecex=1.5,notecol="black",trace='none',key=FALSE,
          lwid=c(.04,.99),lhei=c(.03,.99),margins=c(6,7),cexRow=1.2,cexCol=1.2,cex.lab=1.5,srtCol=33,adjCol=c(1,1))
dev.off()

### compute confidence levels
filepdf="heatmap_gh07_d_conf.pdf"
matrix=matrix(rep(0,length(kk)*length(jj)),nrow=length(kk),ncol=length(jj))
matrix=data.frame(matrix)
rownames(matrix)=rows
colnames(matrix)=cols

for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    matrix[k,j]=(floor(((sum(ratios[k,j,1:runs] > 1,na.rm=TRUE))/(runs-sum(is.na(ratios[k,j,1:runs]))))*100))/100   # floored for graphical purposes
    if (count[k,j] < 0.95 | is.nan(matrix[k,j]) == TRUE) {   # only non-diverging setups considered or get rid of NaNs in graphics
      matrix[k,j] = NA
    }
  }
}

pdf(filepdf)
heatmap.2(as.matrix(matrix),Colv=FALSE,Rowv=FALSE,scale="none",
          xlab=jj.name,ylab=kk.name,
          dendrogram='none',cellnote=round(as.matrix(matrix),2),notecex=1.5,notecol="black",trace='none',key=FALSE,
          lwid=c(.04,.99),lhei=c(.03,.99),margins=c(6,7),cexRow=1.2,cexCol=1.2,cex.lab=1.5,srtCol=33,adjCol=c(1,1))
dev.off()

### compute Welch's tests
filepdf="heatmap_gh07_diff.pdf"
matrix=matrix(rep(0,length(kk)*length(jj)),nrow=length(kk),ncol=length(jj))
matrix=data.frame(matrix)
rownames(matrix)=rows
colnames(matrix)=cols

benchmark_data=ratios[6,6,1:runs]

for (k in 1:length(kk))
{
  for (j in 1:length(jj))
  {
    matrix[k,j]=floor((1-t.test(benchmark_data, ratios[k,j,1:runs], alternative="two.sided", var.equal=FALSE)$p.value)*100)/100   # floored for graphical purposes
    if (is.nan(matrix[k,j]) == TRUE | count[k,j] < 0.95) {   # only non-diverging setups considered or get rid of NaNs in graphics
      matrix[k,j] = NA
    }
    if (is.na(matrix[k,j]) == FALSE && matrix[k,j] == 1) {   # check for exact 1, because nothing is 100% in statistical hypotheses testing...
      matrix[k,j] = 0.99
    }
  }
}

pdf(filepdf)
heatmap.2(as.matrix(matrix),Colv=FALSE,Rowv=FALSE,scale="none",
          xlab=jj.name,ylab=kk.name,
          dendrogram='none',cellnote=round(as.matrix(matrix),2),notecex=1.5,notecol="black",trace='none',key=FALSE,
          lwid=c(.04,.99),lhei=c(.03,.99),margins=c(6,7),cexRow=1.2,cexCol=1.2,cex.lab=1.5,srtCol=33,adjCol=c(1,1))
dev.off()