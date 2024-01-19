#Baypass trial

install.packages("corrplot")
install.packages("ape")
install.packages("WriteXLS")
install.packages("mvtnorm")
source("baypass_utils.R")
library(corrplot)
library(ape)
library(WriteXLS)
library(mtvnorm)
setwd("~/Jenny-Eurytemora_affinis/baypass-corrected")

# Read omega matrix BayPass output
SG.omega=as.matrix(read.table("correct-haplo_mat_omega.out"))

# Identify pool names for plotting
colnames(SG.omega) <- c("2009.start","2009.end","2011.start","2011.end",
                        "2015.start","2015.end","2022.start","2022.end")
rownames(SG.omega) <- c("2009.start","2009.end","2011.start","2011.end",
                        "2015.start","2015.end","2022.start","2022.end")

# Create a correlation matrix of the omega values, which can be used to assess genomic differentiation between pools
cor.mat=cov2cor(SG.omega)
corrplot(cor.mat,method="color",mar=c(2,1,2,2)+0.1,
         main=expression("Correlation map based on"~hat(Omega)))

# We can also assess population differentiation with hierarchical clustering:
bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))
plot(bta14.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

# Read the xtx BayPass output
SG.snp.res=read.table("correct-haplo_summary_pi_xtx.out",h=T)

# Get the Pi Beta distribution for POD generation
SG.pi.beta.coef=read.table("correct-haplo_summary_beta_params.out",h=T)$Mean

# Upload original data to get read counts
SG.data<-geno2YN("genobaypass")

# Simulate POD dataset to use for outlier SNP detection
simu.SG <-simulate.baypass(omega.mat=SG.omega,nsnp=10000,sample.size=SG.data$NN,
                           beta.pi=SG.pi.beta.coef,pi.maf=0,suffix="SG.BP.sim")

#plot omega
plot.omega <- function(omega,PC=c(1,2),pop.names=paste0("Pop",1:nrow(omega)),main=expression("SVD of "*Omega),col=rainbow(nrow(omega)),pch=16,pos=2){
  om.svd=svd(omega)
  eig=om.svd$d
  pcent.var=100*eig/sum(eig)
  plot(om.svd$u[,PC],main=main,pch=pch,col=col,
       xlab=paste0("PC",PC[1]," (",round(pcent.var[PC[1]],2),"%)"),
       ylab=paste0("PC",PC[2]," (",round(pcent.var[PC[2]],2),"%)")
  )
  text(om.svd$u[,PC[1]],om.svd$u[,PC[2]],pop.names,col=col,pos=pos)
  list(PC=om.svd$u,eig=eig,pcent.var=pcent.var)
}

pops=c("2009.start","2009.end","2011.start","2011.end",
       "2015.start","2015.end","2022.start","2022.end")
om.bta.svd=plot.omega(omega=SG.omega,pop.names=pops)

plot.omega(SG.omega,PC=c(1,2),pop.names=paste0("Pop",1:nrow(SG.omega)),
           main=expression("SVD of "*Omega),col=rainbow(nrow(SG.omega)),
           pch=16,pos=2)


#Estimates of the XtX differentiation measures (using the calibrated XtXst estimator)
anacore.snp.res=read.table("correct-haplo_summary_pi_xtx.out",h=T)


#check behavior of the p-values associated to the XtXst estimator
hist(10**(-1*anacore.snp.res$log10.1.pval.),freq=F,breaks=50)
abline(h=1)
layout(matrix(1:2,2,1))
plot(anacore.snp.res$XtXst)
plot(anacore.snp.res$log10.1.pval.,ylab="XtX P-value (-log10 scale)")
abline(h=3,lty=2) #0.001 p--value theshold


# Look at POD xtx values, and identify SNPs where the xtx values are above the 99% significance threshold from the POD. So in the plot, it is those loci (dots) which are above the abline
anacore.snp.res.xtx=read.table("correct-haplo_summary_pi_xtx.out",h=T)$M_XtX
anacore.thresh=quantile(anacore.snp.res.xtx,probs=0.99)
plot(anacore.snp.res$M_XtX)
abline(h=anacore.thresh,lty=2)

anacore.thresh
###
genodi=as.matrix(read.table("correct-haplo_summary_betai_reg.out"))
