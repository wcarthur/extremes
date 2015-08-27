# Title: profileFit.r
# Author: Sanabria A., augusto.sanabria@ga.gov.au
# CreationDate: 2006-10-18
# Description: function to fit a GPD to given data and to calc. 95%
# confidence intervals using a profile likelihood approach.
# The program plots actual retp, retp modelled
# using a GPD (fitted via ML) and conf. interval
# Reference: Function written by AS using some functions from
# package 'extRemes'. 18/10/06:11:58 AM

#SeeAlso: package "ismev" has to be available for loading.
#Version: 1.1
#Modified by: AS
#ModifiedDate: 2006-12-18
#Modification: The function now reads an internal vector called 'dataset'
#intead of reading an external file ('filename')

#Version: 1.2
#Modified by: AS
#ModifiedDate: 2008-04-11
#Modification: The function now calculates the CI using a 'profile
#likelihood method' (more appropriate for non-symmetrical data),
#in the previous version the Delta method was used (ok, for
#symmetrical data).

#Usage:
profileFit = function(dataset,yrs,thres,RP,confInt=0.95,plotRP=FALSE,plotObs=FALSE){

#dataset = vector of clean wind speeds
#yrs = range of years in vector
#thres = threshold for wind speeds
#confInt = confidence interval (default 95%)

#Required:
library("ismev")
#Functions needed:
source(paste(sdir,"plotGPD.r",sep=""))
source(paste(sdir,"paramCI.r",sep=""))
#source(paste(sdir,"RP_xaxis.r",sep="") )

#Define return periods to plot:
#RP <- RP_xaxis(10,7,resolution="high")
#RP <- c(5,10,20,25,50,100,200,250,500,1000,2000,rmcnfInt2500,5000,10000)
mxnp <- length(RP)

#confInt <- 0.95


cov.names<-character(0)
covs <- NULL
sig.cov.cols<-NULL
gam.cov.cols<-NULL
sig.linker<-identity
gam.linker<-identity
xdata <- dataset
nopy <- length(dataset)/yrs


zz <- gpd.fit(xdat=xdata, threshold=thres, npy=nopy,
              ydat=covs, sigl=sig.cov.cols, siglink=sig.linker, 
              shl=gam.cov.cols, shlink=gam.linker, 
              method="BFGS", maxit=10000, show=FALSE)

# Calculate confidence intervals CI using the Profile likelihood method:
CI_low <- c()
CI_up <- c()
zh <- c()
for(i in 1:mxnp){
     confInt <- paramCI(zz, m=RP[i], conf=0.95, nint=100, rl.only=FALSE, xi.only=FALSE,
           rl.xup=NULL, rl.xlow=NULL, xi.xup=NULL,  xi.xlow=NULL, make.plot=FALSE)
     CI_low[i] <- confInt$rl$dn
     CI_up[i] <- confInt$rl$up
     # Calculate corresponding return period for point RP[i]:
     zh[i] <- confInt$rl$mle
}
#Plot RP with confidence interval:
#plotRP <- FALSE
if(plotRP){
    tit1 <- ""
    #tit1 <- try(titulo,silent=T)   #titulo is defined in shell
    ymax<-10*(ceiling(CI_up[mxnp]/10))
    plotGPD_CI(RP, zh, log="x", type="l", ylab="speed (m/s)", 
               xlab="Return Period (yrs)", xlim=c(10,10000),
               ylim=c(10,ymax),main=tit1)
    lines(cbind(RP,CI_low),lty=2,col="blue")
    lines(cbind(RP,CI_up),lty=2,col="blue")
    #text(60,65,"Return Period using a GPD distribution")
    grid(NULL, NULL)
    #If actual return periods in same plot are required change to TRUE:
    include_actualRP <- FALSE
    if(plotObs){
        n <- z$n
        npy <- z$npy
        xdat <- z$xdata
        sdat <- sort(xdat)
        u <- thres
        points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u],col="red")
    }
}

cbind(RP, CI_low, zh,CI_up)
}


