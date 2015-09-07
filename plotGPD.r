#plotGPD(z$mle, z$threshold, z$rate, z$n, z$npy, z$cov, z$data, z$xdata,
#        xlabel="Return period (years)", ylabel="Wind speed (m/s)", plotpoints=TRUE)
library(ggplot2)
plotGPD = function(a, u, la, n, npy, cov, dat, xdat, xlabel="",ylabel="",title="",plotpoints=FALSE)
{
    a <- c(la, a)
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps

    jj <- seq(0.0, 4.5 + log10(npy), by = 0.1)
    m <- c(1/la, 10^jj)
    q <- gpdq2(a[2:3], u, la, m)
    d <- t(gpd.rl.gradient(a = a, m = m))
    mat <- matrix(c((la * (1 - la))/n, 0, 0, 
                    0, cov[1, 1], cov[1, 2], 
                    0, cov[2, 1], cov[2, 2]), nc = 3)
    v <- apply(d, 1, q.form, m = mat)
    
    #qplot(m/npy, q, geom="line", xlim=c(1, max(m)/npy), 
    #      ylim = c(u, max(xdat, q[q > u - 1] + 1.96 * sqrt(v)[q >u - 1])),
    #      xlab = xlabel, ylab = ylabel,main = "Return Level Plot", log='x')
    #plot(m/npy, q, log = "x", type = "n", xlim = c(1, max(m)/npy), 
    #     ylim = c(u, max(xdat, q[q > u - 1] + 1.96 * sqrt(v)[q >u - 1])), 
    #     xlab = xlabel, ylab = ylabel, main = paste("Return Level Plot", title, sep=" - "))
    plot(m/npy, q, log = "x", type = "n", xlim = c(1, 10000), ylim = c(10, 100), 
         xlab = xlabel, ylab = ylabel, main = paste("Return Level Plot", title, sep=" - "))
    grid(equilogs=TRUE)
    lines(m[(m/npy)>1]/npy, q[(m/npy)>1], lwd=2)
    lines(m[(m/npy)>1]/npy, q[(m/npy)>1] + 1.96 * sqrt(v)[(m/npy)>1], col=4)
    lines(m[(m/npy)>1]/npy, q[(m/npy)>1] - 1.96 * sqrt(v)[(m/npy)>1], col=4)
    if(plotpoints==TRUE){
        nl <- n - length(dat) + 1
        sdat <- sort(xdat)
        emprp = (1/(1 - (1:n)/(n + 1))/npy)
        points(emprp[emprp > 1], sdat[emprp > 1], col="red", pch=16)
    }
}

plotGPD_CI = function(dataset,yrs,threshold,RP,confInt=0.95,plotpoints=FALSE)
{
    library("ismev")
    source(paste(sdir,"profileFit.r",sep="") )
    mxnp <- length(RP)
    nyrs <- round(length(dataset)/365,0)
    gpdfit <- profileFit(dataset,nyrs,threshold,RP,confInt)
    ymax<-10*(ceiling(gpdfit[mxnp,4]/10))
    plot(RP,gpdfit[,3],log="x",type="l",ylab="Wind speed (m/s)",
         xlab="Return Period (yrs)",xlim=c(1,10000),ylim=c(10,100), lwd=2)
    lines(cbind(RP,gpdfit[,2]),lty=2,col="blue")
    lines(cbind(RP,gpdfit[,4]),lty=2,col="blue")
    grid(equilogs=FALSE)
    if (plotpoints==TRUE){
      npy=365.25
      n = length(dataset)
      sdat = sort(dataset)
      emprp = (1/(1 - (1:n)/(n + 1))/npy)
      points(emprp[emprp > 1], sdat[emprp > 1], col="red", pch=16)
    }
}
