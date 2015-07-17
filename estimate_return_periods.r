library(evd)
library(ismev)
sdir <<- "C:/WorkSpace/extremes/"
source(paste(sdir,"select_threshold.r",sep="") )
source(paste(sdir,"GPD_distr.r",sep="") )
source(paste(sdir,"GPD_pars.r",sep="") )
source(paste(sdir,"plotGPD.r",sep="") )
source(paste(sdir,"profileFit.r",sep="") )
source(paste(sdir,"bisearch.R",sep=""))
indir = "C:/WorkSpace/data/daily/input/"
outdir = "C:/WorkSpace/data/daily/output/"

setwd(indir)
print(paste("Indir = ",indir,sep="") )

stationFile <- paste(indir,"station_names.txt",sep="")

stations <- read.table(stationFile,sep=",",header=T)
nstns <- length(stations$X)
retper <- c(2,5,10,20,25,50,100,200,250,500,1000,2000,2500,5000,10000)
all_pars_file <- paste(outdir,"station_params.txt",sep="")
all_RPs_file <- paste(outdir,"station_rp.txt",sep="")

table_head<- "Station, thresh, scale, shape, length_excess, length_all"
cat(table_head,fill=T,file=all_pars_file,append=F)

table_head <- paste("Station",toString(retper),sep=" ")
cat(table_head,file=all_RPs_file,fill=T,append=F)
for (i in 1:nstns){
    stnNum <- stations$Number[i]
    stnFile <- paste(indir,stations$File[i],sep="")
    print(paste( "Station file: ",stnFile,sep="") )

    # Load the data file:
    data <- read.table(stnFile,sep=",",header=T)

    # Read out the wind speed data:
    ws <- data$Speed.of.maximum.wind.gust.in.m.s
    #ws <- data$Speed.of.maximum.wind.gust.in.km.h # For 2012 data only

    # Remove missing values:
    ws<- ws[!is.na(ws)]
    
    # Convert to m/s (2012 data only)
    #ws = ws/3.6

    # Set an initial threshold for estimating the most suitable threshold
    init_u <- 0.5*max(ws)

    # Number of years
    nyrs <- round(length(ws)/365,0)

    # Evaluate the most appropriate threshold:
    # C:/WorkSpace/bin/process/extremes/sel_approp_u.r
    u <- select_threshold(ws,nyrs,init_u,0.01)

    # Generate a vector containing the GPD parameters
    # C:/WorkSpace/bin/process/extremes/GPD_pars.r
    par_vector <- c(GPD_pars(ws,nyrs,u[1]),length(ws))

    #RP <- GPD_distr(ws,nyrs,u[1],resolution="std",diag=FALSE,calcCI=T)
    # Calculate the return periods and confidence limits, using a profile-likelihood model
    # C:/WorkSpace/bin/process/extremes/profileFit.r
    RP <- profileFit(ws,nyrs,u[1],retper)

    z <- gpd.fit(ws, u[1], show=F)
    jpeg(paste(outdir,toString(stnNum),"_rl.jpg",sep=""),width=640,height=640,quality=100)
    plot.new()
    plotGPD(z$mle, z$threshold, z$rate, z$n, z$npy, z$cov, z$data, z$xdata,xlabel="Return period (years)",ylabel="Wind speed (m/s)",plotpoints=TRUE)
    dev.off()

    jpeg(paste(outdir,toString(stnNum),".jpg",sep=""),width=640,height=640,quality=100)
    plot.new()
    plotGPD_CI(ws,nyrs,u[1],retper,plotObs=F)
    dev.off()

    jpeg(paste(outdir,toString(stnNum),"_diag.jpg",sep=""),width=640,height=640,quality=100)
    plot.new()
    gpd.diag(z)
    dev.off()

    fname_rp <- paste(outdir,toString(stnNum),"_rp.txt",sep="")
    table_head <- paste("Return period:",toString(RP[,1]),sep=" ")
    cat(table_head,file=fname_rp,fill=T,append=F)
    cat(cbind("Lower CI"," ",t(signif(RP[,2],3)),sep=" " ),fill=T,file=fname_rp,append=T)
    cat(cbind("Mean RPs"," ",t(signif(RP[,3],3)),sep=" " ),fill=T,file=fname_rp,append=T)
    cat(cbind("Upper CI"," ",t(signif(RP[,4],3)),sep=" " ),fill=T,file=fname_rp,append=T)

    fname_pars <- paste(outdir,toString(stnNum),"_pars.txt",sep="")
    cat("Parameters of GPD:",file=fname_pars,fill=T,append=F)
    cat("thresh   scale   shape   length_excess length_all",fill=T,file=fname_pars,append=T)
    cat(cbind(t(signif(par_vector,5)),sep=" " ),fill=T,file=fname_pars,append=T)

    cat(cbind(toString(stnNum),t(signif(RP[,3],3)),sep=" " ),fill=T,file=all_RPs_file,append=T)
    cat(cbind(toString(stnNum),t(signif(par_vector,5)),sep=" " ),fill=T,file=all_pars_file,append=T)
}

