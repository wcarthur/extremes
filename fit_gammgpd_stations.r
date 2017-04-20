suppressPackageStartupMessages(library(evd))
suppressPackageStartupMessages(library(ismev))

sdir <<- "C:/WorkSpace/extremes/"
swcdir <<- "C:/WorkSpace/swc/R/evmix_fit/"
source(paste(sdir, "select_threshold.r", sep = ""))
source(paste(sdir, "GPD_distr.r", sep = ""))
source(paste(sdir, "GPD_pars.r", sep = ""))
source(paste(sdir, "plotGPD.r", sep = ""))
source(paste(sdir, "profileFit.r", sep = ""))
source(paste(sdir, "bisearch.R", sep = ""))
source(paste(swcdir, "evmix_fit.R", sep = ""))

indir = "C:/WorkSpace/data/daily/input/2012/"
outdir = "C:/WorkSpace/data/daily/output/2012/"

setwd(indir)
print(paste("Indir = ", indir ,sep = "") )

stationFile <- paste(indir, "WStn_names.txt", sep = "")

stations <- read.table(stationFile, sep = ",", header = T)
nstns <- length(stations$X)
retper <- c(2,5,10,20,25,30,40,50,75,100,200,250,300,400,500,750,1000,2000,2500,3000,4000,5000,7500,10000)
all_pars_file <- paste(outdir, "station_params_mixture.txt", sep = "")
all_RPs_file <- paste(outdir, "station_ari_bayes.txt", sep = "")

table_head <- "Station, gshape, gscale, u, sigmau, xi, phiu"
cat(table_head, fill = T, file = all_pars_file, append = F)

#table_head <- paste("Station",toString(retper),sep=" ")
#cat(table_head,file=all_RPs_file,fill=T,append=F)
for (i in 1:nstns) {
    stnNum <- stations$Number[i]
    stnName <- stations$Name[i]
    stnFile <- paste(indir,stations$File[i],sep = "")
    print(paste( "Station file: ",stnFile,sep = "") )

    # Load the data file:
    if (file.exists(stnFile)) {
      data <- read.table(stnFile, sep = ",", header = T)
    }else{
      print(paste("Station file ", stnFile, " does not exist", sep = ""))
      next
    }
    # Read out the wind speed data:
    #ws <- data$Speed.of.maximum.wind.gust.in.m.s
    ws <- data$Speed.of.maximum.wind.gust.in.km.h # For 2012 data only

    # Remove missing values:
    ws <- ws[!is.na(ws)]

    # Convert to m/s (2012 data only)
    ws = ws/3.6

    # Create the fit environment for the dataset:
    fitenv = fit_gpd_mixture(data = ws,
                             gpd_threshold_quantile_range = c(0.5, 0.995))

    # Run the MCMC fitting routine:
    mcmc_gpd_mixture(fitenv,
                     annual_event_rate = 365.)

    # Create the plot of ARI values
    print("Plotting return levels with Bayesian credible intervals")
    jpeg(paste(outdir,toString(stnNum), "_mcmc_ari.jpg", sep = ""),
         width = 640, height = 640, quality = 100)
    plot.new()
    mcmc_rl_plot(fitenv, xlim = c(10, 0.0001))
    dev.off()

    # Set an initial threshold for estimating the most suitable threshold
    #init_u <- round(median(ws), 3) # 0.3*max(ws)

    # Number of years
    #nyrs <- round(length(ws)/365,0)

    ## Evaluate the most appropriate threshold:
    #print(paste("Initial threshold: ", init_u, sep = ""))
    #u <- select_threshold(ws, nyrs, init_u, 0.01)
    #if (u[1] == 0) {
    #  print(paste("Fitting routine failed for ",stnNum,sep = ""))
    #  next
    #}
    ## Generate a vector containing the GPD parameters
    #par_vector <- c(GPD_pars(ws,nyrs,u[1]),length(ws))

    ##RP <- GPD_distr(ws,nyrs,u[1],resolution="std",diag=FALSE,calcCI=T)
    ## Calculate the return periods and confidence limits, using a profile-likelihood model
    #RP <- profileFit(ws, nyrs, u[1], retper)

    #z <- gpd.fit(ws, u[1], show = F)
    #print("Plotting return level with observations")
    #jpeg(paste(outdir,toString(stnNum), "_rl.jpg", sep=""),
    #     width=640, height=640, quality=100)
    #plot.new()
    #plotGPD(z$mle, z$threshold, z$rate, z$n, z$npy, z$cov, z$data, z$xdata,
    #        xlabel="Return period (years)", ylabel="Wind speed (m/s)",
    #        title=stnName, plotpoints=TRUE)
    #dev.off()

    #print("Plotting return levels with MLE CIs")
    #jpeg(paste(outdir,toString(stnNum),".jpg",sep=""), width=640, height=640, quality=100)
    #plot.new()
    #plotGPD_CI(ws, nyrs, u[1], retper, plotpoints=T)
    #dev.off()

    #print("Plotting fit diagnostics")
    #jpeg(paste(outdir, toString(stnNum), "_diag.jpg", sep=""),
    #     width=640, height=640, quality=100)
    #plot.new()
    #tryCatch(gpd.diag(z), error=function(e) print("Failed"), finally=print("passed"))
    #dev.off()

    #fname_rp <- paste(outdir,toString(stnNum),"_rp.txt",sep="")
    #table_head <- paste("Return period:",toString(RP[,1]),sep=" ")
    #cat(table_head,file=fname_rp,fill=T,append=F)
    #cat(cbind("Lower CI:"," ",t(round(RP[,2],2)),sep=" " ),fill=T,file=fname_rp,append=T)
    #cat(cbind("Mean RP:"," ",t(round(RP[,3],2)),sep=" " ),fill=T,file=fname_rp,append=T)
    #cat(cbind("Upper CI:"," ",t(round(RP[,4],2)),sep=" " ),fill=T,file=fname_rp,append=T)

    #fname_pars <- paste(outdir,toString(stnNum),"_pars.txt",sep="")
    #cat("Parameters of GPD:",file=fname_pars,fill=T,append=F)
    #cat("thresh   scale   shape   length_excess length_all",fill=T,file=fname_pars,append=T)
    #cat(cbind(t(round(par_vector,5)),sep=" " ),fill=T,file=fname_pars,append=T)

    #cat(cbind(toString(stnNum),t(round(RP[,3],2)),sep=" " ),fill=T,file=all_RPs_file,append=T)
    par_vector = c(fitenv$fit_evmix$gshape, fitenv$fit_evmix$gscale, fitenv$fit_evmix$u,
                   fitenv$fit_evmix$sigmau, fitenv$fit_evmix$xi, fitenv$fit_evmix$phiu)
    cat(cbind(toString(stnNum),t(round(par_vector,5)), sep = " " ),
            fill = T,file = all_pars_file, append = T)
}

