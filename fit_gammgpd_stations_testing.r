swcdir <<- "C:/WorkSpace/swc/R/evmix_fit/"
setwd(swcdir)

source(paste(swcdir, "evmix_fit.R", sep = ""))
source(paste(swcdir, "fgammagpd.r", sep = ""))

indir = "C:/WorkSpace/data/Raw/daily_max_wind_gust/"

stationFile = paste(indir, "DC02D_StnDet_999999999425050.txt", sep="")
stations = read.table(stationFile, sep=",", header=T)
nstns = length(stations$Record.identifier)
nstns
retper <- c(2,5,10,20,25,30,40,50,75,100,200,250,300,400,500,750,1000,2000,2500,3000,4000,5000,7500,10000)
stnNum = stations$Bureau.of.Meteorology.Station.Number[17]

stnNum= sprintf("%06d", stnNum)

stnFile=paste(indir, "DC02D_Data_", stnNum, "_999999999425050.txt", sep="")
stnFile
data <- read.table(stnFile, sep = ",", header = T)
ws <- data$Speed.of.maximum.wind.gust.in.km.h
ws <- ws[!is.na(ws)]
ws = ws/3.6
#source('C:/WorkSpace/swc/R/evmix_fit/fgammagpd.r')
#setwd("C:/WorkSpace/swc/R/evmix_fit")
#source(paste(swcdir, "evmix_fit.R", sep = ""))
fitenv = fit_gpd_mixture(data = ws,
                gpd_threshold_quantile_range = c(0.995, 1),
                continuous = T)

ws_u_upper_limit = sort(ws, decreasing=TRUE)[50]
ws_u_lower_limit = sort(ws, decreasing=FALSE)[50]
par_upper_limits = c(1e+8, 1e+8, ws_u_upper_limit, 1000)
par_lower_limits = c(0, 0, ws_u_lower_limit, -1000)

mcmc_gpd_mixture(fitenv, annual_event_rate = 365.,
  mcmc_nchains = 1,
  mcmc_length = 10000,
  mcmc_burnin = 1000)

mcmc_rl_plot(fitenv, xlim = c(1, 0.0001)); grid()
plot(fitenv$mcmc_chains[[1]][,1])
plot(fitenv$mcmc_chains[[1]][,2])
plot(fitenv$mcmc_chains[[1]][,3])
plot(fitenv$mcmc_chains[[1]][,4])

qqplot(
  ws,
  fitenv$qfun(runif(length(ws))),
  main='QQ-plot of data and a random sample from \n the fitted gamma-GPD model',
  xlab='Data quantile', ylab='Model quantile')
abline(0,1,col='red'); grid()

# Parameters
lapply(fitenv$mcmc_chains, f<-function(x) summary(as.matrix(x)))

# 100-year ARI quantiles
lapply(fitenv$ari_100_chains,
  f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
# 500-year ARI quantiles
lapply(fitenv$ari_500_chains,
  f<-function(x) quantile(x, p=c(0.025, 0.5, 0.975)))
