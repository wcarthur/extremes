library('stringr')

evmix_fit = new.env()
source('evmix_fit.R', local=evmix_fit, chdir=TRUE)


datapath <- "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/2019/Daily"
plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/daily_max_wind/wind/gpd_bayes"

filelist <- list.files(datapath, pattern = "^DC02D_Data.*\\.txt")

for (file in filelist) {
  SKIP <- FALSE
  stnnum <- str_extract(file, "\\d{6}")
  print(sprintf("Processing station %s", stnnum))
  filepath <- paste(datapath, file, sep='/')
  figname <- paste(stnnum, '.png', sep="")
  df <- read.csv(filepath, sep=',', header=TRUE)
  gust <- df$Speed.of.maximum.wind.gust.in.m.s[!is.na(df$Speed.of.maximum.wind.gust.in.m.s)
                                               & df$Quality.of.maximum.gust.speed == 'Y']
  if (length(gust)<100){
    print(sprintf("Insufficient observations - skipping %s", stnnum))
    next
  }
  gust <- gust[gust > 0]
  gust_offset = 0
  gust_mixture_fit = evmix_fit$fit_gpd_mixture(
    data=gust,
    continuous=FALSE,
    data_offset=gust_offset,
    bulk='gamma'
  )
  u_upper_limit = sort(gust, decreasing=TRUE)[10] - gust_offset
  u_lower_limit = sort(gust, decreasing=FALSE)[100] - gust_offset
  
  mcmc_chain_length = 100000 # Might need to increase this for convergence
  mcmc_chain_thin = 1 # No thinning by default, but to save memory can increase this
  mcmc_nchains = 1 # How many chains to run
  mcmc_ncores = 1 # Can run multiple chains in parallel if desired (not on windows)
  # Annual event rate -- estimate it from the data
  start_year = min(df$Year)
  end_year = max(df$Year)
  annual_event_rate = length(df[,1])/(end_year - start_year)
  gust_mixture_fit = evmix_fit$mcmc_gpd_mixture(
    fit_env=gust_mixture_fit, 
    par_lower_limits=c(0, 0, u_lower_limit, 0, 0), 
    par_upper_limits=c(1e+08, 1.0e+08, u_upper_limit, 20, 0.5),
    mcmc_start_perturbation=c(0.4, 0.4, 1., 0.2, 0.01), # If multiple chains, randomly perturb their starting values
    mcmc_length=mcmc_chain_length,
    mcmc_thin=mcmc_chain_thin,
    mcmc_burnin=10000,
    mcmc_nchains=mcmc_nchains,
    mcmc_tune=c(1,1,1,1,1)*1,
    mc_cores=mcmc_ncores,
    annual_event_rate=annual_event_rate,
    verbose=FALSE)
  
  # plot(gust_mixture_fit$mcmc_chains[[1]])
  png(filename = paste(plotpath, figname, sep = "/"))
  evmix_fit$mcmc_rl_plot(gust_mixture_fit)
  dev.off()
  
}