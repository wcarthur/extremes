# Mixture model for extreme value distributions using daily maximum wind speeds extracted from 1-minute observation data
# This one fits a lognormal/GPD distribution
library('evmix')
library('stringr')
library('psych')

#datapath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind"
#plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind/lognormgpd"
#stationfile <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/stationlist.csv"

datapath <- "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/2019/Daily"
plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/daily_max_wind/wind/lognormgpd"
stationfile <- "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/2019/Daily/DC02D_StnDet_999999999632559_updated.txt"

stndf <- read.csv(stationfile)
filelist <- list.files(datapath, pattern = "DC02D_Data.*.txt")
colnames <- c("stnnum", "stnlat", "stnlon", "lnmean", "lnsd", "u", "sigmau", "xi", "phiu", "rsq", "rms")
outdf <- data.frame(stnnum = character(),
                    lnmean = double(),
                    lnsd = double(),
                    u = double(),
                    sigmau = double(),
                    xi = double(),
                    phiu = double(),
                    lat = double(),
                    lon = double(),
                    rsq = double(),
                    rms = double(),
                    stringsAsFactors = FALSE)
names(outdf) <- colnames

for (file in filelist) {
  SKIP <- FALSE
  stnnum <- str_extract(file, "\\d{6}")
  stnrec <- stndf[stndf$Bureau.of.Meteorology.Station.Number==as.integer(stnnum),]
  print(sprintf("Processing station %s", stnnum))
  filepath <- paste(datapath, file, sep = "/")
  figname <- paste(stnnum, ".png", sep = "")
  df <- read.csv(filepath, sep = ",", header = TRUE)
  gust <- df$Speed.of.maximum.wind.gust.in.m.s[df$Quality.of.maximum.gust.speed == "Y"]
  gust <- gust[!is.na(gust)]
  if (length(gust) < 365 * 5) {
    print("Insufficient observations")
    next
  }

  qs <- quantile(gust, seq(0.9, 1, 0.01))

  result <- tryCatch(flognormgpd(gust, useq = qs, fixedu = FALSE),
                     error = function(e) {SKIP <- TRUE})

  if (typeof(result) == "logical") {
    print(paste("Error fitting for station", stnnum))
    next
  } else {
    m.lognormgpd <- result
  }
  if (is.null(m.lognormgpd$cov)){
    rsq <- NA
  } else {
    invcov <- 1 - 1/solve(m.lognormgpd[["cov"]])
    rsq <- invcov[1, 1]
  }
  obscdf <- ecdf(gust)
  mycdf <- with(m.lognormgpd,
                plognormgpd(sort(gust), 
                             lnmean, 
                             lnsd,
                             u, sigmau,
                             xi, 
                             phiu,
                             lower.tail = TRUE))
  
  
  rms <- sqrt(mean((obscdf(sort(gust)) - mycdf)^2))
  row <- c(stnnum,
           stnrec$Latitude.to.4.decimal.places.in.decimal.degrees, 
           stnrec$Longitude.to.4.decimal.places.in.decimal.degrees,
           as.list(with(m.lognormgpd, c(lnmean, lnsd, u, sigmau, xi, phiu))),
           rsq,
           rms
           )
  names(row) <- colnames
  outdf <- rbind(outdf, row, stringsAsFactors = FALSE)

  xp <- seq(0, 200, 1)
  png(filename = paste(plotpath, figname, sep = "/"))

  hist(gust, 100, freq = FALSE, xlab = "Daily maximum wind gust [km/h]",
       main = sprintf("Daily maximum wind gust: %s", stnnum), col = "gray")

  try(with(m.lognormgpd,
           lines(xp,
                 dlognormgpd(xp, lnmean, lnsd, u, sigmau, xi),
                 col = "red", lwd = 2)))

  try(abline(v = m.lognormgpd$u, col = "red", lty = 2))
  dev.off()

  rlfigname <- paste(stnnum, ".rl.png", sep = "")
  png(filename = paste(plotpath, rlfigname, sep = "/"))
  try(rlplot(m.lognormgpd))
  dev.off()
}

write.csv(outdf, paste(plotpath, "params.csv", sep = "/"),
          row.names = FALSE, quote = FALSE, sep = ",")

png(filename = paste(plotpath, "params.png", sep = "/"),
    width = 12, height = 12, units = "in", res = 300)
pairs.panels(outdf[, 4:9],
             method = "pearson",
             density = TRUE,
             ellipses = TRUE)
dev.off()
