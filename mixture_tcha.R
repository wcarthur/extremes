library('evmix')
library('stringr')
library('psych')


fitStn <- function (filename, datapath) {
  stnnum <- str_extract(file, "\\d+")
  filepath <- paste(datapath, filename, sep="/")
  df <- read.csv(filepath, sep=",", header=TRUE)
  gust <- df$wspd
  gust <- gust[gust > 0]
  qs <- quantile(gust, seq(0.95, 1, 0.01))
  result <- tryCatch(fgammagpd(gust, useq = qs, fixedu = TRUE),
                     error = function(e) {SKIP = TRUE})
  if (typeof(result) == "logical") {
    print(paste("Error fitting for station", stnnum))
    next
  } else {
    m.gammagpd <- result
  }
  
  xp <- seq(0, 200, 1)
  
  figname <- paste(stnnum, '.png', sep="")
  png(filename = paste(datapath, figname, sep = "/"))
  
  hist(gust, 100, freq = FALSE, xlab = "Daily maximum wind gust [m/s]",
       main = sprintf("Daily maximum wind gust: %s", stnnum), col = "gray")
  
  try(with(m.gammagpd,
           lines(xp,
                 dgammagpd(xp, gshape, gscale, u, sigmau, xi),
                 col = "red", lwd = 2)))
  
  try(abline(v = m.gammagpd$u, col = "red", lty = 2))
  dev.off()
  
  rlfigname <- paste(stnnum, ".rl.png", sep = "")
  png(filename = paste(datapath, rlfigname, sep = "/"))
  try(rlplot(m.gammagpd, upperfocus = TRUE))
  dev.off()
  
  qqfigname <- paste(stnnum, ".qq.png", sep = "")
  png(filename = paste(datapath, qqfigname, sep = "/"))
  try(qplot(m.gammagpd, upperfocus = TRUE))
  dev.off()
  
}




datapath <- "X:/georisk/HaRIA_B_Wind/projects/qfes_swha/data/derived/hazard/profiles/HISTORICAL_1981-2010/plots/fit/stationdata"
plotpath <- "X:/georisk/HaRIA_B_Wind/projects/qfes_swha/data/derived/hazard/profiles/HISTORICAL_1981-2010/plots/fit/stationdata"
filelist <- list.files(datapath, pattern = "*.csv")
colnames <- c("stnnum", "gshape", "gscale", "u", "sigmau", "xi", "phiu")
outdf <- data.frame(stnnum = character(),
                    gshape = double(),
                    gscale = double(),
                    u = double(),
                    sigmau = double(),
                    xi = double(),
                    phiu = double(),
                    stringsAsFactors = FALSE)
names(outdf) <- colnames

for (file in filelist) {
  SKIP <- FALSE
  stnnum <- str_extract(file, "\\d+")
  print(sprintf("Processing station %s", stnnum))
  filepath <- paste(datapath, file, sep='/')
  figname <- paste(stnnum, '.png', sep="")
  df <- read.csv(filepath, sep=',', header=TRUE)
  #gust <- df$Speed.of.maximum.wind.gust.in.m.s[!is.na(df$Speed.of.maximum.wind.gust.in.m.s)
  #                                             & df$Quality.of.maximum.gust.speed == 'Y']
  gust <- df$wspd
  gust <- gust[gust > 0]
  #if (length(gust) < 365 * 5){
  #  print("Insufficient observations")
  #  next
  #}
  qs <- quantile(gust, seq(0.95, 1, 0.01))

  result <- tryCatch(fgammagpd(gust, useq = qs, fixedu = TRUE),
                     error = function(e) {SKIP = TRUE})
  if (typeof(result) == "logical") {
    print(paste("Error fitting for station", stnnum))
    next
  } else {
    m.gammagpd <- result
  }
  print(m.gammagpd$mle)
  print(m.gammagpd$se)
  row <- c(stnnum, as.list(with(m.gammagpd,
                                c(gshape, gscale, u, sigmau, xi, phiu))))
  names(row) <- colnames
  outdf <- rbind(outdf, row, stringsAsFactors = FALSE)

  xp <- seq(0, 200, 1)
  
  png(filename = paste(plotpath, figname, sep = "/"))

  hist(gust, 100, freq = FALSE, xlab = "Daily maximum wind gust [m/s]",
       main = sprintf("Daily maximum wind gust: %s", stnnum), col = "gray")

  try(with(m.gammagpd,
           lines(xp,
                 dgammagpd(xp, gshape, gscale, u, sigmau, xi),
                 col = "red", lwd = 2)))

  try(abline(v = m.gammagpd$u, col = "red", lty = 2))
  dev.off()

  rlfigname <- paste(stnnum, ".rl.png", sep = "")
  png(filename = paste(plotpath, rlfigname, sep = "/"))
  try(rlplot(m.gammagpd, upperfocus = FALSE))
  dev.off()
  
  qqfigname <- paste(stnnum, ".qq.png", sep = "")
  png(filename = paste(plotpath, qqfigname, sep = "/"))
  try(qplot(m.gammagpd, upperfocus = FALSE))
  dev.off()
}

write.csv(outdf, paste(plotpath, "params.csv", sep = "/"),
          row.names = FALSE, quote = FALSE, sep = ",")

png(filename = paste(plotpath, "params.png", sep = "/"),
    width = 12, height = 12, units = "in", res = 300)
pairs.panels(outdf[, 2:7],
             method = "pearson",
             density = TRUE,
             ellipses = TRUE)
dev.off()
