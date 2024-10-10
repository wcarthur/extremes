# This one fits a gamma/GPD distribution
library('tea')
library('extRemes')
library('ismev')
library('stringr')
library('psych')

datapath <- "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/2019/Daily"
plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/daily_max_wind/wind/seqthreshold"
filelist <- list.files(datapath, pattern = "*.txt")
colnames <- c("stnnum", "u", "sigmau", "xi", "n", "p")
outdf <- data.frame(stnnum = character(),
                    u = double(),
                    sigmau = double(),
                    xi = double(),
                    n = double(),
                    p = double(),
                    stringsAsFactors = FALSE)
names(outdf) <- colnames

for (file in filelist) {
  SKIP <- FALSE
  stnnum <- str_extract(file, "\\d{6}")
  print(sprintf("Processing station %s", stnnum))
  filepath <- paste(datapath, file, sep='/')
  figname <- paste(stnnum, '.png', sep="")
  df <- read.csv(filepath, sep=',', header=TRUE)
  gust <- df$Speed.of.maximum.wind.gust.in.m.s[!is.na(df$Speed.of.maximum.wind.gust.in.m.s)
                                               & df$Quality.of.maximum.gust.speed == 'Y']
  gust <- gust[gust > 0]
  if (length(gust) < 365 * 5){
    print("Insufficient observations")
    next
  }
  qs <- seq(quantile(gust, 0.5), quantile(gust, .995), length.out = 100)

  result <- tryCatch(TH(gust, qs), error = function(e) {SKIP = TRUE})

  if (typeof(result) == "logical") {
    print(paste("Error fitting for station", stnnum))
    next
  } else {
    Adf <- data.frame(tail(result[result[, 'est.shape'] < 0,], n=1))
  
  }
  #print(m.gammagpd$se)
  row <- c(stnnum, as.list(with(Adf,
                                c(threshold, est.scale, est.shape, num.above, p.values))))
  names(row) <- colnames
  outdf <- rbind(outdf, row, stringsAsFactors = FALSE)
  m = gpd.fit(gust, Adf$threshold)
  png(filename = paste(plotpath, figname, sep = "/"))
  try(gpd.diag(m))
  dev.off()

  rlfigname <- paste(stnnum, ".rl.png", sep = "")
  png(filename = paste(plotpath, rlfigname, sep = "/"))
  try(rlplot(m))
  dev.off()
}

write.csv(outdf, paste(plotpath, "params.csv", sep = "/"),
          row.names = FALSE, quote = FALSE, sep = ",")

