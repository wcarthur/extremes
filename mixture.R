library('evmix')
library('stringr')
library('psych')

datapath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind"
plotpath <- "X:/georisk/HaRIA_B_Wind/data/derived/obs/1-minute/wind/gammagpdcon"
filelist <- list.files(datapath, pattern = "*.txt")
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
  stnnum <- str_extract(file, "\\d{6}")
  print(sprintf("Processing station %s", stnnum))
  filepath <- paste(datapath, file, sep='/')
  figname <- paste(stnnum, '.png', sep="")
  df <- read.csv(filepath, sep=',', header=TRUE)
  gust <- df$windgust[df$windgustq == 'Y']
  gust <- gust[gust > 0]
  if (length(gust) < 365 * 5){
    print("Insufficient observations")
    next
  }
  qs <- quantile(gust, seq(0.9, 1, 0.01))

  result <- tryCatch(fgammagpdcon(gust, useq = qs, fixedu = FALSE),
                     error = function(e) {SKIP = TRUE})
  if (typeof(result) == "logical") {
    print(paste("Error fitting for station", stnnum))
    next
  } else {
    m.gammagpd <- result
  }
  row <- c(stnnum, as.list(with(m.gammagpd,
                                c(gshape, gscale, u, sigmau, xi, phiu))))
  names(row) <- colnames
  outdf <- rbind(outdf, row, stringsAsFactors = FALSE)

  xp <- seq(0, 200, 1)
  png(filename = paste(plotpath, figname, sep = "/"))

  hist(gust, 100, freq = FALSE, xlab = "Daily maximum wind gust [km/h]",
       main = sprintf("Daily maximum wind gust: %s", stnnum), col = "gray")

  try(with(m.gammagpd,
           lines(xp,
                 dgammagpdcon(xp, gshape, gscale, u, xi),
                 col = "red", lwd = 2)))

  try(abline(v = m.gammagpd$u, col = "red", lty = 2))

  dev.off()

  rlfigname <- paste(stnnum, ".rl.png", sep = "")
  png(filename = paste(plotpath, rlfigname, sep = "/"))
  try(rlplot(m.gammagpd))
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
