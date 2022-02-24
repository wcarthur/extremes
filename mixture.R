library('evmix')
library('stringr')
datapath = "X:/georisk/HaRIA_B_Wind/data/raw/from_bom/1-minute/wind"
filelist = list.files(datapath)
for (file in filelist){
  stnnum = str_extract(file, "\\d{6}")
  print(sprintf("Processing station %s", stnnum))
  filepath = paste(datapath, file, sep='/')
  df = read.csv(filepath, sep=',', header=TRUE)
  gust = df$windgust[df$windgustq != 'N']
  gust = gust[gust > 0]
  m.lognormgpd = flognormgpdcon(gust, useq=qs, fixedu=TRUE)
  m.gammagpd = fgammagpdcon(gust, useq=qs, fixedu=TRUE)
  m.kdengpdcon = fkdengpdcon(gust, useq=qs, fixedu=TRUE, pvector=pinit)
  m.weibullgpdcon = fweibullgpdcon(gust, useq=qs, fixedu=TRUE)
  xp = seq(0, 200, 1)
  hist(gust, 100, freq=FALSE, xlab="Daily maximum wind gust [km/h]", 
       main=sprintf("Daily maximum wind gust: %d", station))
  with(m.lognormgpd, lines(xp, dlognormgpdcon(xp, lnmean, lnsd, u, sigmau, xi), col='blue'))
  with(m.gammagpd, lines(xp, dgammagpdcon(xp, gshape, gscale, u, xi), col='red'))
  with(m.weibullgpdcon, lines(xp, dweibullgpdcon(xp, wshape, wscale, u, xi), col='brown'))
  legend("topright", c('LN-GPD', 'G-GPD', 'W-GPD'),
         lty = c(rep(1,3),2), col = c('blue', 'red', 'green')
         )
}

# Filter to QC'd observations only
gust = df$windgust[df$windgustq != 'N']

qs = quantile(gust, seq(0.8, 1, 0.01))
pinit = c(1, quantile(gust, 0.9))
# Fit a set of distributions
m.lognormgpd = flognormgpdcon(gust, useq=qs, fixedu=TRUE)
m.gammagpd = fgammagpdcon(gust, useq=qs, fixedu=TRUE)
m.kdengpdcon = fkdengpdcon(gust, useq=qs, fixedu=TRUE, pvector=pinit)
m.weibullgpdcon = fweibullgpdcon(gust, useq=qs, fixedu=TRUE)


hist(gust, 100, freq=FALSE, xlab="Daily maximum wind gust [km/h]", 
     main=sprintf("Daily maximum wind gust: %d", station))
with(m.lognormgpd, lines(xp, dlognormgpdcon(xp, lnmean, lnsd, u, sigmau, xi), col='blue'))
with(m.gammagpd, lines(xp, dgammagpdcon(xp, gshape, gscale, u, xi), col='red'))
#with(m.kdengpdcon, lines(xp, dkdengpdcon(xp, lambda, u, sigmau, xi, bw), col='green'))
with(m.weibullgpdcon, lines(xp, dweibullgpdcon(xp, wshape, wscale, u, xi), col='brown'))
legend("topright", c('LN-GPD', 'G-GPD', 'K-GPD', 'W-GPD'),
       lty = c(rep(1,3),2),
       col=c('blue', 'red', 'green', 'brown'))
