#Title: Function "GPD_distr.r"
#Author: Sanabria A., augusto.sanabria@ga.gov.au
#CreationDate: 2006-08-30
#Description: Function to calc. return periods
#using a GPD distribution. Fitting is done using maximum likelihood - ML.
#The given variables are entered to the funct. via the argument.
#The funct. does plot the RP curve and returns the corresponding
#calculated return periods.
#The user can select a low resolution x-axis (steps of 10) or a
#high resolution x-axis (steps of 5) years.

#Reference: Package "evd"
#SeeAlso: Coles' book: An Introduction to Statistical
#Modeling of Extreme Values.
#Version:

#Required:
library("evd")
#'sdir' is defined globally:
source(paste(sdir,"RP_xaxis.r",sep="")  )


#Usage:
GPD_distr = function(ws,yrs,thres,resolution="low",title="",sub1="",diag = FALSE,calcCI=FALSE){
  
  #ws = clean vector of wind speeds
  #yrs = number of years in dataset (range)
  #thres = Threshold for this data
  #title = main title in plot
  #sub1 = subtitle under x-axis
  #resolution = one of "low","high" or "std" see below.
  #if diag = FALSE generate only GPD plot, otherwise generate diagnostic plots too
  
  #Generate vector of x-axis:
  RP <- RP_xaxis(10,4,"low")        #Return periods, low resolution (normal)
  if(resolution == "high"){RP <- RP_xaxis(10,8,"high")
  }else if(resolution == "std"){RP <- c(5,10,20,25,50,100,200,250,500,1000,2000,2500,5000,10000)
  }else if(resolution == "nwra"){RP <- c(50,100,200,500,1000,20000)
  }else if(resolution == ""){print("**Fatal error: RP resolution not specified!") }
  
  zh <- c()                         #corresponding quantiles
  up <- c()
  lo <- c()
  nopy <- length(ws)/yrs            #number of obs per year
  lastp <- length(RP)
  for(i in 1:lastp) {
    #'fpot' of package evd calculates the parameters of the GPD for a given
    #dateset using ML, it also returns the quantile corresponding to given RP
    q = fpot(ws,thres,mper=RP[i],npp=nopy,std.err=F)
    if(calcCI){
      # Evaluate confidence intervals
      pp = profile(q,which="rlevel",mesh=c(0.05),xmin=c(0))
      ci = confint(pp,parm="rlevel")
      #up[i] = q$estimate[1]+1.96*q$std.err[1]
      #lo[i] = q$estimate[1]-1.96*q$std.err[1]
      up[i] = ci[3]
      lo[i] = ci[1]
    }
    zh[i] = q$estimate[1]                #quantile (speed)
    #Just to check results:
    if(RP[i] == 1000)q1000 <- zh[i]      #quantile at 1000 yrs
    if(RP[i] == 10000)q10000 <- zh[i]    #quantile at 10000 yrs
  }
  
  shape = toString(q$estimate[2],11)
  pars = c(shape)
  #plot.new()
  #Plot wind speeds for given RP, display threshold (top) and shape parameter (bottom).
  #plot(RP,zh,log="x",type="l",ylab="speed (m/s)",xlab="Return Period (yrs)",main=title,sub=sub1 )
  ident <- toString(thres)
  #text(1000,0.98*zh[3],ident)
  #text(RP[2],zh[lastp],"Return Per. via a GPD distr. ")
  #text(2.5*RP[1],zh[lastp-3],"             Parameter sh")
  #text(RP[3],zh[lastp-3],pars[1])
  #grid(NULL, NULL, lwd = 2)
  
  #print(paste("shape =  ",shape,sep="")  )
  #print(paste("scale = ",q$scale,sep="")  )
  
  if(diag == TRUE)plot(q)      #plot diagnostic plots
  
  #To store parameters in external table:
  #source("/mnt/store/Winddata/store_GPDpar.r")
  #store_GPDpar(ws,thres,q$estimate[2],q1000,q10000,-1)
  
  
  #Return results (RP and corresponding quantiles) to do multiple plots:
  cbind(RP,zh,up,lo)
  
}
