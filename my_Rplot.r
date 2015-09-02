
#Title: Function "my_Rplot.r"
#Author: Sanabria A., augusto.sanabria@ga.gov.au
#CreationDate: 2006-06-18
#Description: This function sets up the graph to
#draw multiple plots in it. It also draws the
#first plot given in argument 'Rt'. User can plot subsequent
#plots by using command "lines()".

#x1,x2 = limits for x-axis
#y1,y2 = limits for y-axis
#Rt = 2-D vector of return periods (x,retp)
#if kk == 1 it's a GEV distr. otherwise it's a GPD

#Reference: Coles' book: An Introduction to Statistical 
#Modeling of Extreme Values.
#SeeAlso: 'A Statistical Model of Severe Winds'
#Version:
#Modified by:
#ModifiedDate:
#Modification:

#required ('sdir' must be defined globally):
source(paste(sdir,"RP_xaxis.r",sep="")  )

#Usage:
my_Rplot = function(x1=1,x2=10000,y1=20,y2=60,Rt,t1=titulo,t2="",kk=2){
#my_Rplot = function(x1=50,x2=2000,y1=20,y2=50,Rt,t1="",t2="",kk=2){

#Rt = RP to be plotted
#t1 = main title
#t2 = subtitle
#kk = 2 for GPD or 1 for GEV distributions

 yvec <- seq(y1,y2,by = 5)
 ytxt = c()
 #Generate vector of x-axis:
 resolution <- "high"
 RP <- RP_xaxis(1,5,resolution)        #Return periods, low resolution (normal)
 if(resolution == "high")RP <- RP_xaxis(1,9,"high")

 #Generate string for y-axis:
 for(i in 1:length(yvec) )ytxt[i] <- c(toString(yvec[i]) )
 
 #Generate string for x-axis:
 xtxt = c()
 for(i in 1:length(RP) )xtxt[i] <- c(toString(RP[i]) )
 
 #postscript("hob_comb_obs_only.ps")
 plot.new()
 #Plot a dummy array to force the axis to take given values:
 dummy=cbind(c(x1,x2),c(y1,y2),ncol=2,nrow=2)
 t1 <- try(titulo,silent=T)      #Default title (defined in shell via a global allocation)

 plot(dummy,log="x",ylab="Gust wind speed (m/s)",xlab="Return Period (yrs)",main=t1,sub=t2,axes=F,type="n")
 #plot(dummy,log="x",ylab="Flood (m)",xlab="Return Period (yrs)",main=t1,sub=t2,axes=F,type="n")
 #Type="n", an option to not plot anything
 axis(2,at=yvec,labels=ytxt  )
 axis(1,at=RP,labels=xtxt)
 if(kk == 1)
   text(RP[3],1.5*max(yvec)-0.1*yvec[1],"Return Period using a GEV distr.")
 else
   text(300,max(yvec)-0.2*yvec[1],"Return Period using a GPD distribution")
 grid(NULL, NULL, lwd = 2)
 inbox <- c("Combined","Thunderstorms","Synoptic")
 #legend(RP[5],yvec[3],legend=inbox,lty=1:3,col=1:3)
 #Plot given RP (Rt):
 lines(Rt,col=1)
 #dev.off()
}
