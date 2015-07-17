#Title: Function "RP_axis"
#Author: Sanabria A., augusto.sanabria@ga.gov.au
#CreationDate: 2007-09-17
#Description: This function generates the x-axis for
#return period plots. The user indicates the starting
#point, the number of points to be plotted and the
#axis resolution (either "low" or "high"). 
#"low" = steps of 10, "high" = steps of 5.

#Reference: 'A Statistical Model of Severe Winds'
#SeeAlso:
#Version:
#Modified by:
#ModifiedDate:
#Modification:

#usage:
RP_xaxis = function(first_RP,no_of_RP,resolution="low"){

#first_RP = user given 1st point in RP vector
#no_of_RP = number of elements in RP vector
#resolution = one of "low" or "high"
#For low: RP = 10,100,1000, etc.
#For high: RP = 10,50,100,500, etc.


  retp1 <- log10(first_RP)
  retp_num <- retp1 + no_of_RP - 1
  if(first_RP == 1)retp_num <- retp_num + 1

  RP1 <- 10^seq(retp1,retp_num)
  if(resolution == "high"){RP2 <- 5*10^seq(retp1,retp_num+2)
    xx <- append(RP1,RP2)
    ord <- order(xx,decreasing=F)
    xaxis <- xx[ord]
  }else{
    xaxis <- RP1
  }  
  
  xaxis[1:retp_num]
  
 
}