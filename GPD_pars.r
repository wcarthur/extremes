#Function to find the parameters of a GPD which
#fits a given dataset using package 'evd'
#(GPD fitting via ML)

#Usage:
GPD_pars = function(ws,yrs,thresh){
  
  #ws = clean vector of wind speed
  #yrs = range of years in dataset
  #thresh = threshold to fit the GPD (calculated
  #using function 'sel_approp_u.r')
  
  library("evd")
  
  pars <- fpot(ws,thresh,npp=length(ws)/yrs,std.err = FALSE)
  
  excess <- ws[ws > thresh]
  crate <- length(excess)/length(ws)
  #print(paste("u = ",thresh," mx_exc =",max(excess)," mx_ws = ",max(ws) )  )
  #print(paste("l_exc =",length(excess),"l_all = ",length(ws))  )
  scale <- pars$estimate[1]
  shape <- pars$estimate[2]
  #  cbind(scale,shape,crate)
  cbind(thresh,scale,shape,length(excess))  
}
