# Function to calculate the shape parameter to fit
# a GPD distribution to given data. This parameter
# is calculated using function "fpot" of package 'evd'
# Input:
# u = threshold
# ws = clean vector of wind data
# yrs = No. of years in dataset
# quant = requested quantile (RP in years)
# Returns:
# zh = return value for the given quantile (input variable quant)
# shape = shape parameter for the fitted GPD

calc_qsh = function(u,ws,yrs,RP){

   nopy <- length(ws)/yrs

   #options(warn = -2)    #Supress all warning messages
   #print(paste("Trying to estimate GPD parameters with threshold=",toString(thres)," and RP=",toString(RP)))
   qnt <- try(fpot(ws,u,mper=RP,npp=365.25,std.err=FALSE),silent=TRUE)
   #print(qnt[[1]])
   if(!is.numeric(qnt[[1]])){
        #print(qnt[[1]])
        #print(qnt[[1]])    #Trap errors
        zh <- NA
        shape <- NA
    }else{
        zh <- qnt$estimate[1]
        shape <- qnt$estimate[2]
    }
    cbind(zh,shape)
}


