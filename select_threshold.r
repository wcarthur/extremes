#Title: Function 'sel_approp_u.r'
#Author: Sanabria A., augusto.sanabria@ga.gov.au
#CreationDate: 2006-09-15
#Description: Function to select appropriate threshold to fit
#a GPD distribution to given data.
#Reference:
#Statistical Model of Severe Winds
#SeeAlso: Coles' book.

#Version: 1.1
#Modified by: AS
#ModifiedDate:2006-12-13
#Modification: The parameter epsilon was reduced to -0.01 to
#allow the GPD's shape parameter get closer to 0 to produce
#a higher curve.
#Notice: shape = 0 determines an unbounded (Gumbel)
#distribution, inappropriate to model wind speed,
#a naturally bounded phenomenon.

#Usage:
select_threshold = function(data,yrs,init=0,resol=-1,last=-1){
  
  #ws = a clean Vector of wind speed
  #yrs = number of years covered in the dataset (range)
  #resol = resolution of u (default = 0.25)
  #init = first possible threshold
  #last = last possible threshold
  #If init= -1 and resol = -1, program calculates appropriate
  #start and resolution for process.
  
  #Required:
  source(paste(sdir,"GPD_distr.r",sep=""))
  source(paste(sdir,"my_Rplot.r",sep=""))
  source(paste(sdir,"calc_qsh.r",sep=""))
  source(paste(sdir,"closest_u.r",sep=""))
  
  #eps = -0.03        #min value of shape
  eps <- -0.01
  mx_ws <- max(data)
  count <- 0
  mx_iter <- 50000    #max. steps allowed
  initial_pr <- 1     #Initial value of retp-axis in GPD
  
  u_vec <- c()
  noveru_vec <- c()
  shape_vec <- c()
  mxpnt_vec <- c()
  q1000_vec <- c()
  q10000_vec <- c()
  qdiff_vec <- c()
  
  Retp_vec <- vector("list",mx_iter)    #Vector of lists to store dataframes
  
  start <- ceiling(0.5 + 0.1*max(data) )
  if(init > 0)start <- init     #use user-given "start"
  steps <- 1
  if(resol > 0)steps <- resol   #use user-given "resol"
  u <- start
  
  for (i in 1:mx_iter){
    u <- u + steps
    thresh <- format(u,digits=4)
    noveru <- length(data[data > u])
    #print(paste("N over u = ",toString(noveru),sep=""))
    
    if(noveru < 10) {
      break
    }
    if(last > 0 & u > last){
      #print(paste("Last > 0 & u > last", toString(u), sep=""))
      break
    }
    #Check if this configuration is feasible (use 1000 years RP)
    quan_sh1 <- calc_qsh(u,data,yrs,1000)
    if(!is.numeric(quan_sh1)){
      print("**Error message: threshold ignored, try next one!" )
      next
    }
    q1000 <- quan_sh1[1]
    sh1 <-quan_sh1[2]
    #print(paste("quan_sh1 = ",toString(quan_sh1),sep=""))
    mxpnt1 <- -(mean(data) - sd(data) )/sh1
    #Check if this configuration is feasible (use 10000 years RP)
    quan_sh2 <- calc_qsh(u,data,yrs,10000)
    #print(paste("quan_sh2 = ",toString(quan_sh2),sep=""))
    if(!is.numeric(quan_sh2)){
      print("**Error message: threshold ignored, try next one!" )
      next
    }
    
    q10000 <- quan_sh2[1]
    sh2 <-quan_sh2[2]
    #print(cbind(count," ",thresh," ",noveru," ",sh2)  )
    
    mxpnt2 <- -(mean(ws) - sd(ws) )/sh2
    #'mxpnt1' and 2 should be the same
    #print(cbind(u,sh1,sh2) )
    qdiff <- abs(q10000 - q1000)
    #print(paste("qdiff = ",toString(qdiff),sep=""))
    feasible <- !is.na(q1000) & !is.na(q10000) & sh1 < eps  & sh2 < eps #& qdiff < 0.2*q10000
    if(feasible){
      #print(paste("Feasible (stored): u = ", toString(thresh), sep="")  )
      count <- count + 1
      Retp <- GPD_distr(data, yrs, u, resolution="low", " ", " ", diag=F)
      u_vec[count] <- u
      noveru_vec[count] <- noveru
      shape_vec[count] <- sh1
      mxpnt_vec[count] <- mxpnt1
      q1000_vec[count] <- q1000
      q10000_vec[count] <- q10000
      qdiff_vec[count] <- qdiff
      #Store 2-D object in a vector of lists:
      Retp_vec[[count]] <- Retp
    }else{ 
      #print(paste("Not feasible (eliminated): u = ",toString(thresh),sep="") )
    }
  }
  #print(paste("q1000_vec = ",toString(q1000_vec),sep=""))
  q1000_vec <- q1000_vec[!is.na(q1000_vec)]
  q10000_vec <- q10000_vec[!is.na(q10000_vec)]
  
  if (length(q10000_vec) == 0) {
    return (cbind(c(0), c(0), c(0)))
  }
  
  u_vec <- u_vec[!is.na(u_vec)]
  
  Av1000 <- sum(q1000_vec)/length(q1000_vec)
  Av10000 <- sum(q10000_vec)/length(q10000_vec)
  
  
  #Av1000 <- ceiling(Av1000 + 0.15*Av1000)
  Av1000 <- ceiling(Av1000 + 0.05*Av1000)
  #Av10000 <- ceiling(Av10000 + 0.15*Av10000)
  Av10000 <- ceiling(Av10000 + 0.05*Av10000)
  
  u1000 <- closest_u(u_vec,q1000_vec,Av1000)
  u10000 <- closest_u(u_vec,q10000_vec,Av10000)
  
  
  tit1 <- ""
  #tit1 <- try(titulo,silent=T)                #Title for plot (defined via global allocation in shell)
  filename1 <- paste(outdir, "GPD_table.txt", sep="")
  header1 <- "Summary of GPD fitting for wind speed data"
  header2 <- paste(header1, tit1, sep=" ")
  header3 <- paste(header1, header2, sep=" ")
  write.table(header3, file=filename1, append=F, sep="  ",row.names=F,col.names=F,quote=F)
  tit1 <- "    u         No>u       shape      mx_pnt     q1000     q10000      qdiff"
  write.table(tit1,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  
  for(i in 1:count){
    #Format list to write table of GPD parameters:
    valores <- cbind(u_vec[i],
                     noveru_vec[i],
                     shape_vec[i],
                     mxpnt_vec[i],
                     q1000_vec[i],
                     q10000_vec[i],
                     qdiff_vec[i] )
    val <- format(valores,nsmall=2,digits=2,width=9)
    write.table(val,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  }
  
  tit2 <-  cbind("    max ws =",toString(mx_ws)  )
  write.table(tit2,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  min_sh <- cbind("    acceptable shape < ",toString(eps) )
  write.table(min_sh,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  tit4 <- cbind("    Av1000"," approp u","  Av10000"," approp u")
  write.table(tit4,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  ans <- cbind("   ", Av1000,"   ",u1000,"        ",Av10000,"   ",u10000)
  write.table(ans,file=filename1,append=T,sep="  ",row.names=F,col.names=F,quote=F)
  #print(paste("Summary table 'GPD_table.txt' was written to",outdir,sep=" ")  )
  #
  #Plot all feasible RP curves with different colors:
  #
  #yrange <- c(35,45)      #For Alice Springs temp.
  #yrange <- c(10,70)        #For Wind speed (default)
  #yrange <- c(50,150)        #For Wind speed (default)
  #my_Rplot(1,10000,yrange[1],yrange[2],Retp_vec[[1]],t1=tit1,t2=" ",2)
  #typ <- 1
  #inbox <- c()
  #inbox[1] <- toString(format(u_vec[1],nsmall=2,digits=2)  )
  #for(i in 2:count){
  #     thres <- u_vec[i]
  #     inbox[i] <- toString(format(thres,nsmall=2,digits=2) )
  #     lines(Retp_vec[[i]],lty=i,col=i)
  #}
  #legend(1,60,legend=inbox,lty=1:count,col=1:count)
  #legend(1,130,legend=inbox,lty=1:count,col=1:count)
  #cbind(u_vec[which(q10000_vec == max(q10000_vec))],max(q10000_vec)   )
  
  if(u1000 > u10000){
    sh_mx <- closest_u(shape_vec,q1000_vec,Av1000)
    #print(paste("Returned =",u1000,"Av1000 =",Av1000,"sh = ",sh_mx,sep="  " )  )
    cbind(u1000,Av1000,sh_mx)
  }else{
    sh_mx <- closest_u(shape_vec,q10000_vec,Av10000)
    #print(paste("Returned =",u10000,"Av10000 =",Av10000,"sh = ",sh_mx,sep="  " )  )
    cbind(u10000,Av10000,sh_mx)  }
}
