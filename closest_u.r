#Function to find the closes value to a given parameter

closest_u = function(uu, qvec, qgiven){
  small = 10
  ikeep = 1
  for(i in 1:length(qvec) ){
    if(abs(qgiven - qvec[i]) <= small){
      small <- abs(qgiven-qvec[i])
      ikeep = i 
    }
  }
  uu[ikeep]

}
