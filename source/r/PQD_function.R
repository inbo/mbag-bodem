#  PQ(p,o) prediction quality routine
#  input (p,o) =>  p= predicted value, o = observed value  
#  prediction error = PE = p-o    
#  MPE = Mean prediction error BIAS +/- = MPE? (estimate accuracy)
#  SDPE = SD of prediction error = random variation = SDPE? = estimate precision
#  MSPE = mean squared prediction error (total prediction quality)
#  RMSPE = root mean square prediction error (total prediction quality)
#  Bias= 100*(MPE²)/MSPE
#  Precision= 100*(SDPE²)/MSPE
#  Bruno De Vos - 10/10/2017
### Refined index of Agreement Willmott et al. 2012, 2015

## inputs: p = predicted values / o= observed values, dec= decimals to round

PQD<-function(p,o,dec) 
{  
  dat<-data.frame(na.omit(cbind(p,o)))
  p<-dat[,1]
  o<-dat[,2]
  omean<-mean(o,na.rm=TRUE)
  n<-length(p)
  PE<-p-o 
  APE<-abs(PE)
  AO<-abs(o-omean)
  # see Willmott et al 2015   (wi=1)
  MSE<-(1/n)*sum(PE^2,na.rm=TRUE) # mean squared error
  MAE<-(1/n)*sum(APE,na.rm=TRUE) # mean absolute error = MPE
  MADV<-(1/n)*sum(AO,na.rm=TRUE) #mean absolute deviation
  MPE<-(1/n)*(sum(PE,na.rm=TRUE))
  SDPE<-sqrt(var(PE,na.rm=TRUE))
  MSPE<-(1/n)*(sum(PE^2,na.rm=TRUE))
  RMSPE<-sqrt(MSPE)
  SDSPE<-sqrt(var(PE^2,na.rm=TRUE))
  r<-cor(p,o)
  trueness<-100*((MPE^2)/MSPE)
  precision<-100*((SDPE^2)/MSPE)
  
  MPE<-round(MPE,dec)
  SDPE<-round(SDPE,dec)
  MSPE<-round(MSPE,dec)
  RMSPE<-round(RMSPE,dec)	
  SDSPE<-round(SDSPE,dec)
  COR<-round(r,dec)
  COR2<-round(r^2,dec)
  ERRTRUE<-floor(trueness)
  ERRPREC<-floor(precision)
  
  MSE<-round(MSE,dec)
  MAE<-round(MAE,dec)
  MADV<-round(MADV,dec)
 
  # index of agreement d 
  # The Index of Agreement (d) developed by Willmott (1981) as a standardized measure of the degree of model prediction error and varies between 0 and 1.
  # A value of 1 indicates a perfect match, and 0 indicates no agreement at all (Willmott, 1981). 
  #
  # d = 1 - [ ( sum( (obs - sim)^2 ) ] / sum( ( abs(sim - mean(obs)) + abs(obs - mean(obs)) )^2 ) 
  d<-1-((sum((p-o)^2))/sum((abs(p-mean(o))+abs(o-mean(o)))^2))
  d<-round(d,3)  
  
  # refined index of agreement dr (Willmott et al. 2012,2015)
  # refined index is dimensionless
  # 
  # 2 steps:
  # if MAE <= 2*MAD THEN dr=1-(MAE/(2*MAD))
  # if MAE > 2*MAD THEN dr=((2*MAD)/MAE)-1 
  dr<-ifelse(MAE<=(2*MADV),(1-(MAE/(2*MADV))),((2*MADV)/MAE)-1)
  dr<-round(dr,3)  
  
  print("---- Prediction quality output (p-o)--d decimals--")
  print(paste("n=",n))
  print(paste("MPE=", MPE))
  print(paste("SDPE=",SDPE))
  print(paste("MSPE=",MSPE))
  print(paste("RMSPE=",RMSPE ))
  print(paste("SDSPE=",SDSPE ))
  print(paste("COR=",COR))
  print(paste("R2=",COR2))
  print(paste("trueness (% of total error)=",ERRTRUE))
  print(paste("precision(% of total error)=",ERRPREC))
  print(paste("index of agreement (0-1)=",d))
  print(paste("refined index of agreement (-1 to 1)=",dr))
  print(paste("MSE=", MSE))
  print(paste("MAE=",MAE))
  print(paste("MAD=",MADV))
  print("--------------------------------------------------")
  
  out<-cbind(n, MPE, SDPE, MSPE, RMSPE,COR, COR2,ERRTRUE,ERRPREC,d,dr,MSE,MAE,MADV)
  return(out)
}
