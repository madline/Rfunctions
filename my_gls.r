# gls function for trend and confidence intervals from one timeseries

# input: y=annual-timeseries, time=time in years
# e.g.
#y=erg_ytsurf_B
#time <- year(model_ydate)
# period <- seq(1908,2008,1)

#y=Rsst_y_B
#time=year(r_ydate)
#period=period
#level=.9
#CItype=2

my_gls <- function(y,time,period=seq(1850,2008,1),level=.9,CItype,alpha){ 
  
  
  # alpha determination
  
  id_intvar <- which(time %in% seq(1850,1899,1))
  if (time[1] > 1851 | length((y)<50)) {id_intvar <- 1:50} # falls kurze Periode (Reynolds) -> alpha von Hadley
  #if (is.na(mean(y))) {id_intvar <- 1:(min(which(is.na(y)))-1)}
  #if (length(id_intvar)<10) {id_intvar <- which(y %in% as.numeric(na.contiguous(y)) & period < 1915)}
  if (!is.na(mean(y)) & length(y)>50) {alphaa <- abs(acf(y[id_intvar])$acf[2]) }
  if (is.na(mean(y))) {alphaa <- alpha}
  if (length((y)<50)) {alphaa <- abs(acf(y[hadsst_B])$acf[2]) }
  
  y_period <- y[which(time %in% period)]
  timesteps <- seq(1,length(y_period),1)
  fit <- gls(y_period~timesteps,correlation=corAR1(alphaa,fixed=TRUE),method="ML",na.action=na.omit)
  
  dectrend <- round(fit$coef[[2]]*10,digits=4)
  
  if (CItype==1) { 
  CI <- round((confint(fit,level=level)[2,2]-confint(fit,level=level)[2,1])*10/2,digits=4) 
  return(list("trend"=dectrend,"CI"=CI)) }
  
  if (CItype==2) {CI <-  as.numeric(round(confint(fit,level=level)[2,]*10,digits=4)) 
                  return(list("trend"=dectrend,"CImin"=CI[1],"CImax"=CI[2]))}
  
  if (CItype==3) {pval <-  summary(fit)$tTable[2,4] 
                  return(list("trend"=dectrend,"pval"=pval))} #pvalue for maps
  
  if (CItype==4) { if (is.na(mean(y))) {fitted <- as.numeric();fitted[which(!is.na(y))] <- fitted(fit)}
                   if (!is.na(mean(y))) {fitted <- fitted(fit)}
                   dectrend <- round(fit$coef[[2]]*10,digits=3)
                   CI <- round((confint(fit,level=level)[2,2]-confint(fit,level=level)[2,1])*10/2,digits=3)
                   list("trend"=dectrend,"CI"=CI,"alpha"=alphaa,'slope'=fitted)}

}
