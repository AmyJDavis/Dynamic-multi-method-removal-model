#################################################################################
###
### Processing removal data for use with a dynamic multi-method removal model
###
### January 1, 2022
### Amy J. Davis
###
#################################################################################


### Libraries
library(reshape2)
library(splines)
source("remfuncs.R")
source('RemMultiMethLambdaMCMC.R')

### Read in data
## Should have 7 columns
## Date = date of removal event
## Longitude = longitude of removal event
## Latitude = latitude of removal event
## Take = Number of animals removed 
## Effort = Effort of this removal event
## Methods = Removal method for removal event
## Site = Site number, provide a unique number for each site (area to you want to estimate abundance for) within your study area, this may be only one site
rdata =read.csv("DummyData.csv")

###
### Data set up
###
rdata$Date=as.POSIXct(rdata$Date,format="%m/%d/%Y %H:%M")

### Create a column for the primary period of interest starting at the beginning of the study, must be a discrete number of months, standard is one month
rdata$PPNum=primaryperiodsequence(PPmonths = 2,rdata$Date)

### Make methods as factors and order these based on which would occur earlier in the secondary period (usually day) 
rdata$Methods=factor(rdata$Methods,levels=c("Trap","Aerial","Shooting"),ordered = TRUE)

### Order the data by site then by date, then by method
rdata=rdata[order(rdata$Site,rdata$Date,rdata$Methods),]

### Set the secondary period (SPdays = number of days for secondary period)
rdata$SPNum=secondaryperiod(SPdays=1,rdate=rdata$Date)

### Upload or create an area data frame giving the area for each site in km2, allow for a buffer if removal efforts on the edge of the study area
ardf=data.frame(Site=1:3,Area=c(142.88,127.07,107.17))

### Set up the area of impact for each method, and method names
detdf=data.frame(Method=c("Trap","Aerial","Shooting","0"),Area=c(15,40,5,0))
capnames=detdf$Method[-dim(detdf)[1]]


###
### Processing the data for input in the MCMC model
###
mptake=dcast(data = rdata[,c("Site","PPNum","Take","Methods","SPNum")],Site+PPNum~SPNum+Methods,value.var = "Take",fun.aggregate = sum)
mpeff=dcast(data = rdata[,c("Site","PPNum","Effort","Methods","SPNum")],Site+PPNum~SPNum+Methods,value.var = "Effort",fun.aggregate = sum)

passes=max(apply(mptake,1,function(x)length(which(x>0))))-2
mths=length(unique(mptake$MYNum))
sites=max(rdata$Site)

### Really slow way of processing the data but it works
takeinfo=matrix(0,dim(mptake)[1],passes)
effortinfo=matrix(0,dim(mptake)[1],passes)
metheff=matrix("0",dim(mptake)[1],passes)

for(t in 1:dim(mptake)[1]){
  tkind=which(mptake[t,]>0)[-(1:2)]
  takeinfo[t,seq(1,length(tkind),1)]=as.numeric(mptake[t,tkind])
  effortinfo[t,seq(1,length(tkind),1)]=as.numeric(mpeff[t,tkind])
  metheff[t,seq(1,length(tkind),1)]=gsub(".*_","",names(mptake)[tkind])
}
ymat=cbind(mptake[,1:2],takeinfo)
gbe=(effortinfo)
gtypes=metheff
areatot=ardf$Area
sml=dim(ymat)[1]

### Use the function to calculate the area of impact, this will work with and without lat/long coordinates
Areaper=area.of.impact.finder2(rdata,gtypes,ymat,ardf,GPS=TRUE,detdf,plotareas=FALSE)

### Create a design matrix for the growth rate analysis, Here I used an intercept (necessary) and basis function expansion on time to allow for growth rate to vary with time
Xdat=cbind(1,(ymat$PPNum)/max(ymat$PPNum))
Xd=cbind(1,bs(Xdat[,2],5))

### Set the model inputs
n.tune=500
p.tune=0.15
lam.tune=0.08
sigla=1
siglb=1
sigBa=1
sigBb=1
sigbeta=1
sigPa=1
sigPb=1
pstart=c(0.05,0.2,0.05)
n.mcmc=1000


remres=MultNom.removal.multmeth.lambda.mcmc(ymat,gbe,gtypes,capnames,Areaper,Xd,pstart,sigla,siglb,sigBa,sigBb,
                                            sigPa,sigPb,n.tune,lam.tune,p.tune,n.mcmc)





### 
### Examining the posterior results
###
n.burn=n.mcmc/2

### Posterior abundances
pop.mean=apply(remres$nt.save[n.burn:n.mcmc,],2,mean)
pop.quants=apply(remres$nt.save[n.burn:n.mcmc,],2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=TRUE))

### Posterior growth rates
lamd.mean=apply(remres$lambda.save[n.burn:n.mcmc,],2,mean)
lamd.quants=apply(remres$lambda.save[n.burn:n.mcmc,],2,function(x)quantile(x,c(0.025,0.975)))


### Posterior removal rates
if(length(pstart)==1){
  p.mean=mean(remres$p.save[n.burn:n.mcmc])
  p.quants=quantile(remres$p.save[n.burn:n.mcmc],probs=c(0.025,0.975),na.rm=TRUE)
}else{
  p.mean=colMeans(remres$p.save[n.burn:n.mcmc,])
  p.quants=apply(remres$p.save[n.burn:n.mcmc,],2,function(x)quantile(x,probs=c(0.025,0.975),na.rm=TRUE))
  
}

### Posterior beta estimates
if(dim(Xd)[2]>1){
  beta.mean=colMeans(remres$beta.save[n.burn:n.mcmc,])
  beta.quants=apply(remres$beta.save[n.burn:n.mcmc,],2,function(x)quantile(x,c(0.025,0.975)))
}else{
  beta.mean=mean(remres$beta.save[n.burn:n.mcmc])
  beta.quants=quantile(remres$beta.save[n.burn:n.mcmc],c(0.025,0.975))
}



