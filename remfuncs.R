############################################################################
###
### Data processing functions for multi-method removal data
###
### January 7, 2022
### Amy J. Davis
###
############################################################################






### 
### Function to calculate the area impacted by removal efforts at each time and removal event
###

area.of.impact.finder<-function(rdata,gtypes,ymat,ardf,GPS=TRUE,detdf,plotareas=FALSE){
  require(sp)
  require(rgeos)
  rdata$GroupMth=paste(rdata$Site,rdata$MYNum,sep=".")
  rdata$GroupMthDay=paste(rdata$Site,rdata$MYNum,rdata$Day,sep="-")
  rdata$nm=as.numeric(factor(rdata$GroupMth,ordered=TRUE,levels=unique(rdata$GroupMth)))
  rdata$MthYear=format(rdata$Date,"%b-%Y")
  
  agb=matrix(detdf[match(gtypes,detdf$Method),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2])
  
  
  
  if(GPS==TRUE){
    
    ###
    ### Processing to get the proportion of the area impacted by removal pass
    ###
    ### Creating unique site, month, and event unifiers
    
    ### Getting the number of each event by group and month
    tfv=rdata[,c("GroupMthDay","GroupMth")]
    tfv=tfv[!duplicated(tfv),]
    tfv$nm=as.numeric(factor(tfv$GroupMth,ordered=TRUE,levels=unique(tfv$GroupMth)))
    tfv$Dayn=ave(tfv$nm,tfv$nm,FUN=seq_along)
    rdata$Dayn=tfv[match(rdata$GroupMthDay,tfv$GroupMthDay),"Dayn"]
    
    ardfgm=data.frame(SiteMth=0,Day=0,Dayn=0,nm=0,Area=0)
    
    ### Slow way to get an visualize area of impact by event
    for(k in unique(rdata$GroupMthDay)){
      kind=which(unique(rdata$GroupMthDay)==k)
      kdat=rdata[rdata$GroupMthDay==k,]
      spdf=SpatialPointsDataFrame(coords = rdata[rdata$GroupMthDay==k,c("Longitude","Latitude")],
                                  data = data.frame(1:dim(rdata[rdata$GroupMthDay==k,])[1]))
      sppoly=gConvexHull(spdf,byid = FALSE)
      spbuff=gBuffer(spgeom = sppoly,byid=FALSE,width=0.02)
      ardfgm[kind,"SiteMth"]=kdat$GroupMth[1]
      ardfgm[kind,"Day"]=kdat$Day[1]
      ardfgm[kind,"Dayn"]=kdat$Dayn[1]
      ardfgm[kind,"nm"]=kdat$nm[1]
      ardfgm[kind,"Area"]=gArea(spbuff)/0.00008
      
      if(plotareas==TRUE){
        par(mfrow=c(3,3))
        
        plot(rdata$Longitude[rdata$GroupMthDay==k],rdata$Latitude[rdata$GroupMthDay==k],pch=16,
             xlim=c(min(rdata$Longitude)*1.00019,max(rdata$Longitude)*0.9998),
             ylim=c(min(rdata$Latitude)*0.9992,max(rdata$Latitude)*1.0009),
             xlab="Longitude",ylab="Latitude")
        plot(spbuff,add=TRUE)
        title(paste("Group ",(rdata$MthYear)[which(unique(rdata$GroupMthDay)==k)]))
      }
      
    }
    agbe=agb
    agbe[cbind(ardfgm$nm,ardfgm$Dayn)]=ardfgm$Area
    agbe=agbe/agb
    
    gbea=gbe/agbe
    gbea[is.na(gbea)]=0
    
    Areai=as.matrix(agb)
    Areaper=matrix(pmin(1,Areai/rep(areatot,as.numeric(table(ymat[,1])))),dim(Areai)[1],dim(Areai)[2])
    Areaper[is.na(Areaper)]=0
    
    par(mfrow=c(1,1),mar=c(4,5,1,1)+0.1)
  }else if(GPS==FALSE){
    areatots=matrix(ardf[match(ymat$Site,ardf$Site),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2],byrow = FALSE)
    Areaper=agb^(log(areatots*((gbe*agb)/(areatots+(gbe*agb))))/log(agb))/areatots
    Areaper[is.na(Areaper)]=0
    
  }
  Areaper
  
}




###################
###
### Editing the area finder function to work with primary and secondary period coding instead of month and day coding
###
###################


### 
### Function to calculate the area impacted by removal efforts at each time and removal event
###

area.of.impact.finder2<-function(rdata,gtypes,ymat,ardf,GPS=TRUE,detdf,plotareas=FALSE){
  require(sp)
  require(rgeos)
  rdata$GroupMth=paste(rdata$Site,rdata$PPNum,sep=".")
  rdata$GroupMthDay=paste(rdata$Site,rdata$PPNum,rdata$SPNum,sep="-")
  rdata$nm=as.numeric(factor(rdata$GroupMth,ordered=TRUE,levels=unique(rdata$GroupMth)))
  rdata$MthYear=format(rdata$Date,"%b-%Y")
  
  agb=matrix(detdf[match(gtypes,detdf$Method),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2])
  
  
  
  if(GPS==TRUE){
    
    ###
    ### Processing to get the proportion of the area impacted by removal pass
    ###
    ### Creating unique site, month, and event unifiers
    
    ### Getting the number of each event by group and month
    tfv=rdata[,c("GroupMthDay","GroupMth")]
    tfv=tfv[!duplicated(tfv),]
    tfv$nm=as.numeric(factor(tfv$GroupMth,ordered=TRUE,levels=unique(tfv$GroupMth)))
    tfv$SPNumn=ave(tfv$nm,tfv$nm,FUN=seq_along)
    rdata$SPNumn=tfv[match(rdata$GroupMthDay,tfv$GroupMthDay),"SPNumn"]
    
    ardfgm=data.frame(SiteMth=0,SPNum=0,SPNumn=0,nm=0,Area=0)
    
    ### Slow way to get an visualize area of impact by event
    for(k in unique(rdata$GroupMthDay)){
      kind=which(unique(rdata$GroupMthDay)==k)
      kdat=rdata[rdata$GroupMthDay==k,]
      spdf=SpatialPointsDataFrame(coords = rdata[rdata$GroupMthDay==k,c("Longitude","Latitude")],
                                  data = data.frame(1:dim(rdata[rdata$GroupMthDay==k,])[1]))
      sppoly=gConvexHull(spdf,byid = FALSE)
      spbuff=gBuffer(spgeom = sppoly,byid=FALSE,width=0.02)
      ardfgm[kind,"SiteMth"]=kdat$GroupMth[1]
      ardfgm[kind,"SPNum"]=kdat$SPNum[1]
      ardfgm[kind,"SPNumn"]=kdat$SPNumn[1]
      ardfgm[kind,"nm"]=kdat$nm[1]
      ardfgm[kind,"Area"]=gArea(spbuff)/0.00008
      
      if(plotareas==TRUE){
        par(mfrow=c(3,3))
        
        plot(rdata$Longitude[rdata$GroupMthDay==k],rdata$Latitude[rdata$GroupMthDay==k],pch=16,
             xlim=c(min(rdata$Longitude)*1.00019,max(rdata$Longitude)*0.9998),
             ylim=c(min(rdata$Latitude)*0.9992,max(rdata$Latitude)*1.0009),
             xlab="Longitude",ylab="Latitude")
        plot(spbuff,add=TRUE)
        title(paste("Group ",(rdata$MthYear)[which(unique(rdata$GroupMthDay)==k)]))
      }
      
    }
    agbe=agb
    agbe[cbind(ardfgm$nm,ardfgm$SPNumn)]=ardfgm$Area
    agbe=agbe/agb
    
    gbea=gbe/agbe
    gbea[is.na(gbea)]=0
    
    Areai=as.matrix(agb)
    Areaper=matrix(pmin(1,Areai/rep(areatot,as.numeric(table(ymat[,1])))),dim(Areai)[1],dim(Areai)[2])
    Areaper[is.na(Areaper)]=0
    
    par(mfrow=c(1,1),mar=c(4,5,1,1)+0.1)
  }else if(GPS==FALSE){
    areatots=matrix(ardf[match(ymat$Site,ardf$Site),"Area"],length(unique(rdata$GroupMth)),dim(gtypes)[2],byrow = FALSE)
    Areaper=agb^(log(areatots*((gbe*agb)/(areatots+(gbe*agb))))/log(agb))/areatots
    Areaper[is.na(Areaper)]=0
    
  }
  Areaper
  
}

###################
###
### Creating a flexible primary period sequence based on a set number of months, default is 1 month
###
###################

primaryperiodsequence<-function(PPmonths=1,rdate=rdata$Date){
  ceiling(((as.numeric(format(rdate,"%m"))+(as.numeric(format(rdate,"%Y"))-min(as.numeric(format(rdate,"%Y"))))*12)-min(as.numeric(format(rdate,"%m"))+(as.numeric(format(rdate,"%Y"))-min(as.numeric(format(rdate,"%Y"))))*12)+1)/PPmonths)
}


###################
###
### Creating a flexible secondary period sequence based on a set number of days, default is 1 day
###
###################
secondaryperiod<-function(SPdays=1,rdate=rdata$Date){
  
  ad1=tapply(rdate,rdata$PPNum,function(x)ceiling(as.numeric(format(x,'%j'))/SPdays)) 
  ad2=lapply(ad1,function(x)as.numeric(factor(x,ordered=TRUE,levels=unique(x))))
  unlist(ad2)
}
