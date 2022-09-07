###########################################################################
###
### Removal model using Multinomial 
### Growth rates are linked by study area
### Multiple removal types acceptable (e.g., trap, aerial, ground shooting)
### 
### Sept 26, 2019
### Amy J Davis
###
###########################################################################
MultNom.removal.multmeth.lambda.mcmc<-function(ymat,gbe,gtypes,capnames,Areaper,Xd,pstart,sigla=1,siglb=1,sigBa=1,sigBb=1,
                                               sigPa=1,sigPb=1,n.tune=1000,lam.tune=0.1,p.tune=0.1,n.mcmc){
  
  #### Data inputs
  ##  ymat= matrix where the first column is site ID num, column 2 is month num, and the rest
  ##        of columns are # animals removed by site (rows) and removal pass (columns)
  ##  gbe = the effort, a matrix giving a value for each site, month, and pass (sites*month by passes)
  ##  gtypes = a matrix (of characters I know not a real matrix) denoting the type of removal method being used at each site, month, and pass combination, same dimensions at gbe
  ##  capnames = a character vector naming each of the removal methods to be used
  ##  Areaper = matrix of the same dimensions as gbe denoting the proportion of the site impacted by removal activities at that site, month, and pass (this is calculated using location data an prior information on area of impact)
  ##  Xd   = design matrix of covariates for density, must include at least an intercept, the number of rows should be sites*months 
  ##  pstart = initial starting values for removal rates, a vector of length equal to the number of removal methods (same length as gtypes)
  ##
  
  #### Setting hyperpriors
  ## sigla = the alpha value for the prior distribution on lambda
  ## siglb = the beta value for the prior distribution on lambda
  ## sigBb = the variance value for the prior distribution on betas assuming sigma*I (no correlation structure)
  ## sigPa = the alpha value for the prior distribution on p
  ## sigPb = the beta value for the prior distribution on p
  
  #### Setting inputs and tuning paramters
  ## n.tune = the tuning parameter for the Metropolis-Hastings step for abundance
  ## lam.tune = the tuning parameter for the Metropolis-Hastings step for lambda
  ## p.tune = the tuning parameter for the Metropolis-Hastings step for p
  ## n.mcmc = The number of Marcov Chain Monte Carlo iterations the user would like to use
  
  
  
  ###
  ### Libraries
  ###
  library(MASS)
  library(gplots)
  library(boot)
  library(matlib)
  library(expm)
  library(MCMCpack)
  library(RMKdiscrete)
  
  ###
  ###  Data processing
  ###
  sites=max(ymat[,1])
  mths=max(ymat[,2])
  b.num=dim(Xd)[2]
  passes=dim(gbe)[2]
  sml=dim(ymat)[1]
  
  p=pstart
  pvals=length(p)
  ptype=matrix(p[match(gtypes,capnames)],sml,passes)
  theta=1-(1-ptype)^gbe
  pip=theta*Areaper*cbind(1,t(exp(apply(log((1-Areaper)+Areaper*(1-theta)),1,cumsum)))[,-passes])
  
  
  sigbeta=rinvgamma(1,sigBa,sigBb)
  sigl=rinvgamma(1,sigla,siglb)
  
  mup=rnorm(1,0,1)
  sigP=rinvgamma(1,sigPa,sigPb)
  
  ## normalize the X values
  ymatcum=t(apply(ymat[,-(1:2)],1,cumsum))  
  betas=rnorm(b.num,0,sigbeta)
  
  ## Start lambda as stable
  lamd=rep(1,sml)
  
  
  ntstart=data.frame(Site=1:sites,FirstTime=tapply(ymat[,2],ymat[,1],function(x)min(x)))
  nt=pmax(ymatcum[,passes],rowSums(ymat),rpois(sml,500))
  
  
  ntmat=nt-cbind(0,ymatcum[,-passes])
  ynmat=cbind(ymat[,-(1:2)],nt-rowSums(ymat[,-(1:2)]))
  pipn=cbind(pip,1-rowSums(pip))
  
  ###
  ### Set up matrices to save iterations from the MCMC sample
  ###
  
  nt.save=matrix(0,n.mcmc,sml)   ### Population estimates by site
  lamd.save=matrix(0,n.mcmc,sml)  ### Growth rate by month and not different by site... so far
  psave=matrix(0,n.mcmc,pvals)     ### Population level capture rates
  betaDsave=matrix(0,n.mcmc,b.num) ### Beta values for property level density
  sigBsave=rep(0,n.mcmc)
  siglsave=rep(0,n.mcmc)
  mupsave=rep(0,n.mcmc)
  sigPsave=rep(0,n.mcmc)
  n.burn=round(n.mcmc/2)           ### Create a burn in portion
  
  
  
  ### 
  ### MCMC loop
  ###
  
  for(l in 1:n.mcmc){
    if(l%%1000==0)cat(l," ");flush.console()
    
    
    
    #############################
    ###
    ### Sample from n
    ###
    #############################
    pipn=cbind(pip,1-rowSums(pip))  
    
    #n.star=rpois(sites,nt)
    n.star=rnegbin(sml,mu = nt,theta = n.tune)
    ynmat.star=cbind(ymat[,-(1:2)],n.star-rowSums(ymat[,-(1:2)]))
    
    ### Only update those without negative numbers
    ntmat.star=n.star-cbind(0,ymatcum[,-passes])
    dup=apply(cbind(ynmat.star,ntmat.star),1,function(x)all(x>(-1)))
    dupind=which(dup==TRUE)
    if(length(dupind)>0){
      dupindno1=dupind[!ymat[dupind,2]==1]
      dupindnolast=dupind[!ymat[dupind,2]==mths]
      nMHratio=exp((sum(sapply(dupind,function(x)dmultinom(ynmat.star[x,],prob=pipn[x,],log=TRUE)))+
                      rowSums(cbind(dpois(n.star[dupindno1],lamd[dupindno1-1]*(nt[(dupindno1-1)]-ymatcum[dupindno1-1,passes]),log=TRUE)[match(dupind,dupindno1)],
                                    dpois(nt[dupindnolast+1],lamd[dupindnolast]*(n.star[(dupindnolast)]-ymatcum[dupindnolast,passes]),log=TRUE)[match(dupind,dupindnolast)],
                                    dnegbin(nt[dupind],mu = n.star[dupind],n.tune,log=TRUE)),na.rm = TRUE))-
                     (sum(sapply(dupind,function(x)dmultinom(ynmat[x,],prob=pipn[x,],log=TRUE)))+
                        rowSums(cbind(dpois(nt[dupindno1],lamd[dupindno1-1]*(nt[(dupindno1-1)]-ymatcum[dupindno1-1,passes]),log=TRUE)[match(dupind,dupindno1)],
                                      dpois(nt[dupindnolast+1],lamd[dupindnolast]*(nt[(dupindnolast)]-ymatcum[dupindnolast,passes]),log=TRUE)[match(dupind,dupindnolast)],
                                      dnegbin(n.star[dupind],mu = nt[dupind],n.tune,log=TRUE)),na.rm = TRUE)))
      nMHratio[is.nan(nMHratio)]=0
      tmp.keep=nMHratio>runif(length(dupind))
      nt[dupind][tmp.keep]=n.star[dupind][tmp.keep]
      
      ntmat=nt-cbind(0,ymatcum[,-passes])
      ynmat=cbind(ymat[,-(1:2)],nt-rowSums(ymat[,-(1:2)]))
      
    }
    
    
    
    ### 
    ### Sample from site level growth rate
    ###
    
    lamd.star=exp(rnorm(sml,log(lamd),lam.tune))
    lmhratio=exp((c(dpois(nt[which(ymat[,2]!=1)],lamd.star[which(ymat[,2]!=mths)]*(nt[which(ymat[,2]!=mths)]-ymatcum[which(ymat[,2]!=mths),passes]),log=TRUE)[match(1:sml,which(ymat[,2]!=1))])+
                    dnorm(log(lamd.star),(Xd%*%betas),sigl,log=TRUE))-
                   (c(dpois(nt[which(ymat[,2]!=1)],lamd[which(ymat[,2]!=mths)]*(nt[which(ymat[,2]!=mths)]-ymatcum[which(ymat[,2]!=mths),passes]),log=TRUE)[match(1:sml,which(ymat[,2]!=1))])+
                      dnorm(log(lamd),(Xd%*%betas),sigl,log=TRUE)))
    
    lmhratio[is.na(lmhratio)]=0
    tmp.keep=lmhratio>runif(length(mths-1))
    lamd[tmp.keep]=lamd.star[tmp.keep]
    
    
    
    ###
    ### Sample from betas for growth rate
    ###
    # Cholsky decomposition version
    tmp.1=t(Xd)%*%solve(sigl*diag(sml))
    tmp.chol=chol(tmp.1%*%Xd+solve(sigbeta*diag(b.num)))
    tmp.b=tmp.1%*%(log(lamd))
    betas=backsolve(tmp.chol,backsolve(tmp.chol,tmp.b,transpose = TRUE)+rnorm(b.num))
    
    
    
    ###
    ###  Sample from p for each site, pass, and time
    ###
    p.star=inv.logit(rnorm(pvals,logit(p),p.tune))
    ptype.star=matrix(p.star[match(gtypes,capnames)],sml,passes)
    theta.star=1-(1-ptype.star)^gbe
    pip.star=theta.star*Areaper*cbind(1,t(exp(apply(log((1-Areaper)+Areaper*(1-theta.star)),1,cumsum)))[,-passes])
    pipn.star=cbind(pip.star,1-rowSums(pip.star))
    
    
    pMHratio=exp((sum(sapply(1:(sml),function(x)dmultinom(ynmat[x,],prob=pipn.star[x,],log=TRUE)))+
                    dnorm(logit(p.star),mup,sigP))-
                   (sum(sapply(1:(sml),function(x)dmultinom(ynmat[x,],prob=pipn[x,],log=TRUE)))+
                      dnorm(logit(p),mup,sigP)))
    
    tmp.keep=pMHratio>runif(1)
    tmp.keep=ifelse(is.na(tmp.keep),"FALSE",tmp.keep)
    p[tmp.keep]=p.star[tmp.keep]
    
    ptype=matrix(p[match(gtypes,capnames)],sml,passes)
    theta=1-(1-ptype)^gbe
    pip=theta*Areaper*cbind(1,t(exp(apply(log((1-Areaper)+Areaper*(1-theta)),1,cumsum)))[,-passes])
    
    
    ###
    ### Sample from sigma for Density
    ###
    sigla2=sigla+((sml)/2)
    siglb2=siglb+(sum((log(lamd)-Xd%*%betas)^2))/2
    
    sigl=rinvgamma(1,sigla2,siglb2)
    
    
    ###
    ### Sample from sigma for betas
    ###
    sigba2=sigBa+(b.num/2)
    sigbb2=sigBb+sum((betas)^2)/2
    
    sigbeta=rinvgamma(1,sigba2,sigbb2)
    
    ###
    ### Sample for mean of capture rate
    ###
    mu.1=(sigP/(1+sigP))*(logit(p)/sigP)
    sig.1=(sigP/(1+sigP))
    
    mup=rnorm(1,mu.1,sig.1)
    
    ###
    ### Sample from sigma for capture rate
    ###
    sigpa2=sigPa+pvals/2
    sigpb2=sigPb+sum((logit(p)-mup)^2)/2
    
    sigP=rinvgamma(1,sigpa2,sigpb2)
    
    ###
    ### Save Samples
    ###
    nt.save[l,]=nt
    lamd.save[l,]=lamd
    psave[l,]=p
    betaDsave[l,]=betas
    sigBsave[l]=sigbeta
    siglsave[l]=sigl
    mupsave[l]=mup
    sigPsave[l]=sigP
    
    
  }
  cat("\n")
  
  ###
  ### Make Plots
  ###
  
  ### Plots for population size by site/time
  par(mar=c(2,2,1,1))
  plotsites=min(6,dim(nt.save)[2])
  layout(matrix(c(seq(1,2*min(6,plotsites),by=2),seq(2,2*min(6,plotsites),by=2),seq(2,2*min(6,plotsites),by=2)),min(6,plotsites),3))
  for(i in 1:plotsites){
    plot(nt.save[,i],type="l",main=paste("Trace: Pop Size site ",i),ylab=paste("Nt ",i),xlab="MCMC Iteration")
    plot(density(nt.save[n.burn:n.mcmc,i]),lwd=2,
         main=paste("Posterior and Prior: Nt ",i))
    #curve(dbinom(x,ymatcum[i,passes],pmean),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  ### Plots for population growth rate by site/time
  par(mar=c(2,2,1,1))
  layout(matrix(c(seq(1,2*min(6,plotsites),by=2),seq(2,2*min(6,plotsites),by=2),seq(2,2*min(6,plotsites),by=2)),min(6,plotsites),3))
  for(i in 1:plotsites){
    plot(lamd.save[,i],type="l",main=paste("Trace: Pop growth by time ",i),ylab=paste("lamd ",i),xlab="MCMC Iteration")
    plot(density(lamd.save[n.burn:n.mcmc,i]),lwd=2,
         main=paste("Posterior and Prior: lamd ",i))
    #curve(dbinom(x,ymatcum[i,passes],pmean),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  
  ###
  ### Plots for capture rate (p) by site
  ###
  par(mar=c(2,2,1,1))
  if(pvals==1){
    layout(matrix(c(1,2,2),1,3))
    plot(psave,type="l",main="Trace: p",ylab="Capture rate",xlab="MCMC Iteration")
    plot(density(psave[n.burn:n.mcmc]),lwd=2,main="Posterior for p")
  }else{
    layout(matrix(c(seq(1,2*min(6,pvals),by=2),seq(2,2*min(6,pvals),by=2),seq(2,2*min(6,pvals),by=2)),min(6,pvals),3))
    for(i in 1:pvals){
      plot(psave[,i],type="l",main=paste("Trace: p -",capnames[i]),ylab="Capture rate",xlab="MCMC Iteration")
      plot(density(psave[n.burn:n.mcmc,i]),lwd=2,main=paste("Posterior for p-",capnames[i]))
    }
  }
  
  
  
  ###
  ###  Trace plots for beta values on growth rate
  ###
  par(mar=c(2,2,1,1))
  layout(matrix(c(seq(1,2*min(6,b.num),by=2),seq(2,2*min(6,b.num),by=2),seq(2,2*min(6,b.num),by=2)),min(6,b.num),3))
  for(i in 1:b.num){
    plot(betaDsave[,i],type="l",main=paste("Trace: Beta ",i),ylab=paste("Beta ",i),xlab="MCMC Iteration")
    plot(density(betaDsave[n.burn:n.mcmc,i]),lwd=2,
         main=paste("Posterior and Prior: Beta ",i))
    curve(dnorm(x,0,1),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  
  ### 
  ### Trace plots for sigma values and mean
  ###
  layout(matrix(c(1,2,3,4,5,6,7,8,5,6,7,8),4,3))
  
  plot(sigBsave,type="l",main="Trace: sigma for beta",ylab="Sigma for beta",xlab="MCMC iteration")
  plot(siglsave,type="l",main="Trace: sigma for growth rate",ylab="sigl",xlab="MCMC iteration")
  plot(mupsave,type="l",main="Trace: mean for capture rate",ylab="mu for p",xlab="MCMC iteration")
  plot(sigPsave,type="l",main="Trace: sigma for capture rate",ylab="sig for p",xlab="MCMC iteration")
  
  plot(density(sigBsave[n.burn:n.mcmc]),lwd=2,main="Posterior for sigma for beta")
  plot(density(siglsave[n.burn:n.mcmc]),lwd=2,main="Posterior for sigma for growth rate")
  plot(density(mupsave[n.burn:n.mcmc]),lwd=2,main="Posterior for mean for p")
  plot(density(sigPsave[n.burn:n.mcmc]),lwd=2,main="Posterior for sigma for p")
  
  
  ###
  ### calculate and print posterior means
  ###
  
  pop.mean=apply(nt.save[n.burn:n.mcmc,],2,mean)
  pop.quants=apply(nt.save[n.burn:n.mcmc,],2,function(x) quantile(x,probs=c(0.025,0.975),na.rm=TRUE))
  
  par(mfrow=c(1,1),mar=c(4,5,1,1)+0.1)
  dmax=max(pop.quants)
  dmin=min(pop.quants)
  plot(ymat[,2],pop.mean,xlab="Time",ylab="Study Area Population Size",
       pch=16,ylim=c(dmin,dmax),col=ymat[,1])
  
  
  ####
  ####  Write Output 
  ####
  
  list(nt.save=nt.save,lambda.save=lamd.save,p.save=psave,beta.save=betaDsave)
  
  
  
}
