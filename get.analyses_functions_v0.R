
  

get.analyses<-function(sim,tau.interim,tau.final,Entry.1,D.1,T.1,Entry.0,D.0,T.0,quant,show.simcount=FALSE,direction="GT",
                       get.mean=FALSE,mean.draws=0,show.km=TRUE,get.final=FALSE,show.medians=FALSE,timing.scale=1,
                       titleit.IA=c("KM Interim"),titleit.FA=c("KM Final"),fit1=NULL,fit2=NULL,
                       show.interim.obs=FALSE){
  
  # tau.interim=tau.final --> final look
  if(tau.interim==tau.final) titleit.IA<-titleit.FA
  
  F1.tau<-tau.interim-Entry.1  
  C1.tau<-pmin(D.1,F1.tau)[which(F1.tau>0)] 
  T1.tau<-T.1[which(F1.tau>0)]
  Event1.tau<-ifelse(T1.tau<=C1.tau,1,0)
  Y1.tau<-pmin(T1.tau,C1.tau)
  
  F0.tau<-tau.interim-Entry.0  
  C0.tau<-pmin(D.0,F0.tau)[which(F0.tau>0)] 
  T0.tau<-T.0[which(F0.tau>0)]
  Event0.tau<-ifelse(T0.tau<=C0.tau,1,0) 
  Y0.tau<-pmin(T0.tau,C0.tau)
  
  N1.tau<-length(Y1.tau)
  N0.tau<-length(Y0.tau)
  
  LFU0.tau<-mean(ifelse(D.0[which(F0.tau>0)]<=T0.tau,1,0))
  LFU1.tau<-mean(ifelse(D.1[which(F1.tau>0)]<=T1.tau,1,0))
  LFU.tau<-mean(ifelse(c(D.0[which(F0.tau>0)],D.1[which(F1.tau>0)])<=c(T0.tau,T1.tau),1,0))
  
  Ytau<-c(Y1.tau,Y0.tau)
  Dtau<-c(Event1.tau,Event0.tau)
  Ztau<-c(rep(1,N1.tau),rep(0,N0.tau))  # Treatment indicator
  
  # Control quantile
  temp<-NA.CHR.Weighted(time=Y0.tau,Delta=(Event0.tau==1))
  S0.hat<-temp$S.KM
  m0.interim<-suppressWarnings(min(temp$at.points[which(S0.hat<=quant)]))
  m0.interim<-ifelse(m0.interim<Inf,m0.interim,NA)
  
  # Experimental median
  temp<-NA.CHR.Weighted(time=Y1.tau,Delta=(Event1.tau==1))
  S1.hat<-temp$S.KM
  m1.interim<-suppressWarnings(min(temp$at.points[which(S1.hat<=quant)]))
  m1.interim<-ifelse(m1.interim<Inf,m1.interim,NA)
  
  if(show.simcount) cat("Interim (Control,Exp) quantiles=",c(m0.interim,m1.interim),"\n")
  
  #################################################################################
  #                             Difference in survival estimates
  #################################################################################
  
  oo<-order(Ytau)
  Ytau<-Ytau[oo]; Dtau<-Dtau[oo]; Ztau<-Ztau[oo]
  
  ################
  # Log-rank test
  ################
 
  # Log-rank return 1-sided superiority 
if(direction=="GT") temp<-Logrank(time=Ytau,Delta=Dtau,X=Ztau)
if(direction=="LT") temp<-Logrank(time=Ytau,Delta=Dtau,X=1-Ztau)
  
Z.lr<-temp$Z.lr

pval.lr<-1-pnorm(Z.lr)

  cox.fit<-summary(coxph(Surv(Ytau,Dtau)~Ztau))$coefficients
  bhat<-cox.fit[,1]
  se.bhat<-cox.fit[,3]
  
  #KM.MeanTrunc<-function(Y,Delta,X,L=NULL,draws=0,get.band=FALSE,taus=c(NULL,NULL),tau.seq=0.25,draws.band=0,details=FALSE,plotband=FALSE,alpha=0.025)
  # Truncated mean inference
  mu.diff<-NA
  se.diff<-NA
  if(get.mean){
    temp<-KM.MeanTrunc(Y=Ytau,Delta=Dtau,X=Ztau,draws=mean.draws)
    mu.diff<-temp$dhat
    se.diff<-sqrt(var(temp$dstar))
  }
  
  #colnames(out.interim)<-c("sim","tau.IA","d.IA","Zlr.IA","bhat.IA","se.bhat.IA","m1.IA","m0.IA","drop1.IA","drop0.IA","N1","N0","N")
  out.interim.sim<-c(sim,tau.interim/timing.scale,sum(Dtau),Z.lr,bhat,se.bhat,m1.interim,m0.interim,LFU1.tau,LFU0.tau,mu.diff,se.diff,N1=sum(Ztau),N0=sum(1-Ztau),N=length(Ztau))
  
  if(show.km){
    
    nE<-sum(Dtau)
    
    kmfit<-KM.plot.2sample.weighted(Y=Ytau,E=Dtau,Treat=Ztau,risk.set=TRUE,by.risk=3,tpoints.add=c(-1,0),
                                    xmin=0,xmax=max(Ytau),risk.cex=0.65,risk.add=max(Ytau),
                                    show.logrank=FALSE,show.med=FALSE,show.cox=FALSE,Xlab="Time")
    #title(main=titleit.IA)
    if(show.medians){
      m1<-round(m1.interim,1); m0<-round(m0.interim,1)
      legend("topright",bty = "n",paste("m0=", m0,
                                        ", m1=", m1))
    }
    rp<-vector('expression',1)
    rp[1]=substitute(expression(italic(InterimTime)== MYVALUE), 
                     list(MYVALUE = format(tau.interim,dig=2)))[2]
    legend('top', legend = rp, bty = 'n')              
    
    rp<-vector('expression',1)
    rp[1]=substitute(expression(italic(Events)== MYVALUE), 
                     list(MYVALUE = format(nE,dig=2)))[2]
    legend('top', legend = rp, bty = 'n', inset=0.1)              
    
    rp<-vector('expression',1)
    rp[1]=substitute(expression(italic(p)== MYVALUE), 
                     list(MYVALUE = format(pval.lr,dig=4)))[2]
    legend('top', legend = rp, bty = 'n', inset=0.2)              
    
    show.hr<-FALSE
    if(show.hr){
      hr.ci<-as.matrix(round(summary(coxph(Surv(Ytau,Dtau)~Ztau),conf.int=0.95)$conf.int,4))
      hr.ci<-hr.ci[1,c(1,3,4)]
      rp<-vector('expression',3)
      rp[1]=substitute(expression(italic(HR)== MYVALUE1), 
                       list(MYVALUE1 = format(hr.ci[1],digits = 2)))[2]
      rp[2]=substitute(expression(italic(Lower) == MYVALUE2), 
                       list(MYVALUE2 = format(hr.ci[2], digits = 2)))[2]
      rp[3]=substitute(expression(italic(Upper) == MYVALUE3), 
                       list(MYVALUE3= format(hr.ci[3], digits = 2)))[2]               
      legend('topright', legend = rp, bty = 'n')        
    }

    
    if(show.interim.obs){
      lines(fit1,conf.int=FALSE,type="s",lwd=3,col="grey",lty=1)
      lines(fit2,conf.int=FALSE,type="s",lwd=3,col="brown",lty=1)
    }
  }
  
  
  ####################################################
  # Final analysis?
  # In case want to look at final along with interim
  # But this is highly in-efficient to re-do at every 
  # interim.  Default is to skip.
  ####################################################
  
  # Defunct 
  out.final.sim<-NULL
  
  out<-list(out.final.sim=out.final.sim,out.interim.sim=out.interim.sim)
  return(out)
}
