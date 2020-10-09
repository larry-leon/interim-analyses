



get.analyses.cr<-function(data1,
           data0,
           sim,
           dtaus=c(7,14),
           tau.interim,
           tau.final,
           timing.scale=1,
           cr_event.labels = c("censor", "outcome", "death"),
           quant,
           show.simcount = FALSE,
           direction = "GT",
           show.km = TRUE,
           get.final = FALSE,
           show.medians = FALSE,
           cr.analysis = "FG",
           titleit.IA = c("KM Interim"),
           titleit.FA = c("KM Final"),
           fit1 = NULL,
           fit2 = NULL,
           show.interim.obs = FALSE) {
    
  # tau.interim=tau.final --> label KM plots as final look
  if(tau.interim==tau.final) titleit.IA<-titleit.FA
  
  df<-data1
  Entry<-df$Entry
  Follow<-df$FollowUp
  C<-df$Censoring
  Time<-df$TrueSurvival 
  D<-df$DropOut
  Cause<-df$Cause1
  Event<-df$Event
  
  F.tau<-tau.interim-Entry  
  C.tau<-pmin(D,F.tau) 
  
  Eta.tau<-rep(0,length(Time))
  Eta.tau[which(Time<=C.tau & Cause==1)]<-1
  Eta.tau[which(Time<=C.tau & Cause==0)]<-2
  Event.tau<-ifelse(Eta.tau==1,1,0)  
  
  Event.tau.CR<-factor(Eta.tau,0:2,labels=cr_event.labels)
  
  df1 <-
    data.frame(
      time.cr = pmin(Time, C.tau),
      event.cr = Event.tau.CR,
      treat = 1,
      F.tau = F.tau,
      C.tau = C.tau,
      T.tau = Time,
      Event.tau = Event.tau,
      Y.tau = pmin(Time, C.tau)
    )
  
  df1<-subset(df1,F.tau>0)
  
  
  df<-data0
  Entry<-df$Entry
  Follow<-df$FollowUp
  C<-df$Censoring
  Time<-df$TrueSurvival 
  D<-df$DropOut
  Cause<-df$Cause1
  Event<-df$Event
  
  F.tau<-tau.interim-Entry  
  C.tau<-pmin(D,F.tau) 
  
  Eta.tau<-rep(0,length(Time))
  Eta.tau[which(Time<=C.tau & Cause==1)]<-1
  Eta.tau[which(Time<=C.tau & Cause==0)]<-2
  Event.tau<-ifelse(Eta.tau==1,1,0)  
  
  Event.tau.CR<-factor(Eta.tau,0:2,labels=cr_event.labels)
  
  df0 <-
    data.frame(
      time.cr = pmin(Time, C.tau),
      event.cr = Event.tau.CR,
      treat = 0,
      F.tau = F.tau,
      C.tau = C.tau,
      T.tau = Time,
      Event.tau = Event.tau,
      Y.tau = pmin(Time, C.tau)
    )
  df0<-subset(df0,F.tau>0)
  
  if(show.simcount){
  cat("Experimental arm event types","\n")
  print(table(df1$event.cr))
  
  cat("Control arm event types","\n")
  print(table(df0$event.cr))
    }
  
  df.cr<-as.data.frame(rbind(df1,df0))

  Y1.tau<-df1$Y.tau; D.1<-df1$D.tau; F1.tau<-df1$F.tau; T1.tau<-df1$T.tau; Event1.tau<-df1$Event.tau
  Y0.tau<-df0$Y.tau; D.0<-df0$D.tau; F0.tau<-df0$F.tau; T0.tau<-df0$T.tau; Event0.tau<-df0$Event.tau
  
  N1.tau<-length(Y1.tau)
  N0.tau<-length(Y0.tau)
  
  LFU0.tau<-mean(ifelse(D.0[which(F0.tau>0)]<=T0.tau,1,0))
  LFU1.tau<-mean(ifelse(D.1[which(F1.tau>0)]<=T1.tau,1,0))
  LFU.tau<-mean(ifelse(c(D.0[which(F0.tau>0)],D.1[which(F1.tau>0)])<=c(T0.tau,T1.tau),1,0))
  
  Ytau<-c(Y1.tau,Y0.tau)
  Dtau<-c(Event1.tau,Event0.tau)
  Ztau<-c(rep(1,N1.tau),rep(0,N0.tau))  # Treatment indicator
 
  dth.tau1<-dth.tau2<-dth.tau3<-NULL
  # Fine-Gray analysis
  if(cr.analysis=="FG"){
  # Fine-Gray dataset for cause=outcome
  df.primary<-finegray(Surv(time.cr,event.cr)~.,data=df.cr,etype="outcome")
  fit.outcome<-survfit(Surv(fgstart,fgstop,fgstatus)~treat,data=df.primary,weights=fgwt)
  cox.outcome<-coxph(Surv(fgstart,fgstop,fgstatus)~treat,data=df.primary,weights=fgwt)
  
  
  #df.death<-finegray(Surv(time.cr,event.cr)~.,data=df.cr,etype="death")
  #fit.death<-survfit(Surv(fgstart,fgstop,fgstatus)~treat,data=df.death,weights=fgwt)
  
  # For deaths, censor recoveries at last followup
  # Assume that if recovered, then did not die within 29 days
  
  last.time<-max(Ytau)
  df.cr$time.death<-df.cr$time.cr
  # Define event.FGs where cause-1 are now censored
  # at last follow-up
  df.cr$event.death<-ifelse(df.cr$event.cr=="death",1,0)
  # Set cause-1 event times to last.time
  df.cr$time.death[which(df.cr$event.cr==cr_event.labels[2])]<-c(last.time)
  fit.death<-survfit(Surv(time.death,event.death)~treat,data=df.cr)
  
  # Get difference in mortality estimates at days 7, 14, and 29
  # Day 7 death CIF differences
  CIFd.tau<-1-summary(fit.death,c(dtaus[1]))$surv
  dth.tau1<-CIFd.tau[2]-CIFd.tau[1]
  
  CIFd.tau<-1-summary(fit.death,c(dtaus[2]))$surv
  dth.tau2<-CIFd.tau[2]-CIFd.tau[1]
  
  # Raw proportions
  dth.tau3<-c(with(subset(df.cr,treat==1),mean(event.death))-with(subset(df.cr,treat==0),mean(event.death)))
  
  }
  
  ################################################################
  # FG-simple represents censoring "cause-2" events (eg., deaths)
  # at last follow-up.  We call this simple in the sense 
  # that cause-2 events remain in the risk-set across follow-up
  ################################################################
  if(cr.analysis=="FG-simple"){
  last.time<-max(Ytau)
  df.cr$time.FGs<-df.cr$time.cr
  # Define event.FGs where cause-2 are now censored
  # at last follow-up
  df.cr$event.FGs<-ifelse(df.cr$event.cr=="outcome",1,0)
  # Set cause-2 event times to last.time
  df.cr$time.FGs[which(df.cr$event.cr==cr_event.labels[3])]<-c(last.time)

  fit.outcome<-survfit(Surv(time.FGs,event.FGs)~treat,data=df.cr)
  cox.outcome<-coxph(Surv(time.FGs,event.FGs)~treat,data=df.cr)
    }
  
  table.outcome<-summary(fit.outcome)$table
  m0.interim<-table.outcome[1,7]
  m1.interim<-table.outcome[2,7]
  
  if(show.simcount) cat("Interim (Control,Exp) quantiles=",c(m0.interim,m1.interim),"\n")
  
  cox.summary<-summary(cox.outcome)
  cox.est<-cox.summary$coefficients
  
  bhat<-cox.est[1]
  se.bhat<-cox.est[3]
  
  cox.score<-cox.summary$sctest
  
  Z.lr<-(sign(bhat))*sqrt(cox.score[1])
  
  pval.lr<-1-pnorm(Z.lr)
  
  mu.diff<-table.outcome[2,5]-table.outcome[1,5]
  se.diff<-sqrt(table.outcome[2,6]^2+table.outcome[1,6]^2)
  
  # Set LTFU to represent deaths
  LFU0.tau<-mean(ifelse(df0$event.cr==cr_event.labels[3],1,0))
  LFU1.tau<-mean(ifelse(df1$event.cr==cr_event.labels[3],1,0))
  
  # Interim analyses in terms of timing.scale (tau.interim/timing.scale)
  # E.g., if time is in terms of days then tau.interim/7 represents weeks
  
  out.interim.sim<-c(sim,tau.interim/timing.scale,sum(Dtau),Z.lr,bhat,se.bhat,m1.interim,m0.interim,LFU1.tau,LFU0.tau,mu.diff,se.diff,N1=sum(Ztau),N0=sum(1-Ztau),N=length(Ztau),dth.tau1,dth.tau2,dth.tau3)
  
  
  if(show.km){
    
    nE<-sum(Dtau)
    
    plot(fit.outcome,ylim=c(0.0,1.0),conf.int=F,col=c("red","blue"),lwd=c(1,1),xlab="Time",mark.time=TRUE,fun='event')
    
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
      hr.ci<-as.matrix(round(cox.summary$conf.int,4))
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
  out.final.sim<-NULL
  
  out<-list(out.final.sim=out.final.sim,out.interim.sim=out.interim.sim,df.cr=df.cr)
  return(out)
}





