
# For competing risk data, dgp.j (j=0,1) containg the data-generating-processes for control (j=0) and experimental treatment (j=1)

getdata.looks<-function(sims,N1,N0,med0,hr,d.looks=NULL,d.exact=FALSE,dexact.len=50,seed.start=8316951,get.final=FALSE,
                        S.type="Weibull",cr.analysis="FG",cr_event.labels=c("censor","outcome","death"),
                        timing.scale=1,dtaus=c(7,14,21),
                        nph.cp1=NULL,cp.tau1=NULL,
                        nph.cp0=NULL,cp.tau0=NULL,
                        dgp.1=NULL,dgp.0=NULL,tau.inf=30,
                        Accrue.type='EastProjection',direction="GT",dilution=0.0,AC1=NULL,AC0=NULL,max.follow=Inf,
                                  d.override=TRUE,kmplot.nsims=9,
                                  quant=0.5,sim.status=FALSE,stop.dl.error=FALSE,
                                  interpolate=FALSE,show.curves=FALSE,S1.0=NULL,S0.0=NULL,t1.points=NULL,t0.points=NULL,
                                  show.interim.obs=FALSE,fit1=NULL,fit2=NULL,
                                  dropout.0=0.0,dropout.1=0.0,accrue.time=NULL,tau.min=1,tau.max=100,tau.seq=1,tau.IA.max=NULL,tau.FA.max=NULL,
                                  no.dropout=FALSE,shape.1=1,shape.0=1,tau.drop=tau.max,
                                  drop.type="east",
                                  details=TRUE,
                                  Entry1.obs=NULL,Entry0.obs=NULL,
                                  show.simcount=FALSE,scale.1=NULL,scale.0=NULL,override.ph=FALSE,
                                  get.mean=FALSE,mean.draws=0,
                                  gamma1=NULL,gamma2=NULL,bW=NULL,lamW=NULL,lamT=NULL,beta1=NULL,beta2=NULL,cw1=NULL,cw2=NULL,lamP0=NULL,lamP1=NULL,bP0=NULL,bP1=NULL){


if(S.type=="Sweibull-CR" & (is.null(dgp.1) | is.null(dgp.0)))  stop("For competing-risk simulation, dgp.0 and dgp.1 must be inputed")
if(S.type=="exp_piecewise" & (is.null(nph.cp1) | is.null(cp.tau1)))  stop("For piece-wise exponential parameters nph.cp and cp.tau must be inputed for treatment arm")
if(S.type=="exp_piecewise" & (is.null(nph.cp0) | is.null(cp.tau0)))  stop("For piece-wise exponential parameters nph.cp and cp.tau must be inputed for control arm")
  

t.start <- proc.time()[1]

set.seed(seed.start)

if(is.null(AC1) & is.null(AC0)){
accruals<-get.accrual(N=(N1+N0),rand.ratio=2,Accrue.type="projection",Ramp.up=Ramp.up,nRamp.after=nRamp.after)
AC0<-accruals$AC.0
AC1<-accruals$AC.1
}

# If enrollment pattern is fixed (default)
# The time to complete accrual for both arms is
accrue.time<-ceiling(max(c(AC0,AC1)))
Entry1.obs<-AC1
Entry0.obs<-AC0
N1<-length(Entry1.obs)
N0<-length(Entry0.obs)
# Entry times for patients are then fixed at these times
# (i.e., Entry time are first simulated according to above and then fixed in outcome simulations)
if(details){    
cat("N1=",c(length(Entry1.obs)),"\n")
cat("N0=",c(length(Entry0.obs)),"\n")
}
  
  # Number of looks
  nL<-length(d.looks)
  if(details) cat("Number of looks=",c(nL),"\n")

  
  if(S.type=="Weibull") if(!override.ph & (shape.1 != shape.0)) stop("shapes need to be identical for the PH assumption to hold")
  
  if(S.type=="Weibull"){
    if(!override.ph){
      median.0<-med0
      c.0<-exp(log(-log(0.5))/shape.0)
      scale.0<-median.0/c.0
      scale.1<-scale.0*(hr^{-1/shape.1})
      c.1<-exp(log(-log(0.5))/shape.1)
      median.1<-scale.1*c.1
      # check PH model (at t=3)
      #lam0.3<-shape.0*(scale.0^{-shape.0})*3^{shape.0-1}
      #lam1.3<-shape.1*(scale.1^{-shape.1})*3^{shape.1-1}
      #cat("lam1.3/lam0.3=hr ?",c(lam1.3/lam0.3),"\n")
    }
    if(override.ph){
      c.0<-exp(log(-log(0.5))/shape.0)
      median.0<-scale.0*c.0
      c.1<-exp(log(-log(0.5))/shape.1)
      median.1<-scale.1*c.1
    }
  }
  
  if(S.type=="CrossOver"){
    median.0<-med0
    median.1<-median.0/hr
    c.0<-exp(log(-log(0.5))/shape.0)
    scale.0<-median.0/c.0
    scale.1<-scale.0*(hr^{-1/shape.1})
  }
  
  if(S.type=="Weibull-CR"){
  # Weibull cause-specific times
  # Here shape,scale are specified 
  # So, get medians
  
  shape.1<-dgp.1$est["shape1"]
  scale.1<-dgp.1$est["scale1"]
  
  shape.0<-dgp.0$est["shape1"]
  scale.0<-dgp.0$est["scale1"]
  
# Note: These represent the latent survival time parameters
# But these are not identifiable
# Will use medians defined via CIFs

median.1<-min(dgp.1$tpoints[which(dgp.1$F1.sub>=0.5)])
median.0<-min(dgp.0$tpoints[which(dgp.0$F1.sub>=0.5)])

  }
  
  
  if(drop.type=="uniform"){
  
if(hr==1) scale.1<-scale.0  
    # Uniform censoring parameters
    a.0<-0.0
    b.0<-uniform.b(k=dropout.0,lambda=1/scale.0,tau=tau.drop)
    a.1<-0.0
    b.1<-uniform.b(k=dropout.1,lambda=1/scale.1,tau=tau.drop)
    # If no.dropout --> turn dropout OFF
    if(no.dropout) b.0<-b.1<-999999999
  }
  if(drop.type=="east"){
if(hr==1) dropout.1=dropout.0    
    lambdac.0<--tau.drop/log(1-dropout.0)
    lambdac.1<--tau.drop/log(1-dropout.1)
    if(no.dropout) lambdac.0<-lambdac.1<-999999999
  }
  
  # drop-out
  pcnt.ltfu0<-pcnt.ltfu1<-pcnt.ltfu<-rep(NA,sims)
  pcnt.cens0<-pcnt.cens1<-pcnt.cens<-rep(NA,sims)
  
  d.maxs<-rep(NA,sims)
  Tau.interims<-Tau.finals<-rep(NA,sims)
  
  # CrossOver %'s
  xo1<-xo0<-rep(NA,sims)
  m.xo1<-m.xo0<-rep(NA,sims)
  Tau.times<-seq(tau.min,tau.max,by=tau.seq)
  #######################################################
  # We take data-cuts at times in Tau.times (in months)
  # Cuts are monthly
  #######################################################
  
  if(details) cat("Note: Data-cuts are from tau.first,tau.max=",c(min(Tau.times),max(Tau.times)),"\n")
  
  for(sim in 1:sims){
    if(show.simcount) cat("simulation=",c(sim),"\n")
    
    set.seed(seed.start+1000*sim)
    
    ##############################
    # General "full maximal data"
    ##############################
    
    # Represents data that would be observed if 
    # study was full duration (up to tau.max) 
    
    ####################
    # Experimental data
    ####################

  #cat("S.type=",c(S.type),"\n")
    
  
if(S.type!="Weibull-CR" & S.type!="exp_piecewise"){
  
    if(dilution==0){
    if(drop.type=="uniform"){
      data1<-Data.sim(n=N1,shape=shape.1,scale=scale.1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                      drop.type=drop.type,cmin=a.1,cmax=b.1,Accrue.type=Accrue.type,Entry.obs=Entry1.obs,S.type=S.type,S.0=S1.draw,tpoints.0=t1.draw,
                      R=rep(1,N1),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                      lamP=lamP1,bP=bP1)
    }
    if(drop.type=="east"){
      data1<-Data.sim(n=N1,shape=shape.1,scale=scale.1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                      drop.type=drop.type,scale.drop=lambdac.1,Accrue.type=Accrue.type,Entry.obs=Entry1.obs,S.type=S.type,S.0=S1.draw,tpoints.0=t1.draw,
                      R=rep(1,N1),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                      lamP=lamP1,bP=bP1)
    }
    
    } # No dilution
    
    if(dilution>0){
      if(S.type=="Weibull-CR") stop("Not setup for CR model")
      n1.dilute<-round(dilution*N1,0)
      n1.nodilute<-N1-n1.dilute
      id.treat<-c(1:N1)
      # Randomly select from Entry1.obs
      id.dilute<-sample(id.treat,n1.dilute)
      id.nodilute<-setdiff(id.treat,id.dilute)
      
      Entry1.dilute<-Entry1.obs[id.dilute]
      Entry1.nodilute<-Entry1.obs[id.nodilute]
      
      # For the diluted subpop, n1.dilute, assume these are identical to control
      
      if(drop.type=="uniform"){
        data1a<-Data.sim(n=n1.dilute,shape=shape.0,scale=scale.0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                        drop.type=drop.type,cmin=a.0,cmax=b.0,Accrue.type=Accrue.type,Entry.obs=Entry1.dilute,S.type=S.type,S.0=S0.draw,tpoints.0=t0.draw,
                        R=rep(1,n1.dilute),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                        lamP=lamP1,bP=bP1)
      }
      if(drop.type=="east"){
        data1a<-Data.sim(n=n1.dilute,shape=shape.0,scale=scale.0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                        drop.type=drop.type,scale.drop=lambdac.0,Accrue.type=Accrue.type,Entry.obs=Entry1.dilute,S.type=S.type,S.0=S0.draw,tpoints.0=t0.draw,
                        R=rep(1,n1.dilute),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                        lamP=lamP1,bP=bP1)
      }
      
      # Non-diluted
      
      if(drop.type=="uniform"){
        data1b<-Data.sim(n=n1.nodilute,shape=shape.1,scale=scale.1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                        drop.type=drop.type,cmin=a.1,cmax=b.1,Accrue.type=Accrue.type,Entry.obs=Entry1.nodilute,S.type=S.type,S.0=S1.draw,tpoints.0=t1.draw,
                        R=rep(1,n1.nodilute),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                        lamP=lamP1,bP=bP1)
      }
      if(drop.type=="east"){
        data1b<-Data.sim(n=n1.nodilute,shape=shape.1,scale=scale.1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                        drop.type=drop.type,scale.drop=lambdac.1,Accrue.type=Accrue.type,Entry.obs=Entry1.nodilute,S.type=S.type,S.0=S1.draw,tpoints.0=t1.draw,
                        R=rep(1,n1.nodilute),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                        lamP=lamP1,bP=bP1)
      }
      
      data1<-rbind(data1a,data1b)
      
    } # With dilution
    
    
    # Control data
  if(drop.type=="uniform"){
    data0<-Data.sim(n=N0,shape=shape.0,scale=scale.0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                    drop.type=drop.type,cmin=a.0,cmax=b.0,Accrue.type=Accrue.type,Entry.obs=Entry0.obs,S.type=S.type,S.0=S0.draw,tpoints.0=t0.draw,
                    R=rep(0,N0),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                    lamP=lamP0,bP=bP0)
  }
  
  if(drop.type=="east"){
    data0<-Data.sim(n=N0,shape=shape.0,scale=scale.0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                    drop.type=drop.type,scale.drop=lambdac.0,Accrue.type=Accrue.type,Entry.obs=Entry0.obs,S.type=S.type,S.0=S0.draw,tpoints.0=t0.draw,
                    R=rep(0,N0),gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2,
                    lamP=lamP0,bP=bP0)
  }
  
  # data0,data1 contain outcomes in the *absence* of censoring by deaths;  Only random and administrative censoring
  
  # Output the first simulation
  if(sim==1){
    data1$treat<-1
    data0$treat<-0
    df.sim1<-rbind(data0,data1)
  }
  # For CR simulations, this will be output as the final look dataset
  }
  
  
  # Piecewise exponential
  
  if(S.type=="exp_piecewise"){
    
  
    if(drop.type=="uniform"){
      data1<-Data.sim.piecewise(n=N1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                                nph.cp=nph.cp1,cp.tau=cp.tau1,
                                drop.type=drop.type,cmin=a.1,cmax=b.1,Accrue.type=Accrue.type,Entry.obs=Entry1.obs)
    }
    if(drop.type=="east"){
      data1<-Data.sim.piecewise(n=N1,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                                nph.cp=nph.cp1,cp.tau=cp.tau1,
                                drop.type=drop.type,scale.drop=lambdac.1,Accrue.type=Accrue.type,Entry.obs=Entry1.obs)
    }
    
    
    if(drop.type=="uniform"){
      data0<-Data.sim.piecewise(n=N0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                                nph.cp=nph.cp0,cp.tau=cp.tau0,
                                drop.type=drop.type,cmin=a.0,cmax=b.0,Accrue.type=Accrue.type,Entry.obs=Entry0.obs)
    }
    if(drop.type=="east"){
      data0<-Data.sim.piecewise(n=N0,Accrue=accrue.time,tau=tau.max,max.follow=max.follow,
                                nph.cp=nph.cp0,cp.tau=cp.tau0,
                                drop.type=drop.type,scale.drop=lambdac.0,Accrue.type=Accrue.type,Entry.obs=Entry0.obs)
    }
    
    if(sim==1){
      data1$treat<-1
      data0$treat<-0
      df.sim1<-rbind(data0,data1)
    }
  }
    
  
# Competing-risk DGP

if(S.type=="Weibull-CR"){
  
  data0 <- Data.CR.sim(
    n = N0,
    dgp = dgp.0,
    Accrue = accrue.time,
    tau = tau.max,
    max.follow = max.follow,
    drop.type = drop.type,
    scale.drop = lambdac.0,
    Accrue.type = Accrue.type,
    Entry.obs = Entry0.obs,
    tau.inf = tau.inf
  )
  
  data1 <- Data.CR.sim(
    n = N1,
    dgp = dgp.1,
    Accrue = accrue.time,
    tau = tau.max,
    max.follow = max.follow,
    drop.type = drop.type,
    scale.drop = lambdac.1,
    Accrue.type = Accrue.type,
    Entry.obs = Entry1.obs,
    tau.inf = tau.inf
  )
  
}
  
    n1<-length(data1$Entry)
    n0<-length(data0$Entry)
    
    if(show.simcount) cat("n1,n0=",c(n1,n0),"\n")
    
    ###################################################
    # Find tau.interim such that # events >= d.interim
    ###################################################
    
    if(S.type!="Weibull-CR"){
    Entry1<-data1$Entry
    Follow1<-data1$FollowUp
    C1<-data1$Censoring
    T1<-data1$TrueSurvival 
    D1<-data1$DropOut
    XO1.event<-data1$EventXover
    XO1<-data1$Xover

    Entry0<-data0$Entry
    Follow0<-data0$FollowUp
    C0<-data0$Censoring
    T0<-data0$TrueSurvival 
    D0<-data0$DropOut
    XO0.event<-data0$EventXover
    XO0<-data0$Xover
  
    events.max<-sum(ifelse(c(T0,T1)<=c(C0,C1),1,0))
    }
    
    
    if(S.type=="Weibull-CR"){
      Entry1<-data1$Entry
      Follow1<-data1$FollowUp
      C1<-data1$Censoring
      T1<-data1$TrueSurvival 
      D1<-data1$DropOut
      Cause1<-data1$Cause1
      Event1<-data1$Event
      
      Entry0<-data0$Entry
      Follow0<-data0$FollowUp
      C0<-data0$Censoring
      T0<-data0$TrueSurvival 
      D0<-data0$DropOut
      Cause0<-data0$Cause1
      Event0<-data0$Event
      
      # Cause-1 events
      E0<-ifelse(Event0==1,1,0)
      E1<-ifelse(Event1==1,1,0)
      
      events.max<-sum(c(E0,E1))
    }
    
    d.maxs[sim]<-events.max
    
    if(S.type=="CrossOver"){
      xo1[sim]<-mean(XO1.event)
      xo0[sim]<-mean(XO0.event)
      m.xo1[sim]<-mean(XO1[which(XO1.event==1)])
      m.xo0[sim]<-mean(XO0[which(XO0.event==1)])
    }
    
    Events.tau<-rep(NA,length(Tau.times))
    
    # Find the earliest tau (tau.interim) such that
    # d.interim events are observed.  
    # This is the interim analysis trigger
    
    for(tau.index in 1:length(Tau.times)){
      tau<-Tau.times[tau.index]
      #######################################
      # Experimental data at follow-up = tau
      #######################################
      F1.tau<-tau-Entry1  
      C1.tau<-pmin(D1,F1.tau) 
      
      if(S.type!="Weibull-CR"){
      Event1.tau<-ifelse(T1<=C1.tau,1,0) 
      }
        
      if(S.type=="Weibull-CR"){
      Eta1.tau<-rep(0,n1)
      Eta1.tau[which(T1<=C1.tau & Cause1==1)]<-1
      Eta1.tau[which(T1<=C1.tau & Cause1==0)]<-2
      Event1.tau<-ifelse(Eta1.tau==1,1,0)  
      }

      d1.tau<-sum(Event1.tau[which(F1.tau>0)])

      # Control data at follow-up = tau
      F0.tau<-tau-Entry0  
      C0.tau<-pmin(D0,F0.tau) 
      
      if(S.type!="Weibull-CR"){
        Event0.tau<-ifelse(T0<=C0.tau,1,0) 
      }
      
      if(S.type=="Weibull-CR"){
        Eta0.tau<-rep(0,n0)
        Eta0.tau[which(T0<=C0.tau & Cause0==1)]<-1
        Eta0.tau[which(T0<=C0.tau & Cause0==0)]<-2
        Event0.tau<-ifelse(Eta0.tau==1,1,0)  
      }
      
      d0.tau<-sum(Event0.tau[which(F0.tau>0)])

      Events.tau[tau.index]<-d0.tau+d1.tau
    }
    
    # Find final analysis
    d.final<-max(d.looks)
    tau.final<-ifelse(max(Events.tau)<d.final,max(Tau.times),min(Tau.times[which(Events.tau>=d.final)]))
    
    tau.index<-which(Tau.times==tau.final)
    if(length(tau.index)>1) tau.index<-tau.index[1]
    if(!d.override & Events.tau[tau.index]<d.final) stop("Final number of events not reached --> increase tau.max")
    
    tau1.final<-max(Tau.times[which(Events.tau<=d.final)])
    tau2.final<-min(Tau.times[which(Events.tau>=d.final)])
    
    # If tau2.final=Inf then d.final not met for tau cuts < tau.max
    # Then set d.final to tau.max 
    if(tau2.final==Inf){
      d.final<-tau.max
    }
    
    
    if(d.exact & tau1.final < tau2.final & tau2.final<tau.max){
      #cat("Searching within Final [tl,tu]=",c(tau1.final,tau2.final),"\n")
      # Find intervals containing d.final and then refine search between those
      #tau1.final<-max(Tau.times[which(Events.tau<=d.final)])
      # Search within [tau1,tau2]
      Tau.times.final<-seq(tau1.final,tau2.final,length=dexact.len)
      Events.tau.final<-rep(NA,length(Tau.times.final))
      for(tau.index in 1:length(Tau.times.final)){
        tau<-Tau.times.final[tau.index]
        #######################################
        # Experimental data at follow-up = tau
        #######################################
      
        F1.tau<-tau-Entry1  
        C1.tau<-pmin(D1,F1.tau) 
        
        if(S.type!="Weibull-CR"){
          Event1.tau<-ifelse(T1<=C1.tau,1,0) 
        }
        
        if(S.type=="Weibull-CR"){
          Eta1.tau<-rep(0,n1)
          Eta1.tau[which(T1<=C1.tau & Cause1==1)]<-1
          Eta1.tau[which(T1<=C1.tau & Cause1==0)]<-2
          Event1.tau<-ifelse(Eta1.tau==1,1,0)  
        }
        
        d1.tau<-sum(Event1.tau[which(F1.tau>0)])
        
        # Control data at follow-up = tau
        F0.tau<-tau-Entry0  
        C0.tau<-pmin(D0,F0.tau) 
        
        if(S.type!="Weibull-CR"){
          Event0.tau<-ifelse(T0<=C0.tau,1,0) 
        }
        
        if(S.type=="Weibull-CR"){
          Eta0.tau<-rep(0,n0)
          Eta0.tau[which(T0<=C0.tau & Cause0==1)]<-1
          Eta0.tau[which(T0<=C0.tau & Cause0==0)]<-2
          Event0.tau<-ifelse(Eta0.tau==1,1,0)  
        }
        d0.tau<-sum(Event0.tau[which(F0.tau>0)])
      
        Events.tau.final[tau.index]<-d0.tau+d1.tau
      }
      tau.final<-ifelse(max(Events.tau.final)<d.final,max(Tau.times.final),min(Tau.times.final[which(Events.tau.final>=d.final)]))
      tau.index<-which(Tau.times.final==tau.final)
      if(length(tau.index)>1) tau.index<-tau.index[1]
      if(Events.tau.final[tau.index]<d.final & stop.dl.error) stop("Final number of events not reached for exact matching--> increase tau.max")
    }
    
    
    for(look in 1:length(d.looks)){
      dl<-d.looks[look]
      #tau.interim<-min(min(Tau.times[which(Events.tau>=d.interim)]),max(Tau.times))
      # Will not go beyond max(Tau.times) for interim
      tau.look<-ifelse(max(Events.tau)<dl,max(Tau.times),min(Tau.times[which(Events.tau>=dl)]))
      tau1.look<-max(Tau.times[which(Events.tau<=dl)])
      tau2.look<-min(Tau.times[which(Events.tau>=dl)])
      
      if(tau2.look==Inf){
        tau2.look<-tau.max
      }
      
      if(d.exact & tau1.look < tau2.look){
        # Find intervals containing d.final and then refine search between those
        # Search within [tau1,tau2]
        #cat("Searching within Interim [tl,tu]=",c(tau1.look,tau2.look),"\n")
        
        Tau.times.look<-seq(tau1.look,tau2.look,length=dexact.len)
        Events.tau.look<-rep(NA,length(Tau.times.look))
        for(tau.index in 1:length(Tau.times.look)){
          tau<-Tau.times.look[tau.index]
          #######################################
          # Experimental data at follow-up = tau
          #######################################
          
          F1.tau<-tau-Entry1  
          C1.tau<-pmin(D1,F1.tau) 
          
          if(S.type!="Weibull-CR"){
            Event1.tau<-ifelse(T1<=C1.tau,1,0) 
          }
          
          if(S.type=="Weibull-CR"){
            Eta1.tau<-rep(0,n1)
            Eta1.tau[which(T1<=C1.tau & Cause1==1)]<-1
            Eta1.tau[which(T1<=C1.tau & Cause1==0)]<-2
            Event1.tau<-ifelse(Eta1.tau==1,1,0)  
          }
          
          d1.tau<-sum(Event1.tau[which(F1.tau>0)])
          
          # Control data at follow-up = tau
          F0.tau<-tau-Entry0  
          C0.tau<-pmin(D0,F0.tau) 
          
          if(S.type!="Weibull-CR"){
            Event0.tau<-ifelse(T0<=C0.tau,1,0) 
          }
          
          if(S.type=="Weibull-CR"){
            Eta0.tau<-rep(0,n0)
            Eta0.tau[which(T0<=C0.tau & Cause0==1)]<-1
            Eta0.tau[which(T0<=C0.tau & Cause0==0)]<-2
            Event0.tau<-ifelse(Eta0.tau==1,1,0)  
          }
          d0.tau<-sum(Event0.tau[which(F0.tau>0)])
          
        Events.tau.look[tau.index]<-d0.tau+d1.tau
        }
        tau.look<-ifelse(max(Events.tau.look)<dl,max(Tau.times.look),min(Tau.times.look[which(Events.tau.look>=dl)]))
        
        if(Events.tau.look[which(Tau.times.look==tau.look)] != dl){
          if(stop.dl.error)  stop("Target looks not met")
        }  
      }
      
       ###################
      # Interim analysis
      ###################
  
      if(S.type!="Weibull-CR") {
        Looks <-
          get.analyses(
            sim = sim,
            tau.interim = tau.look,
            tau.final = tau.final,
            timing.scale=timing.scale,
            direction = direction,
            get.final = get.final,
            Entry.1 = Entry1,
            D.1 = D1,
            T.1 = T1,
            Entry.0 = Entry0,
            D.0 = D0,
            T.0 = T0,
            quant = quant,
            show.simcount = show.simcount,
            get.mean = get.mean,
            mean.draws = mean.draws,
            show.km = (sim <= kmplot.nsims),
            titleit.IA = c("KM Interim"),
            titleit.FA = c("KM Final"),
            fit1 = fit1,
            fit2 = fit2,
            show.interim.obs = show.interim.obs
          )
      }
      
      if(S.type=="Weibull-CR") {
        Looks <-
          get.analyses.cr(
            data1 = data1,
            data0 = data0,
            sim = sim,
          dtaus=dtaus,
            tau.interim = tau.look,
            tau.final = tau.final,
            timing.scale=timing.scale,
            cr_event.labels = cr_event.labels,
            direction = direction,
            get.final = get.final,
            cr.analysis = cr.analysis,
            quant = quant,
            show.simcount = show.simcount,
            show.km = (sim <= kmplot.nsims),
            titleit.IA = c("KM Interim"),
            titleit.FA = c("KM Final"),
            fit1 = fit1,
            fit2 = fit2,
            show.interim.obs = show.interim.obs
          )
     
      # For sim=1 and final look, output dataset  
        
if(sim==1 & tau.look==tau.final) df.sim1<-Looks$df.cr           
      }
 
      
      if(exists(paste("InterimLook",look,sep="."),inherits=TRUE)) temp1<-get(paste("InterimLook",look,sep="."))
      if(!exists(paste("InterimLook",look,sep="."),inherits=TRUE)) temp1<-NULL
      
      temp2<-rbind(temp1,Looks$out.interim.sim)
      assign(paste("InterimLook",look,sep="."),temp2)
    }
    rm("temp1","temp2")
    
    
    # Save as datasets
    for(look in 1:length(d.looks)){
      dl<-d.looks[look]
      temp1<-get(paste("InterimLook",look,sep="."))
      
    # Record death rates for CR  
    if(S.type=="Weibull-CR")  colnames(temp1)<-c("sim","tau.IA","d.IA","Zlr.IA","bhat.IA","se.bhat.IA","m1.IA","m0.IA","drop1.IA","drop0.IA","m.diff.IA","se.diff.IA","N1.IA","N0.IA","N.IA","dth1.IA","dth2.IA","dth3.IA")
    if(S.type!="Weibull-CR")  colnames(temp1)<-c("sim","tau.IA","d.IA","Zlr.IA","bhat.IA","se.bhat.IA","m1.IA","m0.IA","drop1.IA","drop0.IA","m.diff.IA","se.diff.IA","N1.IA","N0.IA","N.IA")
    
      assign(paste("InterimLook",look,sep="."),as.data.frame(temp1))
      
      if(dl==max(d.looks)){
        temp<-get(paste("InterimLook",look,sep="."))
        LFU1.tau<-mean(temp[,"drop1.IA"])
        LFU0.tau<-mean(temp[,"drop0.IA"])
        pcnt.ltfu1[sim]<-LFU1.tau
        pcnt.ltfu0[sim]<-LFU0.tau
        pcnt.ltfu[sim]<-(n0*LFU0.tau+n1*LFU1.tau)/(n0+n1)
      }
      
      rm("temp1")
      
    }
    
    
    if(sim.status){
      show.status(ss=sim,t.start=t.start,sims=sims)
    }
  }
  if(details){
    if(S.type=="CrossOver"){
      cat("# Xover Avg percent Control=",c(mean(xo0)),"\n")
      cat("# Among Xovers Avg Control=",c(mean(m.xo0)),"\n")
      cat("# Xover Avg percent Experimental=",c(mean(xo1)),"\n")
      cat("# Among Xovers Avg Experimental=",c(mean(m.xo1)),"\n")
    }
    cat("# Accrual Type",c(Accrue.type),"\n")
    if(Accrue.type=="uniform") cat("# Accrual Time=",c(round(accrue.time,2)),"\n")
    if(Accrue.type=="observed" | Accrue.type=="projection" | Accrue.type=="EastProjection") cat("# Accrual Time=",c(accrue.time),"\n")
    #cat("# Max study duration (final)=",c(tau.max),"\n")
    cat("# Weibull scale parameters (C,T)=",c(scale.0,scale.1),"\n")
    cat("# Weibull shape parameters (C,T)=",c(shape.0,shape.1),"\n")
    if(shape.0==1 & shape.1==1){
      cat("# Shape parameters in terms of exponential model (C,T)=",
          c((1/scale.0),(1/scale.1)),"\n")
    }
    
    if(S.type=="Weibull") cat("# Medians (Cntrl,Exp) =",c(median.0,median.1),"\n")
    if(S.type=="Weibull-CR") cat("# Medians (Cntrl,Exp)=",c(median.0,median.1),"\n")
    
    cat("# HR=,n1,n2=,sims=",c(hr,n1,n0,sims),"\n")
    
    if(!no.dropout){
      cat("# Target dropout rate (Cntrl,Exp)=",c(dropout.0,dropout.1),"\n")
      
      if(drop.type=="uniform"){
        cat("# Control values for U(a,b) drop-out distribution:",c(a.0,b.0),"\n")
        cat("# Treat values for U(a,b) drop-out distribution:",c(a.1,b.1),"\n")
      }
      if(drop.type=="exponential" | drop.type=="east"){
        cat("# Control,Exp drop-out parameters terms of exponential model (C,T)=",
            c((1/lambdac.0),(1/lambdac.1)),"\n")
        cat("# Note: these are analogous to etaC and etaE in gsDesign package","\n")
      }
      cat("# % drop-out (Control,Exp,All): Final Analysis (no stopping)",
          c(mean(pcnt.ltfu0),mean(pcnt.ltfu1),mean(pcnt.ltfu)),"\n")
    }
    
    if(no.dropout & S.type!="Weibull-CR"){
    cat("# % drop-out (Control,Exp,All): Final Analysis (no stopping)",
      c(mean(pcnt.ltfu0),mean(pcnt.ltfu1),mean(pcnt.ltfu)),"\n")
    }  
    
    if(no.dropout & S.type=="Weibull-CR"){
      cat("# % drop-out (Control,Exp,All): Final Analysis (no stopping)",
          c(mean(pcnt.ltfu0),mean(pcnt.ltfu1),mean(pcnt.ltfu)),"\n")
    }  
    cat("# Max events=d(tau.max)","\n")
    cat("# Mean max events=",c(mean(d.maxs)),"\n")

    t.end<-proc.time()[1]
    t.min<-(t.end-t.start)/60
    
    est.1000<-t.min*(1000/sims)
    
    cat("Simulations=",c(sims),"\n")
    cat("Time (min)=",c(t.min),"\n")
    cat("# Estimated minutes per 1,000 simulations=",c(est.1000),"\n")  
    
  }
  
  
  if(details & S.type=="Specified" & interpolate & show.curves){
    win.graph()
    par(mfrow=c(2,2))
    plot(t1.draw,S1.draw,lwd=3,col="grey",type="s")
    lines(t1.points,S1.0,col="blue",lwd=3,type="s",lty=2)
    plot(t0.draw,S0.draw,lwd=3,col="grey",type="s")
    lines(t0.points,S0.0,col="blue",lwd=3,type="s",lty=2)
    
    plot(t1.draw,S1.draw,lwd=3,col="black",type="s")
    lines(t0.draw,S0.draw,lwd=3,col="blue",type="s")
  }
  
  # Output up to 20 looks
  
  Out<-new.env()
  if(exists("InterimLook.1",inherits=TRUE)) Out$Look1<-InterimLook.1
  if(exists("InterimLook.2",inherits=TRUE)) Out$Look2<-InterimLook.2
  if(exists("InterimLook.3",inherits=TRUE)) Out$Look3<-InterimLook.3
  if(exists("InterimLook.4",inherits=TRUE)) Out$Look4<-InterimLook.4
  if(exists("InterimLook.5",inherits=TRUE)) Out$Look5<-InterimLook.5
  if(exists("InterimLook.6",inherits=TRUE)) Out$Look6<-InterimLook.6
  if(exists("InterimLook.7",inherits=TRUE)) Out$Look7<-InterimLook.7
  if(exists("InterimLook.8",inherits=TRUE)) Out$Look8<-InterimLook.8
  if(exists("InterimLook.9",inherits=TRUE)) Out$Look9<-InterimLook.9
  if(exists("InterimLook.10",inherits=TRUE)) Out$Look10<-InterimLook.10
  
  if(exists("InterimLook.11",inherits=TRUE)) Out$Look11<-InterimLook.11
  if(exists("InterimLook.12",inherits=TRUE)) Out$Look12<-InterimLook.12
  if(exists("InterimLook.13",inherits=TRUE)) Out$Look13<-InterimLook.13
  if(exists("InterimLook.14",inherits=TRUE)) Out$Look14<-InterimLook.14
  if(exists("InterimLook.15",inherits=TRUE)) Out$Look15<-InterimLook.15
  if(exists("InterimLook.16",inherits=TRUE)) Out$Look16<-InterimLook.16
  if(exists("InterimLook.17",inherits=TRUE)) Out$Look17<-InterimLook.17
  if(exists("InterimLook.18",inherits=TRUE)) Out$Look18<-InterimLook.18
  if(exists("InterimLook.19",inherits=TRUE)) Out$Look19<-InterimLook.19
  if(exists("InterimLook.20",inherits=TRUE)) Out$Look20<-InterimLook.20
  
  
  Out<-as.list(Out)
  
  return(list(Out=Out,pcnt.ltfu1=pcnt.ltfu1,pcnt.ltfu0=pcnt.ltfu0,pcnt.ltfu=pcnt.ltfu,d.maxs=d.maxs,df.sim1=df.sim1,Out.final=Looks$out.final.sim,minutes=t.min))
}





Data.sim<-function(n,shape,scale,Accrue,tau,cmin=NULL,cmax=NULL,Accrue.type='uniform',Entry.obs,scale.drop=NULL,drop.type='uniform',
                   max.follow=Inf,
                   S.type="Weibull",S.0=NULL,tpoints.0=NULL,
                   R=NULL,gamma1=NULL,gamma2=NULL,bW=NULL,lamW=NULL,lamP=NULL,bP=NULL,lamT=NULL,beta1=NULL,beta2=NULL,cw1=NULL,cw2=NULL){
  
  if(Accrue.type=='uniform') Entry<-runif(n=n,0,Accrue) # Subjects' entry time (here Entry is random)
  if(Accrue.type=='observed' | Accrue.type=='projection' | Accrue.type=='EastProjection') Entry<-Entry.obs # (here entry is fixed at observed)
  
  Entry.sim<-Entry-min(Entry) # Start clock at "time zero" (clock starts from FPFV)
  F.sim<-tau-Entry.sim  # Subject's follow-up time (eg, tau=10, entry=year 1 --> follow-up = 9 years 
  # This will represent an "administrative" censoring variable
  
  # Note: tau can represent an interim look time period prior to the full 
  # study enrollment period (Accrue).  Thus, only subjects with Fsim>0 enter the analysis.
  
  # Drop-out time (from entry).  This is subject-specific relative to entry time
  if(drop.type=='uniform') D.sim<-runif(n=n,min=cmin,max=cmax)  # eg D=5
  if(drop.type=='exponential' | drop.type=='east') D.sim<-rweibull(n=n,shape=1,scale=scale.drop)
  
  
  if(S.type=="Weibull") T.sim<-rweibull(n=n,shape=shape,scale=scale) # eg T=8
  
  if(S.type=="Specified") T.sim<-DrawFromS(n=n,S.draw=S.0,tpoints.draw=tpoints.0)
  
  if(S.type=="CrossOver"){
    temp<-get.ZWsim(R=R,gamma1=gamma1,gamma2=gamma2,bW=bW,lamW=lamW,lamP=lamP,bP=bP,lamT=lamT,beta1=beta1,beta2=beta2,cw1=cw1,cw2=cw2)
    T.sim<-temp$T.obs
    P.sim<-temp$P.true
  }

  # In terms of "Study year"
  # Time "since study" = (Entry+Dsim) whereas Drop-out time (from study entry) is Dsim
  # from a subject-specific analysis time perspective.
  # eg entry = year 1 and drops out 5 years later ---> they drop out in Study year = 6

  # Include censorship at max follow-up

  D.sim<-pmin(D.sim,max.follow)
  
  C.sim<-pmin(D.sim,F.sim) # Censoring time is minimum of Follow-up Time and Drop-out Time
  # eg C=min(5,9) = 5

  Event.sim<-ifelse(T.sim<=C.sim,1,0) 
  # eg T=8, C=5 ---> censored at C=5
  Y.sim<-pmin(T.sim,C.sim)
  # eg Y=5
  Xover.sim<-NULL; Xover.event.sim<-NULL
  if(S.type=="CrossOver"){ 
    Xover.event.sim<-ifelse(P.sim<Y.sim,1,0)
    Xover.sim<-pmin(P.sim,Y.sim)
  }
  if(S.type=="CrossOver"){
    data.sim<-cbind(Entry.sim,F.sim,D.sim,T.sim,C.sim,Event.sim,Y.sim,Xover.sim,Xover.event.sim)
    # Restrict to F.sim>0 
    # So if the analysis is conducted prior to complete accrual then
    # only subjects that entered prior to tau would be available
    data.sim<-data.frame(data.sim[which(F.sim>0),])
    names(data.sim)<-c("Entry","FollowUp","DropOut","TrueSurvival","Censoring","Event","Survival","Xover","EventXover")
    return(data.sim)
  }
  if(S.type!="CrossOver"){
    data.sim<-cbind(Entry.sim,F.sim,D.sim,T.sim,C.sim,Event.sim,Y.sim)
    # Restrict to F.sim>0 
    # So if the analysis is conducted prior to complete accrual then
    # only subjects that entered prior to tau would be available
    data.sim<-data.frame(data.sim[which(F.sim>0),])
    names(data.sim)<-c("Entry","FollowUp","DropOut","TrueSurvival","Censoring","Event","Survival")
    return(data.sim)
  }
}



hweibull<-function (x, shape, scale = 1, log = FALSE) 
{
  if (any(shape <= 0) || any(scale <= 0)) 
    stop("scale and shape must be positive")
  res <- ifelse(x < 0, 0, shape * (x/scale)^(shape - 1)/scale)
  if (log) 
    res <- log(res)
  return(res)
}


Data.CR.sim<-function(n,dgp,Accrue,tau,cmin=NULL,cmax=NULL,Accrue.type='uniform',
Entry.obs,scale.drop=NULL,drop.type='uniform',max.follow=Inf,tau.inf=999){
  
  if(Accrue.type=='uniform') Entry<-runif(n=n,0,Accrue) # Subjects' entry time (here Entry is random)
  if(Accrue.type=='observed' | Accrue.type=='projection' | Accrue.type=='EastProjection') Entry<-Entry.obs # (here entry is fixed at observed)
  
  Entry.sim<-Entry-min(Entry) # Start clock at "time zero" (clock starts from FPFV)
  F.sim<-tau-Entry.sim  # Subject's follow-up time (eg, tau=10, entry=year 1 --> follow-up = 9 years 
  # This will represent an "administrative" censoring variable
  
  # Note: tau can represent an interim look time period prior to the full 
  # study enrollment period (Accrue).  Thus, only subjects with Fsim>0 enter the analysis.
  
  # Drop-out time (from entry).  This is subject-specific relative to entry time
  if(drop.type=='uniform') D.sim<-runif(n=n,min=cmin,max=cmax)  # eg D=5
  if(drop.type=='exponential' | drop.type=='east') D.sim<-rweibull(n=n,shape=1,scale=scale.drop)
  
  shape1<-dgp$est["shape1"]
  scale1<-dgp$est["scale1"]
  
  shape2<-dgp$est["shape2"]
  scale2<-dgp$est["scale2"]
  
  S.0<-dgp$S.overall
  tpoints.0<-dgp$tpoints
  
  T.all.sim<-DrawFromS(n=n,S.draw=S.0,tpoints.draw=tpoints.0)
  # Set inf to tau.inf
  T.all.sim[is.infinite(T.all.sim)]<-c(tau.inf)
  # For these times calculate probability of cause-1
  lam1.cs<-hweibull(T.all.sim,shape=shape1,scale=scale1)
  lam2.cs<-hweibull(T.all.sim,shape=shape2,scale=scale2)
  lam12.cs<-lam1.cs+lam2.cs
  # Probability of event time corresponding to cause-1
  q1<-ifelse(lam12.cs==0,0,lam1.cs/lam12.cs)
  
  Cause1.sim<-rbinom(n=n,size=1,prob=q1)
  
  T.sim<-T.all.sim
  
  # In terms of "Study year"
  # Time "since study" = (Entry+Dsim) whereas Drop-out time (from study entry) is Dsim
  # from a subject-specific analysis time perspective.
  # eg entry = year 1 and drops out 5 years later ---> they drop out in Study year = 6
  
  # Include censorship at max follow-up
  
  D.sim<-pmin(D.sim,max.follow)
  
  C.sim<-pmin(D.sim,F.sim) # Censoring time is minimum of Follow-up Time and Drop-out Time
  # eg C=min(5,9) = 5
  
  # Initiate Event.sim<-0 (censored)
  Event.sim<-rep(0,n)
  Event.sim[which(T.sim<=C.sim & Cause1.sim==1)]<-1
  Event.sim[which(T.sim<=C.sim & Cause1.sim==0)]<-2
  
  # eg T=8, C=5 ---> censored at C=5
  Y.sim<-pmin(T.sim,C.sim)
  # eg Y=5
  data.sim<-cbind(Entry.sim,F.sim,D.sim,T.sim,C.sim,Event.sim,Y.sim,Cause1.sim)
    # Restrict to F.sim>0 
    # So if the analysis is conducted prior to complete accrual then
    # only subjects that entered prior to tau would be available
    data.sim<-data.frame(data.sim[which(F.sim>0),])
    names(data.sim)<-c("Entry","FollowUp","DropOut","TrueSurvival","Censoring","Event","Survival","Cause1")
    return(data.sim)
  }


Data.sim.piecewise<-function(n,Accrue,tau,cmin=NULL,cmax=NULL,Accrue.type='uniform',Entry.obs,scale.drop=NULL,drop.type='uniform',
                   max.follow=Inf,nph.cp=NULL,cp.tau=NULL){
  
  
  if(Accrue.type=='uniform') Entry<-runif(n=n,0,Accrue) # Subjects' entry time (here Entry is random)
  if(Accrue.type=='observed' | Accrue.type=='projection' | Accrue.type=='EastProjection') Entry<-Entry.obs # (here entry is fixed at observed)
  
  Entry.sim<-Entry-min(Entry) # Start clock at "time zero" (clock starts from FPFV)
  F.sim<-tau-Entry.sim  # Subject's follow-up time (eg, tau=10, entry=year 1 --> follow-up = 9 years 
  # This will represent an "administrative" censoring variable
  
  # Note: tau can represent an interim look time period prior to the full 
  # study enrollment period (Accrue).  Thus, only subjects with Fsim>0 enter the analysis.
  
  # Drop-out time (from entry).  This is subject-specific relative to entry time
  if(drop.type=='uniform') D.sim<-runif(n=n,min=cmin,max=cmax)  # eg D=5
  if(drop.type=='exponential' | drop.type=='east') D.sim<-rweibull(n=n,shape=1,scale=scale.drop)
  
  
     # nph.cp are the hazard rates at change-points
    # cp.tau are the change-point times
    if(is.null(nph.cp) | is.null(cp.tau)) stop("Piecewise exponential paramaters are not specified")
   
  # This only allows 4 changepoints 
  #temp<-exp_cdfsim(n=n,theta=nph.cp,tau=cp.tau)
  #T.sim<-temp$time

  T.sim<-rpwexp(n=n,rate=nph.cp,intervals=cp.tau)
  
  # In terms of "Study year"
  # Time "since study" = (Entry+Dsim) whereas Drop-out time (from study entry) is Dsim
  # from a subject-specific analysis time perspective.
  # eg entry = year 1 and drops out 5 years later ---> they drop out in Study year = 6
  
  # Include censorship at max follow-up
  
  D.sim<-pmin(D.sim,max.follow)
  
  C.sim<-pmin(D.sim,F.sim) # Censoring time is minimum of Follow-up Time and Drop-out Time
  # eg C=min(5,9) = 5
  
  Event.sim<-ifelse(T.sim<=C.sim,1,0) 
  # eg T=8, C=5 ---> censored at C=5
  Y.sim<-pmin(T.sim,C.sim)
  # eg Y=5
    data.sim<-cbind(Entry.sim,F.sim,D.sim,T.sim,C.sim,Event.sim,Y.sim)
    # Restrict to F.sim>0 
    # So if the analysis is conducted prior to complete accrual then
    # only subjects that entered prior to tau would be available
    data.sim<-data.frame(data.sim[which(F.sim>0),])
    names(data.sim)<-c("Entry","FollowUp","DropOut","TrueSurvival","Censoring","Event","Survival")
    return(data.sim)
}




