
# alpha.looks represents the efficacy crossing boundary
# alpha.futility.looks represents futility crossing boundary
OC.interims<-function(d.looks,alpha.looks,data.looks,eff.binding=TRUE,details=TRUE,bhat.cut=0.80,d.override=FALSE,direction="GT",
                      alpha.futility.looks=NULL,stop.futility=TRUE){
  # if stop.futility=FALSE, then assume futility crossing is ignored and study is stopped only for efficacy
  
  if(is.null(alpha.futility.looks) | !stop.futility){
  stop.futility=FALSE  
  n.looks<-length(d.looks)
  alpha.futility.looks<-rep(1,n.looks)
  alpha.futility.looks[n.looks]<-alpha.looks[n.looks]
  }
  
  data.stopped<-NULL
  # datasets stored in Out
  data.looks<-data.looks$Out
  # Get first look
  # output<- look, d.target, d.obs, z.boundary, t(cross),t(L),t(U), Avg(hr), p(cross), p(cumulative) 
  d1<-d.looks[1]
  alpha.1<-alpha.looks[1]
  alphaf.1<-alpha.futility.looks[1]
  data.1<-data.looks$Look1
  n.trials<-nrow(data.1)
  #sim tau.IA d.IA    Zlr.IA     bhat.IA se.bhat.IA    m1.IA     m0.IA drop1.IA drop0.IA m.diff.IA se.diff.IA

  
  if(details) cat("Number of trials at first look",c(n.trials),"\n")
  
  pval.1<-1-pnorm(data.1$Zlr.IA)
  
  # Rejecting for efficacy
  reject.1<-c(pval.1<alpha.1)
  
  z.boundary<-qnorm(1-alpha.1)
  
  # Crossing futility: "Accepting Null"
  accept.1<-c(pval.1>alphaf.1)
  
  data.1$reject<-ifelse(pval.1<alpha.1,1,0)
  data.1$accept<-ifelse(pval.1>alphaf.1,1,0)
  
  data.1$look<-1
  data.1$hr<-exp(data.1$bhat.IA)
  data.1$hr.lastlook<-NA
  
  if(!stop.futility) data.stopped<-data.1[which(reject.1==1),]
  if(stop.futility) data.stopped<-data.1[which(reject.1==1 | accept.1==1),]
  
  bhat.boundary<-NULL
  
  #print(summary(data.1[,"bhat.IA"]))
  
  Pbhat.1<-mean(c(exp(data.1$bhat.IA)<=bhat.cut))
  
  bhat.boundary<-c(bhat.boundary,exp(data.1$bhat.IA[reject.1]))
  
if(direction=="GT")  bhat.boundary.1<-max(bhat.boundary,na.rm=TRUE)
if(direction=="LT")  bhat.boundary.1<-min(bhat.boundary,na.rm=TRUE)
  
  
  reject.ByNow<-rep(0.0,n.trials) # initiate Reject.ByNow
  reject.ByNow<-reject.ByNow+reject.1
  
  accept.ByNow<-rep(0.0,n.trials) # initiate Accept.ByNow
  accept.ByNow<-accept.ByNow+accept.1
  
  
  result<-cbind(1,d1,alpha.1,mean(data.1$d.IA),Pbhat.1,bhat.boundary.1,z.boundary,mean(data.1$tau.IA),quantile(data.1$tau.IA,c(0.025)),quantile(data.1$tau.IA,c(0.975)),
                mean(data.1$N1.IA),quantile(data.1$N1.IA,c(0.025)),quantile(data.1$N1.IA,c(0.975)),
                mean(data.1$N0.IA),quantile(data.1$N0.IA,c(0.025)),quantile(data.1$N0.IA,c(0.975)),
                mean(data.1$N.IA),quantile(data.1$N.IA,c(0.025)),quantile(data.1$N.IA,c(0.975)),
                sum(reject.1),sum(reject.1)/n.trials,sum(!reject.ByNow),mean(reject.1),mean(reject.ByNow),
                sum(accept.1),sum(accept.1)/n.trials,sum(!accept.ByNow),mean(accept.1),mean(accept.ByNow))
  
  if(!d.override){
    
    data.last<-data.1 # initialize last-look
    
    for(look in 2:length(d.looks)){
      
      dL<-d.looks[look]
      alpha.L<-alpha.looks[look]
      
      alphaf.L<-alpha.futility.looks[look]
      
      #cat("l,d(l),alpha.l=",c(look,dL,alpha.L),"\n")
      
      z.boundary<-qnorm(1-alpha.L)
      zf.boundary<-qnorm(alphaf.L)
      
      data.L<-switch(look,data.looks$Look1,data.looks$Look2,data.looks$Look3,data.looks$Look4,data.looks$Look5,
                     data.looks$Look6,data.looks$Look7,data.looks$Look8,data.looks$Look9,data.looks$Look10,
                     data.looks$Look11,data.looks$Look12,data.looks$Look13,data.looks$Look14,data.looks$Look15,
                     data.looks$Look16,data.looks$Look17,data.looks$Look18,data.looks$Look19,data.looks$Look20)
      
# Among those that did not cross at last look ("Not Yet Crossed")
      
if(stop.futility)  get.NYC<-c(!reject.ByNow & !accept.ByNow)
if(!stop.futility)  get.NYC<-c(!reject.ByNow)

# Among trials that have not yet crossed for efficacy or futility (if applicable)      
data.L.NYC<-data.L[get.NYC,]

if(details) cat("Number of trials at look=, n=",c(look,nrow(data.L.NYC)),"\n")

pval.L<-1-pnorm(data.L.NYC$Zlr.IA)
reject.L<-c(pval.L<alpha.L)
  
accept.L<-c(pval.L>alphaf.L)

data.L.NYC$reject<-ifelse(pval.L<alpha.L,1,0) # Rejection status for current look
    
data.L.NYC$accept<-ifelse(pval.L>alphaf.L,1,0) # Rejection status for current look

      data.L.NYC$look<-look
      data.L.NYC$hr<-exp(data.L.NYC$bhat.IA)   # hr estimate for current look
      
      sim.match<-match(data.L.NYC$sim,data.last$sim) # Match sim id from current look with last look
      data.L.NYC$hr.lastlook<-data.last$hr[sim.match] # hr estimate from last look for each sim 
      
      if(look < length(d.looks)){
        if(!stop.futility) data.stopped<-rbind(data.stopped,data.L.NYC[which(reject.L==1),]) # concatenate stopped data from last look with current stops
        if(stop.futility) data.stopped<-rbind(data.stopped,data.L.NYC[which(reject.L==1 | accept.L==1),]) 
        }
      
      if(look == length(d.looks)){
      data.stopped<-rbind(data.stopped,data.L.NYC)  # If last look, then stopped data is just the current=last look
      }
      
      Pbhat.L<-mean(c(exp(data.L.NYC$bhat.IA)<=bhat.cut))
      
      # est(HRs) corresponding to efficacy rejections
      bhat.boundary<-c(bhat.boundary,exp(data.L.NYC$bhat.IA[reject.L]))
      
      # max(est(HRs))
 if(direction=="GT") bhat.boundary.L<-max(bhat.boundary,na.rm=TRUE)
 if(direction=="LT") bhat.boundary.L<-min(bhat.boundary,na.rm=TRUE)
      
reject.ByNow[get.NYC]<-reject.ByNow[get.NYC]+reject.L
      
accept.ByNow[get.NYC]<-accept.ByNow[get.NYC]+accept.L
      
# Return summaries for sequences at this look

data.L<-data.L.NYC

      temp<-cbind(look,dL,alpha.L,mean(data.L$d.IA),Pbhat.L,bhat.boundary.L,z.boundary,mean(data.L$tau.IA),quantile(data.L$tau.IA,c(0.025)),quantile(data.L$tau.IA,c(0.975)),
                  mean(data.L$N1.IA),quantile(data.L$N1.IA,c(0.025)),quantile(data.L$N1.IA,c(0.975)),
                  mean(data.L$N0.IA),quantile(data.L$N0.IA,c(0.025)),quantile(data.L$N0.IA,c(0.975)),
                  mean(data.L$N.IA),quantile(data.L$N.IA,c(0.025)),quantile(data.L$N.IA,c(0.975)),
                  sum(reject.L),sum(reject.L)/n.trials,sum(!reject.ByNow),mean(reject.L),mean(reject.ByNow),
                  sum(accept.L),sum(accept.L)/n.trials,sum(!accept.ByNow),mean(accept.L),mean(accept.ByNow))
      
      result<-rbind(result,temp)
      
      rm("temp")
      
      data.last<-data.L.NYC
    }
  }
  
  colnames(result)<-c("Look","d","p(cross)","Avg(d.obs)","P(hr<=cut)","mDD","Z.boundary","t(cross)","t(L)","t(U)",
                      "Avg(N1)","N1(L)","N1(U)",
                      "Avg(N0)","N0(L)","N0(U)",
                      "Avg(N)","N(L)","N(U)",
                      "# Reject H0","% Reject H0","# no X yet","X-NotBefore","X-ByNow",
                      "# Accept H0","% Accept H0","# no Xf yet","Xf-NotBefore","Xf-ByNow")
  
  
  
  if(details) print(result)
  
  return(list(result=result,data.stopped=data.stopped))
  
}







# Details on HR and KM medians

Estimation.properties.multiple<-function(d.looks,alpha.looks,data.looks,binding=FALSE,details=TRUE,d.override=FALSE,direction="GT"){
  # Get first look
  # output<- look, d.target, d.obs, z.boundary, t(cross),t(L),t(U), Avg(hr), p(cross), p(cumulative) 
  d1<-d.looks[1]
  alpha.1<-alpha.looks[1]
  data.1<-data.looks$Look1
  n.trials<-nrow(data.1)
  #sim tau.IA d.IA    Zlr.IA     bhat.IA se.bhat.IA    m1.IA     m0.IA drop1.IA drop0.IA m.diff.IA se.diff.IA
pval.1<-1-pnorm(data.1$Zlr.IA)

  reject.1<-c(pval.1<alpha.1)
  z.boundary<-qnorm(1-alpha.1)
  
  bhat.boundary<-NULL
  
  #print(summary(data.1[,"bhat.IA"]))
  #Pbhat.1<-mean(c(exp(data.1$bhat.IA)<=bhat.cut))
  
  bhat.boundary<-c(bhat.boundary,exp(data.1$bhat.IA[reject.1]))
  
  bhat.boundary.1<-max(bhat.boundary,na.rm=TRUE)
  
  reject.ByNow<-rep(0.0,n.trials) # initiate Reject.ByNow
  reject.ByNow<-reject.ByNow+reject.1
  
  # Are medians estimable?
  estimable.0<-mean(!is.na(data.1$m0.IA))
  estimable.1<-mean(!is.na(data.1$m1.IA))
  
  result<-cbind(1,d1,alpha.1,min(data.1$d.IA),max(data.1$d.IA),min(data.1$tau.IA),max(data.1$tau.IA),mean(exp(data.1$bhat.IA)),quantile(exp(data.1$bhat.IA),c(0.025)),quantile(exp(data.1$bhat.IA),c(0.975)),
                estimable.0,mean(data.1$m0.IA,na.rm=TRUE),quantile(data.1$m0.IA,c(0.025),na.rm=TRUE),quantile(data.1$m0.IA,c(0.975),na.rm=TRUE),
                estimable.1,mean(data.1$m1.IA,na.rm=TRUE),quantile(data.1$m1.IA,c(0.025),na.rm=TRUE),quantile(data.1$m1.IA,c(0.975),na.rm=TRUE))
  
  if(!d.override){
    
    for(look in 2:length(d.looks)){
      
      dL<-d.looks[look]
      alpha.L<-alpha.looks[look]
      
      #cat("l,d(l),alpha.l=",c(look,dL,alpha.L),"\n")
      
      z.boundary<-qnorm(1-alpha.L)
      data.L<-switch(look,data.looks$Look1,data.looks$Look2,data.looks$Look3,data.looks$Look4,data.looks$Look5,
                     data.looks$Look6,data.looks$Look7,data.looks$Look8,data.looks$Look9,data.looks$Look10,
                     data.looks$Look11,data.looks$Look12,data.looks$Look13,data.looks$Look14,data.looks$Look15,
                     data.looks$Look16,data.looks$Look17,data.looks$Look18,data.looks$Look19,data.looks$Look20)
      
      # Among those that did not cross at last look ("Not Yet Crossed")
      
      get.NYC<-!reject.ByNow
      data.L.NYC<-data.L[get.NYC,]
      
      if(binding) data.L<-data.L.NYC # Impose "stopping" if crossed at previous interim
      
    pval.L<-1-pnorm(data.L.NYC$Zlr.IA)
    
      reject.L<-c(pval.L<alpha.L)
      
      #print(summary(data.L[,"bhat.IA"]))
      
      #Pbhat.L<-mean(c(exp(data.L$bhat.IA)<=bhat.cut))
      
      # est(HRs) corresponding to rejections
      bhat.boundary<-c(bhat.boundary,exp(data.L.NYC$bhat.IA[reject.L]))
      
      # max(est(HRs))
      bhat.boundary.L<-max(bhat.boundary,na.rm=TRUE)
      
      reject.ByNow[get.NYC]<-reject.ByNow[get.NYC]+reject.L
      
      temp<-cbind(look,dL,alpha.L,min(data.L$d.IA),max(data.L$d.IA),min(data.L$tau.IA),max(data.L$tau.IA),mean(exp(data.L$bhat.IA)),quantile(exp(data.L$bhat.IA),c(0.025)),quantile(exp(data.L$bhat.IA),c(0.975)),
                  mean(!is.na(data.L$m0.IA)),mean(data.L$m0.IA,na.rm=TRUE),quantile(data.L$m0.IA,c(0.025),na.rm=TRUE),quantile(data.L$m0.IA,c(0.975),na.rm=TRUE),
                  mean(!is.na(data.L$m1.IA)),mean(data.L$m1.IA,na.rm=TRUE),quantile(data.L$m1.IA,c(0.025),na.rm=TRUE),quantile(data.L$m1.IA,c(0.975),na.rm=TRUE))
      result<-rbind(result,temp)
      rm("temp")
    }
  }
  colnames(result)<-c("Look","d","p(cross)","min(d.obs)","max(d.obs)","min(t)","max(t)","Avg(hr)","hr(L)","hr(U)","%m0(est)","Avg(m0)","m0(L)","m0(U)","%m1(est)","Avg(m1)","m1(L)","m1(U)")
  if(details) print(result)
  return(result)
}



est.looks<-function(df,looks,d.looks,z.looks,hr.true){
  # First look
  look<-looks[1]
  df.look<-df$Out[[1]]
  tau.look<-df.look$tau.IA
  hr.look<-exp(df.look$bhat.IA)
  rej.look<-ifelse(df.look$Zlr.IA>=z.looks[1],1,0)
  hr.cross<-hr.look[which(rej.look==1)]
  
  df.hr<-data.frame(hr.cross)
  names(df.hr)<-c("hr")
  df.hr$look<-c(look)
  
  df.est<-data.frame(hr.look,tau.look)
  names(df.est)<-c("hr","tau")
  df.est$look<-c(look)
  
  for(ll in 2:length(looks)){
    look<-looks[ll]
    # Going to next look conditional on NOT crossing at first look
    crossings<-OC.properties.multiple(d.looks=d.looks[c(1:ll)],
                                      alpha.looks=alpha.looks[c(1:ll)],data.looks=df,details=FALSE)
    df.look<-crossings$data.stopped
    # Includes look1 sequences that stopped at look1;
    # Look data --> sequences that did NOT stop at last look
    df.look<-df.look[which(df.look$look==ll),]
    tau.look<-df.look$tau.IA
    hr.look<-exp(df.look$bhat.IA)
    rej.look<-ifelse(df.look$Zlr.IA>=z.looks[ll],1,0)
    hr.cross<-hr.look[which(rej.look==1)]
    
    df.hr.look<-data.frame(hr.cross)
    names(df.hr.look)<-c("hr")
    df.hr.look$look<-c(look)
    
    df.est.look<-data.frame(hr.look,tau.look)
    names(df.est.look)<-c("hr","tau")
    df.est.look$look<-c(look)
    
    df.hr<-rbind(df.hr,df.hr.look)
    df.est<-rbind(df.est,df.est.look)
  }
  df.hr$hr.bias<-df.hr$hr-c(hr.true)
  return(list(df.est=df.est,df.hr=df.hr))
}

