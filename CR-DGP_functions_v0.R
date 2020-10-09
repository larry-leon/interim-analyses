

hweibull<-function (x, shape, scale = 1, log = FALSE) 
{
  if (any(shape <= 0) || any(scale <= 0)) 
    stop("scale and shape must be positive")
  res <- ifelse(x < 0, 0, shape * (x/scale)^(shape - 1)/scale)
  if (log) 
    res <- log(res)
  return(res)
}

Hweibull<-function (x, shape, scale = 1, log.p = FALSE) 
{
  if (any(shape <= 0) || any(scale <= 0)) 
    stop("scale and shape must be positive")
  res <- ifelse(x < 0, 0, (x/scale)^shape)
  if (log.p) 
    res <- log(res)
  return(res)
}


f1.sub<-function(x,shape1,scale1,shape2,scale2){
  lam1.x<-hweibull(x,shape=shape1,scale=scale1)
  lam2.x<-hweibull(x,shape=shape2,scale=scale2)
  G1.x<-Hweibull(x,shape=shape1,scale=scale1)
  G2.x<-Hweibull(x,shape=shape2,scale=scale2)
  S.x<-exp(-G1.x)*exp(-G2.x)
  f1x.sub<-S.x*lam1.x
  return(f1x.sub)
}

f2.sub<-function(x,shape1,scale1,shape2,scale2){
  lam1.x<-hweibull(x,shape=shape1,scale=scale1)
  lam2.x<-hweibull(x,shape=shape2,scale=scale2)
  G1.x<-Hweibull(x,shape=shape1,scale=scale1)
  G2.x<-Hweibull(x,shape=shape2,scale=scale2)
  S.x<-exp(-G1.x)*exp(-G2.x)
  f2x.sub<-S.x*lam2.x
  return(f2x.sub)
}

# A piece-wise CIF target for deaths
# 7 points
# increase to 10 points
cif.death<-function(tpoints,tcuts=seq(1,29,length=15),jcuts=rep(0,15),plotit=FALSE){
cif<-ifelse(tpoints<=tcuts[1],jcuts[1],0)+ifelse(tpoints>tcuts[1] & tpoints<=tcuts[2],jcuts[2],0)+ifelse(tpoints>tcuts[2] & tpoints<=tcuts[3],jcuts[3],0)+ifelse(tpoints>tcuts[3] & tpoints<=tcuts[4],jcuts[4],0)+
  ifelse(tpoints>tcuts[4] & tpoints<=tcuts[5],jcuts[5],0)+ifelse(tpoints>tcuts[5] & tpoints<=tcuts[6],jcuts[6],0)+ifelse(tpoints>tcuts[6] & tpoints<=tcuts[7],jcuts[7],0)+
  ifelse(tpoints>tcuts[7] & tpoints<=tcuts[8],jcuts[8],0)+ifelse(tpoints>tcuts[8] & tpoints<=tcuts[9],jcuts[9],0)+ifelse(tpoints>tcuts[9] & tpoints<=tcuts[10],jcuts[10],0)+
  ifelse(tpoints>tcuts[10] & tpoints<=tcuts[11],jcuts[11],0)+ifelse(tpoints>tcuts[11] & tpoints<=tcuts[12],jcuts[12],0)+ifelse(tpoints>tcuts[12] & tpoints<=tcuts[13],jcuts[13],0)+
  ifelse(tpoints>tcuts[13] & tpoints<=tcuts[14],jcuts[14],0)+ifelse(tpoints>tcuts[14] & tpoints<=tcuts[15],jcuts[15],0)

if(plotit) plot(tpoints,cif,type="s",lty=1,col="grey",lwd=3)
return(cif)
}







obj1.sub_parms<-function(theta,tpoints,target1,target2,wt1=1.0,print=FALSE){
  shape1<-theta[1]
  scale1<-theta[2]
  shape2<-theta[3]
  scale2<-theta[4]
  #atpoints<-c(7,14,21,29)
  atpoints<-tpoints
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    # F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  #error1<-sum(c(F1.sub-target1)^2)
  diff<-F1.sub-target1
  #error1<-diff[which(atpoints==7)]^2+diff[which(atpoints==14)]^2+diff[which(atpoints==21)]^2+diff[which(atpoints==29)]^2
  #error1<-diff[which(atpoints==14)]^2+diff[which(atpoints==21)]^2+diff[which(atpoints==29)]^2
  error1<-0.25*c(diff[which(atpoints==14)]^2)+0.75*c(diff[which(atpoints==29)]^2)
  #error1<-diff[which(atpoints==29)]^2
  #error2<-c(F2.sub[length(tpoints)]-target2)^2
  #lam1<-hweibull(x=c(tpoints),shape=shape1,scale=scale1)
  #lam2<-hweibull(x=c(tpoints),shape=shape2,scale=scale2)
  #p2<-sum(ifelse(lam1+lam2==0,0,lam2/(lam1+lam2)))
  #error2<-(p2-target2)^2
  #error<-error1+error2
  #if(print) print(c(p2,error2))
  return(error1)
}

obj2.sub_parms<-function(theta,tpoints,target1,target2,wt1=1.0,print=FALSE){
  shape1<-theta[1]
  scale1<-theta[2]
  shape2<-theta[3]
  scale2<-theta[4]
  #atpoints<-c(7,14,21,29)
  atpoints<-tpoints
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    # F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  #error1<-sum(c(F1.sub-target1)^2)
  diff<-F1.sub-target1
  #error1<-diff[which(atpoints==7)]^2+diff[which(atpoints==14)]^2+diff[which(atpoints==21)]^2+diff[which(atpoints==29)]^2
  #error1<-diff[which(atpoints==14)]^2+diff[which(atpoints==21)]^2+diff[which(atpoints==29)]^2
  
  error1<-0.10*c(diff[which(atpoints==5)]^2)+0.15*c(diff[which(atpoints==14)]^2)+0.75*c(diff[which(atpoints==29)]^2)
  
  #error1<-diff[which(atpoints==29)]^2
  #error2<-0.1*c(F2.sub[length(tpoints)]-target2)^2
  #lam1<-hweibull(x=c(tpoints),shape=shape1,scale=scale1)
  #lam2<-hweibull(x=c(tpoints),shape=shape2,scale=scale2)
  #p2<-sum(ifelse(lam1+lam2==0,0,lam2/(lam1+lam2)))
  #error2<-(p2-target2)^2
  #error<-error1+error2
  #if(print) print(c(p2,error2))
  return(error1)
}


obj3.sub_parms<-function(theta,tpoints,target1,target2,wt1=1.0,shape1,scale1){
  shape2<-theta[1]
  scale2<-theta[2]
  
  #atpoints<-c(7,14,21,29)
  atpoints<-tpoints
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  
  diff<-F1.sub-target1
  
  #error1<-mean(diff^2)
  
  #error1<-max(abs(diff))
  
  #error1<-diff[which(atpoints==7)]^2+diff[which(atpoints==14)]^2+diff[which(atpoints==21)]^2+diff[which(atpoints==29)]^2
  
  error1<-diff[which(atpoints==14)]^2+diff[which(atpoints==29)]^2
  
  #error1<-0.25*c(diff[which(atpoints==29)]^2)
  
  #error1<-diff[which(atpoints==29)]^2
  
  #error2<-0.75*c((F2.sub[length(tpoints)]-target2)^2)
  
  #lam1<-hweibull(x=c(tpoints),shape=shape1,scale=scale1)
  #lam2<-hweibull(x=c(tpoints),shape=shape2,scale=scale2)
  #p2<-sum(ifelse(lam1+lam2==0,0,lam2/(lam1+lam2)))
  #error3<-(p2-target2)^2
  
  #error<-error1+error2
  #if(print) print(c(p2,error2))
  return(error1)
}



obj4.sub_parms<-function(theta,tpoints,target1,target2,k1=0,k2=0,error.type="CvM"){
  # FH(0,k) type of weighting
  wt1<-target1^k1
  wt2<-target2^k2
  
  wt1<-wt1/sum(wt1)
  wt2<-wt2/sum(wt2)
  
  shape1<-theta[1]
  scale1<-theta[2]
  shape2<-theta[3]
  scale2<-theta[4]
  
  atpoints<-tpoints
  
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  # Give more weight to later time points
  diff1<-c(wt1*(F1.sub-target1))
  diff2<-c(wt2*(F2.sub-target2))
  
  if(error.type=="CvM"){
  error1<-mean(diff1^2)
  error2<-mean(diff2^2)
  }
  
  if(error.type=="KS"){
    error1<-max(abs(diff1))
    error2<-max(abs(diff2))
  }
  
  error<-error1+error2
  return(error)
}


obj5.sub_parms<-function(theta,tpoints,target1,target2,k1=0,k2=0,error.type="CvM"){
  # FH(0,k) type of weighting
  wt1<-target1^k1
  wt2<-target2^k2
  
  wt1<-wt1/sum(wt1)
  wt2<-wt2/sum(wt2)
  # common shape
  shape1<-shape2<-theta[1]
  scale1<-theta[2]
  scale2<-theta[3]
  
  atpoints<-tpoints
  
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  # Give more weight to later time points
  diff1<-c(wt1*(F1.sub-target1))
  diff2<-c(wt2*(F2.sub-target2))
  
  if(error.type=="CvM"){
    error1<-mean(diff1^2)
    error2<-mean(diff2^2)
  }
  
  if(error.type=="KS"){
    error1<-max(abs(diff1))
    error2<-max(abs(diff2))
  }
  
  error<-error1+error2
 
  
  return(error)
  
  
}



obj6.sub_parms<-function(theta,shape.fix=1,tpoints,target1,target2,k1=0,k2=0,error.type="CvM"){
  # FH(0,k) type of weighting
  wt1<-target1^k1
  wt2<-target2^k2
  
  wt1<-wt1/sum(wt1)
  wt2<-wt2/sum(wt2)
  # common shape
  shape1<-shape2<-shape.fix
  scale1<-theta[1]
  scale2<-theta[2]
  
  atpoints<-tpoints
  
  F1.sub<-rep(NA,length(atpoints))
  F2.sub<-rep(NA,length(atpoints))
  for(i in 1:length(atpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=atpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  # Give more weight to later time points
  diff1<-c(wt1*(F1.sub-target1))
  diff2<-c(wt2*(F2.sub-target2))
  
  if(error.type=="CvM"){
    error1<-mean(diff1^2)
    error2<-mean(diff2^2)
  }
  
  if(error.type=="KS"){
    error1<-max(abs(diff1))
    error2<-max(abs(diff2))
  }
  
  error<-error1+error2
  
  
  return(error)
  
  
}



find.CIF.match<-function(tpoints,target1,target2,p.start=c(1,15,2,20),plotit=TRUE,k1=0,k2=0,error.type="CvM"){
  
  temp<-nlm(f=obj4.sub_parms,p=p.start,tpoints=tpoints,target1=target1,target2=target2,k1=k1,k2=k2,error.type=error.type)
  est<-temp$estimate
  shape1<-est[1]
  scale1<-est[2]
  shape2<-est[3]
  scale2<-est[4]
  
  F1.sub<-rep(NA,length(tpoints))
  F2.sub<-rep(NA,length(tpoints))
  for(i in 1:length(tpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  ymin<-min(c(F1.sub,F2.sub,target1,target2))
  ymax<-max(c(F1.sub,F2.sub,target1,target2))
  
  plot(tpoints,F1.sub,type="s",lty=1,col="blue",lwd=3,ylim=c(ymin,ymax))
  lines(tpoints,F2.sub,lty=2,col="grey",lwd=3)
  lines(tpoints,target1,col="red",lty=1,lwd=3)
  lines(tpoints,target2,col="red",lty=2,lwd=3)
  
  names(est)<-c("shape1","scale1","shape2","scale2")
  
  # cause-specific
  lam1.cs<-hweibull(tpoints,shape=shape1,scale=scale1)
  lam2.cs<-hweibull(tpoints,shape=shape2,scale=scale2)
  # overall S(.)
  G1<-Hweibull(tpoints,shape=shape1,scale=scale1)
  G2<-Hweibull(tpoints,shape=shape2,scale=scale2)
  S.overall<-exp(-G1)*exp(-G2)
  
  r1.sub<-ifelse(F1.sub<1,S.overall/(1-F1.sub),0.0)
  lam1.sub<-r1.sub*lam1.cs
  
  r2.sub<-ifelse(F2.sub<1,S.overall/(1-F2.sub),0.0)
  lam2.sub<-r2.sub*lam2.cs
  
  lam12.cs<-lam1.cs+lam2.cs
  # Probability of event time corresponding to cause-1
  q1<-ifelse(lam12.cs==0,0,lam1.cs/lam12.cs)
  q2<-1-q1
  
  
  #cat("Max target 1 and target 2",c(max(target1),max(target2)),"\n")
  #cat("Max F1 and F2 subs",c(max(F1.sub),max(F2.sub)),"\n")
  #cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  merror2<-max(F2.sub)-max(target2)
  if(merror2>=-0.1){
  cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  }
  
  out<-list(est=est,F1.sub=F1.sub,F2.sub=F2.sub,S.overall=S.overall,lam1.sub=lam1.sub,lam2.sub=lam2.sub,tpoints=tpoints,target1=target1,target2=target2,q1=q1,q2=q2)
  return(out)
}


find.CIF.match2<-function(tpoints,target1,target2,p.start=c(1,5,5),plotit=TRUE,k1=0,k2=0,error.type="CvM",opt.type="optimx"){
set.seed(8316951)  

if(opt.type=="optimx"){
get.theta<-optimx(par=p.start,fn=obj5.sub_parms,tpoints=tpoints,target1=target1,target2=target2,k1=k1,k2=k2,error.type=error.type)
# Which solution is best minimizer
locs<-which(get.theta$value==min(get.theta$value))
est<-as.numeric(get.theta[locs,c("p1","p2","p3")])
}

if(opt.type=="nlm"){
get.theta<-nlm(f=obj5.sub_parms,p=p.start,tpoints=tpoints,target1=target1,target2=target2,k1=k1,k2=k2,error.type=error.type)
est<-get.theta$estimate
}  
  
# common shape
shape1<-shape2<-c(est[1])
scale1<-c(est[2])
scale2<-c(est[3])
  
  F1.sub<-rep(NA,length(tpoints))
  F2.sub<-rep(NA,length(tpoints))
  for(i in 1:length(tpoints)){
  F1.sub[i]<-integrate(f1.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  F2.sub[i]<-integrate(f2.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  ymin<-min(c(F1.sub,F2.sub,target1,target2))
  ymax<-max(c(F1.sub,F2.sub,target1,target2))
  
  plot(tpoints,F1.sub,type="s",lty=1,col="blue",lwd=3,ylim=c(ymin,ymax))
  lines(tpoints,F2.sub,lty=2,col="grey",lwd=3)
  lines(tpoints,target1,col="red",lty=1,lwd=3)
  lines(tpoints,target2,col="red",lty=2,lwd=3)
  
  est.out<-c(shape1,scale1,shape2,scale2)
  
  names(est.out)<-c("shape1","scale1","shape2","scale2")
  
  # cause-specific
  lam1.cs<-hweibull(tpoints,shape=shape1,scale=scale1)
  lam2.cs<-hweibull(tpoints,shape=shape2,scale=scale2)
  # overall S(.)
  G1<-Hweibull(tpoints,shape=shape1,scale=scale1)
  G2<-Hweibull(tpoints,shape=shape2,scale=scale2)
  S.overall<-exp(-G1)*exp(-G2)
  
  r1.sub<-ifelse(F1.sub<1,S.overall/(1-F1.sub),0.0)
  lam1.sub<-r1.sub*lam1.cs
  
  r2.sub<-ifelse(F2.sub<1,S.overall/(1-F2.sub),0.0)
  lam2.sub<-r2.sub*lam2.cs
  
  lam12.cs<-lam1.cs+lam2.cs
  # Probability of event time corresponding to cause-1
  q1<-ifelse(lam12.cs==0,0,lam1.cs/lam12.cs)
  q2<-1-q1
  
  
  #cat("Max target 1 and target 2",c(max(target1),max(target2)),"\n")
  #cat("Max F1 and F2 subs",c(max(F1.sub),max(F2.sub)),"\n")
  #cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  merror2<-max(F2.sub)-max(target2)
  if(merror2>=-0.1){
    cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  }
  
  out<-list(est=est.out,F1.sub=F1.sub,F2.sub=F2.sub,S.overall=S.overall,lam1.sub=lam1.sub,lam2.sub=lam2.sub,tpoints=tpoints,target1=target1,target2=target2,q1=q1,q2=q2)
  return(out)
}




find.CIF.match3<-function(tpoints,target1,target2,p.start=c(5,5),plotit=TRUE,k1=0,k2=0,error.type="CvM",opt.type="optimx",shape.fix=1){
  set.seed(8316951)  
  
  if(opt.type=="optimx"){
    get.theta<-optimx(par=p.start,fn=obj6.sub_parms,tpoints=tpoints,target1=target1,target2=target2,k1=k1,k2=k2,error.type=error.type,shape.fix=shape.fix)
    # Which solution is best minimizer
    locs<-which(get.theta$value==min(get.theta$value))
    est<-as.numeric(get.theta[locs,c("p1","p2")])
  }
  
  if(opt.type=="nlm"){
    get.theta<-nlm(f=obj6.sub_parms,p=p.start,tpoints=tpoints,target1=target1,target2=target2,k1=k1,k2=k2,error.type=error.type,shape.fix=shape.fix)
    est<-get.theta$estimate
  }  
  
  # common shape
  shape1<-shape2<-shape.fix
  scale1<-c(est[1])
  scale2<-c(est[2])
  
  F1.sub<-rep(NA,length(tpoints))
  F2.sub<-rep(NA,length(tpoints))
  for(i in 1:length(tpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2)$value
  }
  
  ymin<-min(c(F1.sub,F2.sub,target1,target2))
  ymax<-max(c(F1.sub,F2.sub,target1,target2))
  
  plot(tpoints,F1.sub,type="s",lty=1,col="blue",lwd=3,ylim=c(ymin,ymax))
  lines(tpoints,F2.sub,lty=2,col="grey",lwd=3)
  lines(tpoints,target1,col="red",lty=1,lwd=3)
  lines(tpoints,target2,col="red",lty=2,lwd=3)
  
  est.out<-c(shape1,scale1,shape2,scale2)
  
  names(est.out)<-c("shape1","scale1","shape2","scale2")
  
  # cause-specific
  lam1.cs<-hweibull(tpoints,shape=shape1,scale=scale1)
  lam2.cs<-hweibull(tpoints,shape=shape2,scale=scale2)
  # overall S(.)
  G1<-Hweibull(tpoints,shape=shape1,scale=scale1)
  G2<-Hweibull(tpoints,shape=shape2,scale=scale2)
  S.overall<-exp(-G1)*exp(-G2)
  
  r1.sub<-ifelse(F1.sub<1,S.overall/(1-F1.sub),0.0)
  lam1.sub<-r1.sub*lam1.cs
  
  r2.sub<-ifelse(F2.sub<1,S.overall/(1-F2.sub),0.0)
  lam2.sub<-r2.sub*lam2.cs
  
  lam12.cs<-lam1.cs+lam2.cs
  # Probability of event time corresponding to cause-1
  q1<-ifelse(lam12.cs==0,0,lam1.cs/lam12.cs)
  q2<-1-q1
  
  
  #cat("Max target 1 and target 2",c(max(target1),max(target2)),"\n")
  #cat("Max F1 and F2 subs",c(max(F1.sub),max(F2.sub)),"\n")
  #cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  merror2<-max(F2.sub)-max(target2)
  if(merror2>=-0.1){
    cat("parameters=",c(est.out),"\n")
    cat("F2-target",c(max(F2.sub)-max(target2)),"\n")
  }
  
  out<-list(est=est.out,F1.sub=F1.sub,F2.sub=F2.sub,S.overall=S.overall,lam1.sub=lam1.sub,lam2.sub=lam2.sub,tpoints=tpoints,target1=target1,target2=target2,q1=q1,q2=q2)
  return(out)
}



get.CIF<-function(thetas,tpoints,plotit=FALSE){
  
  shape1<-thetas[1]
  scale1<-thetas[2]
  shape2<-thetas[3]
  scale2<-thetas[4]
  
  F1.sub<-rep(NA,length(tpoints))
  F2.sub<-rep(NA,length(tpoints))
  for(i in 1:length(tpoints)){
    F1.sub[i]<-integrate(f1.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2,stop.on.error=FALSE)$value
    F2.sub[i]<-integrate(f2.sub,lower=0,upper=tpoints[i],shape1=shape1,scale1=scale1,shape2=shape2,scale2=scale2,stop.on.error=FALSE)$value
  }
  
  if(plotit){
    ymin<-min(c(F1.sub,F2.sub))
    ymax<-max(c(F1.sub,F2.sub))
    
    plot(tpoints,F1.sub,type="s",lty=1,col="blue",lwd=3,ylim=c(ymin,ymax))
    lines(tpoints,F2.sub,lty=2,col="grey",lwd=3)
    lines(tpoints,cdf0,col="red",lty=1,lwd=3)
    lines(tpoints,cdf2,col="red",lty=2,lwd=3)
  }
  
  # cause-specific
  lam1.cs<-hweibull(tpoints,shape=shape1,scale=scale1)
  lam2.cs<-hweibull(tpoints,shape=shape2,scale=scale2)
  # overall S(.)
  G1<-Hweibull(tpoints,shape=shape1,scale=scale1)
  G2<-Hweibull(tpoints,shape=shape2,scale=scale2)
  S.overall<-exp(-G1)*exp(-G2)
  
  r1.sub<-ifelse(F1.sub<1,S.overall/(1-F1.sub),0.0)
  lam1.sub<-r1.sub*lam1.cs
  
  r2.sub<-ifelse(F2.sub<1,S.overall/(1-F2.sub),0.0)
  lam2.sub<-r2.sub*lam2.cs
  
  lam12.cs<-lam1.cs+lam2.cs
  # Probability of event time corresponding to cause-1
  q1<-ifelse(lam12.cs==0,0,lam1.cs/lam12.cs)
  q2<-1-q1
  
  out<-list(est=thetas,F1.sub=F1.sub,F2.sub=F2.sub,S.overall=S.overall,lam1.sub=lam1.sub,lam2.sub=lam2.sub,tpoints=tpoints,q1=q1,q2=q2,G1=G1,G2=G2)
  return(out)
}

