
get.accrual<-function(N,rand.ratio=1,Accrue.type="EastProjection",Ramp.up,nRamp.after){
  N0<-round(N*(1/(1+rand.ratio)))
  N1<-N-N0
  
  Ramp.up<-as.integer(Ramp.up) # Convert to integer
  N<-N1+N0
  n.left<-N-sum(Ramp.up)
  N.ramp<-c(Ramp.up,rep(nRamp.after,floor(n.left/nRamp.after)))
  
  N.k<-c(N.ramp,N-sum(N.ramp))
  
  N0.k<-floor(N.k/(rand.ratio+1))
  N1.k<-N.k-N0.k
  
  # Intervals where N1>N0
  dk<-which(N1.k-N0.k==rand.ratio)
  # Take from N1.k to N0.k for every other dk
  for(nn in 1:length(dk)){
    if(nn%%(rand.ratio+1)==1){
      N1.k[dk[nn]]<-N1.k[dk[nn]]-1
      N0.k[dk[nn]]<-N0.k[dk[nn]]+1
    }
  }
  
  AC.obs<-NULL
  nk<-N1.k
  tau.k<-c(0,seq(1,length(nk)))
  for(k in 1:(length(tau.k)-1)){
    #cat("# index,n.k=",c(tau.k[k],tau.k[k+1],nk[k]),"\n")
    if(Accrue.type=="projection") AC<-runif(nk[k],tau.k[k],tau.k[k+1])
    if(Accrue.type=="EastProjection") AC<-rep((tau.k[k]+(tau.k[k+1]-tau.k[k])/2),nk[k])
    # If time scale is in months
    # Entry time in months:  nk subjects enter the study (each arm) uniformly between (tau[k],tau[k+1])
    # That is, for nk there will be nk entry times simulated from U(tau[k],tau[k+1])
    AC.obs<-c(AC.obs,AC)
  }
  AC1.obs<-AC.obs
  
  
  AC.obs<-NULL
  nk<-N0.k
  tau.k<-c(0,seq(1,length(nk)))
  for(k in 1:(length(tau.k)-1)){
    #cat("# index,n.k=",c(tau.k[k],tau.k[k+1],nk[k]),"\n")
    if(Accrue.type=="projection") AC<-runif(nk[k],tau.k[k],tau.k[k+1])
    if (Accrue.type=="EastProjection") AC<-rep((tau.k[k]+(tau.k[k+1]-tau.k[k])/2),nk[k])
    # Entry time in months:  nk subjects enter the study (each arm) uniformly between (tau[k],tau[k+1])
    # That is, for nk there will be nk entry times simulated from U(tau[k],tau[k+1])
    #cat("AC =",c(AC),"\n")
    #print(length(AC))
    AC.obs<-c(AC.obs,AC)
  }
  AC0.obs<-AC.obs
  
  #if(length(AC1.obs)!=N1) stop("Not all experimental subjects assigned accrual time")
  #if(length(AC0.obs)!=N0) stop("Not all control subjects assigned accrual time")
  
  if(length(AC1.obs)!=N1){
    n1.left<-N1-length(AC1.obs)
    k.last<-length(tau.k)-1
    k<-k.last
    if(Accrue.type=="projection") AC1.left<-runif(n1.left,tau.k[k],tau.k[k+1])
    if(Accrue.type=="EastProjection") AC1.left<-rep((tau.k[k]+(tau.k[k+1]-tau.k[k])/2),n1.left)
    AC1.obs<-c(AC1.obs,AC1.left)
    
    # Force randomization ratio
    n0.over<-length(AC0.obs)-N0
    if(n0.over>0){
      # Randomly exclude n0.over entry times
      id.exclude<-sample(c(1:length(AC0.obs)),n0.over,replace=FALSE)
      AC0.obs<-AC0.obs[-c(id.exclude)]
    }
  }   
  
  
  return(list(AC.0=AC0.obs,AC.1=AC1.obs))
  
}




show.status<-function(ss,t.start,sims){
  if(ss==10){
    # clean this up
    cat("First 10 sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
  }
  
  if(round(ss/sims,digits=4)==0.10){
    cat("10% sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
    t.left<-(est.final-t.sofar)
    cat("# Estimated mins left=",c(t.left),"\n")
  }
  
  if(round(ss/sims,digits=4)==0.25){
    cat("25% sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
    t.left<-(est.final-t.sofar)
    cat("# Estimated mins left=",c(t.left),"\n")
  }
  
  
  if(round(ss/sims,digits=4)==0.50){
    cat("50% sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
    t.left<-(est.final-t.sofar)
    cat("# Estimated mins left=",c(t.left),"\n")
  }
  
  
  if(round(ss/sims,digits=4)==0.75){
    cat("75% sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
    t.left<-(est.final-t.sofar)
    cat("# Estimated mins left=",c(t.left),"\n")
  }
  
  
  if(round(ss/sims,digits=4)==0.90){
    cat("90% sims done","\n")
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    cat("Sims and time (mins) so far=",c(ss,t.sofar),"\n")
    
    est.final<-t.sofar*(sims/ss)
    cat("# Estimated mins=",c(est.final),"\n")
    t.left<-(est.final-t.sofar)
    cat("# Estimated mins left=",c(t.left),"\n")
  }
  if(ss==sims){
    t.now<-proc.time()[1]
    t.sofar<-(t.now-t.start)/60
    est.final<-t.sofar*(1000/sims)
    cat("# Estimated mins per 1000=",c(est.final),"\n")
  }
  }


get.ZWsim<-function(R,gamma1,gamma2,bW,lamW,lamP,bP,lamT,beta1,beta2,cw1=5,cw2=2){

# W ~ lamW0*exp(R*bW)
# Pj ~lamPj*exp(gamma1*x1+gamma2*x2+bPj*a(t,W)), j=0,1
# T ~ lamT*exp(beta1*R+beta2*a(t,P))  (P--> corresponding Pj)

n<-length(R)

#cat("N=",c(n),"\n")

    # x1 is baseline covariate;
    x1 <- rbinom(n, 1, 0.5)
    # x2 is baseline covariate;
    x2 <- rnorm(n)
uW<-runif(n)
W<--log(uW)/(lamW*exp(R*bW))

    uP <- runif(n)
    # P is the time to initiate post-study treatment, which depends on treatment, baseline covariates
    # x1, x2 and time-dependent confounder x3temp=W (occurence of x3temp=W delay initiation of P);    
    P <- NULL
    for (i in 1:n) {
            P[i]<- -log(uP[i])/(lamP*exp(gamma1*x1[i]+gamma2*x2[i]))
            if (P[i]>=W[i]){ 
                    if(W[i]>6) P[i]<-W[i]+(P[i]-W[i])/exp(bP) 
                    if(W[i]>=3 & W[i]<=6) P[i]<-W[i]+(P[i]-W[i])/(cw2*exp(bP))
                    if(W[i]<3) P[i]<-W[i]+(P[i]-W[i])/(cw1*exp(bP))                   
                           }
                else W[i]=NA; 
                  }
uT<-runif(n)
T <- -log(uT)/(lamT*exp(beta1*R))   
T.true<-T
    for (i in 1:n) {
        # post-study treatment prolongs survival
        if (T[i]>P[i]) T[i]<- P[i]+(T[i]-P[i])/exp(beta2)
    }
    
  # plot(T)    
return(list(T.true=T.true,T.obs=T,P.true=P))
}


sumZ.toleft<-function(x,Y,Z){
sum(Z[Y<=x])
}
sumZ.toright<-function(x,Y,Z){
sum(Z[Y>=x])
}
py.info<-function(x,Y){
sum(pmin(Y,x))
}
py.interim<-function(x,Y,A){
id<-(A<=x)
sum(pmin(Y[id],x))
}
events.interim<-function(x,D,Y,A){
id<-(A<=x & Y<=x)
sum(D[id])
}



##############################################################################
#               Functions for determining censoring distribution             #
# Want to find U(0,b) censoring distribution to yield (1-k)% censoring rate  #
##############################################################################
uniform.b<-function(k,lambda,tau){
F.tau<-1-exp(-lambda*tau)
c.tau<-(1/lambda)*((tau*lambda+1)*exp(-lambda*tau)-1)
b<-c.tau/((1-k)-F.tau)
return(b)
}

dH.uniform.ab<-function(s,a,b,lambda){
dFs<-lambda*exp(-lambda*s)
#G1<-0
#if(s>=a & s<=b) G1<-(b-s)/(b-a)
#if(s<a) G1<-1
#if(s>b) G1<-0.0
G1<-ifelse(s>=a & s<=b,(b-s)/(b-a),ifelse(s<a,1,0))
dH<-G1*dFs
return(dH)
}

uniform.ab<-function(parms,k,lambda,tau){
a<-parms[1]; b<-parms[2]
H<-integrate(dH.uniform.ab,a=a,b=b,lambda=lambda,lower=0,upper=tau)$value
error<-abs(H-(1-k))
return(error)
}

# Exponential censoring
dH.exp<-function(s,lambda.c,lambda){
dFs<-lambda*exp(-lambda*s)
G1<-1-exp(-lambda.c*s)
dH<-G1*dFs
return(dH)
}

exp.censoring<-function(lambda.c,k,lambda,tau){
H<-integrate(dH.exp,lambda.c=lambda.c,lower=0,upper=tau,lambda=lambda)$value
error<-(H-(1-k))
return(error)
}

#get.exp<-uniroot(f=exp.censoring,interval=c(0,10),lambda=1,k=0.2,tau=10)
#print(get.exp)
#lambdac.hat<-get.exp$root
#exp.censoring(lambda.c=lambdac.hat,k=0.2,lambda=1,tau=10)

#get.exp<-uniroot(f=exp.censoring,interval=c(-10,10),lambda=1/1.646912,k=0.1,tau=30/12)
#print(get.exp)
#lambdac.hat<-get.exp$root
#exp.censoring(lambda.c=lambdac.hat,k=0.1,lambda=1/1.646912,tau=30/12)


# test dH
#check<-dH.uniform.ab(s=3,a=0,b=7,lambda=1)
checking<-FALSE
if(checking){
H<-integrate(dH.uniform.ab,a=0,b=7,lambda=1,lower=0,upper=10)
Ks<-seq(0.05,0.5,by=0.05)
lambda<-2.5
for(kk in 1:length(Ks)){
# test
b.u<-uniform.b(k=Ks[kk],lambda=lambda,tau=10)
#[1] 9.999546
H<-integrate(dH.uniform.ab,a=0,b=b.u,lambda=lambda,lower=0,upper=10)$value
print(c(Ks[kk],1-H,b.u))
}
# check derived "uniform.b" algorithm
uniform.b.check<-function(b,a,k,lambda,tau){
H<-integrate(dH.uniform.ab,a=a,b=b,lambda=lambda,lower=0,upper=tau)$value
error<-H-(1-k)
#print(error)
return(error)
}
#get.uniform.ab<-optim(par=c(-10,10),fn=uniform.ab,lambda=1,k=0.2,tau=10)
#print(get.uniform.ab)
lambdas<-seq(0.5,5,length=25)
bhat.deriv<-bhat.num<-rep(NA,length(lambdas))
for(ll in 1:length(lambdas)){
lambda<-lambdas[ll]
get.uniform.b<-uniroot(f=uniform.b.check,interval=c(0.01,1000),lambda=lambda,a=0,k=0.3,tau=5)
bhat.num[ll]<-get.uniform.b$root
bhat.deriv[ll]<-uniform.b(k=0.3,lambda=lambda,tau=5)
}
plot(lambdas,bhat.deriv,type="l",lty=1,col="grey",lwd=2)
lines(lambdas,bhat.num,lty=2,col="black",lwd=0.25)
}

##############################################################################
#             End functions for determining censoring distribution           #
##############################################################################


Logrank<-function(time,Delta,X){
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; Delta<-Delta[id]; X<-X[id]
}
at.points<-sort(unique(c(time[Delta==1])))

U0<-time[which(X==0)]
D0<-Delta[which(X==0)]

# Control group
# Risk and Counting processes
risk.z0<-unlist(lapply(as.list(at.points),R.Weighted,error=U0))
counting<-unlist(lapply(as.list(at.points),N.Weighted,error=U0[D0==1]))
N.z0<-counting
dN.z0<- diff(c(0, counting))

U1<-time[which(X==1)]
D1<-Delta[which(X==1)]

# Control group
# Risk and Counting processes
risk.z1<-unlist(lapply(as.list(at.points),R.Weighted,error=U1))
counting<-unlist(lapply(as.list(at.points),N.Weighted,error=U1[D1==1]))
N.z1<-counting
dN.z1<- diff(c(0, counting))

dN.pooled<-dN.z0+dN.z1
risk.pooled<-risk.z0+risk.z1

K<-(risk.z0*risk.z1)/(risk.pooled)

term0<-sum(ifelse(risk.z0>0,(K/risk.z0)*dN.z0,0.0))
term1<-sum(ifelse(risk.z1>0,(K/risk.z1)*dN.z1,0.0))

lr<-term0-term1

# variance
h0<-ifelse(risk.z0==0,0,(K^2/risk.z0))
h1<-ifelse(risk.z1==0,0,(K^2/risk.z1))
dJ<-ifelse(risk.pooled==1,0,(dN.pooled-1)/(risk.pooled-1))
dL<-ifelse(risk.pooled==0,0,dN.pooled/risk.pooled)
sig2s<-(h0+h1)*(1-dJ)*dL
sig2<-sum(sig2s)

Z.lr<-lr/sqrt(sig2)
result<-list(Z.lr=Z.lr,lr=lr,pval=1-pchisq(Z.lr^2,1))
return(result)
}


get.Tdraw<-function(U,S,tpoints){
if(S[length(S)]>0){
S<-c(S,0)
t.mod<-c(tpoints,Inf)
}
else{
t.mod<-tpoints
}
min(t.mod[(1-S)>=U])
}

DrawFromS<-function(n,S.draw,tpoints.draw){
u.draw<-runif(n)
T.draw<-unlist(lapply(u.draw,get.Tdraw,S=S.draw,tpoints=tpoints.draw))
return(T.draw)
}


# Piecewise-exponential

require(Hmisc)
# From cpsurvsim
exp_cdfsim<-function (n, endtime=Inf, theta, tau = NA) 
{
  if (is.numeric(n) == FALSE | n <= 0) {
    stop("n must be an integer greater than 0.")
  }
  if (Hmisc::all.is.numeric(theta, "test") == FALSE | 
      all(theta > 0) == FALSE | is.numeric(endtime) == FALSE | 
      endtime <= 0) {
    stop("Endtime and theta must be numeric and > 0.")
  }
  if (length(tau) > 4) {
    stop("This function only allows for up to 4 change-points.")
  }
  n <- as.integer(n)
  x <- stats::rexp(n)
  if (is.na(tau[1]) == TRUE) {
    t <- x/theta
  }
  if (is.na(tau[1]) == FALSE) {
    if (Hmisc::all.is.numeric(tau, "test") == FALSE | 
        all(tau > 0) == FALSE) {
      stop("Tau must be numeric and > 0.")
    }
    if (endtime < tau[length(tau)]) {
      warning("Warning: Change-points occur after endtime.")
    }
    if (length(theta) != (length(tau) + 1)) {
      stop("Length of theta and tau not compatible.")
    }
  }
  if (length(tau) == 1 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau
    cdfcp1 <- function(v) {
      ifelse(v < first, v/theta[1], ((v - first)/theta[2]) + 
               tau)
    }
    t <- cdfcp1(x)
  }
  if (length(tau) == 2 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    cdfcp2 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ((v - second)/theta[3]) + 
                                             tau[2]))
    }
    t <- cdfcp2(x)
  }
  if (length(tau) == 3 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    cdfcp3 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ifelse(v < third, 
                                                                                   ((v - second)/theta[3]) + tau[2], ((v - third)/theta[4]) + 
                                                                                     tau[3])))
    }
    t <- cdfcp3(x)
  }
  if (length(tau) == 4 & is.na(tau[1]) == FALSE) {
    first <- theta[1] * tau[1]
    second <- first + theta[2] * (tau[2] - tau[1])
    third <- second + theta[3] * (tau[3] - tau[2])
    fourth <- third + theta[4] * (tau[4] - tau[3])
    cdfcp4 <- function(v) {
      ifelse(v < first, v/theta[1], ifelse(v < second, 
                                           ((v - first)/theta[2]) + tau[1], ifelse(v < third, 
                                                                                   ((v - second)/theta[3]) + tau[2], ifelse(v < 
                                                                                                                              fourth, ((v - third)/theta[4]) + tau[3], 
                                                                                                                            ((v - fourth)/theta[5]) + tau[4]))))
    }
    t <- cdfcp4(x)
  }
  endtime <- as.numeric(endtime)
  C <- rep(endtime, length(x))
  time <- pmin(t, C)
  censor <- as.numeric(time != endtime)
  dta <- data.frame(time = time, censor = censor)
  return(dta)
}


#df1<-exp_cdfsim(n=1000,theta=c(0.1,0.4),tau=c(6))
#df0<-exp_cdfsim(n=1000,theta=c(0.1,1.2),tau=c(6))
#df0$treat<-0
#df1$treat<-1
#dfa<-rbind(df1,df0)
#dfa$event<-1
#kma<-survfit(Surv(time,event)~treat,data=dfa)
#plot(kma)

rpwexp <- function(n, rate=1, intervals=NULL, cumulative=FALSE){
  if(is.null(intervals)){
    if (cumulative){return(cumsum(rexp(n,rate[1])))}else
      return(rexp(n,rate[1]))}
  k <- length(rate)
  if (k==1){
    if(cumulative){return(cumsum(rexp(n,rate)))}else
      return(rexp(n,rate))
  }
  if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
  tx <- 0
  j <- 1
  times <- array(0,n)
  timex <- cumsum(intervals)
  indx <- array(TRUE,n)
  for(i in 1:k){
    nindx <- sum(indx)
    if (nindx==0) break
    increment <- rexp(nindx,rate[i])
    if (cumulative) times[indx] <- tx + cumsum(increment)
    else times[indx] <- tx + increment
    if (i<k){
      tx <- timex[i]
      indx <- (times > timex[i])
    }
  }
  return(times)
}

#t1<-rpwexp(n=1000,rate=c(0.1,0.4),intervals=c(6))
#df1<-data.frame(time=t1)
#t0<-rpwexp(n=1000,rate=c(0.1,1.2),intervals=c(6))
#df0<-data.frame(time=t0)
#df0$treat<-0
#df1$treat<-1
#dfa<-rbind(df1,df0)
#dfa$event<-1
#kma<-survfit(Surv(time,event)~treat,data=dfa)
#plot(kma)




