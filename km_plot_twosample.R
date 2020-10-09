


###################################################################################################################################################
# K-M is weighted via Weight (Weight==1 is standard KM)
# Cox model is either via Weighting OR Covariates (if Covariates = NULL, then Cox model is weighted;  Otherwise, Cox model included Covariates)
###################################################################################################################################################

KM.plot.2sample.weighted<-function(Y,E,Treat,Covariates=NULL,Weight=rep(1,length(Y)),show.cox=FALSE,col.cox="red",cox.cex=1.0,show.logrank=FALSE,stop.onerror=FALSE,check.KM=FALSE,
details=FALSE,put.legend.cox='topright',put.legend.lr='top',del.med=0.075,lr.digits=2,
direction="any",
tpoints.add=c(-1),by.risk=NULL,Xlab="time",Ylab="proportion surviving",col.1="black",col.2="blue",show.med=FALSE,med.cex=1.25,risk.cex=1,
quant=0.5,
censor.mark.all=TRUE,show.ticks=TRUE,risk.set=TRUE,ymin=-0.125,ymax=1,
add.segment=FALSE,risk.add=NULL,xmin=0,xmax=NULL,x.truncate=NULL,time.zero=0.0,prob.points=c(seq(0.1,1.0,by=0.1))){

#if(censor.mark.all==TRUE) warning("All censorings will be marked regardless of *tying to an event* [DOES NOT REPLICATE R Survfit()]")

Wgt<-Weight
tpoints.add.plot<-tpoints.add # This will be for plotting individual curves
# Common points to estimate both arms to include all event times (for differences)
tpoints.add<-sort(unique(c(tpoints.add,sort(unique(Y[which(E==1)])))))

Y1<-Y[which(Treat==1)]
E1<-E[which(Treat==1)]
W1<-Weight[which(Treat==1)]

Y0<-Y[which(Treat==0)]
E0<-E[which(Treat==0)]
W0<-Weight[which(Treat==0)]

Time<-c(Y1,Y0)
Event<-c(E1,E0)
Strata<-c(rep(2,length(Y1)),rep(1,length(Y0)))
Weight<-c(W1,W0)


if(!is.null(Covariates)){ 
coxfit<-coxph(Surv(Y,E)~Treat+Covariates,robust=TRUE)  
hr.ci<-as.matrix(round(summary(coxfit)$conf.int,4))
pval.cox<-summary(coxfit)$coefficients[1,6]
}
if(is.null(Covariates)){ 
coxfit<-coxph(Surv(Y,E)~Treat,robust=TRUE)  
hr.ci<-as.matrix(round(summary(coxfit)$conf.int,4))
pval.cox<-summary(coxfit)$coefficients[6]
}
# ipw weighting
if(is.null(Covariates) & !all(Wgt==1)){ 
coxfit<-coxph(Surv(Y,E)~Treat,weights=Wgt,robust=TRUE)
hr.ci<-as.matrix(round(summary(coxfit)$conf.int,4))
pval.cox<-summary(coxfit)$coefficients[6]
}

# Get HR for Treat
hr.ci<-hr.ci[1,c(1,3,4)]

if(details){
if(!is.null(Covariates)){
cat("Cox adjusted HR=",c(paste(hr.ci[1])),"\n")
cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}


if(is.null(Covariates) & all(Wgt==1)){
cat("Cox un-adjusted HR=",c(paste(hr.ci[1])),"\n")
cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}

if(is.null(Covariates) & !all(Wgt==1)){
cat("Cox ipw-adjusted HR=",c(paste(hr.ci[1])),"\n")
cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")}
}


if(is.null(by.risk)){
tt<-sort(unique(Time))
by.risk<-quantile(tt-c(0,tt[-length(tt)]),c(0.75),na.rm=TRUE)
}

stratums<-c(sort(unique(Strata),decreasing=TRUE))
if(any(is.na(Strata))) stop("Missing strata information")

x1<-Time[which(Strata==stratums[1])]
to.get<-!is.na(x1)
x1<-x1[to.get]
eta1<-Event[which(Strata==stratums[1])]
eta1<-eta1[to.get]
w1<-Weight[which(Strata==stratums[1])]
w1<-w1[to.get]

x2<-Time[which(Strata==stratums[2])]
to.get<-!is.na(x2)
x2<-x2[to.get]
eta2<-Event[which(Strata==stratums[2])]
eta2<-eta2[to.get]
w2<-Weight[which(Strata==stratums[2])]
w2<-w2[to.get]

xx<-c(x1,x2)
ee<-c(eta1,eta2)

at.points<-sort(c(unique(c(x1,x2,tpoints.add))))
at.points1<-sort(c(unique(c(x1,tpoints.add))))
at.points2<-sort(c(unique(c(x2,tpoints.add))))

at.points1.plot<-sort(c(unique(c(x1,tpoints.add.plot))))
at.points2.plot<-sort(c(unique(c(x2,tpoints.add.plot))))


# Experimental
fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1)
S1.KM<-fit1$S.KM
SE1.KM<-fit1$se.KM

# Control
fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2)
S2.KM<-fit2$S.KM
SE2.KM<-fit2$se.KM

# LR() denotes the log-rank scale
# Integral of "LR differences" = log-rank test (non-standardized)
ZLR.dpoints<-LR.dpoints<-S1.dpoints<-SE1.dpoints<-S2.dpoints<-SE2.dpoints<-NULL

tp1.match<-match(tpoints.add,at.points1)
tp2.match<-match(tpoints.add,at.points2)

S1.dpoints<-S1.KM[tp1.match]
SE1.dpoints<-SE1.KM[tp1.match]
Y1.dpoints<-fit1$n.risk[tp1.match]
dN1.dpoints<-fit1$dN[tp1.match]

S2.dpoints<-S2.KM[tp2.match]
SE2.dpoints<-SE2.KM[tp2.match]
Y2.dpoints<-fit2$n.risk[tp2.match]
dN2.dpoints<-fit2$dN[tp2.match]

K<-(Y1.dpoints*Y2.dpoints)/(Y1.dpoints+Y2.dpoints)

Y.dpoints<-Y1.dpoints+Y2.dpoints
dN.dpoints<-dN1.dpoints+dN2.dpoints

term1<-cumsum(ifelse(Y1.dpoints>0,(K/Y1.dpoints)*dN1.dpoints,0.0))
term2<-cumsum(ifelse(Y2.dpoints>0,(K/Y2.dpoints)*dN2.dpoints,0.0))

# variance
h0<-ifelse(Y1.dpoints==0,0,(K^2/Y1.dpoints))
h1<-ifelse(Y2.dpoints==0,0,(K^2/Y2.dpoints))
dJ<-ifelse(Y.dpoints==1,0,(dN.dpoints-1)/(Y.dpoints-1))
dL<-ifelse(Y.dpoints==0,0,dN.dpoints/Y.dpoints)
sig2s<-(h0+h1)*(1-dJ)*dL
sig2<-cumsum(sig2s)

LR.dpoints<-term2-term1
ZLR.dpoints<-LR.dpoints/sqrt(sig2)
# Zlr() is weighted log-rank process
# Zlr(tau) for tau end of followup is 
# the "final" z statistic (1-sided)
z.lr<-ZLR.dpoints[length(ZLR.dpoints)]
if(direction=="superiority") p.lr<-1-pnorm(z.lr)
if(direction=="any") p.lr<-2*(1-pnorm(abs(z.lr)))


if(check.KM & all(Wgt==1)){
p.2side<-2*p.lr
logrnk<-survdiff(Surv(Y,E)~Treat)
p.val<-1-pchisq(logrnk$chisq,length(logrnk$n) - 1)
#print(c(p.2side,p.val))
if(round(p.2side-p.val,6)>0) warning("Discrepancy with log-rank process and survdiff")
}


#########################################
# Note: want to confirm with R function
#########################################
#           Check KM fits               #
if(check.KM){

# only at observed unique observations to align with survfit
# Experimental
S1.check<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=sort(unique(x1)))$S.KM
S2.check<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=sort(unique(x2)))$S.KM

KM.fit<-survfit(Surv(Y1,E1)~1,weights=W1)

max.error<-round(max(abs(KM.fit$surv-S1.check)),8)
print(max.error)
if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM1 fit")

KM.fit<-survfit(Surv(Y0,E0)~1,weights=W0)
max.error<-round(max(abs(KM.fit$surv-S2.check)),8)
print(max.error)
if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM2 fit")
#        End Check KM fits               #
}

events1<-sort(unique(x1[eta1==1]))
events2<-sort(unique(x2[eta2==1]))

if(!censor.mark.all){
# Censored points
# Curves are marked at each censoring time which is not also a death time
x1.cens<-x1[eta1==0] 
x2.cens<-x2[eta2==0]
x1.cens<-setdiff(x1.cens,events1)  # Censorings that are NOT an event (this is the R version which follows estimation convention)
x2.cens<-setdiff(x2.cens,events2)
x1.match<-match(x1.cens,at.points1)
x2.match<-match(x2.cens,at.points2)
}

if(censor.mark.all){
# Censored points
x1.cens<-x1[eta1==0] 
x2.cens<-x2[eta2==0]
x1.match<-match(x1.cens,at.points1)
x2.match<-match(x2.cens,at.points2)
}

risk.points<-round(c(seq(0,max(at.points),by=by.risk)))
risk.points<-sort(unique(c(risk.points,risk.add)))

risk.2<-unlist(lapply(as.list(risk.points),R.Weighted,error=x2,W=w2))
risk.1<-unlist(lapply(as.list(risk.points),R.Weighted,error=x1,W=w1))

risk.2<-round(risk.2)
risk.1<-round(risk.1)

# Note: "med" denotes general quantile, but median is default
med.1<-suppressWarnings(min(at.points1[S1.KM<=quant]))
# log transform
KM1.lower<-exp(log(S1.KM)-1.96*SE1.KM/S1.KM)
KM1.upper<-exp(log(S1.KM)+1.96*SE1.KM/S1.KM)
med.1.lower<-suppressWarnings(min(at.points1[KM1.lower<=quant]))
med.1.upper<-suppressWarnings(min(at.points1[KM1.upper<=quant]))

if(details){
cat("Stratum,n,events=",c(stratums[1],length(x1),sum(eta1)),"\n")
cat("Median, Lower, Upper=",c(med.1,med.1.lower,med.1.upper),"\n")
}

med.2<-suppressWarnings(min(at.points2[S2.KM<=quant]))
KM2.lower<-exp(log(S2.KM)-1.96*SE2.KM/S2.KM)
KM2.upper<-exp(log(S2.KM)+1.96*SE2.KM/S2.KM)
med.2.lower<-suppressWarnings(min(at.points2[KM2.lower<=quant]))
med.2.upper<-suppressWarnings(min(at.points2[KM2.upper<=quant]))

if(details){
cat("Stratum,n,events=",c(stratums[2],length(x2),sum(eta2)),"\n")
cat("Median, Lower, Upper=",c(med.2,med.2.lower,med.2.upper),"\n")
}


# Experimental
fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1.plot)
S1.KM.plot<-fit1$S.KM

# Control
fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2.plot)
S2.KM.plot<-fit2$S.KM

Y1.plot<-S1.KM.plot 
Y2.plot<-S2.KM.plot

Y1<-S1.KM
Y2<-S2.KM

if(is.null(x.truncate)) xmax<-max(c(at.points,xmax))
if(!is.null(x.truncate)) xmax<-x.truncate

plot(at.points1.plot,Y1.plot,type="s",ylim=c(ymin,ymax),xlim=c(xmin,xmax),lty=1,col=col.1,lwd=2,xlab=Xlab,ylab=Ylab,axes=FALSE)
lines(at.points2.plot,Y2.plot,lty=1,type="s",col=col.2,lwd=2)

if(show.ticks==TRUE){
points(x1.cens,Y1[x1.match],pch=3,col=col.1)
points(x2.cens,Y2[x2.match],pch=3,col=col.2)
}

#if(what.toplot=="KM") abline(v=0,col="grey",lwd=3,lty=2)

abline(h=0,lty=1,col=1)
#axis(2,at=c(0.0,0.2,0.4,0.5,0.6,0.8,1.0))
axis(2,at=c(prob.points))
axis(1,at=risk.points,labels=c(risk.points+time.zero))
#axTicks(side=1, axp =c(risk.points+5.25))
 
box()
if(risk.set){
text(c(risk.points),-0.075,c(risk.2),col=col.2,cex=risk.cex)
text(c(risk.points),ymin,c(risk.1),col=col.1,cex=risk.cex)
}
if(add.segment) segments(-10,1,0,1)  

if(show.med){
if(med.1!=Inf & med.2!=Inf){
text(xmax-4,ymax-0.25,paste(round(med.1,2)),col=col.1,cex=med.cex)
text(xmax-4,ymax-(0.25+del.med),paste(round(med.2,2)),col=col.2,cex=med.cex)
}
if(med.1!=Inf & med.2==Inf){
text(xmax-2,ymax-0.2,paste(round(med.1,2)),col=col.1,cex=med.cex)
text(xmax-2,ymax-(0.2+del.med),paste("median NA"),col=col.2,cex=1)
}
if(med.1==Inf & med.2!=Inf){
text(xmax-2,ymax-0.2,paste("median NA"),col=col.1,cex=1)
text(xmax-2,ymax-(0.2+del.med),paste(round(med.2,2)),col=col.2,cex=med.cex)
}
if(med.1==Inf & med.2==Inf){
text(xmax-4,ymax-0.2,paste("median NA"),col=col.1,cex=1)
text(xmax-4,ymax-(0.2+del.med),paste("median NA"),col=col.2,cex=1)
}

}
if(show.cox){
rp<-vector('expression',4)
rp[1]=substitute(expression(italic(HR)== MYVALUE1), 
        list(MYVALUE1 = format(hr.ci[1],digits = 2)))[2]
rp[2]=substitute(expression(italic(Lower) == MYVALUE2), 
        list(MYVALUE2 = format(hr.ci[2], digits = 2)))[2]
rp[3]=substitute(expression(italic(Upper) == MYVALUE3), 
        list(MYVALUE3= format(hr.ci[3], digits = 2)))[2]               
rp[4]=substitute(expression(italic(pval.cox) == MYVALUE3), 
                 list(MYVALUE3= format(pval.cox, digits = 2)))[2]               
legend(put.legend.cox, legend = rp, bty = 'n')
}

if(show.logrank){
rp<-vector('expression',1)
rp[1]=substitute(expression(italic(p[superiority])== MYVALUE), 
        list(MYVALUE = format(p.lr,digits=lr.digits)))[2]
legend(put.legend.lr, legend = rp, bty = 'n')
}


return(list(S1.KM=S1.KM,S2.KM=S2.KM,KM1.lower=KM1.lower,KM2.lower=KM2.lower,KM1.upper=KM1.upper,KM2.upper=KM2.upper,
            pval.cox=pval.cox,pval.lr=p.lr,cpoints=tpoints.add,
S1.dpoints=S1.dpoints,SE1.dpoints=SE1.dpoints,
S2.dpoints=S2.dpoints,SE2.dpoints=SE2.dpoints,
LR.dpoints=LR.dpoints,ZLR.dpoints=ZLR.dpoints,
at.points1=at.points1,at.points2=at.points2,
med.1=med.1,med.2=med.2,events.1=sum(eta1),events.2=sum(eta2),risk.1=risk.1,risk.2=risk.2,risk.points=risk.points))
}



