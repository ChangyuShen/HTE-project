est.main<-function(subpop,n, prior, n.sam=100000,n.iter=1)
{
  relu<-0.001
  relu.0<-0.0001
  relu.1<-0.0001
  close.0<-seq(relu.0,relu-relu.0,relu.0)
  close.1<-seq(1-relu+relu.1,1-relu.1,relu.1)
  pct<-c(close.0,seq(relu,1-relu,relu),close.1)
  pct.target<-c(0.005,0.025,0.5,0.955,0.975)
  dif.pos<-9+length(close.0)+pct.target/relu
  or.pos<-9+2*length(close.0)+length(close.1)+(1-relu)/relu+pct.target/relu
  out.style<-c(1:9,dif.pos,or.pos)

  sel<-which(prior[,3]>0 & prior[,6]>0)
  true.par<-cbind(ps=prior[sel,9]/sum(prior[sel,9]),pc=prior[sel,2]/prior[sel,3],pt=prior[sel,5]/prior[sel,6])

  Mu<-0.5*c(1,sum(true.par[,1]*true.par[,2]),sum(true.par[,1]*true.par[,3]))
  Sigma<-0.25*xpnd(c(sum(true.par[,1]^2),sum(true.par[,1]^2*true.par[,2]),sum(true.par[,1]^2*true.par[,3]),
       sum(true.par[,1]^2*true.par[,2]^2),sum(true.par[,1]^2*true.par[,2]*true.par[,3]),
       sum(true.par[,1]^2*true.par[,3]^2)))

  temp<-foreach(i=1:n.iter,.combine=c,.options.RNG=123) %dorng% const.cal(subpop,mu=Mu,sigma=Sigma,n=n,n.sample=n.sam)
  temp<-log(n.sam)-temp
  max.c<-max(temp)
  c.large<-log(n.iter)+log(n.sam)-max.c-log(sum(exp(temp-max.c)))
  temp1<-2*max.c+log(sum(exp(2*(temp-max.c))))-log(n.iter)-2*log(n.sam)
  temp1<-0.5*(temp1+log(1-exp(-2*c.large-temp1))-log(n.iter))
  relative.err<-exp(c.large+temp1)

  set.seed(312)
  post1d<-post.seminormal(subpop,mu=Mu,sigma=Sigma,ci=T,n=n,const=c.large,probs=pct,n.sample=n.sam)
  return(post1d[out.style])
}

# sample.subpop samples sub-populations and calculate sub-population proportion, the event rate of
# the control and intervention.
# subpop: sampled sub-populations
# ind: a vector of strings storing the "name" of the sampled subpopulations
# cordi: subpopulation proportion, event rates of the control and intervention

sample.subpop<-function(b=1000)
{
  n<-sum(all[,9])
  n.cell<-nrow(all)
  all.subpop<-as.data.frame(matrix(NA,b,ncol(all)))
  colnames(all.subpop)<-colnames(all)
  bin.ind<-NULL
  for (i in 1:b)
   {
     pos<-rbinom(n.cell,1,0.5)
     temp<-apply(all[which(pos==1),],2,sum)
     if (temp[1]>0 & temp[2]>0 & temp[4]>0 & temp[5]>0)
       {
         all.subpop[i,]<-temp
         bin.ind<-c(bin.ind,paste(pos,collapse=""))
       }
   }
  cordi<-cbind(ps=all.subpop$total/n,pc=all.subpop$event.ctl/all.subpop$total.ctl,pt=all.subpop$event.trt/all.subpop$total.trt)
  return(list(subpop=all.subpop,ind=bin.ind,cordi=cordi))
}

# sample.subpop1 essentialy is the same as sample.subpop, except it takes "dat" as input and the return object is different. 
# dat: cell level data

sample.subpop1<-function(dat,b=1000)
{
  n<-sum(dat[,9])
  n.cell<-nrow(dat)
  all.subpop<-as.data.frame(matrix(NA,b,ncol(dat)))
  colnames(all.subpop)<-colnames(dat)
  bin.ind<-NULL
  for (i in 1:b)
   {
     pos<-rbinom(n.cell,1,0.5)
     temp<-apply(dat[which(pos==1),],2,sum)
     if (temp[1]>0 & temp[2]>0 & temp[4]>0 & temp[5]>0)
       {
         all.subpop[i,]<-temp
         bin.ind<-c(bin.ind,paste(pos,collapse=""))
       }
   }
  rownames(all.subpop)<-bin.ind
  return(all.subpop)
}



# simu.cordi does bootstrap and subpopulation sampling using cell-level format in simulation
# dat: a data frame at cell level. Rows correspond to cells and columns correspond to the nine statistics

simu.cordi<-function(dat,b=1000)

{
 n<-sum(dat[,9])
 out<-matrix(NA,b,3)
 colnames(out)<-c("ps","pc","pt")
 rmulti<-function(x) return(rmultinom(n=1,size=x[1],prob=x[2:length(x)])) 
 for (i in 1:b)
   {
     bs.sample<-as.vector(rmultinom(1,size=n,prob=dat[,9]))
     names(bs.sample)<-rownames(dat)
     pos<-rbinom(length(bs.sample),1,0.5)
     bs.subpop<-bs.sample[pos==1]
     bs.subpop<-bs.subpop[bs.subpop>0]
     all.subpop<-apply(cbind(as.numeric(bs.subpop),dat[rownames(dat)%in%(names(bs.subpop)),c(1,2,4,5)]),1,rmulti)
     subpop.n<-apply(all.subpop,1,sum)
     if (subpop.n[1]>0 & subpop.n[2]>0 & subpop.n[3]>0 & subpop.n[4]>0)
          out[i,]<-c(sum(subpop.n)/n,subpop.n[2]/sum(subpop.n[1:2]),subpop.n[4]/sum(subpop.n[3:4]))
   }
  return(out)
}



  
# post.moment caculates the poterior mean,variance and covariance of the log odds for control and treatment arms
# and return the probability that the difference of the log odds (control minus treatment) is greater than 0
# y: vector that stores the nonevent,event and total count for control, treatment and combined of a sub-population
# n: total number of subjects in the original data sets
# mu.pc: posteiror mean of the log odds of the control
# mu.pt: posteiror mean of the log odds of the intervention
# mu.dif: posteiror log odds ratio
# sd.dif: large sample standard error of the posterior log odds ratio
# z: z value
# probgt0: prob(z>0)

  
post.moment<-function(y,n=4000)
{
  phat<-c(y[9]/n,y[2]/y[3],y[5]/y[6])
  grid.den<-myderiv(cordi)
  f.df<-nbr_fun(grid.den,phat)$est
  if (f.df[1]==0) return("NA")
  else
    {
       mu.pc<-f.df[2]/f.df[1]/y[3]+digamma(y[2]+1)-digamma(y[1]+1)
       mu.pt<-f.df[3]/f.df[1]/y[6]+digamma(y[5]+1)-digamma(y[4]+1)
       var.pc<-(f.df[4]/f.df[1]-f.df[2]^2/f.df[1]^2)/y[3]^2+trigamma(y[2]+1)+trigamma(y[1]+1)
       var.pt<-(f.df[6]/f.df[1]-f.df[3]^2/f.df[1]^2)/y[6]^2+trigamma(y[5]+1)+trigamma(y[4]+1)
       cov.pc.pt<-(f.df[5]/f.df[1]-f.df[2]*f.df[3]/f.df[1]^2)/(y[3]*y[6])
       mu.dif<-mu.pc-mu.pt
       #var.dif<-1/y[1]+1/y[2]+1/y[4]+1/y[5]
       var.dif<-var.pc+var.pt-2*cov.pc.pt
       if (var.dif<=0) 
           return (NA) else
       {
         sd.dif<-sqrt(var.dif)
         z<-mu.dif/sd.dif
         return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,probgt0=pnorm(z)))
       }
     }
}
  
# Arguments and outputs same as post.moment. The method is based on 1d density and derivative estimation.
# n.grid is the number of nearest data points around the given data point.

post1d.moment<-function(y,n.grid=5000,n=4000)
{
    phat<-as.numeric(c(y[9]/n,y[2]/y[3],y[5]/y[6]))
    tem<-cordi[,1]-phat[1]
    d.pc<-cbind(tem,cordi[,3]-phat[3])
    d.pc<-d.pc%*%t(U.pc)
    d.pc<-d.pc[,1]^2+d.pc[,2]^2
    x.pc<-cordi[order(d.pc)[1:n.grid],2]
    d.pt<-cbind(tem,cordi[,2]-phat[2])
    d.pt<-d.pt%*%t(U.pt)
    d.pt<-d.pt[,1]^2+d.pt[,2]^2
    x.pt<-cordi[order(d.pt)[1:n.grid],3]   
    df.pc<-myderiv(x.pc)
    loc.pc<-sum(phat[2]-df.pc$x.grid[[1]]>=0)
    if (loc.pc==0) loc.pc<-1
    df.pt<-myderiv(x.pt)
    loc.pt<-sum(phat[3]-df.pt$x.grid[[1]]>=0) 
    if (loc.pt==0) loc.pt<-1
    if (df.pc$fhat[loc.pc]==0 | df.pt$fhat[loc.pt]==0) 
          return(c(rep(NA,6)))  else
       {
          mu.pc<-as.numeric(df.pc$fhat_grad[loc.pc]/df.pc$fhat[loc.pc]/y[3]+digamma(y[2]+1)-digamma(y[1]+1))
          mu.pt<-as.numeric(df.pt$fhat_grad[loc.pt]/df.pt$fhat[loc.pt]/y[6]+digamma(y[5]+1)-digamma(y[4]+1))
          mu.dif<-mu.pc-mu.pt
          var.dif<-1/y[1]+1/y[2]+1/y[4]+1/y[5]
          sd.dif<-as.numeric(sqrt(var.dif))
          z<-mu.dif/sd.dif
          return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,probgt0=pnorm(z)))
       }
}


# Arguments and outputs same as post.moment. Calculate standard statistics.

  
naive.moment<-function(y)
{
  phat<-c(y[2]/y[3],y[5]/y[6])
  odds.pc<-phat[1]/(1-phat[1])
  odds.pt<-phat[2]/(1-phat[2])       
  odds.ratio<-odds.pc/odds.pt
  var.logor<-1/y[1]+1/y[2]+1/y[4]+1/y[5]
  if (var.logor<=0) 
      return (NA) else
  {
      sd.logor<-sqrt(var.logor)
      sd.or<-odds.ratio*sd.logor
      z<-log(odds.ratio)/sd.logor
      return(c(mu.pc=phat[1],mu.pt=phat[2],mu.dif=phat[1]-phat[2],odds.pc=odds.pc,odds.pt=odds.pt,odds.ratio=odds.ratio,sd.or=sd.or,z=z,probgt0=pnorm(z)))
  }
}


post.cellsimu<-function(y,truep,thre=0,ci=F,n.sample=10000,n=4000,r=0.5,prob.thre=c(0.025,0.05,0.975))
{
ff<-function(k)
{
  a<-truep[rbinom(nrow(truep),1,0.5)==1,]
  a<-apply(cbind(a[,1],a[,1]*a[,2],a[,1]*a[,3]),2,sum)
  return(a)
} 
condi.disX<-function(par,y0)
{
  d1<-dbinom(y0[9],n,par[1],log=T)
  d2<-dbinom(y0[3],y0[9],r,log=T)
  d3<-dbinom(y0[2],y0[3],par[2],log=T)
  d4<-dbinom(y0[5],y0[6],par[3],log=T)
  return(d1+d2+d3+d4)
}

theta<-t(sapply(1:n.sample,ff))
logodds.pc<-log(theta[,2]/(theta[,1]-theta[,2]))
logodds.pt<-log(theta[,3]/(theta[,1]-theta[,3]))
logor<-logodds.pc-logodds.pt

theta[,2]<-theta[,2]/theta[,1]
theta[,3]<-theta[,3]/theta[,1]
if (is.vector(y))
{
wt<-apply(theta,1,condi.disX,y0=y)
wt<-exp(wt-max(wt))
wt<-wt/sum(wt)

mu.pc<-sum(wt*logodds.pc)
mu.pt<-sum(wt*logodds.pt)
mu.dif<-sum(wt*logor)
sd.dif<-sqrt(sum(wt*logor^2)-mu.dif^2)
prb<-sum(wt*(logor>thre))
if (ci==F)
  return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,prb=prb)) else
  {
    ci<-Hmisc:::wtd.quantile(logor,weights=wt*n.sample,probs=prob.thre,normwt=T)
    names(ci)<-paste("q",prob.thre*100,"%",sep="")
    return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,prb=prb,ci))
  }
} else
{
  wt<-NULL
  for (i in 1:nrow(y))
     wt<-cbind(wt,apply(theta,1,condi.disX,y0=as.numeric(y[i,])))
  wt<-exp(t(wt)-apply(wt,2,max))
  wt<-t(wt/apply(wt,1,sum))

  mu.pc<-apply(wt*logodds.pc,2,sum)
  mu.pt<-apply(wt*logodds.pt,2,sum)
  mu.dif<-apply(wt*logor,2,sum)
  sd.dif<-sqrt(apply(wt*logor^2,2,sum)-mu.dif^2)
  prb<-apply(wt*(logor>thre),2,sum)
  if (ci==F)
    return(cbind(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,prb=prb)) else
  {
    ci<-apply(wt*n.sample,2,Hmisc:::wtd.quantile,x=logor,probs=prob.thre,normwt=T)
    colnames(ci)<-paste("q",prob.thre*100,"%",sep="")
    return(cbind(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,prb=prb,ci))
  }
}
}

true.post<-function(par,y0,mu,sigma,n,log=T)
{
  d1<-dbinom(y0[9],n,par[1],log=T)
  d2<-dbinom(y0[2],y0[3],par[2],log=T)
  d3<-dbinom(y0[5],y0[6],par[3],log=T)
  p1<-dnorm(par[1],mean=mu[1],sd=sqrt(sigma[1,1]),log=T)
  p2<-mvtnorm:::dmvnorm(par[2:3],mean=(mu[2:3]+((par[1]-mu[1])/sigma[1,1])*as.numeric(sigma[2:3,1]))/par[1],
      sigma=(sigma[2:3,2:3]-sigma[2:3,1]%*%t(sigma[1,2:3])/sigma[1,1])/par[1]^2,log=T)
  if (log)
     return(d1+d2+d3+p1+p2) else
  return(exp(d1+d2+d3+p1+p2))
}

true.post.matrix<-function(par,y0,mu,sigma,n,log=T)
{
  d1<-dbinom(y0[9],n,par[,1],log=T)
  d2<-dbinom(y0[2],y0[3],par[,2],log=T)
  d3<-dbinom(y0[5],y0[6],par[,3],log=T)
  p1<-dnorm(par[,1],mean=mu[1],sd=sqrt(sigma[1,1]),log=T)
  p2<-mvtnorm:::dmvnorm(par[,2:3],mean=(mu[2:3]+((par[1]-mu[1])/sigma[1,1])*as.numeric(sigma[2:3,1]))/par[1],
      sigma=(sigma[2:3,2:3]-sigma[2:3,1]%*%t(sigma[1,2:3])/sigma[1,1])/par[1]^2,log=T)
  if (log)
     return(d1+d2+d3+p1+p2) else
  return(exp(d1+d2+d3+p1+p2))
}

const.cal<-function(y0,mu,sigma,n.sample=100000,n)
{
  theta<-mvtnorm:::rmvnorm(n.sample,mean=mu,sigma=sigma)
  theta[,2:3]<-theta[,2:3]/theta[,1]
  theta<-(theta>=1)*0.9999999+(theta<=0)*0.0000001+(theta>0 & theta<1)*theta
  d1<-dbinom(y0[9],n,theta[,1],log=T)
  d2<-dbinom(y0[2],y0[3],theta[,2],log=T)
  d3<-dbinom(y0[5],y0[6],theta[,3],log=T)
  d<-d1+d2+d3
  max.d<-max(d)
  out<-log(n.sample)-max.d-log(sum(exp(d-max.d))) 
  return(out)
}


post.seminormal<-function(y,mu,sigma,thre=1,ci=F,n.sample=50000,n=4000,probs=c(seq(0.005,0.05,0.005),seq(0.95,0.995,0.005)),const)
{

theta.hat<-c(y[9]/n,y[2]/y[3],y[5]/y[6])
omega.hat.inv<-diag(c(n,y[3],y[6])/theta.hat/(1-theta.hat))
delta.inv<-diag(c(1,mu[1],mu[1]))
sigma.inv<-solve(sigma)

v<-solve(delta.inv%*%sigma.inv%*%delta.inv+omega.hat.inv)
u<-v%*%(delta.inv%*%sigma.inv%*%mu+omega.hat.inv%*%theta.hat)

theta.sample<-mvtnorm:::rmvnorm(n.sample,mean=u,sigma=v)
theta.sample<-(theta.sample>=1)*0.9999999+(theta.sample<=0)*0.0000001+(theta.sample>0 & theta.sample<1)*theta.sample

log.true.density<-true.post.matrix(theta.sample,y0=y,mu=mu,sigma=sigma,n=n)
log.imp.density<-mvtnorm:::dmvnorm(theta.sample,mean=u,sigma=v,log=T)
wt<-exp(const+log.true.density-log.imp.density)

risk.dif<-theta.sample[,2]-theta.sample[,3]
odds.pc<-theta.sample[,2]/(1-theta.sample[,2])
odds.pt<-theta.sample[,3]/(1-theta.sample[,3])
or<-odds.pc/odds.pt

mu.pc<-mean(wt*theta.sample[,2])
mu.pt<-mean(wt*theta.sample[,3])
mu.dif<-mean(wt*risk.dif)
sd.dif<-sqrt(mean(wt*(risk.dif-mu.dif)^2))

odds.pc<-mean(wt*odds.pc)
odds.pt<-mean(wt*odds.pt)
odds.ratio<-mean(wt*or)
sd.ratio<-sqrt(mean(wt*(or-odds.ratio)^2))
prb<-min(mean(wt*(or>thre)),1)

#return(cbind(theta.sample,wt))
if (ci==F)
  return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,odds.pc=odds.pc,odds.pt=odds.pt,odds.ratio=odds.ratio,sd.ratio=sd.ratio,prb=prb)) else
  {
    ci.dif<-Hmisc:::wtd.quantile(risk.dif,weights=wt,type="i/n",probs=probs,normwt=T)
    ci.ratio<-Hmisc:::wtd.quantile(or,weights=wt,type="i/n",probs=probs,normwt=T)
    return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,odds.pc=odds.pc,odds.pt=odds.pt,odds.ratio=odds.ratio,sd.ratio=sd.ratio,
           prb=prb,ci.dif,ci.ratio))
  }
}


bs.est<-function(dat)
  {
    n<-sum(dat[,9])
    rmulti<-function(x) return(rmultinom(n=1,size=x[1],prob=x[2:length(x)])) 
    bs<-as.vector(rmultinom(1,size=n,prob=dat[,9]))
    use<-which(bs>0)
    bs<-t(apply(cbind(as.numeric(bs[use]),dat[use,c(1,2,4,5)]),1,rmulti))
    n.pc<-bs[,1]+bs[,2]
    n.pt<-bs[,3]+bs[,4]
    use<-which(n.pc>0 & n.pt>0)
    bs.par<-cbind(ps=(n.pc[use]+n.pt[use])/sum(n.pc[use]+n.pt[use]),pc=bs[use,2]/n.pc[use],pt=bs[use,4]/n.pt[use])
    return(bs.par)
  }

bs.est1<-function(dat)
  {
    n<-sum(dat[,9])
    rmulti<-function(x) return(rmultinom(n=1,size=x[1],prob=x[2:length(x)]))
    bs<-as.vector(rmultinom(1,size=n,prob=dat[,9]))
    use<-which(bs>0)
    bs<-t(apply(cbind(as.numeric(bs[use]),dat[use,c(1,2,4,5)]),1,rmulti))
    n.pc<-bs[,1]+bs[,2]
    n.pt<-bs[,3]+bs[,4]
    use<-which(n.pc>0 & n.pt>0)
    bs.par<-cbind(ps=(n.pc[use]+n.pt[use])/sum(n.pc[use]+n.pt[use]),pc=bs[use,2]/n.pc[use],pt=bs[use,4]/n.pt[use])
    mu<-0.5*c(1,sum(bs.par[,1]*bs.par[,2]),sum(bs.par[,1]*bs.par[,3]))
    sigma<-0.25*c(sum(bs.par[,1]^2),sum(bs.par[,1]^2*bs.par[,2]),sum(bs.par[,1]^2*bs.par[,3]),
       sum(bs.par[,1]^2*bs.par[,2]^2),sum(bs.par[,1]^2*bs.par[,2]*bs.par[,3]),
       sum(bs.par[,1]^2*bs.par[,3]^2))
    return(c(mu,sigma))
  }


stocha.search<-function(dat,method=naive.moment,n.iter=15000,sgn=1)
{
  n.total<-sum(dat[,9])
  n.cell<-nrow(dat)
  temp<-rep(0,9)
  dif<-NA
  while (temp[3]==0 | temp[6]==0 | is.na(dif))
   {
      pop.old<-rbinom(n.cell,1,0.5)
      temp<-apply(dat[pop.old==1,],2,sum)
      if (temp[3]>0 & temp[6]>0) 
          dif<-method(temp,n.total)[5]
   } 
  trace.dif<-dif
  trace.pop<-paste(pop.old,collapse="")
  for (iter in 1:n.iter)
  {
    tmp<-(iter<=7000)*0.3+(iter>7000 & iter<=13000)*0.05+(iter>13000)*0.005
    prb<-(iter<=7000)*0.2+(iter>7000 & iter<=13000)*0.3+(iter>13000)*0.5
    temp<-rep(0,9)
    dif.new<-NA
    while (temp[3]==0 | temp[6]==0 | is.na(dif.new))
     {
       pop.new<-perturb(pop.old,prb.stay=prb)
       temp<-apply(dat[pop.new==1,],2,sum)
       if (temp[3]>0 & temp[6]>0) 
          dif.new<-method(temp,n.total)[5]
     } 
  
    if (log(runif(1))*tmp<(dif.new-trace.dif[iter])*sgn)
     {   
       trace.dif<-c(trace.dif,dif.new)
       trace.pop<-c(trace.pop,paste(pop.new,collapse=""))
     }  else
     { 
       trace.dif<-c(trace.dif,trace.dif[iter])
       trace.pop<-c(trace.pop,trace.pop[iter])
     }
     pop.old<-pop.new
   }
   return(list(sub.pop=trace.pop,dif=trace.dif))
}


# perturb generates new sub population based on existing one based on a bernoulli distribution
# sub: a binary vector indicating current sub-population (
# level: a vector indicating the number of category (L+1) for each component
# prb.stay: the probability of staying in current category. 

perturb<-function(sub,prb.stay=0.5)
{
   n.sub<-length(sub)
   sub.new<-sub
   to.change<-which(rbinom(n.sub,1,1-prb.stay)==1)
   sub.new[to.change]<-1-sub[to.change]
   return(sub.new)
}


post1d.refine<-function(y,n.grid=5000,n=4000,epi=0.0001)
{
    phat<-as.numeric(c(y[9]/n,y[2]/y[3],y[5]/y[6]))
    tem<-cordi[,1]-phat[1]
    d.pc<-cbind(tem,cordi[,3]-phat[3])
    d.pc<-d.pc%*%t(U.pc)
    d.pc<-d.pc[,1]^2+d.pc[,2]^2
    x.pc<-cordi[order(d.pc)[1:n.grid],2]
    d.pt<-cbind(tem,cordi[,2]-phat[2])
    d.pt<-d.pt%*%t(U.pt)
    d.pt<-d.pt[,1]^2+d.pt[,2]^2
    x.pt<-cordi[order(d.pt)[1:n.grid],3]
    df.pc<-myderiv(x.pc)
    loc.pc<-sum(phat[2]-df.pc$x.grid[[1]]>=0)
    if (loc.pc==0) loc.pc<-1
    df.pt<-myderiv(x.pt)
    loc.pt<-sum(phat[3]-df.pt$x.grid[[1]]>=0)
    if (loc.pt==0) loc.pt<-1
    d.pc.star<-cbind(tem,cordi[,3]-phat[3]-epi)
    d.pc.star<-d.pc.star%*%t(U.pc)
    d.pc.star<-d.pc.star[,1]^2+d.pc.star[,2]^2
    x.pc.star<-cordi[order(d.pc.star)[1:n.grid],2] 
    df.pc.star<-myderiv(x.pc.star)
    loc.pc.star<-sum(phat[2]-df.pc.star$x.grid[[1]]>=0)
    if (loc.pc.star==0) loc.pc.star<-1
    if (df.pc$fhat[loc.pc]==0 | df.pt$fhat[loc.pt]==0)
          return(c(rep(NA,6)))  else
       {
          mu.pc<-as.numeric(df.pc$fhat_grad[loc.pc]/df.pc$fhat[loc.pc]/y[3]+digamma(y[2]+1)-digamma(y[1]+1))
          mu.pt<-as.numeric(df.pt$fhat_grad[loc.pt]/df.pt$fhat[loc.pt]/y[6]+digamma(y[5]+1)-digamma(y[4]+1))
          mu.dif<-mu.pc-mu.pt
          #var.dif<-1/y[1]+1/y[2]+1/y[4]+1/y[5]
          var.pc<-(df.pc$fhat_curv[loc.pc]/df.pc$fhat[loc.pc]-df.pc$fhat_grad[loc.pc]^2/df.pc$fhat[loc.pc]^2)/y[3]^2+trigamma(y[2]+1)+trigamma(y[1]+1) 
          var.pt<-(df.pt$fhat_curv[loc.pt]/df.pt$fhat[loc.pt]-df.pt$fhat_grad[loc.pt]^2/df.pt$fhat[loc.pt]^2)/y[6]^2+trigamma(y[5]+1)+trigamma(y[4]+1)
          if (df.pc.star$fhat[loc.pc.star]>0)
            {
              cov.pcpt<-(df.pc.star$fhat_grad[loc.pc.star]/df.pc.star$fhat[loc.pc.star]-df.pc$fhat_grad[loc.pc]/df.pc$fhat[loc.pc])/epi
              cov.pcpt<-cov.pcpt/y[3]/y[6]
            }  else
          cov.pcpt<-0
          var.dif<-var.pc+var.pt-2*cov.pcpt
          sd.dif<-as.numeric(sqrt(var.dif))
          z<-mu.dif/sd.dif
          return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,sd.dif=sd.dif,z=mu.dif/sd.dif,probgt0=pnorm(z)))
       }
}

find.truth<-function(s)
{
  a<-as.numeric(unlist(strsplit(s,"")))
  temp<-true.par[a==1,]
  temp<-apply(cbind(temp[,1],temp[,1]*temp[,2],temp[,1]*temp[,3]),2,sum)
  mu.pc<-temp[2]/temp[1]
  mu.pt<-temp[3]/temp[1]
  mu.dif<-mu.pc-mu.pt
  odds.pc<-mu.pc/(1-mu.pc)
  odds.pt<-mu.pt/(1-mu.pt)
  odds.ratio<-odds.pc/odds.pt
  return(c(mu.pc=mu.pc,mu.pt=mu.pt,mu.dif=mu.dif,odds.pc=odds.pc,odds.pt=odds.pt,odds.ratio=odds.ratio))
}

# clap.subpop return the same subpopulation with some cells merged. The merged cell will
# have a label in the form of multiple-digit number for certain variables, e.g. 
# a cell with 1 210 10 means level 1 for first variable, levels 2, 1 and 0 for second
# variable and levels 1 and 0 for the third variable 
# sub is a vector with each element (an element is a cell) being a character string that records the level of the corresponding
# variable

clap.subpop<-function(sub)
{ 
  n<-length(sub)
  sub<-sapply(sub,strsplit,split="")
  sub<-lapply(sub,as.numeric)
  n.str<-length(sub[[1]])
  keep.going<-1
  while (keep.going>0)
  {
     keep.going<-0
     i<-1
     while (i<length(sub))
      {
        j<-i+1
        while (j<=length(sub))
          {
            temp<-sub[[i]]-sub[[j]]
            zero<-sum(temp==0)  
            if (zero==n.str-1)
               {
                  pos<-which(temp!=0)
                  a<-c(as.numeric(unlist(strsplit(as.character(sub[[i]][pos]),split=""))),
                       as.numeric(unlist(strsplit(as.character(sub[[j]][pos]),split=""))))
                  a<-sort(a,decreasing=T)
                  sub[[i]][pos]<-as.numeric(paste(a,collapse=""))
                  sub<-sub[-j]
                  keep.going<-keep.going+1
               }  else
               j<-j+1
          }     
        i<-i+1
      } 
  }
  return(sub)
} 
