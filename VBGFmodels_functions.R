#### Defining the five von Bertalannfy growth model parameterisations ####
## to be used by the general purpose optimisation function optim()
# must run entire script before running unif.growth.R

# 1: EXPONENTIAL
predict.Aeq1<- function(k, m0, the.data){
  # Calculate the predictions for the model with A=B=1
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  prediction<- c()
  for (i in 1:length(individuals)){
    subset<- the.data$time[the.data$individual == individuals[i]]
    prediction<- c(prediction, m0[i]*exp(k[i]*(subset - subset[1])))
  }  
  return(prediction)
}


m0.start<- function(the.data){
  # Sets the parameter m0 to equal the first measured mass value for each individual
  individuals<- unique(the.data$individual)
  m0<- rep(0, length(individuals))
  for (i in 1:length(individuals)){
    subset<- the.data$mass[the.data$individual == individuals[i]]
    m0[i]<- subset[1]
  }
  return(m0)
}

# A wrapper so that optim can be used to fit the same models
sumsq.optim<- function(par, the.data, predict.fun, m0=NULL, log=F){
  # 
  prediction<- predict.fun(par, the.data, m0)
  if (!log){
    return(sum((prediction - the.data$mass)^2))
  }else{
    return(sum((log(prediction) - log(the.data$mass))^2))
  }
  
}
pred.optim.Ae1<- function(par, the.data, m0){
  
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  
  k<-par[1:n]
  if (is.null(m0)){
    # if m0 is null, then that means it is contained within par
    m0<- par[n+(1:n)]
  }
  return(predict.Aeq1(k, m0, the.data))
}


# 2: GOMPERTZ
predict.gomp<- function(b, k, m0, the.data){
  # Calculate the predictions for the Gompertz model
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  prediction<- c()
  
  for (i in 1:length(individuals)){
    subset<- the.data$time[the.data$individual == individuals[i]]
    prediction<- c(prediction, m0[i]*exp(-b[i]*(exp(-k[i]*(subset-subset[1]))-1)))
  }  
  return(prediction)
}
predict.gomp.nls<- function(b.t, k.t, m0, the.data){
  # wrapper that takes unconstrained parameters used in the fitting and 
  # converts to the parameters used by predict.Alt1
  new.params.gomp<- untransform.params.gomp(list(b.t=b.t, k.t=k.t))
  return(with(new.params.gomp, predict.gomp(b, k, m0, the.data)))
}
pred.optim.gomp<- function(par, the.data, m0){
  
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  
  b.t<-par[1:n]
  k.t<-par[n+(1:n)]
  if (is.null(m0)){
    # if m0 is null, then that means it is contained within par
    m0<- par[2*n+(1:n)]
  }
  return(predict.gomp.nls(b.t, k.t, m0, the.data))
}
library(boot)
sigmoid<- inv.logit
inv.sigmoid<- logit
transform.params.gomp<- function(par.u){
  par.t<-list(b.t= log(par.u$b))
  par.t$k.t<- log(par.u$k)
  return(par.t)
}
untransform.params.gomp<- function(par.t){
  par.u<-list(b=exp(par.t$b.t))
  par.u$k<- exp(par.t$k.t)
  return(par.u)
}
b.start<- function(I.start, m0){
  return(log(I.start/m0))
}

# 3: GENERALISED-VBGF
predict.Alt1.third<- function(A, f, k, m0, the.data){
  # Calculate the predictions for the Gompertz model
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  prediction<- c()
  
  for (i in 1:length(individuals)){
    subset<- the.data$time[the.data$individual == individuals[i]]
    kk<- exp(-k[i]*(subset-subset[1]))
    prediction<- c(prediction,  m0[i]*((1-f*kk)/(1-f))^(-1/(A-1)))
  }  
  return(prediction)
}
predict.Alt1.third.nls<- function(A.t, f.t, k.t, m0, the.data){
  # wrapper that takes unconstrained parameters used in the fitting and 
  # converts to the parameters used by predict.AltB
  
  new.params.Alt1<- untransform.params.Alt1(list(A.t=A.t, f.t=f.t, k.t=k.t))
  return(with(new.params.Alt1, predict.Alt1.third(A, f, k, m0, the.data)))
}

transform.params.Alt1<- function(par.u){
  par.t<-list(A.t= inv.sigmoid(par.u$A))
  par.t$f.t<- inv.sigmoid(par.u$f)
  par.t$k.t<- log(par.u$k)
  return(par.t)
}
untransform.params.Alt1<- function(par.t){
  par.u<- list(A=sigmoid(par.t$A)) # change to A.t????
  par.u$f<- sigmoid(par.t$f.t)
  par.u$k<- exp(par.t$k.t)
  return(par.u)
}

# Code to choose starting values for A and f=1-Z from the GOmpertz parameter b, such that these are both in the range (0,1) 
start.af<- function(b, amax=0.1, fmax=0.9){
  bmax<- max(b)
  a<- min(amax, fmax/bmax)
  return(list(A=1-a, f=b*a))
}
pred.optim.Alt1.third.nls<- function(par, the.data, m0){
  
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  
  A.t<- par[1]
  f.t<- par[1+(1:n)]
  k.t<- par[1+n+(1:n)]
  if (is.null(m0)){
    # if m0 is null, then that means it is contained within par
    m0<- par[1+2*n+(1:n)]
  }
  return(predict.Alt1.third.nls(A.t, f.t, k.t, m0, the.data))
}
# 4: PURE ISOMORPHY

pred.optim.A23.third.nls<- function(par, the.data, m0){
  
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  
  A.t=inv.sigmoid(2/3)
  f.t<- par[(1:n)]
  k.t<- par[n+(1:n)]
  if (is.null(m0)){
    # if m0 is null, then that means it is contained within par
    m0<- par[2*n+(1:n)]
  }
  return(predict.Alt1.third.nls(A.t, f.t, k.t, m0, the.data))
}


# 5: SUPRA-EXPONENTIAL 
predict.Agt1<- function(A, Z, K, m0, the.data){
  # Calculate the predictions for the Gompertz model
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  prediction<- c()
  
  for (i in 1:length(individuals)){
    subset<- the.data$time[the.data$individual == individuals[i]]
    ff<- 1-Z[i]
    kk<- exp(K[i]*(A-1)*(subset-subset[1]))
    prediction<- c(prediction, m0[i]*((1-ff*kk)/(Z[i]))^(-1/(A-1)))
  }  
  
  return(prediction)
}
predict.Agt1.nls<- function(alpha.t, Z.t, s.t, m0, the.data){
  # wrapper that takes unconstrained parameters used in the fitting and 
  # converts to the parameters used by predict.Agt1
  #print(c(alpha.t, Z.t, s.t, m0))
  # First transform from unconstrained scale to 
  new.params<- untransform.params.Agt1(list(alpha.t=alpha.t, Z.t=Z.t, s.t=s.t))
  # Now calculate biological parameters
  
  A<- 1/new.params$alpha
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  t0<- rep(0, n)
  tmax<- rep(0, n)
  for (i in 1:n){
    subset<- the.data$time[the.data$individual == individuals[i]]
    t0[i]<- min(subset)
    tmax[i]<- max(subset)
  }  
  tstar.t0<- -log(1-new.params$Z)/(A-1)
  K<- new.params$s*tstar.t0/(tmax-t0)
  
  
  return(predict.Agt1(A, new.params$Z, K, m0, the.data))
}
# Code for transforming betweeen unconstrained parameters to be optimised by 
# NLS and bioligically meaningful params
library(boot)
sigmoid<- inv.logit
inv.sigmoid<- logit
transform.params.Agt1<- function(par.u){
  par.t<-list(alpha.t= inv.sigmoid(par.u$alpha))
  par.t$Z.t<-inv.sigmoid(par.u$Z)
  par.t$s.t<- inv.sigmoid(par.u$s)
  return(par.t)
}
untransform.params.Agt1<- function(par.t){
  par.u<- list(alpha=sigmoid(par.t$alpha.t))
  par.u$Z<- sigmoid(par.t$Z.t)
  par.u$s<- sigmoid(par.t$s.t)
  return(par.u)
}


start.params.Agt1.from.exp<- function(k, the.data, Z=1e-3, A=1.01){
  # Gives candidate starting parameters for an extended von Bertalanffy 
  # model with A>1, B=1, to be close to an exponential exp(kt)
  
  alpha<- 1/A
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  t0<- rep(0, n)
  tmax<- rep(0, n)
  for (i in 1:n){
    subset<- the.data$time[the.data$individual == individuals[i]]
    t0[i]<- min(subset)
    tmax[i]<- max(subset)
  }  
  
  # The following code works provided Z and k are either scalars or vectors
  # of length n
  if(length(Z) == 1){Z=rep(Z, n)}
  
  K<- k*Z
  tstar.t0<- -log(1-Z)/(K*(A-1))
  
  s<- (tmax-t0)/tstar.t0
  
  par.u<- list(alpha=alpha, Z=Z, s=s)
  
  #print(par.u)
  #print(n)
  return(transform.params.Agt1(par.u))
}

pred.optim.Agt1.nls<- function(par, the.data, m0){
  
  individuals<- unique(the.data$individual)
  n<- length(individuals)
  
  alpha.t<- par[1]
  Z.t<- par[1+(1:n)]
  s.t<- par[1+n+(1:n)]
  if (is.null(m0)){
    # if m0 is null, then that means it is contained within par
    m0<- par[1+2*n+(1:n)]
  }
  return(predict.Agt1.nls(alpha.t, Z.t, s.t, m0, the.data))
}

#### END ####
