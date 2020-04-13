### profile likelihood code for calculating confidence intervals
# for parameter A for the Generalised VBGF ###
########################
## load the parameter estimates

#######################
library(boot)
sigmoid<- inv.logit
inv.sigmoid<- logit
A<-sigmoid(params.Alt1[1]) # untransformed
f<-sigmoid(params.Alt1[2]) # untransformed
k<-exp(params.Alt1[3]) # untransformed
m0<-params.Alt1[4]
par.ml<-params.Alt1 # transformed parameters

#####################
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
#####################
# negative log likelihood function 
nll1<- function(pars, the.data, m0, pred.func){
  pp1<- pred.optim.Alt1.third.nls(pars, the.data, m0)
  S<-sumsq.optim(pars, the.data, pred.func, m0, log=T)
  #print(S)
  N<- length(pp1)
  LL<-((N/2)*(log(2*pi)+log(S/N)+1))
  
  return(LL)
  
}

nll1(par.ml, growthdata, params.Alt1[4], pred.optim.Alt1.third.nls) 
#####################
# Now let's compute confidence limits for the coefficient for A using a profile likelihood method

# First define a wrapper suitable for feeding to optim.
# The first argument is the set of parameters over which we optimise;
# the second argument is the one for A
nll1.wrap<- function(parwrap, A, the.data, pred.func){
  par<- rep(0, 4)
  par[2:4]<- parwrap[1:3] # k, f and mo
  par[1]<- A # parameter A 
  return(nll1(par, the.data, params.Alt1[4], pred.func))
  
}

# Now let's define a function for the profile likelihood:
nll.prof<- function(A, the.data, pred.func, ml.params, return.pars=F){
  par.start<- ml.params[2:4] # I added this as an argument to this function because otherwise the code is hard to follow
  oo<- optim(par.start, nll1.wrap, gr=NULL, A, the.data, pred.func)
  if (return.pars){ # I added the option to return all the details of the fit in case these might be useful
    return(oo)      # Thanks to R's design this is backward compatible with earlier code that does not set the param return.pars
  }else{
    return(oo$value)
  }
}

par.start<- c(par.ml[2:4]) # f, k, m0


######################

# Now for some code to find the confidence interval, i.e. the range of values where the nll is within a particular chi-squared statistic of the
# max lik values.
nll.ml<- nll1(par.ml[1:3], growthdata, par.ml[4], pred.optim.Alt1.third.nls)

nll.target<- nll.ml + 0.5*qchisq(0.95, 1) # the confidence limits are where the nll equals this value

# Set a range of values for (transformed) A
A.vals<- seq(par.ml[1], 0, length=1000)

# Taken from the profile likelihood script.  Each time A is changed, the optimisation starts from the best fitting values form the previous value of A
nll.vec<- c()
par.start<- par.ml
for (A in A.vals){
  nll.temp<-  nll.prof(A, growthdata, pred.optim.Alt1.third.nls, par.start, return.pars=T)
  nll.vec<- c(nll.vec, nll.temp$value)
  par.start<- c(0, nll.temp$par)
}
# We will find values for A where (nll - nll.target) equals zero, using the standard "uniroot" function.
# To do this, we need to find a pair of values of A that bracket the zero.
# Let's use another wrapper to define a function that is zero at
# these confidence limits
nll.cl<- function(A, the.data, pred.fun, ml.params){
  return(nll.prof(A, the.data, pred.fun, ml.params) - nll.target)
}

# First, the upper limit
# Need first to find values A.low and A.up that bracket each confidence interval.
# The first obvious choice is (0, A.ml)

A.low<- 0
A.up<- par.ml[1]

# If this window of parameters does not bracket the zero ,then keep moving it until it does
# Note that par.ml[1] might actually be negative, but this does not actually matter for the code
# because A.low is shifted by (-par.ml[1])


while (nll.cl(A.low, growthdata, pred.optim.Alt1.third.nls, par.ml) <  0){
  A.up<- A.low
  A.low<- A.low - par.ml[1]
}

uu<- uniroot(nll.cl, c(A.low, A.up), growthdata, pred.optim.Alt1.third.nls, par.ml) 
A.ci1<- uu$root # One confidence limit - the upper one is par.ml[1] is positive, otherwise the lower one

# Now the other limit

A.low<- par.ml[1]
A.up<- 2*par.ml[1]

# If this window of parameters does not bracket the zero ,then keep moving it until it does

while (nll.cl(A.up, growthdata, pred.optim.Alt1.third.nls, par.ml) <  0){
  A.low<- A.up
  A.up<- A.up  + par.ml[1]
}

uu<- uniroot(nll.cl, c(A.low, A.up), growthdata, pred.optim.Alt1.third.nls, par.ml) 
A.ci2<- uu$root # lower confidence interval


############################################
### visualisation 
# Lets look at these on a graph. x-axis is *untransformed* A
# Re-generate the nll vector to make sure it covers the confidence intervals.

xmin<- 1.2*A.ci1 -0.2*A.ci2
xmax<- 1.2*A.ci2 - 0.2*A.ci1
A.vals<- seq(xmin, xmax, length=100)


nll.vec<- c()
par.start<- par.ml
for (A in A.vals){
  nll.temp<-  nll.prof(A, growthdata, pred.optim.Alt1.third.nls, par.start, return.pars=T)
  nll.vec<- c(nll.vec, nll.temp$value)
  par.start<- c(0, nll.temp$par)
}

plot(sigmoid(A.vals), nll.vec, xlab="A", ylab="neg log likelihood", t="l")
abline(h=nll.target, col="green")
abline(v=sigmoid(A.ci1), col="green")
abline(v=sigmoid(A.ci2), col="green")

#### END ####
