### Code for fitting the five VBGF parameterisations to growth data


##### FUNCTION FOR LOOPING THE 5 MODELS OVER INDIVIDUALS IN A SPECIES ####
## to fit for a single species with one or more individual(s)
# this uses log least squares method of estimation 

unif.growth.mod<-function(the.data){ 
  # index the data
  individuals<-unique(the.data$individual)
  n<-length(individuals)
  # create lists to store results of each model 
  mod.exp.results<- list() 
  mod.gomp.results<- list()
  mod.Alt1.results<- list()
  mod.A23.results<- list()
  mod.Agt1.results<- list()
  A.estimates.Alt1<-list()
  A.estimates.Agt1<- list()
  sumsq.exp<- list()
  sumsq.gomp<- list()
  sumsq.Alt1<- list()
  sumsq.A23<- list()
  sumsq.Agt1<- list()
  # to loop models over individuals of a species
  for(i in 1:length(individuals)){
    individuals<-unique(the.data$individual)
    subset<- the.data[the.data$individual == individuals[i], ] 
    # now state conditions
    m0<-m0.start(subset) # start value for m0
    N<-length(subset$mass) 
    n<-length(individuals)
    
    # now the optimisation itself for each of the 5 growth models
    ### 1: EXPONENTIAL 
    # 1a: optimise with m0 fixed
    
    par.optim.exp<-rep(0.1, 1)
    pp.exp.fix<- optim(par.optim.exp, sumsq.optim, gr=NULL, subset, pred.optim.Ae1, m0, T, control=list(maxit=100000), method="Brent", lower=0, upper=1) # T for log least squares
    repeat{ 
      par.optim.exp<-pp.exp.fix$par
      pp.exp.fix<-optim(par.optim.exp, sumsq.optim, gr=NULL, subset, pred.optim.Ae1, m0, T, control=list(maxit=100000), method="Brent", lower=0, upper=1)
      par.optim.exp1<-pp.exp.fix$par
      identical<- sum((par.optim.exp-par.optim.exp1)^2) == 0 
      if(identical) break 
    }
    
    # 1b: optimise with m0 fit
    
    par.optim.exp.fit<-c(pp.exp.fix$par[1], m0)
    pp.exp.fit<- optim(par.optim.exp.fit, sumsq.optim, gr=NULL, subset, pred.optim.Ae1, NULL,  T, control=list(maxit=100000))
    repeat{ 
      par.optim.exp.fit<-pp.exp.fit$par
      pp.exp.fit<-optim(par.optim.exp.fit, sumsq.optim, gr=NULL, subset, pred.optim.Ae1, NULL, T, control=list(maxit=100000))
      par.optim.exp.fit1<-pp.exp.fit$par
      identical<- sum((par.optim.exp.fit - par.optim.exp.fit1)^2) == 0 
      if(identical) break
      
    }
    pred.exp<- predict.Aeq1(pp.exp.fit$par[1], pp.exp.fit$par[2], subset)
    
    mod.exp.results<-c(mod.exp.results, list(list(k.estimate=pp.exp.fit$par[1], m0.estimate=pp.exp.fit$par[2], sumsq=pp.exp.fit$value, 
                                                  counts=pp.exp.fit$counts, convergence=pp.exp.fit$convergence, prediction.exp=pred.exp))) # save m0 fit results for exp model
    sumsq.exp<-c(sumsq.exp, list(list(sumsq.exp=pp.exp.fit$value)))
    
    ### 2: GOMPERTZ
    # 2a: m0 fixed
    k0.t<- pp.exp.fit$par[1] # start value
    
    b0.t<-b.start(c(2*(max(subset$mass))), m0) # start value for b 
    
    # Transform parameters into values optim() understands
    par.start.gomp<- list(b0=b0.t, k0=k0.t)
    par.t.gomp<- transform.params.gomp(par.start.gomp)
    par.optim.gomp.fix<- unlist(par.t.gomp) # unlist so optim() understands
    pp.gomp.fix<- optim(par.optim.gomp.fix, sumsq.optim, gr=NULL, subset, pred.optim.gomp, m0, T, control=list(maxit=100000))
    repeat{
      par.optim.gomp.fix<-pp.gomp.fix$par
      pp.gomp.fix<-optim(par.optim.gomp.fix, sumsq.optim, gr=NULL, subset, pred.optim.gomp, m0, T, control=list(maxit=100000)) 
      par.optim.gomp.fix1<-pp.gomp.fix$par
      identical<- sum((par.optim.gomp.fix - par.optim.gomp.fix1)^2) ==0 
      if(identical) break
    }
    
    # 2b: m0 fit
    par.optim.gomp.fit<- c(pp.gomp.fix$par[1], pp.gomp.fix$par[2], m0) # b, k, m0
    pp.gomp.fit<- optim(par.optim.gomp.fit, sumsq.optim, gr=NULL, subset, pred.optim.gomp, NULL, T, control=list(maxit=100000))
    # add repeat
    repeat{
      par.optim.gomp.fit<-pp.gomp.fit$par
      pp.gomp.fit<-optim(par.optim.gomp.fit, sumsq.optim, gr=NULL, subset, pred.optim.gomp, NULL, T, control=list(maxit=100000))
      par.optim.gomp.fit1<-pp.gomp.fit$par
      identical<- sum((par.optim.gomp.fit - par.optim.gomp.fit1)^2) ==0
      if(identical) break
      
    }
    pred.gomp<-predict.gomp.nls(pp.gomp.fit$par[1], pp.gomp.fit$par[2], pp.gomp.fit$par[3], subset)
    
    mod.gomp.results<-c(mod.gomp.results, list(list(b.estimate=pp.gomp.fit$par[1], k.estimate=pp.gomp.fit$par[2], 
                                                    m0.estimate=pp.gomp.fit$par[3], sumsq=pp.gomp.fit$value, counts=pp.gomp.fit$counts, 
                                                    convergence=pp.gomp.fit$convergence, prediction.gomp=pred.gomp))) # save m0 fit results for Gomp model 
    
    sumsq.gomp<-c(sumsq.gomp, list(list(sumsq.gomp=pp.gomp.fit$value)))
    
    ### 3: A<1 MODEL
    # 3a: m0 fixed
    # starting values, obtained from the fitted Gompertz model
    b0.t.Alt1<- pp.gomp.fit$par[1]
    # untransform
    b0.t.Alt1<- exp(b0.t.Alt1[1]) 
    startp<- start.af(b0.t.Alt1) # starting A and b value(s)
    startp$k<- exp(pp.gomp.fit$par[2]) # starting k value
    startp.t<- transform.params.Alt1(startp) # starting A, f(b) and k values
    par.optim.Alt1.fix<- with(startp.t, c(A.t, f.t, k.t))
    print(par.optim.Alt1.fix)
    
    pp.Alt1.fix<- optim(par.optim.Alt1.fix, sumsq.optim, gr=NULL, subset, pred.optim.Alt1.third.nls, m0, T, control=list(maxit=100000))
    print(pp.Alt1.fix) 
    
    repeat{ 
      par.optim.Alt1.fix<-pp.Alt1.fix$par
      pp.Alt1.fix<- optim(par.optim.Alt1.fix, sumsq.optim, gr=NULL, subset, pred.optim.Alt1.third.nls, m0, T, control=list(maxit=100)) # reduced maxit as error appears with species 69 ind 1
      par.optim.Alt1.fix1<-pp.Alt1.fix$par
      identical<-sum((par.optim.Alt1.fix - par.optim.Alt1.fix1)^2) ==0 
      if(identical) break
      
    }
    print(pp.Alt1.fix) 
    # 3b: m0 fit
    par.optim.Alt1.fit<-c(pp.Alt1.fix$par[1], pp.Alt1.fix$par[2], pp.Alt1.fix$par[3], m0) # (A, f, k)
    pp.Alt1.fit<- optim(par.optim.Alt1.fit, sumsq.optim, gr=NULL, subset, pred.optim.Alt1.third.nls, NULL, T, control=list(maxit=100000))
    print(pp.Alt1.fit)
    repeat{
      par.optim.Alt1.fit<-pp.Alt1.fit$par
      pp.Alt1.fit<-optim(par.optim.Alt1.fit, sumsq.optim, gr=NULL, subset, pred.optim.Alt1.third.nls, NULL, T, control=list(maxit=100000))
      par.optim.Alt1.fit1<-pp.Alt1.fit$par
      identical<-sum((par.optim.Alt1.fit - par.optim.Alt1.fit1)^2) ==0
      if(identical) break
    }
    print(pp.Alt1.fit)
    # save estimates for m0 fit
    pred.Alt1<-predict.Alt1.third.nls(pp.Alt1.fit$par[1], pp.Alt1.fit$par[2], pp.Alt1.fit$par[3], pp.Alt1.fit$par[4], subset)
    
    mod.Alt1.results<-c(mod.Alt1.results, list(list(A.estimate=sigmoid(pp.Alt1.fit$par[1]), f.estimate=(pp.Alt1.fit$par[2]), k.estimate=(pp.Alt1.fit$par[3]),
                                                    m0.estimate=pp.Alt1.fit$par[4], sumsq=pp.Alt1.fit$value, counts=pp.Alt1.fit$counts, 
                                                    convergence=pp.Alt1.fit$convergence, prediction.Alt1=pred.Alt1)))
    
    A.estimates.Alt1<-c(A.estimates.Alt1, list(list(A.estimate=sigmoid(pp.Alt1.fit$par[1]))))
    
    sumsq.Alt1<-c(sumsq.Alt1, list(list(sumsq.Alt1=pp.Alt1.fit$value)))
    
    ### 4: A=2/3 MODEL
    # 4a: m0 fixed
    par.optim.23.fix<-c(pp.Alt1.fit$par[2], pp.Alt1.fit$par[3]) 
    pp.23.fix<- optim(par.optim.23.fix, sumsq.optim, gr=NULL, subset, pred.optim.A23.third.nls, m0, T, control=list(maxit=10000000)) # Laura
    
    repeat{
      par.optim.23.fix<-pp.23.fix$par
      pp.23.fix<-optim(par.optim.23.fix, sumsq.optim, gr=NULL, subset, pred.optim.A23.third.nls, m0, T, control=list(maxit=10000000))
      par.optim.23.fix1<-pp.23.fix$par
      identical<-sum((par.optim.23.fix - par.optim.23.fix1)^2) ==0
      if(identical) break
    }
    # 4b: m0 fit
    par.optim.23.fit<- c(pp.23.fix$par[1], pp.23.fix$par[2], m0) # (f, k, m0) 
    pp.23.fit<- optim(par.optim.23.fit, sumsq.optim, gr=NULL, subset, pred.optim.A23.third.nls, NULL, T, control=list(maxit=100000))
    repeat{
      par.optim.23.fit<-pp.23.fit$par
      pp.23.fit<-optim(par.optim.23.fit, sumsq.optim, gr=NULL, subset, pred.optim.A23.third.nls, NULL, T, control=list(maxit=100000))
      par.optim.23.fit1<-pp.23.fit$par
      identical<- sum((par.optim.23.fit - par.optim.23.fit1)^2) ==0
      if(identical) break
    }
    print(pp.23.fit)
    pred.A23<-predict.Alt1.third.nls((inv.sigmoid(2/3)), pp.23.fit$par[1], pp.23.fit$par[2], pp.23.fit$par[3], subset)
    
    mod.A23.results<-c(mod.A23.results, list(list(f.estimate=(pp.23.fit$par[1]), k.estimate=(pp.23.fit$par[2]), m0.estimate=pp.23.fit$par[3],
                                                  sumsq=pp.23.fit$value, counts=pp.23.fit$counts, convergence=pp.23.fit$convergence,prediction.A23=pred.A23)))
    sumsq.A23<-c(sumsq.A23, list(list(sumsq.A23=pp.23.fit$value)))
    
    ### 5: A>1 MODEL
    # 5a: m0 fixed
    par.t.Agt1<- start.params.Agt1.from.exp(0.05, subset) 
    par.optim.Agt1.fix<- unlist(par.t.Agt1)
    pp.Agt1.fix<- optim(par.optim.Agt1.fix, sumsq.optim, gr=NULL, subset, pred.optim.Agt1.nls, m0, T, control=list(maxit=100000))
    repeat{
      par.optim.Agt1.fix<-pp.Agt1.fix$par
      pp.Agt1.fix<-optim(par.optim.Agt1.fix, sumsq.optim, gr=NULL, subset, pred.optim.Agt1.nls, m0, T, control=list(maxit=100000))
      par.optim.Agt1.fix1<-pp.Agt1.fix$par
      identical<-sum((par.optim.Agt1.fix - par.optim.Agt1.fix1)^2) <=1e-20
      if(identical) break
    }
    # 5b: m0 fit
    par.optim.Agt1.fit<- c(pp.Agt1.fix$par[1],pp.Agt1.fix$par[2],pp.Agt1.fix$par[3], m0) # (alpha, Z, s, m0)
    pp.Agt1.fit<- optim(par.optim.Agt1.fit, sumsq.optim, gr=NULL, subset, pred.optim.Agt1.nls, NULL, T, control=list(maxit=100000))
    repeat{
      par.optim.Agt1.fit<-pp.Agt1.fit$par
      pp.Agt1.fit<-optim(par.optim.Agt1.fit, sumsq.optim, gr=NULL, subset, pred.optim.Agt1.nls, NULL, T, control=list(maxit=100000))
      par.optim.Agt1.fit1<-pp.Agt1.fit$par
      identical<- sum((par.optim.Agt1.fit - par.optim.Agt1.fit1)^2) <=1e-20
      if(identical) break
    }
    
    pred.Agt1<-predict.Agt1.nls(pp.Agt1.fit$par[1], pp.Agt1.fit$par[2], pp.Agt1.fit$par[3], pp.Agt1.fit$par[4], subset)
    
    mod.Agt1.results<-c(mod.Agt1.results, list(list(Alpha.estimate=(pp.Agt1.fit$par[1]), Z.estimate=pp.Agt1.fit$par[2], s.estimate=pp.Agt1.fit$par[3],
                                                    m0.estimate=pp.Agt1.fit$par[4], sumsq=pp.Agt1.fit$value, counts=pp.Agt1.fit$counts,
                                                    convergence=pp.Agt1.fit$convergence, predict.Agt1=pred.Agt1)))
    
    A.estimates.Agt1<-c(A.estimates.Agt1, list(list(A.estimate=(1/sigmoid(pp.Agt1.fit$par[1])))))
    
    sumsq.Agt1<-c(sumsq.Agt1, list(list(sumsq.Agt1=pp.Agt1.fit$value)))
  }
  
  return(list(Exp=mod.exp.results, Gomp=mod.gomp.results, Alt1=mod.Alt1.results, A23=mod.A23.results, Agt1=mod.Agt1.results, A.estimates.Alt1=A.estimates.Alt1, A.estimates.Agt1=A.estimates.Agt1, 
              sumsq.exp=sumsq.exp, sumsq.gomp=sumsq.gomp, sumsq.Alt1=sumsq.Alt1, sumsq.A23=sumsq.A23,sumsq.Agt1, sumsq.Agt1=sumsq.Agt1))
  
}


##### FUNCTION FOR LOOPING THE 5 MODELS OVER SEVERAL SPECIES ####
### the following is code for an 'all-in-one' function to run all models for 
# several species by looping over species 

unif.growth.mod.spec<- function(the.data){
  species<-the.data$species
  species.loop<-unique(the.data$species)
  growth.species<-list() # create empty list to store results 
  
  for(i in 1:length(species.loop)){ # looping over species
    species.loop<-unique(the.data$species)
    subset<-the.data[the.data$species == species.loop[i], ]
    
    growth.mod.spec<- unif.growth.mod(subset) 
    
    
    
    growth.species<- c(growth.species, list(growth.mod.spec)) # store all results in a list
    
  }
  
  return(growth.species) 
}


#### END #####
