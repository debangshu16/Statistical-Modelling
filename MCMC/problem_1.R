################## Problem 1############
cat ("Problem 1")

library(datasets)

plot(faithful,pch=20)

hist(faithful$waiting, nclass=20,probability = T)

summary(faithful$waiting)

length(faithful$waiting[(faithful$waiting<88) & (faithful$waiting > 72)])/length(faithful$waiting)
#so there is about 50.7% of data which is between 72 and 88

### Model waiting as a mixture of two gamma distributions

## X : waiting time
## X ~ f(x) = p* Gamma(alpha1,beta1) + (1-p)*Gamma(alpha2,beta2)

## theta =(p,alpha1,beta1,alpha2,beta2)


dmixGamma = function(x,p,alpha1,beta1,alpha2,beta2){
  pdf = p*dgamma(x,shape = alpha1, scale = beta1)+(1-p)*dgamma(x,shape = alpha2, scale = beta2)
  return(pdf)
}

dHalfCauchy = function(theta,log=FALSE){
  if(log==FALSE){
    return(2/(pi*(1+theta^2)))
  }else{
    return(log(2)-log(pi)-log(1+theta^2))
  }
}

x=faithful$waiting

logLikeMixGamma = function(x,theta){
  p = theta[1]
  alpha1 = theta[2]
  beta1 = theta[3]
  alpha2 = theta[4]
  beta2 = theta[5]
  
  ll = sum(log(dmixGamma(x=x,p=p,
                         alpha1 = alpha1,beta1 = beta1,
                         alpha2 = alpha2,beta2 = beta2)))
  return(ll)
  
}

theta.init = c(0.5,50,1,70,1)

logLikeMixGamma(x=faithful$waiting,theta = theta.init)


### define prior

### p ~ unif(0,1)
### alpha1,alpha2 ~ Half Cauchy(0,1)
### beta1, beta2 ~ Half Cauchy(0,1)




logPrior = function(theta)
{
  p = theta[1]
  alpha1 = theta[2]
  beta1 = theta[3]
  alpha2 = theta[4]
  beta2 = theta[5]
  
  lp = dunif(p,log = TRUE)+ dHalfCauchy(alpha1,log = TRUE)+dHalfCauchy(alpha2, log = TRUE)+dHalfCauchy(beta1, log = T) + dHalfCauchy(beta2,log = TRUE)
  return(lp)
}


logPrior(theta = theta.init)

logPosterior = function(x,theta){
  lp = logLikeMixGamma(x=x,theta=theta)+logPrior(theta=theta)
  return(lp)
}

logPosterior(x=faithful$waiting,theta = theta.init)


################## MCMC
set.seed(100)

sim.size = 10000
burn.in =2000

theta.star = matrix(NA,nrow = (burn.in+sim.size), ncol = 5)
colnames(theta.star) =  c('p','alpha1','beta1','alpha2','beta2')


theta.star[1,] = theta.init

#for convergence of parameters, I played around with the propsal intervals
# Burnin of 2000 is required for these values

for( i in 2:(burn.in+sim.size)){
  if( i %% 100 ==0) cat('i = ',i,'\n')
  
  ### Metropolis-Hasting for p
  min_p = theta.star[i-1,'p']-0.05
  max_p = theta.star[i-1,'p']+0.05
  
  p.proposal = runif(1,min = min_p,max = max_p)
  
  proposal = theta.star[i-1,]
  proposal['p'] = p.proposal
  r =exp(logPosterior(x=x,theta=proposal) - logPosterior(x=x,theta=theta.star[i-1,]))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'p'] = p.proposal
  }else{
    theta.star[i,'p'] = theta.star[i-1,'p']
  }
  
  ### Metropolis-Hasting for alpha1
  min_alpha1 = theta.star[i-1,'alpha1']-2
  max_alpha1 = theta.star[i-1,'alpha1']+2
  
  alpha1.proposal = runif(1,min = min_alpha1,max = max_alpha1)
  
  proposal = theta_hat = c(theta.star[i,'p'],theta.star[i-1,c('alpha1','beta1','alpha2','beta2')])
  proposal['alpha1'] = alpha1.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'alpha1'] = alpha1.proposal
  }else{
    theta.star[i,'alpha1'] = theta_hat['alpha1']
  }
  
  ### Metropolis-Hasting for beta1
  min_beta1 = theta.star[i-1,'beta1']-0.05
  max_beta1 = theta.star[i-1,'beta1']+0.05
  
  beta1.proposal = runif(1,min = min_beta1,max = max_beta1)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha1')],theta.star[i-1,c('beta1','alpha2','beta2')])
  proposal['beta1'] = beta1.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'beta1'] = beta1.proposal
  }else{
    theta.star[i,'beta1'] = theta_hat['beta1']
  }
  
  
  ### Metropolis-Hasting for alpha2
  min_alpha2 = theta.star[i-1,'alpha2']-2
  max_alpha2 = theta.star[i-1,'alpha2']+2
  
  alpha2.proposal = runif(1,min = min_alpha2,max = max_alpha2)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha1','beta1')],theta.star[i-1,c('alpha2','beta2')])
  proposal['alpha2'] = alpha2.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'alpha2'] = alpha2.proposal
  }else{
    theta.star[i,'alpha2'] = theta_hat['alpha2']
  }
  
  ### Metropolis-Hasting for beta2
  min_beta2 = theta.star[i-1,'beta2']-0.05
  max_beta2 = theta.star[i-1,'beta2']+0.05
  
  beta2.proposal = runif(1,min = min_beta2,max = max_beta2)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha1','beta1','alpha2')],theta.star[i-1,c('beta2')])
  proposal['beta2'] = beta2.proposal
  
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'beta2'] = beta2.proposal
  }
  else{
    theta.star[i,'beta2'] = theta_hat['beta2']
  }
  
  
}


theta.star = data.frame(theta.star)

theta_hat = apply(theta.star[(burn.in+1):(burn.in+sim.size),],2,mean)

plot(ts(theta.star$p))
hist(theta.star$p)

plot(ts(theta.star$alpha1))
hist(theta.star$alpha1)

plot(ts(theta.star$beta1))
hist(theta.star$beta1)

plot(ts(theta.star$alpha2))
hist(theta.star$alpha2)

plot(ts(theta.star$beta2))
hist(theta.star$beta2)

######### Bayesian Density Estimate ##########

waiting = sort(faithful$waiting)
hist(waiting, nclass=20,probability = T,ylim=c(0,0.08))
lines(density(waiting),col='blue',lwd=2)

y_hats = matrix(0, nrow = 100, ncol = length(waiting))
for(j in 1:100){
  
  j1 = j+burn.in
  pdf = dmixGamma(x=waiting,p=theta.star$p[j1]
                  ,alpha1=theta.star$alpha1[j1]
                  ,beta1=theta.star$beta1[j1]
                  ,alpha2=theta.star$alpha2[j1]
                  ,beta2=theta.star$beta2[j1])
  lines(waiting,pdf,lwd=2,col='red')
  y_hats[j,] = pdf 
}

density_estimate = apply(y_hats, 2, mean)
lines(waiting, density_estimate, col = 'black', lwd = 2)

########  P(72<waiting<88 | data) #############

sum_prob = 0
for (j in 1:sim.size)
{
  j1 = burn.in + j
  p = theta.star$p[j1]
  alpha1 = theta.star$alpha1[j1]
  beta1 = theta.star$beta1[j1]
  alpha2 = theta.star$alpha2[j1]
  beta2 = theta.star$beta2[j1]
  
  
  prob = integrate(dmixGamma, lower = 72, upper = 88, p = p, alpha1 = alpha1, beta1 = beta1,
                   alpha2 = alpha2, beta2 = beta2)
  
  sum_prob = sum_prob + prob$value
}

sum_prob = sum_prob/sim.size
sum_prob
#So our probability value of y between 72 and 88 given data is very close to the value of the data.
