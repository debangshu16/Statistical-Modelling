################## Problem 2 #######################
cat ("Problem 2")

library(datasets)

plot(faithful,pch=20)

hist(faithful$eruptions, nclass=20,probability = T)

summary(faithful$eruptions)

length(faithful$eruptions[faithful$eruptions<3])/length(faithful$eruptions)
#So, about 35.6% of the data is less than 3

### Model eruptions as a mixture of one gamma and one normal distribution

## X : eruptions time
## X ~ f(x) = p* Gamma(alpha,beta) + (1-p)*Normal(miu, sigma)

## theta =(p,alpha, beta, miu, sigma)

dHalfCauchy = function(theta,log=FALSE){
  if(log==FALSE){
    return(2/(pi*(1+theta^2)))
  }else{
    return(log(2)-log(pi)-log(1+theta^2))
  }
}

dmix = function(x,p,alpha,beta,miu, sigma){
  pdf = p*dgamma(x,shape = alpha, scale = beta)+(1-p)*dnorm(x,mean = miu, sd = sigma)
  return(pdf)
}

dmix(4, 0.5, 2, 1, 4.5, 2)

logLikeMix = function(x,theta){
  p = theta[1]
  alpha = theta[2]
  beta = theta[3]
  miu = theta[4]
  sigma = theta[5]
  
  ll = sum(log(dmix(x=x,p=p,
                    alpha = alpha,beta = beta,
                    miu = miu,sigma = sigma)))
  return(ll)
  
}

theta.init = c(0.5, 3, 1, 4, 1)

logLikeMix(x=faithful$eruptions,theta = theta.init)


### define prior

### p ~ unif(0,1)
### alpha, miu ~ Cauchy(0,1)
### beta, sigma ~  Half Cauchy(0,1)




logPrior = function(theta)
{
  p = theta[1]
  alpha = theta[2]
  beta = theta[3]
  miu = theta[4]
  sigma = theta[5]
  
  lp = dunif(p,log = TRUE)+ dcauchy(alpha,log = TRUE)+dcauchy(miu, log = TRUE)+ dHalfCauchy(beta, log = TRUE) + dHalfCauchy(sigma, log = TRUE)
  return(lp)
}


logPrior(theta = theta.init)

logPosterior = function(x,theta){
  lp = logLikeMix(x=x,theta=theta)+logPrior(theta=theta)
  return(lp)
}

logPosterior(x=faithful$eruptions,theta = theta.init)


##################
set.seed(100)
sim.size = 10000
burn.in =2000
x=faithful$eruptions

theta.star = matrix(NA,nrow = (burn.in+sim.size), ncol = 5)
colnames(theta.star) =  c('p','alpha','beta','miu','sigma')


theta.star[1,] = theta.init


for( i in 2:(burn.in+sim.size)){
  if( i %% 100 ==0) cat('i = ',i,'\n')
  
  ### Metropolis-Hasting for p
  min_p = theta.star[i-1,'p']-0.05
  max_p = theta.star[i-1,'p']+0.05
  
  p.proposal = runif(1,min = min_p,max = max_p)
  
  proposal = theta.star[i-1,]
  proposal['p'] = p.proposal
  r =exp(logPosterior(x=x,theta=proposal) - logPosterior(x = x,theta=theta.star[i-1,]))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'p'] = p.proposal
  }else{
    theta.star[i,'p'] = theta.star[i-1,'p']
  }
  
  ### Metropolis-Hasting for alpha
  min_alpha = theta.star[i-1,'alpha']-2
  max_alpha = theta.star[i-1,'alpha']+2
  
  alpha.proposal = runif(1,min = min_alpha,max = max_alpha)
  
  proposal = theta_hat = c(theta.star[i,'p'],theta.star[i-1,c('alpha','beta','miu','sigma')])
  proposal['alpha'] = alpha.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'alpha'] = alpha.proposal
  }else{
    theta.star[i,'alpha'] = theta_hat['alpha']
  }
  
  ### Metropolis-Hasting for beta
  min_beta = theta.star[i-1,'beta']-0.01
  max_beta = theta.star[i-1,'beta']+0.01
  
  beta.proposal = runif(1,min = min_beta,max = max_beta)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha')],theta.star[i-1,c('beta','miu','sigma')])
  proposal['beta'] = beta.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'beta'] = beta.proposal
  }else{
    theta.star[i,'beta'] = theta_hat['beta']
  }
  
  
  ### Metropolis-Hasting for miu
  min_miu = theta.star[i-1,'miu']-0.5
  max_miu = theta.star[i-1,'miu']+0.5
  
  miu.proposal = runif(1,min = min_miu,max = max_miu)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha','beta')],theta.star[i-1,c('miu','sigma')])
  proposal['miu'] = miu.proposal
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'miu'] = miu.proposal
  }else{
    theta.star[i,'miu'] = theta_hat['miu']
  }
  
  ### Metropolis-Hasting for sigma
  min_sigma = theta.star[i-1,'sigma']-0.01
  max_sigma = theta.star[i-1,'sigma']+0.01
  
  sigma.proposal = runif(1,min = min_sigma,max = max_sigma)
  
  proposal = theta_hat = c(theta.star[i,c('p','alpha','beta','miu')],theta.star[i-1,c('sigma')])
  proposal['sigma'] = sigma.proposal
  
  r =exp(logPosterior(theta=proposal,x=x) - logPosterior(theta=theta_hat,x=x))
  prob = min(r,1)
  if (runif(1) < prob){
    theta.star[i,'sigma'] = sigma.proposal
  }
  else{
    theta.star[i,'sigma'] = theta_hat['sigma']
  }
  
  
}


theta.star = data.frame(theta.star)

theta_hat = apply(theta.star[(burn.in+1):(burn.in+sim.size),],2,mean)

plot(ts(theta.star$p))
hist(theta.star$p)

plot(ts(theta.star$alpha))
hist(theta.star$alpha)

plot(ts(theta.star$beta))
hist(theta.star$beta)

plot(ts(theta.star$miu))
hist(theta.star$miu)

plot(ts(theta.star$sigma))
hist(theta.star$sigma)

######### Bayesian Density Estimate ##########

eruptions = sort(faithful$eruptions)
hist(eruptions, nclass=20,probability = T)
lines(density(eruptions),col='blue',lwd=2)

y_hats = matrix(0, nrow = 100, ncol = length(eruptions))
for(j in 1:100){
  
  j1 = j+burn.in
  pdf = dmix(x=eruptions,p=theta.star$p[j1]
             ,alpha=theta.star$alpha[j1]
             ,beta=theta.star$beta[j1]
             ,miu=theta.star$miu[j1]
             ,sigma=theta.star$sigma[j1])
  lines(eruptions,pdf,lwd=2,col='red')
  y_hats[j, ] = pdf
}

density_estimate = apply(y_hats, 2, mean)
lines(eruptions, density_estimate, col = 'black', lwd = 2)

########  P(0<eruptions<3 | data) #############

sum_prob = 0
for (j in 1:sim.size)
{
  j1 = burn.in + j
  p = theta.star$p[j1]
  alpha = theta.star$alpha[j1]
  beta = theta.star$beta[j1]
  miu = theta.star$miu[j1]
  sigma = theta.star$sigma[j1]
  
  
  prob = integrate(dmix, lower = 0, upper = 3, p = p, alpha = alpha, beta = beta,
                   miu = miu, sigma = sigma)
  
  sum_prob = sum_prob + prob$value
}

sum_prob = sum_prob/sim.size
sum_prob


#####Posterior Mode using optim routine#######

param = theta.init
data<- faithful$eruptions

negative_log_posterior <- function(param)
{
  return (-logPosterior(data, param))
}
set.seed(100)
options(warn = -1)
negativeLogPosterior_optimization = optim(theta.init, negative_log_posterior)
posterior_mode = negativeLogPosterior_optimization$par

posterior_mode