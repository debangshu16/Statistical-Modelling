################ Problem 3  ###################################
cat("Problem 3")

#Read data
#par(mfrow = c(1,1))
data_train=AirPassengers[1:96]
data_test = AirPassengers[97:144]
tme = time(AirPassengers)
tme_min = min(tme)
tme_max =max(tme)
plot(NULL,xlim=c(tme_min,tme_max),ylim=c(100,650)
     ,xlab='',ylab='Air Passengers')
lines(tme[1:96],data_train,lwd=2,col='purple')
lines(tme[97:144],data_test,lwd=2,col='blue')

#Represent data as dataframe
data.train = cbind.data.frame(time=tme[1:96],AirPassengers=data_train)
data.test = cbind.data.frame(time=tme[97:144],AirPassengers=data_test)
data.train$t = data.train$time -1949
data.test$t = data.test$t - 1949
head(data.train)

tail(data.train)
head(data.test)
tail(data.test)

#fitting model y = exp(alpha0 + alpha1*t + beta1*sin(omega*t) + beta2*cos(omega*t) + epsilon)
# where omega = 2*pi and epsilon follows N(0, sigma^2)

omega=2*pi
fit.lm = lm(log(AirPassengers)~t+sin(omega*t)+cos(omega*t),data=data.train)
summary(fit.lm)

#predict test data

data.test$y_hat = exp(predict(fit.lm,newdata = data.test))
head(data.test)

#Plotting prediction
plot(NULL,xlim=c(tme_min,tme_max),ylim=c(100,650)
     ,xlab='',ylab='Air Passengers')
lines(tme[1:96],data_train,lwd=2,col='purple')
lines(tme[97:144],data_test,lwd=2,col='blue')
lines(tme[97:144],data.test$y_hat,lwd=2,col='red')

############Problem 3#################

#a) Assume Cauchy(0,1) prior on alpha0, alpha1, beta1, beta2 and Half-Cauchy(0,1) prior
#on sigma and find posterior mode using optim routine.

dHalfCauchy = function(theta,log=FALSE){
  if(log==FALSE){
    return(2/(pi*(1+theta^2)))
  }else{
    return(log(2)-log(pi)-log(1+theta^2))
  }
}
logLikelihood<- function(param, t, y)
{
  #y = exp(alpha0 + alpha1*t + beta1*sin(omega*t) + beta2*cos(omega*t) + epsilon
  omega = 2*pi
  alpha0  = param[1]
  alpha1 = param[2]
  beta1 = param[3]
  beta2 = param[4]
  sigma = param[5]
  
  pred = alpha0 + alpha1 * t + beta1* sin(omega*t) + beta2 * cos(omega*t)
  res = log(y) - pred
  
  #res follows N(0, sigma^2) : pdf = (1/root(2*pi*sigma))*exp(-(x^2)/2*(sigma^2))
  
  log_likelihoods = -(0.5*(log(2) + log(pi)))  - log(sigma) - (1/2*((res/sigma)^2)) 
  #og_likelihoods = dnorm(res, 0, sigma, log = T)
  return (sum(log_likelihoods))
}

logPrior<-function(param)
{
  alpha0  = param[1]
  alpha1 = param[2]
  beta1 = param[3]
  beta2 = param[4]
  sigma = param[5]
  
  alpha_priors = dcauchy(alpha0,0,1, log = T) + dcauchy(alpha1,0,1, log = T)
  beta_priors = dcauchy(beta1,0,1, log = T) + dcauchy(beta2,0,1,log=T)
  sigma_prior = dHalfCauchy(sigma, log = T)
  
  return (alpha_priors + beta_priors + sigma_prior)
}

logPosterior<- function(param, x, y)
{
  like = logLikelihood(param=param,t = x, y = y)
  prior = logPrior(param=param)
  post  = like + prior
  return ( post )
}

param.init = c(0, 1, 1, 1, 1)
logLikelihood(param.init, t = data.train$t, y = data.train$AirPassengers)
logPrior(param.init)
logPosterior(param.init, x = data.train$t, y = data.train$AirPassengers)

negativeLogPosterior<- function(param)
{
  return (- logPosterior(param, x = data.train$t, y = data.train$AirPassengers))
}

posteriorMode = optim(param.init, negativeLogPosterior)$par

posteriorMode
OLS_params = c(fit.lm$coefficients, sigma(fit.lm))
OLS_params
## The MLE estimate of linear model fit and the posterior mode estimates almost match upto the decimal

## b) MCMC sampling

#install.packages('mvtnorm')
library('mvtnorm')

proposalfunction = function(param,x){
  X=cbind(rep(1,length(x)),x, sin(omega*x), cos(omega*x))
  
  S=(param[5]^2)*solve(t(X)%*%X)
  prop = c(rmvnorm(1
                   ,mean = param[1:4]
                   ,sigma = S)
           ,rgamma(1,param[5]*5,5))
  return(prop)
}

run_metropolis = function(startvalue, N.sim, burnin){
  iterations = N.sim + burnin
  chain = array(dim = c(iterations+1,5))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,],x=x)
    
    probab = exp(logPosterior(param=proposal
                              ,y=y,x=x) 
                 - logPosterior(param=chain[i,]
                                ,y=y,x=x))
    
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}
omega = 2*pi
options(warn=-1)
y=data.train$AirPassengers
x=data.train$t
startvalue = c(1, 1, 1, 1, 1)
N.sim=10000
burnin=1000

set.seed(1)
chain = run_metropolis(startvalue=startvalue
                       ,N.sim=N.sim
                       ,burnin=burnin)
colnames(chain)=c("alpha0","alpha1","beta1", "beta2", "sigma")
converged_chain=chain[(burnin+1):nrow(chain),]
options(warn=0)

### c) Trace and Density plot
par(mfrow=c(5,2))
plot(ts(chain[,"alpha0"]),lwd=2
     ,ylab = "alpha0")
plot(density(chain[,"alpha0"]),lwd=2
     ,main = "",col='red'
     ,xlab = "alpha0")

plot(ts(chain[,"alpha1"]),lwd=2
     ,ylab = "alpha1")
plot(density(chain[,"alpha1"]),lwd=2
     ,main = "",col='red'
     ,xlab = "alpha1")

plot(ts(chain[,"beta1"]),lwd=2
     ,ylab = "beta1")
plot(density(chain[,"beta1"]),lwd=2
     ,main = "",col='red'
     ,xlab = "beta")

plot(ts(chain[,"beta2"]),lwd=2
     ,ylab = "beta2")
plot(density(chain[,"beta2"]),lwd=2
     ,main = "",col='red'
     ,xlab = "beta2")

plot(ts(chain[,"sigma"]),lwd=2
     ,ylab = "sigma")
plot(density(chain[,"sigma"]),lwd=2
     ,main = "",col='red'
     ,xlab = "sigma")

theta_hat = apply(converged_chain, 2, mean)

## d) First 100 predicted path
par(mfrow = c(1,1))
plot(NULL,xlim=c(tme_min,tme_max),ylim=c(100,650)
     ,xlab='',ylab='Air Passengers')
lines(tme[1:96],data_train,lwd=2,col='purple')
#lines(tme[97:144],data_test,lwd=2,col='black')

n = nrow(data.test)
y_hats = matrix(data = 0,nrow = n, ncol = N.sim)


for (i in 1:100)
{
  j1 = burnin + i
  x = data.test$t
  
  theta_i = chain[j1,]
  
  alpha0 = theta_i[1]
  alpha1 = theta_i[2]
  beta1 = theta_i[3]
  beta2 = theta_i[4]
  sigma = theta_i[5]
  
  eps_t = rnorm(length(x), mean = 0, sd = sigma)
  y = exp(alpha0 + (alpha1*x) + (beta1*sin(omega*x)) + (beta2*cos(omega*x)) + eps_t)
  lines(tme[97:144],y, type = 'l', lwd = 2, col = 'cyan')
}

for (i in 1:N.sim)
{
  j1 = burnin + i
  x = data.test$t
  
  theta_i = chain[j1,]
  
  alpha0 = theta_i[1]
  alpha1 = theta_i[2]
  beta1 = theta_i[3]
  beta2 = theta_i[4]
  sigma = theta_i[5]
  
  eps_t = rnorm(length(x), mean = 0, sd = sigma)
  y = exp(alpha0 + (alpha1*x) + (beta1*sin(omega*x)) + (beta2*cos(omega*x)) + eps_t)
  y_hats[,i] = y
}

## e) 95% Bayesian Interval
y_hat_lb = apply(y_hats, 1, quantile, 0.025)
y_hat_ub = apply(y_hats, 1, quantile, 0.975)
y_hat = apply(y_hats, 1, mean)

lines(tme[97:144],y_hat_lb, type = 'l', lwd = 2, col = 'black')
lines(tme[97:144],y_hat_ub, type = 'l', lwd = 2, col = 'black')
lines(tme[97:144],y_hat, type = 'l', lwd = 2, col ='green')
lines(tme[97:144],data_test,lwd=2,col='red')
legend('topleft',legend = c('Train data','Simulated Predicted paths',
                            'Confidence bounds','Mean Prediction path',
                            'Test data'),
       col = c('purple','cyan','black','green','red'), lty = 1,
       cex = 0.8
)

##3.f) If we increase the number of seasonality components,
#i.e, instead of only sin(omega*t) + cos(omega*t) we now add
#sin (omega*t) + cos(omega*t) + sin(2*omega*t) + cos(2*omega*t) and so on.

## Showing difference for OLS fit

#current model summary
summary(fit.lm)

rmse_test = sqrt(sum((data.test$AirPassengers - data.test$y_hat)**2)/nrow(data.test)) 
rmse_test

#current linear model's test data RMSE = 63.93219

# Increasing model seasonal components:
fit.lm = lm(log(AirPassengers)~t+sin(omega*t)+cos(omega*t) + sin(2*omega*t) + cos(2*omega*t),data=data.train)
summary(fit.lm)

#We can see that this model's adjusted R2(0.97) is better than the previous model's adjusted R2(0.94).
data.test$y_hat = exp(predict(fit.lm,newdata = data.test))
rmse_test = sqrt(sum((data.test$AirPassengers - data.test$y_hat)**2)/nrow(data.test)) 
rmse_test

#Also the RMSE has decreased for this new model. Similarly by increasing the seasonality components sin cos terms with increasing omega's (omega, 2*omega, 3*omega and so on) , the model gets better.

#new model prediction plot for OLS
plot(NULL,xlim=c(tme_min,tme_max),ylim=c(100,650)
     ,xlab='',ylab='Air Passengers')
lines(tme[1:96],data_train,lwd=2,col='purple')
lines(tme[97:144],data_test,lwd=2,col='blue')
lines(tme[97:144],data.test$y_hat,lwd=2,col='red')

#by adding one more seasonal component it gets even better

fit.lm = lm(log(AirPassengers)~t+sin(omega*t)+cos(omega*t) + 
              sin(2*omega*t) + cos(2*omega*t)
            + sin(3*omega*t) + cos(3*omega*t),data=data.train)
summary(fit.lm)

#We can see that this model's adjusted R2 is even better.
data.test$y_hat = exp(predict(fit.lm,newdata = data.test))
rmse_test = sqrt(sum((data.test$AirPassengers - data.test$y_hat)**2)/nrow(data.test)) 
rmse_test

#Also the RMSE has decreased even further.
#new model prediction plot for OLS
plot(NULL,xlim=c(tme_min,tme_max),ylim=c(100,650)
     ,xlab='',ylab='Air Passengers')
lines(tme[1:96],data_train,lwd=2,col='purple')
lines(tme[97:144],data_test,lwd=2,col='blue')
lines(tme[97:144],data.test$y_hat,lwd=2,col='red')
