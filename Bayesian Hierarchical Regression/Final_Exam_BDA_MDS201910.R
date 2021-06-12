### Bayesian Data Analysis 
### Final Big Assignment
### Instructor : Sourish Das

### Student Name: Debangshu Bhattacharya
### Roll Number: MDS201910
### Email: debangshub@cmi.ac.in

### Please save the file as Final_Exam_BDA_ROLL_Number.R
### For exam if my roll number is MDS202000
### then I will save the file as Final_Exam_BDA_MDS202000.R

### Last Date of Submission : 30th April, 2021

### Following data is downloading the global COVID-19 data repository 
### from the John Hopkins University
### You have two tasks to complete in this assignment

### Task 1: Consider the data since 1st Dec 2020. 

rm(list=ls())

## loading data
state_wise=read.csv(file="https://api.covid19india.org/csv/latest/state_wise_daily.csv",header=TRUE,stringsAsFactors = F)
dim(state_wise)
tail(state_wise,n=10)

test_data = read.csv(file = 'https://api.covid19india.org/csv/latest/statewise_tested_numbers_data.csv')

## data processing prepare the data so that we can implement
## the Bayesian hierarchical regression model

## "Dadra and Nagar Haveli and Daman and Diu" are considered
##  as single UT in test data
##  In state-wise data it is considered as two sepearate UT
##  We decide to drop the "Dadra and Nagar Haveli and Daman and Diu"
##  from the analysis

state_nm_test_data = unique(test_data$State)
state_nm_test_data = state_nm_test_data[-8]

state_nm_state_data = colnames(state_wise)[5:ncol(state_wise)]
state_nm_state_data = state_nm_state_data[-c(8,9,38)]

state_nm = cbind.data.frame(test_data=state_nm_test_data
                            ,state_data=state_nm_state_data)


## prepare the data with confirmed cases and 
## test positivity rate

data_confirmed<-state_wise[state_wise$Status=='Confirmed',]



data_confirmed$Date_YMD = as.Date(data_confirmed$Date_YMD,format='%Y-%m-%d')
data_confirmed=data_confirmed[data_confirmed$Date_YMD>=as.Date("2020-12-01",format="%Y-%m-%d"),]




test_data$Updated.On=as.Date(test_data$Updated.On,format = '%d/%m/%Y')

test_data=test_data[test_data$Updated.On>=as.Date('2020-12-01',format='%Y-%m-%d'),]


i = 1
test_data_sub=test_data[test_data$State==state_nm[i,'test_data']
                        ,c('Updated.On','State','Total.Tested')]

#tail(test_data_sub)
col_nm = c('Date_YMD',state_nm[i,'state_data'])
data_confirmed_sub=data_confirmed[,col_nm]
data_sub = merge(test_data_sub,data_confirmed_sub
                 ,by.x="Updated.On",by.y='Date_YMD')
colnames(data_sub)[4] = 'new_cases'
data_sub$test_postivity_rate=  ((data_sub$new_cases+1)/(data_sub$Total.Tested+1))*100
data_final = data_sub

for(i in 2:nrow(state_nm)){
  test_data_sub=test_data[test_data$State==state_nm[i,'test_data']
                          ,c('Updated.On','State','Total.Tested')]
  
  #tail(test_data_sub)
  col_nm = c('Date_YMD',state_nm[i,'state_data'])
  data_confirmed_sub=data_confirmed[,col_nm]
  data_sub = merge(test_data_sub,data_confirmed_sub
                   ,by.x="Updated.On",by.y='Date_YMD')
  colnames(data_sub)[4] = 'new_cases'
  data_sub$test_postivity_rate= ((data_sub$new_cases+1)/(data_sub$Total.Tested+1))*100
  data_final = rbind.data.frame(data_final,data_sub)
}

### clean data where new cases < 0
n = nrow(data_final)
for (i in 1:n)
{
  if (data_final$new_cases[i] < 0)
  {
     data_final$new_cases[i] =  - data_final$new_cases[i]
     data_final$test_postivity_rate[i] = ((data_final[i, "new_cases"]+1)/(data_final[i, "Total.Tested"]+1))*100
  }
}
#data_final$test_postivity_rate= ((data_final$new_cases+1)/(data_final$Total.Tested+1))*100

# Remove missing data
data_final[rowSums(is.na(data_final))>0,]
data_final<- na.omit(data_final)

########################################
plot(data_final$test_postivity_rate
     ,data_final$new_cases,pch=20
     ,xlab='test postivity rate'
     ,ylab='new cases')

###########

plot(log(data_final$test_postivity_rate)
     ,log(data_final$new_cases),pch=20
     ,xlab='test postivity rate in log scale'
     ,ylab='new cases in log scale'
     ,col='purple')

fit = lm(log(new_cases+1)~log(test_postivity_rate)
         ,data=data_final)## 1 added to new case to avoid log(0)
abline(fit,col='black',lwd=2,lty=2)

points(log(data_final$test_postivity_rate[data_final$State=="Maharashtra"])
       ,log(data_final$new_cases[data_final$State=="Maharashtra"])
       ,col='orange',pch=20)

points(log(data_final$test_postivity_rate[data_final$State=="Meghalaya"])
       ,log(data_final$new_cases[data_final$State=="Meghalaya"])
       ,col='red',pch=20)

## test positivity rate indicates high new cases
## in order to reduce the new cases; if we increase the 
## total number of cases then test positivity rate will drop
## As a result we expect number of new cases will drop as well


## fit the following model:
## log(new_cases) = a + b*log(test_postivity_rate) + error

## a) fit the model at state level
## b) create the following region of the states: (1) North, (2) South,
##    (3) West, (4) East and (5) North-East
## c) Fit the model at the regional levels
## d) What is the interpretation of high b?
## e) What is the interpretation of low b?
## f) Find the state where estimates b is maximum. 
## g) Find the state where estimates b is minimum.


## a) Model at state level
library(bayesm)
library(MCMCpack)

# creating data for bayesian hierarchical regression with one level (only states)
regdata<-NULL
nreg = length(state_nm_test_data)
i = 1
for (i in 1:nreg) { 
  filter <- data_final$State==state_nm_test_data[i] 
  
  y <- log(data_final[filter, "new_cases"] + 1) 
  X <- cbind(1, log(data_final[filter, "test_postivity_rate"])) 
  regdata[[i]] <- list(y=y, X=X) 
}

# fitting bayesian hierarchical regression model
Data <- list(regdata=regdata) 
sim.size<-21000
burn.in<-1000
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

# sampling beta values after burn.in
beta_draw <- out$betadraw[,1 , (burn.in+1):sim.size]
a <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 ,  (burn.in+1):sim.size]
b <- round(apply(beta_draw,1,mean),2)

beta_states <-cbind.data.frame(state_nm_test_data, a, b)

#displaying states and their a and b coefficient values for the model
head(beta_states)

###plotting a vs b coefficient values
plot(beta_states[, 2:3], pch = 20)

beta_states1 = beta_states
## create the following region of the states: (1) North, (2) South,
##    (3) West, (4) East and (5) North-East
state_nm_test_data

north = c("Chandigarh", "Delhi", "Haryana", "Himachal Pradesh", "Jammu and Kashmir", 
          "Ladakh", "Punjab", "Uttar Pradesh", "Uttarakhand")

south = c("Andaman and Nicobar Islands", "Andhra Pradesh", "Karnataka" , "Kerala",
          "Lakshadweep", "Puducherry", "Tamil Nadu", "Telangana"
          )

west = c("Goa", "Gujarat", "Madhya Pradesh", "Maharashtra", "Rajasthan")

east = c("Bihar", "Chhattisgarh" , "Jharkhand", "Odisha", "Tripura", "West Bengal")

northeast = c("Arunachal Pradesh", "Assam" ,
              "Manipur", "Meghalaya", "Mizoram", "Nagaland", "Sikkim")




### c) Fit the model at the regional levels


# creating data with 2 level hierarchy
states = c(north, south, west, east, northeast)
n0 = length(states)

regdata<-NULL
nreg = length(states)

i = 1
for (i in 1:nreg) { 
  filter <- data_final$State==states[i] 
  
  y <- log(data_final[filter, "new_cases"] + 1) 
  X <- cbind(1, log(data_final[filter, "test_postivity_rate"])) 
  
  
  regdata[[i]] <- list(y=y, X=X) 
}

# creating indicator matrix for which region state falls in
n1 = length(north)
n2 = length(south)
n3 = length(west)
n4 = length(east)
n5 = length(northeast)


Z = matrix(c(rep(1, n1), rep(0, (n0-n1)), rep(0, n1), rep(1, n2), rep(0, (n0-n1-n2)),
             rep(0,n1), rep(0,n2), rep(1,n3), rep(0, (n0-n1-n2-n3)), rep(0,n1), rep(0,n2),
             rep(0,n3), rep(1, n4), rep(0, (n0-n1-n2-n3-n4)), rep(0, (n0-n5)), rep(1, n5)), ncol = 5)
colnames(Z) = c("North", "South", "West", "East", "Northeast")      
nz = ncol(Z)

Z_ext <- cbind.data.frame(states,Z)

# fitting bayesian hierarchical regression model at 2 level hierarchy (regions and states)
sim.size<-2000
burn.in<-floor(sim.size*0.1)
Data <- list(regdata=regdata, Z=Z) 
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

#averaging over beta values after burn.in
beta_draw <- out$betadraw[,1 , (burn.in + 1):sim.size]
a <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 , (burn.in + 1):sim.size]
b <- round(apply(beta_draw,1,mean),2)


beta_states <-cbind.data.frame(states,a,b)
colnames(beta_states)<-c("States","a","b")
head(beta_states)


### Regional level analysis

Deltadraw <- out$Deltadraw
Deltadraw_star <- matrix(Deltadraw[200,],nrow = nz,byrow = F)
Delta_hat<-apply(Deltadraw[(burn.in + 1):sim.size,],2,mean)
Delta_hat<-matrix(Delta_hat,nrow = nz,byrow = F)
colnames(Delta_hat)<-c("a","b")
rownames(Delta_hat)<-colnames(Z)
Delta_hat

#We see from the regions, b value is maximum at West and is minimum in Northeast.

# the b value indicates the slope of the fitted regression line between 
#log (new cases) and log (test positivity rate)

### d) Interpretation of high b: 
## So, higher b value indicates that the rate of growth in new cases is higher.
## For two different states which have the same test positivity rate, 
## if the b values differ, the state with the higher b
## will have higher growth in new cases.

### e) Interpretation of low b:
## A low b value indicates that the rate of growth in new cases is low.
## For two different states which have the same test positivity rate, 
## if the b values differ, the state with the lower b
## will have lower growth in new cases.



### f) State where b estimate is maximum
state_bs = beta_states$b
max_ind = which.max(state_bs)

beta_states[max_ind, ]
# So, we see that in Maharashtra the b value is maximum

### g) State where b estimate is minimum
min_ind = which.min(state_bs)
beta_states[min_ind, ]

#So, we see that in Lakshadweep the b value is minimum

## Task 2: try to fit the following model:
## log(new_cases) = a + b*log(test_postivity_rate) + c*log(time)+d*log(time)^2 + error
## fit the model at different hierarchy
## a) Find the state where estimates b is maximum. 
## b) Find the state where estimates b is minimum. 

#creating the time variable
data_final$time = as.integer(as.factor(data_final$Updated.On))


#creating the data for bayesian hierarchical regression model
regdata<-NULL
nreg = length(states)

i = 1
for (i in 1:nreg) { 
  filter <- data_final$State==states[i] 
  
  y <- log(data_final[filter, "new_cases"] + 1) 
  X <- cbind(1, log(data_final[filter, "test_postivity_rate"]),
             log(data_final[filter, "time"]), (log(data_final[filter, "time"])^2)) 
  
  
  regdata[[i]] <- list(y=y, X=X) 
}

# fitting the hierarchical regression model
sim.size<-2000
burn.in<-floor(sim.size*0.1)
Data <- list(regdata=regdata, Z=Z) 
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

# averaging beta values after burn.in
beta_draw <- out$betadraw[,1 , (burn.in + 1):sim.size]
a <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 , (burn.in + 1):sim.size]
b <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,3 , (burn.in + 1):sim.size]
c <- round(apply(beta_draw,1,mean),2)

beta_draw <- out$betadraw[,4 , (burn.in + 1):sim.size]
d <- round(apply(beta_draw,1,mean),2)


beta_states <-cbind.data.frame(states,a, b, c, d)

#displaying some states along with the coefficient values
colnames(beta_states)<-c("States","a","b","c","d")
head(beta_states)

# regional level
Deltadraw <- out$Deltadraw
Deltadraw_star <- matrix(Deltadraw[200,],nrow = nz,byrow = F)
Delta_hat<-apply(Deltadraw[(burn.in + 1):sim.size,],2,mean)
Delta_hat<-matrix(Delta_hat,nrow = nz,byrow = F)
colnames(Delta_hat)<-c("a","b","c","d")
rownames(Delta_hat)<-colnames(Z)
Delta_hat
# Here, the b value is maximum in West and minimum in south

## a) state where b estimate is maximum
state_bs = beta_states$b
max_ind = which.max(state_bs)
beta_states[max_ind,]

#Maharashtra has the highest b value

## b) state where b estimate is minimum
min_ind = which.min(state_bs)
beta_states[min_ind,]

#Lakshadweep has the lowest b value

## Task 3: Consider the cheese data set
## fit the hierarchical regression model at retailer level
## log(Volume) = a + b*log(Price) + error
## 
## a) If you decide to provide a 10% discount on the expected price then identify 
## the top 5 retailers where the profit will be maximum
## b) If you decide to provide a 10% discount on the expected price then identify 
## the retailers where the you should not provide the discount at all.

# removing environment variables from previous problems
rm(list = ls())
#reading cheese dataset
data(cheese)
head(cheese)

#model log(volume) = a + b* log(Price) + error

#creating dataset based on retailer as level
retailer = levels(cheese$RETAILER)
nreg = length(retailer)

regdata<-NULL 
for (i in 1:nreg) { 
  filter <- cheese$RETAILER==retailer[i] 
  y <- log(cheese$VOLUME[filter]) 
  X <- cbind(1,      # intercept placeholder 
             log(cheese$PRICE[filter])) 
  
  regdata[[i]] <- list(y=y, X=X) 
}

# fitting bayesian hierarchical regression model
Data <- list(regdata=regdata) 
sim.size<-21000
burn.in<-1000
Mcmc <- list(R=sim.size)
set.seed(7831)
system.time(out <- bayesm::rhierLinearModel(
  Data=Data,
  Mcmc=Mcmc))

#averaging beta values after burn.in
beta_draw <- out$betadraw[,1 , (burn.in+1):sim.size]
int <- round(apply(beta_draw,1,mean),2)
beta_draw <- out$betadraw[,2 ,  (burn.in+1):sim.size]
b_price <- round(apply(beta_draw,1,mean),2)

beta_retailer <-cbind.data.frame(retailer,int,b_price)
head(beta_retailer)


#for each retailer checking profit ater 10% discount is given
retailer_no = 1
profits = rep(0, length(retailer))

for (retailer_no in 1: length(retailer))
{
  retailer_summary <- round(apply(cheese[cheese$RETAILER==retailer[retailer_no],c(2,4)],2,mean,na.rm=T),3)
  retailer_summary
  
  current_price<-retailer_summary[2]
  current_offer <- c(1,log(current_price))
  
  # applying 10% discount
  alternate_price<-0.9 * current_price
  alternate_offer <- c(1,log(alternate_price))
  
  #averaging after burn.in
  beta_draw <- out$betadraw[retailer_no, , (burn.in+1):sim.size]
  beta_mean <-apply(beta_draw,1,mean)
  
  #calculating change in volume
  vol_change<-rep(NA,ncol(beta_draw))
  for(i in 1:ncol(beta_draw)){
    volume1 <- exp(beta_draw[,i]%*%current_offer)
    volume2 <- exp(beta_draw[,i]%*%alternate_offer)
    vol_change[i]<-((volume2-volume1)/volume1)*100
  }
  sum_change <-summary(vol_change)
  sum_change
  
  #calculating expected revenue without discount
  expected_revenue = retailer_summary["VOLUME"]*retailer_summary["PRICE"]
  names(expected_revenue)<-"current expected revenue"
  
  #calculating changed revenue with discount and changed volume
  changed_exp_revenue<-rep(NA,length(vol_change))
  for(i in 1:length(vol_change))
    changed_exp_revenue[i] = retailer_summary["VOLUME"]*(1+vol_change[i]/100)*alternate_price

  new_expected_revenue = mean(changed_exp_revenue)
  names(new_expected_revenue)<-"new expected revenue"
  
  #calculating profit as a result of discount
  profit = new_expected_revenue - expected_revenue
  names(profit) = "profit"
  profits[retailer_no] = profit
}

retailer_profits<- data.frame("Retailer" = retailer, 
           "Profits" = profits)

#sorting retailer profits
sorted_idx = sort(retailer_profits$Profits, decreasing = T, index.return = T)$ix
retailer_profits = retailer_profits[sorted_idx,]

#displaying top 5 retailers which will profit the most with giving 10% discount
head(retailer_profits, 5)

# the top 5 retailers and their respective profits:
#                         Retailer   Profits
#                  CHICAGO - JEWEL 13734.578
# BUFFALO/ROCHESTER - TOPS MARKETS 10077.553
#               CHICAGO - DOMINICK  9840.472
#        NEW YORK (NEW) - PATHMARK  9337.370
#      BUFFALO/ROCHESTER - WEGMANS  8765.090

#b) Retailers which should not provide discount at all

#displaying retailer which will incur a loss after giving discount
retailer_profits[retailer_profits$Profits < 0,]

# retailers which should not provide discount:
#                   Retailer   Profits
#JACKSONVILLE,FL - FOOD LION -131.2593
#           CHARLOTTE - BI LO -140.6583