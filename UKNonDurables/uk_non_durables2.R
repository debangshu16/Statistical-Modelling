library(tseries)
library(forecast)
library(dplyr)
library(moments)

data = read.table('UKNonDurables.csv',sep=',',header=T,stringsAsFactors = F)
head(data)

#Only keeping two columns : Time and Value
data = data[, c('time','value')]
head(data)

original_values = data$value
data$value = log(data$value)

#Ordering data based on time
data = data[order(data$time),]

#Check for missing data 
paste("Missing values in dataset is",sum(is.na(data)))
#So there are no missing values in the dataset

#Time series plot of the whole data
start_date = data$time[1]
value = ts(data$value, start = start_date,freq =4)

ts.plot(value,main = "log of Consumption of Non Durables in UK")

#Cleaning data
value = tsclean(value)
data$value = value

tm = data$time

#Split data into train data and test data with train data = 95% of the data
period<- 8
train_lim = length(data$value) - period

train_data = data[1:train_lim,]
test_data = data[(train_lim+1):length(data$value),]

n = length(data$value)
m = length(train_data$value)

value = ts(train_data$value,start = start_date, freq = 4)
#Descriptive statistics of the data
descriptive_stats = list('Mean' = mean(value),'Median' = median(value),'Minimum' = min(value),
                         'Maximum' = max(value), 'Standard Deviation' = sd(value) )
descriptive_stats$skewness = skewness(value)
descriptive_stats$kurtosis = kurtosis(value)

f <-function(x,mean1,sd1)
{
  #print (paste (mean1,sd1))
  count = 0
  for (i in c(1:length(x)))
  {
    if ((x[i]>=(mean1-3*sd1)) && (x[i]<=(mean1+3*sd1)))
      count = count+1
  }
  
  return ((count*100)/(length(x)))
}
descriptive_stats$perct_within_3sd = f(value,descriptive_stats$Mean,descriptive_stats$`Standard Deviation`)

str(descriptive_stats)

decomposed = decompose(value,type = "multiplicative")
plot(decomposed)

seasonal = decomposed$seasonal
seasonal.components = seasonal[1:4]

value = value - seasonal
ts.plot(value)

adf = adf.test(value)
print (adf)
kpss = kpss.test(value)
print (kpss)
pp = pp.test(value)
print (pp)


ts.plot(diff(value),main = "Difference of log of Consumption of Non Durables in UK")

adf = adf.test(diff(value))
print (adf)
kpss = kpss.test(diff(value))
print (kpss)
pp = pp.test(value)
print (pp)

d = ndiffs(value)

d

aic_model = auto.arima(value,start.p =0, start.q = 0,seasonal = F,max.p = 5, max.q = 5,d =d, ic = "aic")
aic_model

bic_model = auto.arima(value,start.p =0, start.q = 0,max.p = 5,seasonal = F, max.q = 5,d =d, ic = "bic")
bic_model


#plot(forecast(bic_model,h=period),type='l')
#lines(c(130:136),test_data$value)

##forecast_values = forecast(bic_model,h = period)
#forecast_values = data.frame(forecast_values)

forecasted_values = forecast(aic_model, h= period)
forecasted_values = as.numeric(forecasted_values$mean)

starting_test_index = data$time[(m+1)]
starting_quarter = ((starting_test_index - floor(starting_test_index))*100/25)+1

quarter = starting_quarter
for (i in c(1:period))
{
  forecasted_values[i] = forecasted_values[i] + seasonal.components[quarter]
  quarter = ((quarter + 1)%%4) +1
}

forecasted_values = lapply(forecasted_values,exp)
forecasted_values = as.numeric(forecasted_values)

true_values = original_values[(m+1):(m+period)]
mse = mean((true_values-forecasted_values)^2)

sqrt(mse)



plot(NULL,xlim = c(tm[1],tm[n]),ylim = c(20000,70000),ylab = 'Non Durables consumption',xlab='')
points(train_data$time, exp(train_data$value),type='l',lwd=2,col='purple')
points(test_data$time, exp(test_data$value),type='l',lwd=2,col='blue')
points(test_data$time, forecasted_values,type='l',lwd=2,col='black')
legend("bottomright",
      legend =  c("train_values","test_values","forecasted_values"),
      col = c("purple","blue","black"),lty = 1,
      cex = 0.8)




