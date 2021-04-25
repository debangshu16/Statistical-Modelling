state_wise=read.csv(file="https://api.covid19india.org/csv/latest/state_wise_daily.csv",header=TRUE,stringsAsFactors = F)

test_data = read.csv(file = 'https://api.covid19india.org/csv/latest/statewise_tested_numbers_data.csv')


start_date = as.Date("01/05/2020", format = "%d/%m/%Y")
end_date = as.Date("24/04/2021", format = "%d/%m/%Y")

tm = seq.Date(from = start_date, to = end_date, by = "day")
n = length(tm)

get_daily_tested<-function(test_data, state = "West Bengal", start_date = as.Date("01/05/2020", format = "%d/%m/%Y"),
              end_date = as.Date("24/04/2021", format = "%d/%m/%Y"))
{
  tm = seq.Date(from = start_date, to = end_date, by = "day")
  n = length(tm)
  
  state_test_data = test_data[(test_data$State==state),]
  tested_before = state_test_data[(as.Date(state_test_data$Updated.On, format = "%d/%m/%Y") == (start_date-1)),]$Total.Tested
  state_test_data = state_test_data[(as.Date(state_test_data$Updated.On, format = "%d/%m/%Y")>=start_date),]
  state_test_data = state_test_data[(1:n),]
  
  daily_tested = rep(0,n)

  total_tested = state_test_data$Total.Tested
  for (i in 1:n)
  {
    if (i==1)
      daily_tested[i] = total_tested[i] - tested_before
    else
      daily_tested[i] = total_tested[i] - total_tested[i-1]
  }
  
  return (daily_tested)
}

get_daily_confirmed<- function(state_wise, state = "WB", start_date = as.Date("01/05/2020", format = "%d/%m/%Y"),
                               end_date = as.Date("24/04/2021", format = "%d/%m/%Y"))
{
  state_data = state_wise[(state_wise$Date_YMD>=start_date),c("Date_YMD","Status",state)]
  state_data_confirmed = state_data[(state_data$Status=="Confirmed"),state]
  return (state_data_confirmed)
}

confirmed_per_tested<- function(confirmed, tested)
{
  n = length(confirmed)
  res = rep(0, n)
  
  for (i in 1:n)
  {
    res[i] = confirmed[i]/tested[i]
  }
  return (res)
}



# test data state names: "West Bengal", "Delhi", "Maharashtra", "Tamil Nadu"
# state data state names: "WB", "DL", "MH", "TN"


## West bengal
daily_tested_wb = get_daily_tested(test_data, "West Bengal")
plot(tm, daily_tested_wb, type = 'l')

daily_confirmed_wb = get_daily_confirmed(state_wise, state = "WB")
plot(tm, daily_confirmed_wb, type = 'l')


tpr_wb = confirmed_per_tested(daily_confirmed_wb, daily_tested_wb) * 100

plot(tm, tpr_wb, type = 'l')

## Delhi
daily_tested_dl = get_daily_tested(test_data, "Delhi")
plot(tm, daily_tested_dl, type = 'l')

daily_confirmed_dl = get_daily_confirmed(state_wise, state = "DL")
plot(tm, daily_confirmed_dl, type = 'l')


tpr_dl = confirmed_per_tested(daily_confirmed_dl, daily_tested_dl) * 100
plot(tm, tpr_dl, type = 'l')

## Maharashtra
daily_tested_mh = get_daily_tested(test_data, "Maharashtra")
plot(tm, daily_tested_mh, type = 'l')

daily_confirmed_mh = get_daily_confirmed(state_wise, state = "MH")
plot(tm, daily_confirmed_mh, type = 'l')


tpr_mh = confirmed_per_tested(daily_confirmed_mh, daily_tested_mh) * 100
plot(tm, tpr_mh, type = 'l')


## Tamil Nadu
daily_tested_tn = get_daily_tested(test_data, "Tamil Nadu")
plot(tm, daily_tested_tn, type = 'l')

daily_confirmed_tn = get_daily_confirmed(state_wise, state = "TN")
plot(tm, daily_confirmed_tn, type = 'l')


tpr_tn = confirmed_per_tested(daily_confirmed_tn, daily_tested_tn) * 100
plot(tm, tpr_tn, type = 'l')


plot(x = tm, y = tpr_wb, type = 'l', col = 'red',   ylim = c(0,50),
     xlab = "Time", ylab = "True Positive Rate", main = "Tpr for various states")
lines(x = tm, y = tpr_mh, type = 'l', col = "green")
lines(x = tm, y = tpr_dl, type = 'l', col = "blue")     
lines(x = tm, y = tpr_tn, type = 'l', col = "brown")
legend('topleft',
       legend =  c("West Bengal", "Maharashtra", "Delhi", "Tamil Nadu"),
       col = c("red", "green" , "blue", "brown"),
       lty = 1,
       cex = 0.6)


### Fit ARIMA model on true positive rate data for West Bengal
library(tseries)
library(forecast)

n
period = 7
train_tpr_wb = tpr_wb[1:(n-period)]
test_tpr_wb = tpr_wb[(n-period+1):n]

kpss.test(train_tpr_wb)
# Data is not stationary

kpss.test(diff(train_tpr_wb))
# 1st difference Data stationary

d = ndiffs(train_tpr_wb)
d

arima_model = auto.arima(train_tpr_wb,start.p =0, start.q = 0,seasonal = F,
                         max.p = 10, max.q = 10, ic = "aicc")
arima_model

forecasted_values = forecast(arima_model, h= period)
forecasted_values = as.numeric(forecasted_values$mean)

mse = sum((test_tpr_wb - forecasted_values)^2)
mse

prediction = data.frame("Actual" = test_tpr_wb, "Predicted" = forecasted_values)
rownames(prediction) = c(as.Date(tm[(n-period+1):n], format = "%Y-%m-%d"))

prediction

# plotting Arima prediction
plot(x = tm, y = tpr_wb, col = 'blue', type = 'l')
lines(x = c(as.Date(tm[(n-period+1):n], format = "%Y-%m-%d")),
      y = forecasted_values, type = 'l', col = 'red')

legend('topleft',
legend =  c("Data", "Forecasted last 7 days"),
col = c("blue", "red" ),
lty = 1,
cex = 0.6)

# Arima captures the trend but cannot capture the huge spike

# Model diagnostics
acf(arima_model$residuals)

#looking at residuals the lag at 3 is significant. Hence model is not capturing perfectly
Box.test(arima_model$residuals,lag = 10, type = "Box-Pierce")
# null hypothesis of no autocorrelation of residuals is rejected. So, model has lack of fit

qqnorm(arima_model$residuals)
shapiro.test(arima_model$residuals)

# all model diagnostics of ARIMA model fails. So it is not the correct model in capturing this pattern


######## GARCH model