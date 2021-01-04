library(tseries)
library(zoo)
library(forecast)
library(rugarch)

start_date <- as.Date("2019-11-01")
end_date <- as.Date("2020-10-31")
google_stock<-get.hist.quote(instrument = "GOOGL"
                                      ,start=start_date,end=end_date
                                      ,quote="Close",provider = "yahoo") 


plot(google_stock, xlab = "Time", ylab ="Close", main = "Google Stock Close vs Time
     ")


data = ts(google_stock)
#Stationarity test of time series
adf.test(data)
kpss.test(data)
pp.test(data)

#Hence the stock price series is not stationary

lrt = diff(log(data))
plot(lrt, xlab = "Time", ylab = "Log Return",
     main = "Google stock Log return Time Series")

adf.test(lrt)
kpss.test(lrt)
pp.test(lrt)

# The log return series is stationary

period = 6
n<- length(lrt)
train_lim = n - period
train_data = lrt[c(1:train_lim)]
test_data = lrt[c((train_lim+1):n)]

length(train_data)
length(test_data)

d<-ndiffs(train_data)

aic_model = auto.arima(train_data,max.p = 5, max.q = 5,d = d, ic = "aic")
aic_model

bic_model = auto.arima(train_data,max.p = 5, max.q = 5,d =d, ic = "bic")
bic_model

pred_arima<- forecast(aic_model, h = period)
arima_forecast<- pred_arima$mean


plot(forecast(aic_model,h=period),type='l')
lines(c((train_lim+1):n),test_data, col = 'black')

residuals = aic_model$residuals
residuals2 = residuals ^ 2


Box.test(residuals, lag = 5, type = c("Ljung-Box"))
Box.test(residuals2,lag = 5, type = c("Ljung-Box"))

acf(residuals)

acf(residuals2)
pacf(residuals2)

spec = ugarchspec(variance.model = list(model = "sGARCH", 
                                        garchOrder = c(1,1)),
                  mean.model     = list(armaOrder = c(1, 0)),
                  distribution.model = 'norm' ) 
                                        
garch_lrt <- ugarchfit(spec = spec, data = train_data)
garch_lrt
pred_garch = ugarchforecast(garch_lrt, data = train_data, n.ahead = period)

forecast_garch<- pred_garch@forecast$seriesFor

sigma_forecast<- pred_garch@forecast$sigmaFor


plot(NULL,xlim = c(1,n),ylim = c(-0.2,0.2),ylab = 'Market Log Return',xlab='Time',
     main = 'Actual and Forecasted values of market log return')
points(c(1:train_lim), train_data,type='l',lwd=2,col='grey')
points(c(1:n), lrt, type = 'l', col = 'grey')
points(c((train_lim+1):n), test_data,type='l',lwd=2,col='blue')
points(c((train_lim+1):n), arima_forecast, type='l', lwd = 2, col = 'red')
points(c((train_lim+1):n), forecast_garch,type='l',lwd=2,col='black')
abline(v = (train_lim+1), col = 'grey', lwd = 2)
legend("bottomright",
       legend =  c("train_values","test_values","GARCH Forecast", "ARIMA Forecast"),
       col = c("grey","blue","black", "red"),lty = 1,
       cex = 0.8)

mse_garch = sum((test_data - forecast_garch)^2)
mse_garch

mse_arima = sum((test_data - arima_forecast)^2)
mse_arima


