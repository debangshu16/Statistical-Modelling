This folder contains the R codes for the time series analaysis performed on the UKNon Durables dataset.
The dataset contains quarterly time series data of the Non Durables Consumption in UK from 1955 to 1989.  
As we proceed with the course, we try out various methodologies to get a good model to forecast the non durables consumption in UK.

The first folder "Non_Seasonal_Basic" ignores the seasonality in the data and fits an ARIMA model with the lowest AIC/BIC score. Obviously, not considering seasonality at all is not a good thing to do but this is the ground zero from where we start. We also perform model diagnostics to show that the model selected is not a good one.
