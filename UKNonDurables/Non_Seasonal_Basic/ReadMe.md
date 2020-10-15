This folder contains building a simple ARIMA model to forecast Non Durables Consumption in UK.
We blindly assume that the data is not seasonal (unreasonable assumption), but we just want to quick build a simple model as a starter and then proceed from there.
We split the data into train and test sets with test sets covering the last 6 quarters. We fit the ARIMA model on the train set and compare the forecasted values with the original test set values.
We also perform model diagnostics to see how good/bad the model fitted actually is.
