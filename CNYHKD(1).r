#install.packages("TSA")
library(ggplot2)
library(tseries)

setwd("E:/ио©н/STAT4601/GP/GP")
# Set the file path of the CSV file
file_path <- "E:/ио©н/STAT4601/GP/CNYHKD=X.csv"

# Read the CSV file
data <- read.csv(file_path)

# Select the first and fifth column and exclude the last 5 rows
CNYHKD <- data[1:(nrow(data)-5), c(1, 5)]

# Print the data
print(data)

# Convert the Date column to a date format
CNYHKD$Date <- as.Date(CNYHKD$Date)

##############################################################################
# PART I. Stationary check, transformation, and model identification

# Draw a time plot using ggplot2
ggplot(CNYHKD, aes(x = Date, y = Close)) +
  geom_line() +
  labs(x = "Date", y = "Close", title = "CNY/HKD Exchange Rate Time Plot")

# Check ACF and PACF plots
acf(CNYHKD$Close)
pacf(CNYHKD$Close)

# Perform the Dickey-Fuller test
result <- tseries::adf.test(CNYHKD$Close)

# Print the test results
print(result)

diff_CNYHKD <- diff(CNYHKD$Close)

ggplot(data.frame(Date = CNYHKD$Date[-1], Close = diff_CNYHKD), aes(x = Date, y = Close)) +
  geom_line() +
  labs(x = "Date", y = "Close", title = "Differenced CNY/HKD Exchange Rate Time Plot")

log_CNYHKD <- log(CNYHKD$Close[-1])

ggplot(data.frame(Date = CNYHKD$Date[-1], Close = log_CNYHKD), aes(x = Date, y = Close)) +
  geom_line() +
  labs(x = "Date", y = "Close", title = "Log-transformed CNY/HKD Exchange Rate Time Plot")

diff_log_CNYHKD <- diff(log_CNYHKD)
ggplot(data.frame(Date = CNYHKD$Date[-c(1,2)], Close = diff_log_CNYHKD), aes(x = Date, y = Close)) +
  geom_line() +
  labs(x = "Date", y = "Close", title = "Diff_Log-transformed CNY/HKD Exchange Rate Time Plot")


acf(diff_CNYHKD)
#pacf(diff_CNYHKD)

acf(log_CNYHKD)
#pacf(log_CNYHKD)

# Perform the Dickey-Fuller test
result <- tseries::adf.test(log_CNYHKD)
print(result)

acf(diff_log_CNYHKD)
pacf(diff_log_CNYHKD)

# Perform the Dickey-Fuller test
result <- tseries::adf.test(diff_log_CNYHKD)
print(result)

par(mfrow=c(1,2),mai=c(0.8,0.8,0.8,0.1))
acf(diff_log_CNYHKD); pacf(diff_log_CNYHKD) 

#par(mfrow=c(1,3),mai=c(0.8,0.8,0.8,0.1))

###############################################################################################################

#  Part III Parameter estimation;
# Fit MA model with q = 2
ma_model <- arima(diff_log_CNYHKD, order = c(0, 0, 2), method='CSS-ML')
ma_model
summary(ma_model)
# Fit AR model with p = 2
ar_model <- arima(diff_log_CNYHKD, order = c(2, 0, 0), method='CSS-ML')
ar_model
summary(ar_model)

#  Part IV Model diagnostics;
# IV.1 Time plot of residuals
par(mfrow=c(1,1),mai=c(0.8,0.8,0.8,0.1))
plot(ma_model$residuals, ylab='Residuals'); abline(h=0)
plot(ar_model$residuals, ylab='Residuals'); abline(h=0)

# IV.2 Normality checking
qqnorm(ma_model$residuals); qqline(ma_model$residuals);
qqnorm(ar_model$residuals); qqline(ar_model$residuals); 
hist(ma_model$residuals)
hist(ar_model$residuals)

shapiro.test(ma_model$residuals)
shapiro.test(ar_model$residuals)

# IV.3 Autocorrelation checking to check residual autocorrelations individually
par(mfrow=c(1,2),mai=c(0.8,0.8,0.8,0.1))
acf(ma_model$residuals); pacf(ma_model$residuals)
acf(ar_model$residuals); pacf(ar_model$residuals)

par(mfrow=c(1,1),mai=c(0.8,0.8,0.8,0.1))

# IV.4 Box-Pierce or Ljung-Box test to check residual autocorrelations jointly
Box.test(ma_model$residuals, lag = 10, type = c("Ljung-Box"), fitdf = 1)
Box.test(ar_model$residuals, lag = 10, type = c("Ljung-Box"), fitdf = 1)

# overparameter
model_3<-arima(diff_log_CNYHKD,order = c(0,0,3), method='CSS-ML')
model_12<-arima(diff_log_CNYHKD, order = c(1,0,2), method='CSS-ML')
model_4<-arima(diff_log_CNYHKD, order = c(0,0,4), method='CSS-ML')


AIC(ma_model)
AIC(model_3)
AIC(model_12)
AIC(model_4)

BIC(ma_model)
BIC(model_3)
BIC(model_12)
BIC(model_4)




