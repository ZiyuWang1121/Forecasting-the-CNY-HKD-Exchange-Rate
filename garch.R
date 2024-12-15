library(TSA)
library(ggplot2)
library(lmtest)
library(tidyverse)
library(FinTS)
library(texreg)
library(rugarch)
library(forecast)
library(fGarch)

setwd("C:\\Users\\87840\\Desktop\\4601 gp")
data <- read.csv("CNYHKD=X.csv")

# Select the first and fifth column and exclude the last 5 rows
CNYHKD <- data[1:(nrow(data)-5), c(1, 5)]

CNYHKD$Date <- as.Date(CNYHKD$Date)

log_CNYHKD <- log(CNYHKD$Close[-1])
diff_log_CNYHKD <- diff(log_CNYHKD)

# Fit AR(2) model
ar2_model1 <- arima(diff_log_CNYHKD, order = c(2, 0, 0), method = "ML")
ar2_model1

ar2_model2 <- arima(diff_log_CNYHKD, order = c(2, 0, 0), method = "CSS")
ar2_model2
arima(diff_log_CNYHKD, order = c(2, 0, 0), method = "CSS-ML")

ma2_model <- arima(diff_log_CNYHKD,order = c(0,0,2))

plot(ar2_model1$residuals, ylab='Residuals',main="AR(2) Residual Plot"); abline(h=0)
plot(ma2_model$residuals, ylab='Residuals',main="MA(2) Residual Plot"); abline(h=0)

resid <- ar2_model1$residuals

ArchTest(ar2_model1$residuals)
ArchTest(ma2_model$residuals)

acf(residuals(ar2_model1)^2)
pacf(residuals(ar2_model1)^2)


#fit the garch model
gfit <- garchFit(~garch(1,1), data=diff_log_CNYHKD, trace=FALSE, cond.dist = c("norm"),
                 include.mean = FALSE)
summary(gfit)

# Predict the next 5 values
predict(gfit, n.ahead = 5, plot=TRUE, conf=.95, nx=100)

#Estimate Garch model

# OPTIONS
ar_lag <- 2 # lag used for ar term in mean equation (0 in paper)
ma_lag <- 0 # lag used for ma term in mean equation (0 in paper)
arch_lag <- 1 # lag in arch effect (1 in paper)
garch_lag <- 1 # lag in garch effect (1 in paper)
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch manual for more
distribution_to_estimate <- c('norm','sstd') # distribution used in all models
my_html_file <- 'tabs/tab04-estimation_garch.html'

# END OPTIONS

# get all combinations of models
df_grid <- expand_grid(ar_lag,
                       ma_lag,
                       arch_lag,
                       garch_lag,
                       models_to_estimate,
                       distribution_to_estimate)

#define function for estimation
estimate_garch <- function(ar_lag,
                           ma_lag,
                           arch_lag,
                           garch_lag,
                           models_to_estimate,
                           distribution_to_estimate) {
  
  message('Estimating ARMA(',ar_lag,',', ma_lag, ')', '-',
          models_to_estimate, '(', arch_lag, ',', garch_lag, ') ', 
          'dist = ', distribution_to_estimate)
  
  # estimate model
  my_spec <- ugarchspec(variance.model = list(model = models_to_estimate,
                                              garchOrder = c(arch_lag, 
                                                             garch_lag)),
                        mean.model = list(armaOrder = c(ar_lag,
                                                        ma_lag)), 
                        distribution.model = distribution_to_estimate)
  
  my_garch <- ugarchfit(spec = my_spec, data = diff_log_CNYHKD)
  
  return(my_garch)
}

# estimate all models
l_args <- as.list(df_grid)
l_models <- pmap(.l = l_args, .f = estimate_garch)
l_models

find_best_arch_model <- function(x, type_models, dist_to_use, 
                                 max_lag_AR, max_lag_MA, max_lag_ARCH, max_lag_GARCH) {
  
  require(rugarch)
  
  # Initialize empty data frame to store results
  df_results <- data.frame()
  
  # Loop over each model type and distribution
  for(model in type_models) {
    for(dist in dist_to_use) {
      
      # Loop over each possible order combination
      for(ar_lag in 0:max_lag_AR) {
        for(ma_lag in 0:max_lag_MA) {
          for(arch_lag in 0:max_lag_ARCH) {
            for(garch_lag in 0:max_lag_GARCH) {
              
              # Define specification of the GARCH model
              spec <- ugarchspec(variance.model = list(model = model,
                                                       garchOrder = c(arch_lag, garch_lag)), 
                                 mean.model = list(armaOrder = c(ar_lag, ma_lag), 
                                                   include.mean = TRUE), 
                                 distribution.model = dist)
              
              # Try to fit the GARCH model and skip to next iteration if there is an error
              tryCatch({
                fit <- ugarchfit(spec, data = x)
                
                # Create a data frame with the results
                df_temp <- data.frame(
                  model_name = paste(model, dist, sep = "_"),
                  type_model = model,
                  type_dist = dist,
                  AR_lag = ar_lag,
                  MA_lag = ma_lag,
                  ARCH_lag = arch_lag,
                  GARCH_lag = garch_lag,
                  AIC = infocriteria(fit)[1],
                  BIC = infocriteria(fit)[2],
                  stringsAsFactors = FALSE
                )
                
                # Bind the temporary data frame to the main data frame
                df_results <- rbind(df_results, df_temp)
              }, error = function(e) {
                message(paste("Error fitting model with parameters: ",
                              "model =", model,
                              "dist =", dist,
                              "AR_lag =", ar_lag,
                              "MA_lag =", ma_lag,
                              "ARCH_lag =", arch_lag,
                              "GARCH_lag =", garch_lag))
              })
            }
          }
        }
      }
    }
  }
  
  return(df_results)
}
max_lag_AR <- 1 # used 1 in paper
max_lag_MA <- 1 # used 1 in paper
max_lag_ARCH <- 2 # used 2 in paper
max_lag_GARCH <- 1 # used 1 in paper
dist_to_use <- c('norm', 'std') # see rugarch::ugarchspec help for more
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch::rugarchspec help for more
out <- find_best_arch_model(x = diff_log_CNYHKD, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)



# 1. Reshape data from wide to long
df_long <- out %>%
  tidyr::pivot_longer(
    cols = c("AIC", "BIC"),
    names_to = "criteria",
    values_to = "value"
  )

# 2. Identify the best model (with the lowest AIC and BIC)
best_model_AIC <- df_long %>%
  filter(criteria == "AIC") %>%
  slice(which.min(value))

best_model_BIC <- df_long %>%
  filter(criteria == "BIC") %>%
  slice(which.min(value))

best_models <- rbind(best_model_AIC, best_model_BIC)

# 3. Create the plot
p1 <- ggplot(df_long, 
             aes(x = reorder(model_name, value), 
                 y = value, 
                 shape = type_dist, 
                 color = type_model)) + 
  geom_point(size = 3.5, alpha = 0.65) + 
  coord_flip() + 
  theme_bw(base_family = "TT Times New Roman") + 
  facet_wrap(~criteria, scales = 'free_x') + 
  geom_point(data = best_models, 
             aes(x = reorder(model_name, value), 
                 y = value), 
             color = 'blue', size = 5, shape = 8) +
  labs(title = 'Selecting GARCH Models by Fitness Criteria', 
       subtitle = 'The best model is the one with lowest AIC or BIC (with star)',
       x = '',
       y = 'Value of Fitness Criteria',
       shape = 'Type of Dist.',
       color = 'Type of Model') + 
  theme(legend.position = "right")

p1


spec = ugarchspec(variance.model = list(model = "gjrGARCH", 
                                        garchOrder = c(1, 0)), 
                  mean.model = list(armaOrder = c(1, 0), 
                                    include.mean = TRUE), 
                  distribution.model = "std")

# Fit the model to your data
fit = ugarchfit(spec, data = diff_log_CNYHKD)
fit

garch_resid <- residuals(fit)

plot(garch_resid, type='l', ylab='Log Return', xlab='Time')

hist(garch_resid,breaks = 10)

qqnorm(garch_resid)
qqline(garch_resid)

# ACF plot of residuals
acf(garch_resid)

# PACF plot of residuals
pacf(garch_resid)

Box.test(residuals(fit)^2, type = "Ljung-Box",lag = 10)
Box.test(residuals(fit)^2, type = "Ljung-Box",lag = 15)
Box.test(residuals(fit)^2, type = "Ljung-Box",lag = 20)

shapiro.test(garch_resid)

volatility = sigma(fit)
plot(volatility, main = "Conditional Volatility")

# Make the forecast
forecast = ugarchforecast(fit, n.ahead = 5)

# Print the forecast
print(forecast$mean)
# Print the forecasted mean values
print(forecast@forecast$seriesFor)

# Print the forecasted volatility
print(forecast@forecast$sigmaFor)

plot(forecast)
