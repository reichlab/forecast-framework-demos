library(R6)
library(ForecastFramework)
library(dplyr)
library(forecast)
library(ggplot2)
library(gridExtra)
setwd("C://Users/house/Documents/GitHub/Reich_Lab/ForecastFramework Project/forecast-framework-demos/models")

source('GamModel.R')
source('SARIMATD1Model.R')
source('ContestModel.R')

# Upload data
dat <- read.csv('../data/province-biweek-counts.csv')
dat$date_sick <- as.Date(strptime(dat$date_sick,"%m/%d/%Y")) # convert to date values
training_data <- dat %>% filter(year != 2013)
testing_data <- dat %>% filter(year == 2013)

# Convert Training and Testing data to incidence matrices
# time.in.year is a requirement for the GAM model
preprocess_inc <- function(dat){
  dat$time.in.year = dat$biweek
  dat$t = dat$year + dat$biweek/26
  inc = ObservationList$new(dat)
  inc$formArray('province','t',val='cases',
                dimData = list(NULL,list('biweek','year','time.in.year','t')),
                metaData = list(t.step = 1/26,max.year.time = 26))
  return(inc)
}

# Preprocess Data and convert to Inc Mat
training_inc <- preprocess_inc(training_data)
testing_inc <- preprocess_inc(testing_data)

# Defining the Models
nsim <- 100 # Number of simulations 
prov_nums <- 1 # Number of Provinces
sarima_model <- SARIMATD1Model$new(period = 26, nsim = nsim)
gam_model <- GamModel$new(numProvinces = prov_nums,nSims = nsim)

# Fit the Models 
sarima_model$fit(training_inc)
gam_model$fit(training_inc)

# Forecasting the Models
steps <- 26 # forecast ahead 26 biweeks
forecast_gam <- gam_model$forecast(step = steps)
forecast_sarimaTD <- sarima_model$forecast(step = steps)

# Function to create time series with Visualization
visualize_model <- function(forecast, rib_color, main_color){
  data_X_3years <- dat %>% filter(year > 2011)
  preds_df <- data.frame(as.table(t(forecast$data$mat)))
  preds_df[["date_sick"]] <-as.Date(testing_data$date_sick, format = "%d/%m/%Y")
  preds_df <- preds_df %>%
    mutate( 
      pred_total_cases = as.vector(forecast$median()$mat),
      pred_95_lb = as.vector(forecast$quantile(0.025)$mat),
      pred_95_ub = as.vector(forecast$quantile(0.975)$mat),
      pred_80_lb = as.vector(forecast$quantile(0.05)$mat),
      pred_80_ub = as.vector(forecast$quantile(0.95)$mat),
      pred_50_lb = as.vector(forecast$quantile(0.25)$mat),
      pred_50_ub = as.vector(forecast$quantile(0.75)$mat)
    )
  plot <- ggplot() +
    geom_ribbon(
      mapping = aes(x = date_sick, ymin = pred_95_lb, ymax = pred_95_ub),
      fill = rib_color,
      alpha = 0.2,
      data = preds_df) +
    geom_ribbon(
      mapping = aes(x = date_sick, ymin = pred_80_lb, ymax = pred_80_ub),
      fill = rib_color,
      alpha = 0.2,
      data = preds_df) +
    geom_ribbon(
      mapping = aes(x = date_sick, ymin = pred_50_lb, ymax = pred_50_ub),
      fill = rib_color,
      alpha = 0.2,
      data = preds_df) +
    geom_line(
      mapping = aes(x = date_sick, y = pred_total_cases),
      color = main_color,
      size = 1,
      data = preds_df) + 
    geom_line(mapping = aes(x = date_sick, y = cases),
              size=0.7,
              data = data_X_3years) +
    xlab("") + ylab("Number of Cases") +
    coord_cartesian(ylim = c(0, 1000))
  return(plot)
}

# Create each plot
plot1 <- visualize_model(forecast_sarimaTD, "coral1", "red") + ggtitle("SARIMATD")
plot2 <- visualize_model(forecast_gam, "cornflowerblue" , "blue")+ ggtitle("GAM")

# Display side-by-side plots
grid.arrange(plot1, plot2, ncol=2)

# Easily Calculate Errors with ForecastFramework
MAE <- function(forecast_inc, test_inc){
  err <- test_inc$mat-forecast_inc$mean()$mat
  abs_err <- abs(err)
  return(mean(abs_err))
}

MSE <- function(forecast_inc, test_inc){
  err <- test_inc$mat-forecast_inc$mean()$mat
  err_sq <- err^2
  return(mean(err_sq))
}

MAPE <- function(forecast_inc, test_inc){
  err <- test_inc$mat-forecast_inc$mean()$mat
  abs_err <- abs(err)
  abs_err_percent <- abs_err/test_inc$mat
  return(mean(abs_err_percent)*100)
}

RMSE <- function(forecast_inc, test_inc){
  err <- test_inc$mat-forecast_inc$mean()$mat
  err_sq <- err^2
  mse <-mean(err_sq)
  return(sqrt(mse))
}

# Create Arrays with each Metric
evaluate <- function(forecast){
  MAE <- MAE(forecast, testing_inc)
  MSE <- MSE(forecast, testing_inc)
  MAPE <- MAPE(forecast, testing_inc)
  RMSE <- RMSE(forecast, testing_inc)
  return(c(MAE,MSE,MAPE,RMSE))
}

gam_eval <- evaluate(forecast_gam)
sarimaTD_eval <- evaluate(forecast_sarimaTD)


df <- data.frame("MAE" = c(gam_eval[1],sarimaTD_eval[1]),
                 "MSE" = c(gam_eval[2],sarimaTD_eval[2]),
                 "MAPE" = c(gam_eval[3],sarimaTD_eval[3]),
                 "RMSE" = c(gam_eval[4],sarimaTD_eval[4]),
                 stringsAsFactors = FALSE)
row.names(df) <- c('GAM','sarimaTD')


