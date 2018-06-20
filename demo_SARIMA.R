# Import all Libraries
require(ForecastFramework)
require(R6)
require(forecast)
require(dplyr)
require(ggplot2)

# Source R6 Files
source('models/ContestModel.R')
source('models/SARIMAModel.R')

###########################
# ~~~~~~ Model Fit ~~~~~~
###########################

nsim <- 1000 # Number of SARIMA simulations 
sarima_model <- SARIMAModel$new(period = 26, nsim=nsim)

# training data for province 10, years 2006 - 2012 
dat <- read.csv('data/province-biweek-counts-training.csv')
dat$date_sick <- as.Date(strptime(dat$date_sick,"%m/%d/%Y")) # convert to date values
inc <- IncidenceMatrix$new(1+reshape2::acast(dat,province~date_sick,value.var='cases'))

# define how many provinces to model
# this demo only has data for one province: Province 10
prov_nums <- 1
nmodels <- length(prov_nums)
sarima_model$fit(inc$subset(rows = prov_nums, mutate = FALSE))

##########################
# ~~~~ Model Forecast ~~~~
###########################

steps <- 26 # forecast ahead 26 biweeks
forecast_X <- sarima_model$forecast(steps = steps)

# add date_sick to forecast_X$data$mat by first
# converting predictions to a dataframe and using dplyr
preds_df <- data.frame(as.table(t(forecast_X$data$mat)))

# import testing dataset which has 2013 data
data_X <- read.csv('data/province-biweek-counts-testing.csv')
data_X$date_sick <- strptime(data_X$date_sick,"%m/%d/%Y") # convert to date values
data_X$date_sick <- as.Date(data_X$date_sick)
preds_dates <- data_X %>%
          filter(year == 2013)
# add prediction dates to original forecast
preds_df[["date_sick"]] <-as.Date(preds_dates$date_sick, format = "%d-%m-%Y")

# add prediction columns to original forecast
preds_df <- preds_df %>%
   mutate( 
      pred_total_cases = as.vector(forecast_X$median()$mat),
      pred_95_lb = as.vector(forecast_X$quantile(0.025)$mat),
      pred_95_ub = as.vector(forecast_X$quantile(0.975)$mat),
      pred_80_lb = as.vector(forecast_X$quantile(0.05)$mat),
      pred_80_ub = as.vector(forecast_X$quantile(0.95)$mat),
      pred_50_lb = as.vector(forecast_X$quantile(0.25)$mat),
      pred_50_ub = as.vector(forecast_X$quantile(0.75)$mat)
   )

##########################
# ~~~~ Model Plot ~~~~
###########################
ggplot() +
  geom_ribbon(
    mapping = aes(x = date_sick, ymin = pred_95_lb, ymax = pred_95_ub),
    fill = "cornflowerblue",
    alpha = 0.2,
    data = preds_df) +
  geom_ribbon(
    mapping = aes(x = date_sick, ymin = pred_80_lb, ymax = pred_80_ub),
    fill = "cornflowerblue",
    alpha = 0.2,
    data = preds_df) +
  geom_ribbon(
    mapping = aes(x = date_sick, ymin = pred_50_lb, ymax = pred_50_ub),
    alpha = 0.2,
    data = preds_df) +
  geom_line(
    mapping = aes(x = date_sick, y = pred_total_cases),
    color = "royalblue",
    size = 2,
    data = preds_df) +
  geom_line(mapping = aes(x = date_sick, y = cases),
            size=0.7,
            data = data_X) +
  xlab("") + ylab("Number of Cases") +
  coord_cartesian(ylim = c(0, 1000))