require(ForecastFramework)
require(R6)
require(forecast)


SARIMATDModel <- R6Class(
  inherit = ContestModel,
  private = list(
    .data = NULL,        ## every model should have this
    .models = list(),    ## specific to models that are fit separately for each location
    .nsim = 1000,         ## models that are simulating forecasts need this
    .lambda = list(),    ## specific to ARIMA models
    .period = integer(0) ## specific to SARIMA models
  ),
  public = list(
    ## data will be MatrixData
    fit = function(data) {
      if("fit" %in% private$.debug){browser()}
      ## stores data for easy access and checks to make sure it's the right class
      private$.data <- IncidenceMatrix$new(data)
      
      ## for each location/row
      for (row_idx in 1:private$.data$nrow) {
        ### need to create a y vector with incidence at time t
        y <- private$.data$subset(rows = row_idx, mutate = FALSE)
        
        ## private$.models[[row_idx]] <- something
        y_ts <- ts(as.vector(y$mat), frequency = private$.period)
        private$.lambda[[row_idx]] <- BoxCox.lambda(y_ts)
        private$.models[[row_idx]] <- fit_sarima(y = y_ts,
                                                 transformation = "box-cox",
                                                 seasonal_difference = TRUE)
      }
      
    },
    forecast = function(newdata = private$.data, steps) {
      if("forecast" %in% private$.debug){browser()}
      
      nmodels <- length(private$.models)
      # sim_forecasts <- SimulatedIncidenceMatrix$new(array(dim = c(nmodels, steps, private$.nsim)))
      sim_forecasts <- array(dim = c(nmodels, steps, private$.nsim))
      dimnames(sim_forecasts) <- list(newdata$rnames, 1:steps, NULL)
      
      for(model_idx in 1:length(private$.models)) {
        tmp_arima <-  simulate(object = private$.models[[model_idx]],
                               nsim = 1000,
                               seed = 1,
                               newdata = as.vector(newdata$mat),
                               h = steps
        )
        tmp_arima <- t(tmp_arima)
        sim_forecasts[model_idx, , ] <- tmp_arima
      }
      
      private$output <- SimulatedIncidenceMatrix$new(sim_forecasts)
      return(IncidenceForecast$new(private$output, forecastTimes = rep(TRUE, steps)))
    },
    initialize = function(period = 26, nsim=100) { 
      ## this code is run during SARIMAModel$new()
      ## need to store these arguments within the model object
      private$.nsim <- nsim
      private$.period <- period
    },
    predict = function(newdata) {
      stop("predict method has not been written.")
    }
  ),
  active = list(
    ## This list determines how you can access pieces of your model object
    data = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the data is not allowed.")
      return(private$.data)
    },
    models = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the models is not allowed.")
      return(private$.models)
    },
    nsim = function(value) {
      ## use this form when you want to be able to change this parameter
      private$defaultActive(type="private", ".nsim", val=value)
    },
    lambda = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to lambda is not allowed.")
      return(private$.lambda)
    },
    period = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the model period is not allowed.")
      return(private$.period)
    }
  )
)