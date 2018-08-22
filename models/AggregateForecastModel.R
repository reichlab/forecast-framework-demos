library(ForecastFramework)
############################AggregateForecastModel##########################
#' @title AggregateForecastModel
#' @description A collection of models for recursively predicting multiple
#'   timesteps into the future.  Instead of using the same model to forecast
#'   every row, this model uses a different model for each row.
#' @docType class
#' @importFrom R6 R6Class
#' @export AggregateForecastModel
#' @keywords model forecast
#' @family Models ForecastModels
#' @example AggregateForecastModel.R
AggregateForecastModel <- R6Class(
  classname="AggregateForecastModel",
  inherit=RecursiveForecastModel,
  private = list(
    .models = list(),
    .nrow = 0
  ),
  public = list(
    #' @method fit_ Fit the model for predicting each of the rows. Assumes the data has been  put into place.
    #' @param private$.data The data to use in fitting the model is here. private$newdata is an AbstractIncidenceMatrix object, which is mainly a matrix.  The dimensions correspond to Variable - Named row with the variable name as rowname. Instance - Each instance represents a slice of data to be predicted.
    fit_ = function(){
      for(row in 1:private$.nrow){
        self$fitRow_(row)
      }
    },
    #' @method fitRow Fit the model to a particular row.
    #' @param data The data to use in fitting the model is here. private$newdata is an AbstractIncidenceMatrix object, which is mainly a matrix.  The dimensions correspond to Variable - Named row with the variable name as rowname. Instance - Each instance represents a slice of data to be predicted.
    #' @param row The row to fit a model for.
    fitRow = function(data,row){
      self$prepareFitData(data)
      self$fitRow_(row)
    },
    #' @method fitRow_ This method \bold{must} be extended.  Fit a model for a particular row.  Assumes the data has been put into place
    #' @reference private$newdata The data to use in fitting the model is here. private$newdata is an IncidenceMatrix object, which is mainly a 3D array.  The three dimensions correspond to Variable - Named row with the variable name as rowname. Instance - Each instance represents a slice of data to be predicted
    #' @param row The row to predict.
    fitRow_ = function(row){
      private$.defaultAbstract()
    }
  ),
  active = list(
    #' @field models  A list of models, over which to aggregate 
    models = function(value){
      if(missing(value)){
        return(private$.models)
      }
      stop("Use fit to generate the list of models.")
    }
  )
)
