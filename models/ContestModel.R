require(ForecastFramework)
require(R6)

ContestModel <- R6Class(
  inherit = ForecastModel,
  private = list(
    .nsim = integer(1),
    forecastRun = FALSE,
    newdata = IncidenceMatrix$new(),
    output = SimulatedIncidenceMatrix$new()
  ),
  public = list(
    #' @method peakTime
    #' @param newdata
    #' @param start Time Step relative to end of the data to start.  May be negative.  For example, -3 would mean 3 time steps before the last one in the data.  Note that the start is included in the range.
    #' @param end Time Step relative to the end of the data to stop.  May be negative.  For example 6 would mean 6 time steps after the last time step in the data.  Note that the end is included in the range.
    #' @return A SimulatedForecast containing the time step containing the highest incidence.  These should be stored as ordinal time steps beginning with start.  So if start is -3, and the answer is 1, then this value should be 5.
    peakTime = function(newdata,start,end){
      if("peakTime" %in% private$.debug){browser()}
      private$newdata = IncidenceMatrix$new(newdata)
      dnames = list(
        private$newdata$rnames,
        paste('[',start,',',end,')'),
        NULL
      )
      if(start <= 0){
        private$newdata$tail(1-start)
        prev.max = apply(private$newdata$mat,1,max)
        prev.time = apply(private$newdata$mat,1,which.max)
      } else {
        prev.max = 0
        prev.time = 1
      }
      if(end > 0){
        self$forecast(newdata,end)
        cur.max = apply(private$output$arr,c(1,3),max)
        cur.time = apply(private$output$arr,c(1,3),which.max)
      } else {
        cur.max = matrix(0,private$newdata$nrow,self$nsim)
        cur.time = matrix(1,private$newdata$nrow,self$nsim)
      }
      data <- SimulatedIncidenceMatrix$new(
        array(
          (matrix(prev.max,nrow(cur.max),ncol(cur.max)) - cur.max > 0) * prev.time +
            (matrix(prev.max,nrow(cur.max),ncol(cur.max)) - cur.max <= 0) * cur.time,
          c(private$newdata$nrow,1,self$nsim),
          list(
            private$newdata$rnames,
            'peak',
            NULL
          )
        )
      )
      data$dnames <- dnames
      
      return(
        IncidenceForecast$new(
          data,
          forecastTimes = TRUE
        )
      )
      private$defaultAbstract()
    },
    #' @method peakIncidence 
    #' @param newdata
    #' @param start Time Step relative to end of the data to start.  May be negative.  For example, -3 would mean 3 time steps before the last one in the data.  Note that the start is included in the range.
    #' @param end Time Step relative to the end of the data to stop.  May be negative.  For example 6 would mean 6 time steps after the last time step in the data.  Note that the end is included in the range.
    #' @return A SimulatedForecast containing the projected incidence at the time step containing the highest incidence.
    peakIncidence = function(newdata,start,end){
      if("peakIncidence" %in% private$.debug){browser()}
      private$newdata = IncidenceMatrix$new(newdata)
      dnames = list(
        private$newdata$rnames,
        paste('[',start,',',end,')'),
        NULL
      )
      if(start <= 0){
        private$newdata$tail(1-start)
        prev.max = apply(private$newdata$mat,1,max)
      } else {
        prev.max = 0
      }
      if(end > 0){
        self$forecast(newdata,end)
        cur.max = apply(private$output$arr,c(1,3),max)
      } else {
        cur.max = matrix(0,private$newdata$nrow,self$nsim)
      }
      data <- SimulatedIncidenceMatrix$new(
        array(
          (matrix(prev.max,nrow(cur.max),ncol(cur.max)) - cur.max > 0) * prev.max +
            (matrix(prev.max,nrow(cur.max),ncol(cur.max)) - cur.max <= 0) * cur.max,
          c(private$newdata$nrow,1,self$nsim),
          list(
            private$newdata$rnames,
            'peak',
            NULL
          )
        )
      )
      data$dnames <- dnames
      
      return(
        IncidenceForecast$new(
          data,
          forecastTimes = TRUE
        )
      )
      private$defaultAbstract()
    },
    totalIncidence = function(newdata,start,end){
      if("totalIncidence" %in% private$.debug){browser()}
      private$newdata = IncidenceMatrix$new(newdata)
      dnames = list(
        private$newdata$rnames,
        paste('[',start,',',end,')'),
        NULL
      )
      if(start <= 0){
        private$newdata$tail(1-start)
        prev.sum = apply(private$newdata$mat,1,sum,na.rm=T)
      } else {
        prev.sum = 0
      }
      if(end > 0){
        self$forecast(newdata,end)
        cur.sum = apply(private$output$arr,c(1,3),sum)
      } else {
        cur.sum = matrix(0,private$newdata$nrow,self$nsim)
      }
      data <- SimulatedIncidenceMatrix$new(
        array( prev.sum + cur.sum,
          c(private$newdata$nrow,1,self$nsim),
          list(
            private$newdata$rnames,
            'total',
            NULL
          )
        )
      )
      data$dnames <- dnames
      
      return(
        IncidenceForecast$new(
          data,
          forecastTimes = TRUE
        )
      )
      private$defaultAbstract()
    }
  ),
  active = list(
      nsim = function(value){
        if(missing(value)){return(private$.nsim)}
        private$defaultActive('.nsim','private',as.integer(value))
      }
  )
)

#' @description This class is for testing a model.  If a shifted version of a model scores better than the unshifted version, it likely means there's a problem in the code somewhere.
ShiftedModel <- R6Class(
  inherit = ForecastModel,
  private = list(model = ForecastModel$new(),shift=numeric(0)),
  public = list(
    initialize = function(model,shift){
      self$model = model
      self$shift = shift
    },
    forecast = function(newdata,steps){
      return(self$model$forecast(newdata=newdata,steps=steps+self$shift))
    },
    peakTime = function(newdata,start,end){
      return(self$model$peakTime(newdata=newdata,start=start+self$shift,end=end+self$shift))
    },
    peakIncidence = function(newdata,start,end){
      return(self$model$totalIncidence(newdata=newdata,start=start+self$shift,end=end+self$shift))
    },
    totalIncidence = function(newdata,start,end){
      return(self$model$totalIncidence(newdata=newdata,start=start+self$shift,end=end+self$shift))
    }
  )
)
