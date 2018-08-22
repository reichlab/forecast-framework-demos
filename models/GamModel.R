library(mgcv)
library(stats)
library(dplyr)
library(reshape2)
source('../models/AggregateForecastModel.R')
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##             Spatial Prediction Model Object Oriented Framework             ##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
#' @include Models.R
#' @include IncidenceMatrix.R SimulatedIncidenceMatrix.R
####################################GamModel###################################
#' @title GamModel
#' @description A compartmental autoregressive model, which fits lagged versions
#'   of rows to the data.  Instead of attempting to use all of the lagged rows,
#'   at each time step, it selects the rows with the highest correlation.
#' @docType class
#' @importFrom R6 R6Class
#' @export GamModel
#' @exportClass GamModel
#' @keywords model forecast
#' @family Models ForecastModels
GamModel <- R6Class(
  classname="GamModel",
  inherit=AggregateForecastModel,
  private = list(
    .data = IncidenceMatrix$new(),
    newdata = IncidenceMatrix$new(),
    transdata = IncidenceMatrix$new(),
    output = SimulatedIncidenceMatrix$new(),
    .numProvinces = as.integer(3),
    .nsim = as.integer(3),
    .knots = as.integer(3),
    .topProvinceIndex = list(),
    .topProvinceNames = list(),
    .models = list(),
    .scale = "log",
    .ignoreSelf = F,
    .useDeltas = F,
    ##Switch to using these instead of that loop
    transform = NA,
    untransform = NA,
    .useOldCorCalculation=T
  ),
  public = list(
    initialize = function(
      predCols = c(1),
      numProvinces=3,
      nSims=3,
      scale='log',
      useDeltas=F,
      ignoreSelf=F,
      useOldCorCalculation=T,
      splineKnots = 10,
      stochastic='Poisson'
    ){
      "Create a new GamModel."
      "@param predCols A vector containing the lags which are used to
        predict a particular column.  0 represents the column itself, 1 the
        previous  column etc."
      "@param numProvinces The number of provinces to use as part of the
        model for each time lag."
      "@param nSims The number of simulations to use as part of forecasting."
      "@param scale The transformation necessary for model fitting."
      "@param ignoreSelf Whether or not to allow each province to be
        correlated with itself"
      "@param useOldCorCalculation Whether or not to calculate the top
        provinces using the old method or the new method."
      "@param splineKnots.  Passed to gam.  Defines the number of knots in the
         cyclic spline."
      ##Update scale description
      self$predCols = as.integer(predCols)
      self$numProvinces= as.integer(numProvinces)
      self$nsim = as.integer(nSims)
      self$scale = scale
      self$ignoreSelf = ignoreSelf
      self$useDeltas = useDeltas
      self$useOldCorCalculation=T
      self$knots = as.integer(splineKnots)
      self$stochastic = stochastic
    },
    fitRow_ = function(row){
      "Fit a model for a particular row.  Assumes the data has been put into
         place "
      "@reference private$transdata The data to use in fitting the model is
         here. private$transdata is an IncidenceMatrix object, which is mainly
         a 3D array."
      if('fitRow_' %in% private$.debug){
        browser()
      }
      ##For now, assuming ignore.self = F
      if(private$.ignoreSelf){
        ##stop("ignoreSelf=T is not currently implemented")
      }
      ##There's a problem with this
      ##Replace with: private$transdata$rowData$topProvinceIndex[[row]]
      ## We need to preallocate the topProvinceIndex somewhere...  Probably
      ##   prepareFitData
      ##private$.topProvinceIndex = array(NA,c(
      ##    private$.data$nrow,
      ##    private$transdata$ncol,
      ##    self$numProvinces
      ##))
      ## This only 
      if(private$.numProvinces > 0){
        for(i in 1:length(private$.predCols)){
          ##Extract the indices of the row names
          #' @importFrom dplyr starts_with
          tempRnames = private$transdata$rnames %>%
            starts_with(
              match=paste("L",private$.predCols[[i]],sep=''),
              ignore.case=TRUE
            )
          ##Find the top correlated provinces for each lag
          #' @importFrom stats cor
          ##This is for doing the cor call the way that spatialpred does it.
          ##I DON'T like the idea of initializing something in the fitRow_ method
          ##  like this.  This is a hack and should be fixed.
          tmpTopProvinces= array(NA,private$.numProvinces+private$.useDeltas)
          if(private$.useOldCorCalculation){
            ##We need to deal with useDeltas separately in this case
            if(private$.useDeltas){
              ##Using predCols[[i]] instead of i here wastes space.  This should
              ##  be changed later
              tmpTopProvinces =
                tempRnames[(
                  (
                    -cor(
                      ##Because of this subtraction
                      private$transform(private$.data$mat[row,]) -
                        private$transform(
                          private$.data$lag(indices=1,mutate=F,na.rm=F)$mat[row,]
                        ),
                      t(private$transdata$mat[tempRnames,]),
                      use='pairwise.complete.obs'
                    ) %>%
                      sort(index.return=T)
                  )$ix
                )[1:(private$.numProvinces+private$.ignoreSelf)]
                ]
            } else{
              ##This is the old case still, but no useDeltas
              tmpTopProvinces =
                tempRnames[(
                  (
                    -cor(
                      ##So just transforming is enough here.
                      private$transform(private$.data$mat[row,]),
                      t(private$transdata$mat[tempRnames,]),
                      use='pairwise.complete.obs'
                    ) %>%
                      sort(index.return=T)
                  )$ix
                )[1:(private$.numProvinces+private$.ignoreSelf)]
                ]
            }
          } else {
            ##This is the new case.  Here we don't transform the output variable
            tmpTopProvinces =
              tempRnames[(
                (
                  -cor(
                    ##See, no private$transform call
                    private$.data$mat[row,],
                    t(private$transdata$mat[tempRnames,]),
                    use='pairwise.complete.obs'
                  ) %>%
                    sort(index.return=T)
                )$ix
              )[1:(private$.numProvinces+private$.ignoreSelf)]
              ]
          }
          if(private$.ignoreSelf){
            if(tempRnames[row] %in% tmpTopProvinces){
              private$.topProvinceIndex[row,private$.predCols[[i]],] =
                tmpTopProvinces[tmpTopProvinces != tempRnames[row]]
            } else{
              private$.topProvinceIndex[row,private$.predCols[[i]],] =
                tmpTopProvinces[-length(tmpTopProvinces)]
            }
          } else{
            private$.topProvinceIndex[row,private$.predCols[[i]],] =
              tmpTopProvinces[]
          }
        }
        ##private$.topProvinceIndex[[row,lag]] = unique(
        ##  1+((-cor(
        ##        private$.data$mat[row,],
        ##        t(private$transdata$mat),
        ##        use='complete.obs'
        ##  ) %>% sort(index.return=T)
        ##)$ix %% private$.nrow))[1:private$.numProvinces]
        private$.topProvinceNames[[row]] =
          private$transdata$rnames[private$.topProvinceIndex[row,,][]]
      } else {
        private$.topProvinceNames = lapply(1:self$data$nrow,function(x){c()})
      }
      # eqn = paste("y~s(time.in.year,bs=\"cc\",k=",private$.knots,")",sep='')
      eqn = paste("y~s(time.in.year",")",sep='')
      rownames = paste(private$.topProvinceNames[[row]],collapse=" + ")
      paste(eqn,rownames,sep= " + ") -> eqn
      ##replicate(private$.numProvinces,private$.predCols) %>%
      ##  aperm(c(2,1)) ->
      ##  lagNames
      ##rownames =paste(
      ##  "L",lagNames,"R",private$.topProvinceNames[[row]],
      ##  sep="",
      ##  collapse= " + "
      ##)
      ##paste(eqn," + ",rownames) -> eqn
      #' @importFrom stats offset
      ##This line means we *must* create an offset column separately.
      ##eqn = paste(eqn," + ","offset(log(1+off))",sep='')
      eqn = paste(eqn," + ","offset(off)",sep='')
      #' @importFrom mgcv gam
      private$.models[[row]] <- gam(
        #' @importFrom stats formula
        formula(eqn),
        data=data.frame(
          t(private$transdata$mat),
          y=self$data$mat[row,],
          ##This next one is a hack,
          off=log(1+c(NA,self$data$mat[row,-self$data$ncol])),
          year=private$.data$colData$year,
          time.in.year=private$.data$colData$time.in.year
        ),
        #' @importFrom stats poisson
        family=poisson
        ##family=quasipoisson
      )
    },
    predict_ = function(col){
      "Using a model previously fit with \\code{fit} to predict each row of the
         next column. This function assumes that all of the data preprocessing
         has already been taken care of."
      "@reference private$transdata The data to use in fitting the model is
         here. private$transdata is an SimulatedIncidenceMatrix object, which is
         mainly a 3D array."
      if('predict_' %in% private$.debug){
        browser()
      }
      if(missing(col)){
        col = 1:private$transdata$ncol
      }
      rows=1:private$.nrow
      ##tmp should be given a home in the GamModel, and this should be done in
      ##  prepareFitData
      private$transdata$simulations %>%
        #' @importFrom reshape2 melt
        melt() %>%
        #' @importFrom reshape2 acast
        #' @importFrom stats formula
        acast(formula=formula('Var1 ~ ...')) %>%
        t %>%
        as.data.frame ->
        tmp
      private$newdata$simulations %>%
        #' @importFrom reshape2 melt
        melt() %>%
        #' @importFrom reshape2 acast
        #' @importFrom stats formula
        acast(formula=formula('Var1 ~ ...')) %>%
        t %>%
        as.data.frame ->
        tmp2
      ##There is a bug when ncol, and nsim are both 1.  The resulting array is
      ##  demoted to a matrix and so the aperm call doesn't work.  Maybe figure
      ##  out something else to do in that case, or use an lapply call...
      if(private$.nsim * private$transdata$ncol > 1){
        vapply(
          1:private$.nrow,
          function(row){
            predict(
              object = private$.models[[row]],
              newdata = data.frame(tmp,off=log(1+tmp2[,row])),
              newdata.garaunteed=TRUE,
              type='response'
            )
          },
          FUN.VALUE = array(
            1.,
            c(
              private$transdata$ncol,
              private$.nsim
            )
          ),
          USE.NAMES = FALSE
        ) %>% aperm(c(3,1,2)) %>%
          private$output$mutate(cols=col,rows=rows,sims=1:private$.nsim)
      } else{
        vapply(
          1:private$.nrow,
          function(row){
            predict(
              object = private$.models[[row]],
              newdata = data.frame(tmp,off=log(1+tmp2[,row])),
              newdata.garaunteed=TRUE,
              type='response'
            )
          },
          FUN.VALUE = array(
            1.,
            c(
              private$transdata$ncol,
              private$.nsim
            )
          ),
          USE.NAMES = FALSE
        ) %>%
          private$output$mutate(cols=col,rows=rows)
      }
    },
    predictRow_ = function(row,col){
      "Using a model previously fit with \\code{fitRow} to predict the
       \\code{row}th row of the next column. This function assumes that all of
       the data preprocessing has already been taken care of."
      "@param \\code{row} The row to predict the value of."
      "fit the model for predicting all of the rows. Assumes the data has been
         put into place."
      "@reference private$newdata The data to use in fitting the model is here.
         private$newdata is an SimulatedIncidenceMatrix object, which is mainly
         a 3D array."
      "@mutate private$output The output data.  The dimensions of the output
         matrix should always be m by n by k where
         m is private$.nrow (the number of locations in the data),
         n is the number of instances in private$newdata
         k is the number of simulations."
      self$predict_(col)
      return()
      stop("predictRow is broken right now.")
      ##This is the old code for when I come back to it.  Its never worked, but
      ##  the idea is that we take the function from predict_ and run it here
      ##  directly.  However, we need to generate tmp and tmp2 for the code to
      ##  run.
      predict(
        object = private$.models[[row]],
        newdata = data.frame(tmp,off=log(1+tmp2[,row])),
        newdata.garaunteed=TRUE,
        type='response'
      ) %>%
        aperm(c(3,1,2)) %>%
        private$output$mutate(cols=col,rows=row)
    },
    prepareFitData = function(data){
      "Prepares the internal \\code{data} within the object in order to allow
         for fitting and predicting."
      "@param \\code{data} A IncidenceMatrix of data to prepare"
      "@reference \\code{predCols} Which columns to use for prediction.  See
         \\code{IncidenceMatrix$lag} for details."
      "@mutate \\code{self$data} The data used to fit the model will be stored
         here."
      "@mutate \\code{private$.transdata} The data with predCols accounted for
         will be stored here."
      "@mutate \\code{private$.nrow} We update this so the model knows how many
         rows the data should have."
      if('prepareFitData' %in% private$.debug){
        browser()
      }
      ##Consider adding an option to not store the fit data.
      private$.data = IncidenceMatrix$new(data)
      private$transdata = self$data$clone(TRUE)
      ## Error checking
      if(!all(c('t','time.in.year','year') %in% names(self$data$colData))){
        stop("This model requires input data to provide a t, a year, and time.in.year for each column")
      }
      if(!all(c('t.step','max.year.time') %in% names(self$data$metaData))){
        stop("This model requires input data to provide a t.step, and max.year.time")
      }
      if(self$data$nrow <= self$numProvinces){
        warning(paste("There are not enough rows in the data.  Reducing number of correlated provinces from ",self$numProvinces,"to",self$data$nrow-1))
        self$numProvinces <- as.integer(self$data$nrow - 1)
      }
      ##This next line replaces the commented out lines below
      ##See active::scale for details of how it works
      private$transdata$scale(private$transform)
      ##if(private$.scale == 'log'){
      ##  private$transdata$scale(function(x){log(x+1)})
      ##} else if(private$.scale == 'sqrt'){
      ##  private$transdata$scale(function(x){sqrt(x)})
      ##}else if(private$.scale != 'identity'){
      ##  stop("This scale is not currently supported.")
      ##}
      if(private$.useDeltas){
        private$transdata$diff()
      }
      private$.nrow = private$.data$nrow
      private$transdata$lag(indices=self$predCols,na.rm=FALSE)
      private$.topProvinceIndex = array(
        NA,
        c(
          private$.data$nrow,
          length(private$.predCols),
          self$numProvinces
        )
      )
    },
    preparePredictData = function(data){
      "Prepares the internal \\code{data} within the object in order to allow
         for fitting and predicting."
      "@param \\code{data} A IncidenceMatrix of data to prepare"
      "@param \\code{simulations} The number of Monte Carlo simulations to
         prepare for."
      "@reference \\code{predCols} Which columns to use for prediction.  See
         \\code{IncidenceMatrix$lag} for details."
      "@mutate \\code{private$.newdata} The data with predCols accounted for
         will be stored here."
      if('preparePredictData' %in% private$.debug){
        browser()
      }
      ##Consider checking to see if useDeltas is false and scale is 'identity'
      ##  here...
      if(!missing(data)){
        private$newdata = SimulatedIncidenceMatrix$new(data,private$.nsim)
      }
      if(!all(c('t','time.in.year','year') %in% names(private$newdata$colData))){
        stop("This model requires input data to provide a year and time.in.year for each column")
      }
      if(!all(c('t.step','max.year.time') %in% names(private$newdata$metaData))){
        stop("This model requires input data to provide a year and time.in.year for each column")
      }
      private$transdata = private$newdata$clone(TRUE)
      if(private$.nrow != private$transdata$nrow){
        stop("The data to predict musthave the same format as the fit data.")
      }
      ##We lag the newdata matrix so it is easier to find the offset term
      ##The following code is reused.  It should be cleaned up.  cols is
      ##  always 1, so most of the formulae below should simplify easily
      cols = 1
      private$newdata$addColumns(cols)
      private$newdata$colData$t[private$newdata$ncol-cols:1+1] =
        private$newdata$colData$t[private$newdata$ncol-cols] +
        private$newdata$metaData$t.step * 1:cols
      private$newdata$colData$year[private$newdata$ncol-cols:1+1] =
        floor(private$newdata$colData$t[private$newdata$ncol-cols:1+1])
      private$newdata$colData$time.in.year[private$newdata$ncol-cols:1+1] =
        (
          (private$newdata$colData$time.in.year[private$newdata$ncol-cols] +
             1:cols - 1
          ) %% private$newdata$metaData$max.year.time) + 1
      private$newdata$lag(1,na.rm=FALSE)
      ##Add dummy columns to the end of the data so that when we lag we get
      ##  of the data.
      if(min(self$predCols) > 0){
        cols = min(self$predCols)
        private$transdata$addColumns(cols)
        private$transdata$colData$t[private$transdata$ncol-cols:1+1] =
          private$transdata$colData$t[private$transdata$ncol-cols] +
          private$transdata$metaData$t.step * 1:cols
        private$transdata$colData$year[private$transdata$ncol-cols:1+1] =
          floor(private$transdata$colData$t[private$transdata$ncol-cols:1+1])
        private$transdata$colData$time.in.year[private$transdata$ncol-cols:1+1] =
          (
            (private$transdata$colData$time.in.year[private$transdata$ncol-cols] +
               1:cols - 1
            ) %% private$transdata$metaData$max.year.time) + 1
      }
      ##This next line replaces the commented out lines below
      ##See active::scale for details of how it works
      private$transdata$scale(private$transform)
      ##if(private$.scale == 'log'){
      ##  private$transdata$scale(function(x){log(x+1)})
      ##} else if(private$.scale == 'sqrt'){
      ##  private$transdata$scale(function(x){sqrt(x)})
      ##} else if(private$.scale != 'identity'){
      ##  stop("This scale is not currently supported.")
      ##}
      if(private$.useDeltas){
        private$transdata$diff(1)
      }
      ##We have to add the offset now, which is the raw case counts.
      private$transdata$lag(indices=private$.predCols,na.rm=FALSE)
      ##Rows to add: time.in.year, year, and off
      private$transdata$addRows(3)
      ##For the first two rows, we put the time.in.year and year
      private$transdata$mutate(
        rows=private$transdata$nrow-1:2,
        data=matrix(
          c(
            private$transdata$colData$time.in.year,
            private$transdata$colData$year
          ),
          2,
          private$transdata$ncol,
          dimnames=list(c('time.in.year','year'),private$.cnames),
          byrow=TRUE
        )
      )
      private$output = SimulatedIncidenceMatrix$new(
        IncidenceMatrix$new(
          matrix(NA,private$.data$nrow+2,private$newdata$ncol)
        ),
        private$.nsim
      )
    },
    prepareForecastData = function(data){
      "Prepares the internal \\code{data} within the object in order to allow
         for fitting and predicting."
      "@param \\code{data} A IncidenceMatrix of data to prepare"
      "@reference \\code{predCols} Which columns to use for prediction.  See
         \\code{IncidenceMatrix$lag} for details."
      "@mutate \\code{private$.transdata} The data with predCols accounted for
         will be stored here."
      if('prepareForecastData' %in% private$.debug){
        browser()
      }
      ##This should be removed...
      ##if(!(1 %in% private$.predCols)){
      ##  private$.predCols = c(private$.predCols,1)
      ##}
      private$newdata = SimulatedIncidenceMatrix$new(data,private$.nsim)
      if(!all(c('t','time.in.year','year') %in% names(private$newdata$colData))){
        stop("This model requires input data to provide a year and time.in.year for each column")
      }
      if(!all(c('t.step','max.year.time') %in% names(private$newdata$metaData))){
        stop("This model requires input data to provide a year and time.in.year for each column")
      }
      private$transdata = SimulatedIncidenceMatrix$new(data,private$.nsim)
      cols = 1
      private$newdata$addColumns(cols)
      private$newdata$colData$t[private$newdata$ncol-cols:1+1] =
        private$newdata$colData$t[private$newdata$ncol-cols] +
        private$newdata$metaData$t.step * 1:cols
      private$newdata$colData$year[private$newdata$ncol-cols:1+1] =
        floor(
          private$newdata$colData$t[private$newdata$ncol-cols:1+1]
        )
      private$newdata$colData$time.in.year[private$newdata$ncol-cols:1+1] =
        (
          (
            private$newdata$colData$time.in.year[private$newdata$ncol-cols] +
              1:cols-1
          ) %% private$newdata$metaData$max.year.time
        ) + 1
      private$newdata$lag(1,na.rm=FALSE)
      private$newdata$tail(direction=2,k=1)
      if(min(self$predCols) > 0){
        cols = min(self$predCols)
        private$transdata$addColumns(cols)
        private$transdata$colData$t[private$transdata$ncol-cols:1+1] =
          private$transdata$colData$t[private$transdata$ncol-cols] +
          private$transdata$metaData$t.step * 1:cols
        private$transdata$colData$year[private$transdata$ncol-cols:1+1] =
          floor(
            private$transdata$colData$t[private$transdata$ncol-cols:1+1]
          )
        private$transdata$colData$time.in.year[private$transdata$ncol-cols:1+1] =
          (
            (
              private$transdata$colData$time.in.year[private$transdata$ncol-cols]+
                1:cols - 1
            ) %% private$transdata$metaData$max.year.time
          ) + 1
      }
      ##This next line replaces the commented out lines below
      ##See active::scale for details of how it works
      private$transdata$scale(private$transform)
      ##if(private$.scale == 'log'){
      ##  private$transdata$scale(function(x){log(x+1)})
      ##} else if(private$.scale == 'sqrt'){
      ##  private$transdata$scale(function(x){sqrt(x)})
      ##} else if(private$.scale != 'identity'){
      ##  stop("This scale is not currently supported.")
      ##}
      if(private$.useDeltas){
        private$transdata$diff(1)
      }
      private$transdata$lag(indices=unique(c(1,self$predCols)),na.rm=TRUE)
                                        #time.in.year,year and off
      private$transdata$addRows(3)
      private$transdata$mutate(
        rows=private$transdata$nrow-1:2,
        data=matrix(
          c(
            private$transdata$colData$time.in.year,
            private$transdata$colData$year
          ),
          2,
          private$transdata$ncol,
          dimnames=list(c('time.in.year','year'),private$.cnames),
          byrow=TRUE
        )
      )
      ##private$transdata$mutate(
      ##rows=private$transdata$nrow,
      ##data=NA,
      ##dimnames = list('off')
      ##)
    },
    prepareOutputData = function(inputData,steps){
      if(missing(steps)){
        steps = min(private$.predCols)
      }
      private$output = SimulatedIncidenceMatrix$new(inputData,private$.nsim)
      private$output$addColumns(steps)
      if(steps > 0){
        private$output$colData$t[private$output$ncol-steps:1+1] =
          private$output$colData$t[private$output$ncol-steps] +
          private$output$metaData$t.step * 1:steps
        private$output$colData$year[private$output$ncol-steps:1+1] =
          floor(
            private$output$colData$t[private$output$ncol-steps:1+1]
          )
        private$output$colData$time.in.year[private$output$ncol-steps:1+1] =
          (
            (
              private$output$colData$time.in.year[private$output$ncol-steps] +
                1:steps - 1
            ) %% private$output$metaData$max.year.time
          ) + 1
      }
    },
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
      private$defaultActive('.nsim','private',value)
    },
    topProvinceNames = function(value){
      if(missing(value)){
        return(private$.topProvinceNames)
      }
      stop("Do not write directly to the top provinces")
    },
    topProvinceIndex = function(value){
      if(missing(value)){
        return(private$.topProvinceIndex)
      }
      stop("Do not write directly to the top provinces")
    },
    models = function(value){
      if(missing(value)){
        return(private$.models)
      }
      stop("Do not write directly to the models")
    },
    numProvinces = function(value){
      private$defaultActive('.numProvinces','private',value)
    },
    scale = function(value){
      if(missing(value)){
        return(private$.scale)
      }
      if(!(value %in% c('log','sqrt','identity'))){
        stop("Currently only log, sqrt, and identity are supported.")
      }
      if(value == 'identity'){
        private$transform = function(x){x}
        private$untransform = function(y){y}
      }
      if(value == 'sqrt'){
        private$transform = function(x){sqrt(x)}
        private$untransform = function(y){y^2}
      }
      if(value == 'log'){
        private$transform = function(x){log(1+x)}
        private$untransform = function(y){exp(y)-1}
      }
      private$defaultActive('.scale','private',value)
    },
    ignoreSelf = function(value){
      private$defaultActive('.ignoreSelf','private',value)
    },
    useDeltas = function(value){
      private$defaultActive('.useDeltas','private',value)
    },
    maxPredCol = function(value){
      "The farthest back column from the output value used in prediction."
      if(missing(value)){
        return(
          private$defaultActive('.maxPredCol','private',value)+
            private$.useDeltas
        )
      }
      stop("Do not write directly to the maxPredCol")
    },
    useOldCorCalculation = function(value){
      private$defaultActive('.useOldCorCalculation','private',value)
    },
    knots = function(value){
      private$defaultActive(".knots","private",value)
    }
  )
)
