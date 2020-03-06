#'Predict sales for a new data
#'
#'Predict sales for a new data.
#'
#'@param object A FAIRforecast object obtained with \code{\link{FAIR_train}}.
#'@param new_data A dataframe object that contains observations with
#'the following variables:
#'store identifier (factor), category identifier (factor), time id (numeric),
#'seasonality and calendar variables (factor or numeric),
#'and marketing variables (numeric).
#'@param parallel Logical. Should parallel computing be used? The default is
#'  \code{F}
#'
#'@details All components of the data is supplied in one data.frame and the
#'  column names should correspond to the ones used to obtain \code{object}.
#'@details \code{FAIR_predict} can forecast for periods which are after the
#'  last period of the training data and are in the forecasting horizon.
#'
#'@return A data.frame of the final forecast and its components for each
#'  \code{store}-\code{category}-\code{time_id} combination in \code{new_data}.
#'  The columns are:
#'  \item{Seasonal_and_trend_pattern}{The forecast based on established
#'  seasonality and trend in sales and marketing. Obtained with Stage 1 and
#'  adjusted for bias.}
#'  \item{Marketing_deviation_multiplier}{Multiplier based on the marketing
#'  scenario deviation from the established pattern.}
#'  \item{Disturbance_multiplier}{Multiplier reflecting the changes in base
#'  sales.}
#'  \item{Forecast}{The final forecast obtained by multiplying the components
#'  and subtracting 1. If more complex parts of the forecast (multipliers) are
#'  set to NA, they are ignored.}
#'
#' @references Gür Ali, Ö. and Gürlek, R. (2020) Automatic Interpretable Retail
#' Forecasting (FAIR) with Promotional Scenarios. International Journal of
#' Forecasting, forthcoming.
#'
#'@seealso \code{\link{FAIR_train}}
#'
#'@export
FAIR_predict <- function(object, new_data, parallel = F) {
  pars <- object$pars
  list2env(pars, environment())
  lag_data <- object$lag_data
  list2env(object$models, environment())
  remove(object)

  temp <-
    new_data[, c(marketing)]
  if (sum(!sapply(temp, is.numeric)) > 0){
    stop(
      "Sales and marketing columns should be numeric!"
    )
  }
  if (sum(temp < 0) > 0) {
    stop(
      "Sales and marketing columns should be non-negative!"
    )
  }

  if (max(new_data[, time_id]) > (last_week + horizon) |
      min(new_data[, time_id]) <= last_week) {
    stop("Any 'time_id' in 'new_data' should be greater than 'max(train_data[, time_id])' and no greater than 'max(train_data[, time_id]) + horizon'")
  }

  new_data[, sales] <- NA
  new_data <-
    new_data[, c(store, category, time_id,
                 seasonality[seasonality != time_id],
                 marketing, sales)]
  if (!setequal(colnames(lag_data), colnames(new_data))) {
    stop("Colunm names of new_data and the training data do not match!")
  }
  temp <- rbind(lag_data, new_data)
  remove(lag_data)
  category_list <- unique(temp[, category])

  deseason_new <- function(my_data) {
    my_category <- my_data[, category][1]
    my_store <- my_data[, store][1]
    var_list <- c(sales, marketing)
    fit <- function(my_var) {
      my_model <- Stage1[[paste0(my_store, ".", my_category)]][[my_var]]
      if (is.numeric(my_model))
        return(rep(my_model, nrow(my_data)))
      if (is.null(my_model))
        return(rep(NA, nrow(my_data)))
      return(suppressWarnings(predict.lm(my_model, my_data)))
    }

    new_vars <- lapply(var_list, fit)
    new_vars <- do.call(data.frame, new_vars)
    my_data[, var_list] <- new_vars
    colnames(my_data)[colnames(my_data) == sales] <- "prediction_1"

    return(my_data)
  }

  #deseasonalize data
  temp <- plyr::ddply(temp, c(category, store), deseason_new,
                      .parallel = parallel)

  #append marketing variables of other categories and required lags as columns
  temp <-
    var_cre(
      temp,
      category = category,
      category_list = category_list,
      lag = lag,
      time_id = time_id,
      marketing = marketing,
      cc_marketing = cc_marketing,
      store = store,
      parallel = parallel
    )
  temp <- temp[temp[, time_id] > last_week,]

  my_index <-
    match(paste0(new_data[, store], new_data[, category], new_data[, time_id]),
          paste0(temp[, store], temp[, category], temp[, time_id]))
  new_data$prediction_1 <- temp[my_index, "prediction_1"]

  #step2
  ext_marketing <-
    c(marketing, ext_marketing) #marketing variables are excluded
  #in FAIR_train function
  step2_new <- function(my_data) {
    my_category <- my_data[1, category]

    for (i in cc_marketing) {
      #exclude focal marketing dublicates
      ext_marketing <-
        ext_marketing[ext_marketing != paste0(my_category, "_", i)]
    }

    my_model <- Stage2[[my_category]]
    my_data$prediction_2 <-
      as.numeric(glmnet::predict.glmnet(
        my_model, as.matrix(my_data[, ext_marketing])))
    my_data$residual_2 <-
      my_data[, "prediction_1"] - my_data$prediction_2
    return(my_data)
  }

  temp <-
    plyr::ddply(temp, category, step2_new, .parallel = parallel)

  #step3
  ext_marketing <- ext_marketing[!(ext_marketing %in% marketing)]
  temp <- temp[,!(colnames(temp) %in% ext_marketing)]

  step3_new <- function(x) {
    my_store <- x[1, store]
    my_category <- x[1, category]
    my_model <- Stage3[[paste0(my_store, ".", my_category)]]
    if (is.null(my_model))
      return(NULL)
    predictions <-
      as.numeric(forecast::forecast(my_model, h = horizon)$mean)
    predictions <-
      data.frame((last_week + 1):(last_week + horizon), predictions)
    colnames(predictions) <- c(time_id, "prediction_3")
    x <- merge(x, predictions, by = time_id, all.x = T)
    return(x)
  }

  temp <- plyr::ddply(temp, c(store, category), step3_new)
  temp$residual_3 <- temp$residual_2 - temp$prediction_3

  #create outputs
  new_data <- merge(new_data[, c(store, category, time_id, "prediction_1")],
                    temp[, c(store,
                             category,
                             time_id,
                             "prediction_2",
                             "prediction_3")],
                    all.x = T)
  colnames(new_data)[4:6] <-  c("Seasonal_and_trend_pattern",
                                "Marketing_deviation_multiplier",
                                "Disturbance_multiplier")
  new_data <-
    merge(new_data,
          in_sample_bias,
          by = c(store, category),
          all.x = T)
  new_data$Seasonal_and_trend_pattern <-
    exp(new_data$Seasonal_and_trend_pattern) * new_data$Bias_cor
  new_data$Marketing_deviation_multiplier <-
    exp(new_data$Marketing_deviation_multiplier)
  new_data$Disturbance_multiplier <- exp(new_data$Disturbance_multiplier)
  new_data$Forecast <- (new_data$Seasonal_and_trend_pattern *
                        ifelse(is.na(new_data$Marketing_deviation_multiplier), 1,
                               new_data$Marketing_deviation_multiplier) *
                        ifelse(is.na(new_data$Disturbance_multiplier), 1,
                               new_data$Disturbance_multiplier) - 1)
  new_data <- new_data[, !colnames(new_data) %in% "Bias_cor"]

  return(new_data)
}
