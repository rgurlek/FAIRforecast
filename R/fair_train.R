#'Fit a FAIR model
#'
#'\code{FAIR_train} fits a FAIR model to the provided data.
#'
#'@param train_data a dataframe object that contains observations with the following
#'  variables: store identifier (factor), category identifier (factor), time id
#'  (numeric), seasonality and calendar variables (factor or numeric), marketing
#'  variables (numeric), and sales (numeric).
#'@param store a string. The name of the store identifier
#'  column.
#'@param category a string. The name of the category
#'  identifier column.
#'@param time_id a string. The name of the time
#'  identifier column.
#'@param seasonality a character vector. The names of the seasonality and
#'  calendar columns.
#'@param marketing a character vector. The names of the marketing columns.
#'  The corresponding columns are expected to represent continuous variables and
#'  to be non-negative.
#'@param cc_marketing a character vector. The names of the
#'  cross-category marketing columns. Cross-category marketing variables are a
#'  subset of the marketing variables that are assumed to have a cross-series
#'  impact. The default for this is \code{NULL}.
#'@param sales a string. The name of the sales column. The corresponding column
#'  is expected to be non-negative.
#'@param lag numeric. Number of lags of the marketing variables to consider.
#'  The default is \code{2}
#'@param K numeric. Number of folds in CV. The default is \code{6}
#'@param alpha numeric vector of alpha values. These are used as alpha parameter
#'  of \code{\link{glmnet}} function employed in Stage 2. The best value is
#'  chosen with a CV. All values should be in [0,1]. The default is
#'  \code{c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)}
#'@param frequency numeric. Frequency of the timeseries. The default is
#'  \code{365.25 / 7}
#'@param horizon numeric. Leadtime of the forecast. The default is \code{1}
#'@param parallel Logical. Should parallel computing be used? The default is
#'  \code{F}
#'
#'@details This function fits a FAIR model to forecast store-category-level
#'  retail sales time series. The data is assumed to be of panel nature where
#'  the series in the same store affects each other.
#'@details All components of \code{train_data} is supplied in one
#'  data.frame. The user specifies the names of the individual components with
#'  the corresponding arguments.
#'@details  The series do not have to be of the same size and they may have gaps,
#'  i.e., the input data may lack rows for some
#'  \code{time_id}. The algorithm can work with missing values. However, if the proportion
#'  of the missing values in a series exceeds a certain threshold the series is
#'  ignored and the user is informed.
#'
#'@return The function outputs a list with the following objects:
#'  \item{fit}{A data.frame of the sales estimate and its components for each
#'  \code{store}-\code{category}-\code{time_id} combination in \code{train_data}.
#'  See the value section of \code{\link{FAIR_predict}} for details.}
#'  \item{models}{the model for each stage and a data.frame for in-sample bias}
#'  \item{pars}{parameters used for the training}
#'  \item{lag_data}{data.frame of the last periods, used by
#'  \code{FAIR_predict}}
#'
#'@author Özden Gür Ali and Ragıp Gürlek
#'
#'Maintainer: Ragıp Gürlek \email{rgurlek@@yahoo.com.tr}
#'
#' @references Gür Ali, Ö. and Gürlek, R. (2020) Automatic Interpretable Retail
#' Forecasting (FAIR) with Promotional Scenarios. International Journal of
#' Forecasting, forthcoming.
#'
#'@seealso \code{\link{FAIR_predict}}
#'
#'@export
FAIR_train <- function(train_data,
                       store,
                       category,
                       time_id,
                       seasonality,
                       marketing,
                       cc_marketing = NULL,
                       sales,
                       lag = 2,
                       K = 6,
                       alpha = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1),
                       frequency = 365.25 / 7,
                       horizon = 1,
                       parallel = F,
                       pool_seasonality = F,
                       NA_threshold = 0.3,
                       trim_Stage1 = T,
                       regularized_seasonality = F) {
  #drop unnecessary columns
  train_data <- train_data[, c(store, category, time_id, seasonality, marketing, sales)]
  temp <- train_data[, c(marketing, sales)]
  if (sum(!sapply(temp, is.numeric)) > 0){
    stop("Sales and marketing columns should be numeric!")
  }
  if (sum(temp < 0) > 0) {
    stop("Sales and marketing columns should be non-negative!")
  }
  last_week <- max(train_data[, time_id])
  t <- seq(to = last_week + 1 - horizon,
           by = horizon,
           length.out = K)
  #t is the vector of first time_ids of each block's validation set
  len_first_block <- t[1] - 1 - min(train_data[, time_id]) + 1
  if(len_first_block / frequency < 2){
    stop("FAIR requires the first CV block to cover at least two time-series cycles! Try a smaller 'K' or input a longer data.")
  }
  #append time_id to seasonality vector
  if(! time_id %in% seasonality) seasonality <- c(seasonality, time_id)

  category_list <- unique(train_data[, category])
  train_data[, category] <- as.factor(train_data[, category])

  par_mat <- data.frame()
  lambdas <- list()#used by step2 function
  for (t1 in t) {
    #block loop
    #deseasonalize data
    is_val <- t1 != utils::tail(t, 1)
    if(pool_seasonality){
      pool <- store
    } else {
      pool <- c(store,category)
    }
    if (is_val) {
      agg_fun <- plyr::ddply
    } else {
      agg_fun <- plyr::dlply
    }
    temp <- agg_fun(
      train_data,
      pool,
      deseason,
      t1 = t1,
      horizon = horizon,
      sales = sales,
      time_id = time_id,
      marketing = marketing,
      category = category,
      seasonality = seasonality,
      is_val = is_val,
      pool_seasonality = pool_seasonality,
      trim_model = trim_Stage1,
      parallel = parallel,
      .parallel = parallel,
      regularized_seasonality = regularized_seasonality
    )
    if(!is_val){
      s1_models <- extract(temp, 2)
      temp <- extract(temp, 1)
      temp <- do.call(rbind, temp)
    }
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
    #extended marketing vector
    ext_marketing <- marketing
    for (i in marketing) {
      for (j in 0:lag) {
        if (j == 0)
          next()
        ext_marketing <- append(ext_marketing, paste0(i, "_lag", j))
      }
    }
    for (i in cc_marketing) {
      for (j in category_list) {
        ext_marketing <- append(ext_marketing, paste0(j, "_" , i))
      }
    }

    #do step 2 for all categories and get par_mat
    par_mat <- rbind(
      par_mat,
      plyr::ddply(
        temp,
        category,
        step2,
        t1 = t1,
        alpha = alpha,
        marketing = ext_marketing,
        category = category,
        lag = lag,
        time_id = time_id,
        sales = sales,
        category_list = category_list,
        cc_marketing = cc_marketing,
        is_val = T,
        .parallel = parallel
      )
    )
  }

  best_pars <- function(my_data) {
    my_data$alpha <- as.numeric(as.character(my_data$alpha))
    my_data$lambda <- as.numeric(as.character(my_data$lambda))
    parameter_ev1 <- stats::aggregate(
      my_data$RMSE,
      by = list(my_data$alpha, my_data$lambda),
      FUN = "mean",
      na.rm = TRUE,
      drop = TRUE
    )
    names(parameter_ev1)[ncol(parameter_ev1)] <- "mean"
    parameter_ev2 <- stats::aggregate(
      my_data$RMSE,
      by = list(my_data$alpha, my_data$lambda),
      FUN = "sd",
      na.rm = TRUE,
      drop = TRUE
    )
    names(parameter_ev2)[ncol(parameter_ev2)] <- "sd"
    parameter_ev <- merge(parameter_ev1, parameter_ev2)
    colnames(parameter_ev)[1:2] <- c("alpha", "lambda")

    #choose the parameters with least mean+SD RMSE
    my_index <- which.min(parameter_ev$mean + parameter_ev$sd)
    par_mat <- data.frame()
    par_mat[1, 1] <- parameter_ev[my_index, "alpha"]
    par_mat[1, 2] <- parameter_ev[my_index, "lambda"]
    return(par_mat)
  }
  par_mat_best <- plyr::ddply(par_mat, category, best_pars)
  colnames(par_mat_best) <- c("category", "alpha", "lambda")

  my_index <-
    match(paste0(train_data[, store], train_data[, category], train_data[, time_id]),
          paste0(temp[, store], temp[, category], temp[, time_id]))
  prediction_1 <- (log(temp[my_index, "Original_Sales"] + 1) -
                          temp[my_index, sales])

  #step2
  temp <-
    plyr::dlply(
      temp,
      category,
      step2,
      t1 = last_week + 1,
      alpha = alpha,
      marketing = ext_marketing,
      category = category,
      lag = lag,
      time_id = time_id,
      sales = sales,
      category_list = category_list,
      cc_marketing = cc_marketing,
      is_val = F,
      par_mat_best = par_mat_best,
      .parallel = parallel
    )

  s2_models <- extract(temp, 2)
  temp <- extract(temp, 1)
  temp <- do.call(rbind, temp)

  #step3
  ext_marketing <- ext_marketing[!(ext_marketing %in% marketing)]
  temp <- temp[, !(colnames(temp) %in% ext_marketing)]
  series <- time_series(temp, store, time_id, category, NA_threshold)
  predictions_3 <- plyr::dlply(series,
                               c(store, category),
                               step3,
                               time_id = time_id,
                               frequency = frequency)

  s3_models <- extract(predictions_3, 2)
  predictions_3 <- extract(predictions_3, 1)
  predictions_3 <- do.call(rbind, predictions_3)
  match_vector <-
    match(
      paste0(temp[, store], temp[, category], temp[, time_id]),
      paste0(predictions_3[, store], predictions_3[, category],
             predictions_3[, time_id])
    )

  temp$prediction_3 <- predictions_3[match_vector, "STL"]
  temp$residual_3 <- temp$residual_2 - temp$prediction_3
  temp$Forecast <-
    exp(log(temp$Original_Sales + 1) - temp[, sales] +
          ifelse(is.na(temp$prediction_2), 0, temp$prediction_2) +
          ifelse(is.na(temp$prediction_3), 0, temp$prediction_3)) - 1

  #bias correction Snowdon, P. (1991). ratio estimator
  in_sample_bias <- plyr::ddply(temp, c(store, category),
                                function(x) {
                                  x <- x[!is.na(x$Forecast), ]
                                  return(sum(x$Original_Sales) /
                                           sum(x$Forecast))
                                })
  colnames(in_sample_bias)[3] <- "Bias_cor"

  #create outputs
  lag_data <- train_data[train_data[, time_id] > (last_week - lag), ]
  my_index <-
    match(paste0(train_data[, store], train_data[, category], train_data[, time_id]),
          paste0(temp[, store], temp[, category], temp[, time_id]))
  train_data <- data.frame(train_data[, c(store, category, time_id)],
                     prediction_1,
                     temp[my_index, c("prediction_2", "prediction_3")])
  colnames(train_data)[4:6] <-  c("Base_Sales",
                            "Marketing_deviation_multiplier",
                            "Disturbance_multiplier")
  train_data <-
    merge(train_data,
          in_sample_bias,
          by = c(store, category),
          all.x = T)
  train_data$Base_Sales <-
    exp(train_data$Base_Sales) * train_data$Bias_cor
  train_data$Marketing_deviation_multiplier <- exp(train_data$Marketing_deviation_multiplier)
  train_data$Disturbance_multiplier <- exp(train_data$Disturbance_multiplier)
  train_data$Sales_estimate <- (train_data$Base_Sales *
                    ifelse(is.na(train_data$Marketing_deviation_multiplier), 1,
                            train_data$Marketing_deviation_multiplier) *
                    ifelse(is.na(train_data$Disturbance_multiplier), 1,
                            train_data$Disturbance_multiplier) - 1)
  train_data <- train_data[, !colnames(train_data) %in% "Bias_cor"]

  pars <- list(
    store,
    category,
    time_id,
    seasonality,
    marketing,
    cc_marketing,
    sales,
    lag,
    frequency,
    horizon,
    last_week,
    ext_marketing,
    pool_seasonality,
    regularized_seasonality
  )
  names(pars) <- c(
    "store",
    "category",
    "time_id",
    "seasonality",
    "marketing",
    "cc_marketing",
    "sales",
    "lag",
    "frequency",
    "horizon",
    "last_week",
    "ext_marketing",
    "pool_seasonality",
    "regularized_seasonality"
  )


  return(list(
    fit = train_data,
    models = list(
      Stage1 = s1_models,
      Stage2 = s2_models,
      Stage3 = s3_models,
      in_sample_bias = in_sample_bias
    ),
    #CV = par_mat,
    pars = pars,
    lag_data = lag_data
  ))
}
