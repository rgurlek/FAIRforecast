#' FAIR: A package for FAIR forecast
#'
#' The main functionality of this package is to implement the FAIR
#' method proposed by Gür Ali and Gürlek (2019) to any panel data of
#' store-category-level sales time series of a retail. It utilizes
#' orthogonalization, regularization via Elastic net, and extrapolation
#' to account for seasonality, focal- and cross-category marketing activity
#' as well as random disturbances to base sales. It can handle large
#' dimensionality and deal with regularization induced confounding.
#'
#' @details Package:  FAIR
#' @details Type:     Package
#' @details Version:  0.1.0
#' @details Release Date: 09-19-19
#' @details License:  GPL-2
#'
#'@author Özden Gür Ali and Ragıp Gürlek
#'
#'Maintainer: Ragıp Gürlek \email{rgurlek@@yahoo.com.tr}
#'
#' @references Gür Ali, Ö. and Gürlek, R. (2019) Automatic Interpretable Retail
#' Forecasting with Promotional Scenarios
#' @references Gür Ali, Ö. and Pınar, E. (2016). Multi-period-ahead forecasting
#'   with residual extrapolation and information sharing - Utilizing a multitude
#'   of retail series. International Journal of Forecasting, 32(2):502–517.
#' @references Hahn, P. R., Carvalho, C. M., Puelz, D., and He, J. (2018).
#' Regularization and confounding in linear regression for treatment effect
#' estimation. Bayesian Analysis, 13(1):163– 182.
#' @references Zou, H. and Hastie, T. (2005), Regularization and variable
#' selection via the elastic net. Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology), 67: 301-320.
#'
#' @docType package
#' @name FAIR
NULL

#' A sample retail sales dataset
#'
#' A dataset containing the weekly category-store-level sales of a retail chain.
#' It provides the promotional activity as well as the seasoanlity and calendar
#' variables.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#' \item{price}{price, in US dollars}
#' \item{carat}{weight of the diamond, in carats}
#'
#' }
#' Note that distributional variables do not add up to 1 because they exclude
#' non-discounted SKUs.
"sample_data"

extract <-
  function(my_list, index) {
    # when a list of lists created, this function
    #creates a new list with a specific element of sublists
    return(plyr::llply(my_list, function(my_element)
      my_element[[index]]))
  }

deseason <- function(my_data,
                     t1,
                     horizon,
                     sales,
                     time_id,
                     marketing,
                     store,
                     seasonality,
                     is_val,
                     parallel) {
  #used in the cross validation. Does the first step.
  #t1 is the first time_id of the validation set.
  s1_models <- list()
  my_data$Original_Sales <- my_data[, sales]
  my_data <- my_data[my_data[, time_id] < t1 + horizon, ]
  var_list <- c(sales, marketing)
  fit <- function(my_var) {
    train <- my_data[my_data[, time_id] < t1, c(my_var, seasonality)]
    if (sd(train[, my_var]) == 0){
      s1_models[[my_var]] <<- train[1, my_var]
      return(my_data[, my_var])
      #If there is no variance in the training data, return
      #training + validation as it is
    } else{
      myformula <- as.formula(paste0("log(", my_var, "+1) ~ ."))
      my_model <- lm(formula = myformula,
                     data = train,
                     model = F)
      s1_models[[my_var]] <<- my_model
      #above line ensures obs until t1+horizon-1 are predicted.
      #If the data is not validation (training), it is same as t1-1 since
      #all points are less than t1
      return(log(my_data[, my_var] + 1) - predict.lm(my_model, my_data))
    }
  }
  new_vars <- lapply(var_list, fit)
  new_vars <- do.call(data.frame, new_vars)
  my_data[, var_list] <- new_vars
  if (is_val)
    return(my_data)
  my_list <- list(my_data, s1_models)
  return(my_list)
}
#'@export
var_cre <- function(my_data,
                    category,
                    category_list,
                    lag,
                    time_id,
                    marketing,
                    cc_marketing,
                    store,
                    parallel) {
  #create stage 2 variables (marketing of the other categories and lag of
  #focal category)
  lag_func <- function(temp, l) {
    my_week <- temp[, time_id]
    var_loop <- function(my_var){
      a <- temp[match(my_week - l, temp[, time_id], nomatch = NA), my_var]
    }
    new_vars <- lapply(marketing, var_loop)
    new_vars <- do.call(data.frame, new_vars)
    colnames(new_vars) <- paste0(marketing, "_lag", l)
    return(cbind(temp, new_vars))
  }
  for (i in 0:lag) {
    if (i == 0)
      next()
    my_data <-
      plyr::ddply(my_data,
                  c(store, category),
                  lag_func,
                  l = i,
                  .parallel = parallel)
  }

  cross_fun <- function(temp) {
    for (i in cc_marketing) {
      for (k in category_list) {
        a <- temp[match(k, temp[, category], nomatch = NA), i]
        temp[, paste0(k , "_", i)] <- a
      }
    }
    return(temp)
  }
  my_data <-
    plyr::ddply(my_data, c(store, time_id), cross_fun, .parallel = parallel)
  return(my_data)
}

step2 <-
  function(my_data,
           t1,
           alpha,
           marketing,
           category,
           lag,
           time_id,
           sales,
           category_list,
           cc_marketing,
           is_val,
           par_mat_best = NULL) {
    my_category <- my_data[1, category]

    for (i in cc_marketing) {
      #exclude focal marketing dublicates
      marketing <- marketing[marketing != paste0(my_category, "_", i)]
    }

    train <- my_data[my_data[, time_id] < t1,]
    train <- train[complete.cases(train[, marketing]),]

    if (is_val) {
      test <- my_data[!my_data[, time_id] < t1,]

      par_ev <- data.frame()
      for (i in alpha) {
        nam <- paste0("lambda_", my_category)
        lambdas <- get("lambdas", envir = parent.frame(5))
        if (t1 == get("t", envir = parent.frame(5))[1]) {
          my_model <-
            glmnet::glmnet(
              as.matrix(train[, marketing]),
              train[, sales],
              family = "gaussian",
              alpha = i,
              nlambda = 50
            )
          lambdas[[paste0(my_category, "_", i)]] <- my_model$lambda
          assign("lambdas", lambdas, envir = parent.frame(5))
        } else {
          my_lambda <- lambdas[[paste0(my_category, "_", i)]]
          my_model <-
            glmnet::glmnet(
              as.matrix(train[, marketing]),
              train[, sales],
              family = "gaussian",
              alpha = i,
              lambda = my_lambda
            )
        }
        predictions <-
          glmnet::predict.glmnet(my_model, as.matrix(test[, marketing]))
        predictions <- exp(log(test$Original_Sales + 1) - test[, sales]
                           + predictions) - 1
        sq_errors <- (test$Original_Sales - predictions) ^ 2
        rmse <- sqrt(colMeans(sq_errors, na.rm = TRUE))
        temp <- data.frame(rep(t1, length(rmse)),
                           rep(i, length(rmse)),
                           as.factor(my_model$lambda),
                           rmse)
        colnames(temp) <- c("t", "alpha", "lambda", "RMSE")
        par_ev <- rbind(par_ev, temp)
      }
      return(par_ev)
    } else{
      pars <- par_mat_best[par_mat_best$category == my_category,]
      a <- as.numeric(as.character(pars$alpha[1]))
      l <- as.numeric(as.character(pars$lambda[1]))

      my_model <-
        glmnet::glmnet(
          as.matrix(train[, marketing]),
          train[, sales],
          family = "gaussian",
          alpha = a,
          lambda = l
        )
      train$prediction_2 <-
        as.numeric(glmnet::predict.glmnet(
          my_model, as.matrix(train[, marketing])))
      train$residual_2 <- train[, sales] - train$prediction_2
      return(list(train, my_model))
    }
  }

time_series <- function(my_data, store, time_id, category) {
  sum_func <- function(x) if (all(is.na(x))) NA_integer_  else sum(x)
  series <- aggregate(
    my_data$residual_2,
    by = list(my_data[, store], my_data[, time_id], my_data[, category]),
    FUN = sum_func,
    drop = FALSE
  )
  colnames(series) <- c(store, time_id, category, "residual_2")

  f_week <- min(series[, time_id])
  l_week <- max(series[, time_id])
  total_week <- (l_week - f_week + 1)
  mis_fun <- function(x) {
    total_mis <- sum(is.na(x$residual_2))
    is_valid <- (total_mis / total_week) <= 0.3
    if (is_valid) {
      return(x)
    } else {
      warning(
        "Series for store ", x[1, store], " and category ", x[1, category],
        " has missing values that are more than 30% of the lenth of the series.",
        " Thus, it is excluded and its forecasts are set to NA."
      )
      return(NULL)
    }
  }

  series <- plyr::ddply(series, c(store, category), mis_fun)

  return(series)
}

step3 <- function(x, time_id, frequency) {
  x <- x[order(x[, time_id]),]
  my_series <- ts(x$residual_2, frequency = frequency)
  my_model <- forecast::stlm(forecast::na.interp(my_series))
  x$STL <- as.numeric(my_model$fitted)
  return(list(x, my_model))
}

#'@export
variable_importance <- function(object, n_vars = 10, categories = NULL,
                                only_own = T){
  if(is.null(categories)){
    models <- object$models$Stage2
  } else {
    models <- object$models$Stage2[categories]
  }
  models <- lapply(models, function(category){
    abs(as.matrix(category$beta))
  })
  vars <- Reduce('+', models) / length(models)
  vars <- data.frame(variable = row.names(vars),
                     magnitude = as.numeric(vars[, 1]))
  if(only_own){
    cross <- expand.grid(names(object$models$Stage2), object$pars$cc_marketing)
    cross <- paste(cross$Var1, cross$Var2, sep = "_")
    vars <- vars[!vars$variable %in% cross, ]
  }
  vars <- vars[order(vars[, "magnitude"], decreasing = T), ]
  n_vars <- min(n_vars, nrow(vars))
  #vars matrix has the own-category marketing variables sorted by their
  #average absolute coefficients in Stage 2 calculated over categories

  p <- plotly::layout(plotly::plot_ly(
    x = factor(vars[1:n_vars, "variable"],
               levels = vars[1:n_vars, "variable"]),
    y = vars[1:n_vars, "magnitude"],
    type = "bar"),
    title = "Variable Importance")
  return(p)
}

#'@export
cross_category <- function(object){
  models <- object$models$Stage2
  category_list <- names(object$models$Stage2)
  if(is.null(object$pars$cc_marketing))
    stop("The model is not trained with cross-category effects")
  var_list <- object$pars$cc_marketing
  cc_list <- expand.grid(category_list, var_list)
  cc_list <- paste(cc_list[, 1], cc_list[, 2], sep = '_')
  remove(object)
  cross_mat <- matrix(NA, nrow = length(category_list), ncol = length(category_list),
                      dimnames = list(category_list, category_list))
  for(i in 1:length(models)){
    my_mat <- abs(as.matrix(models[[i]]$beta))
    my_mat <- data.frame(variable = row.names(my_mat),
                         magnitude = as.numeric(my_mat[, 1]))
    my_mat <- my_mat[my_mat[, 'variable'] %in% cc_list, ]
    my_mat$category <- substr(my_mat$variable, 1,
                              regexpr("_", my_mat$variable) - 1)
    my_mat <- aggregate(my_mat$magnitude, list(my_mat$category), FUN = sum)
    colnames(my_mat) <- c('Category', 'magnitude')
    for(j in 1:nrow(my_mat)){
      cross_mat[names(models)[i], my_mat[j, 'Category']] <-
        my_mat[j, 'magnitude']
    }
  }
  cross_mat <- cross_mat / length(var_list)
  cross_mat <- ifelse(is.na(cross_mat), 0, cross_mat)
  p <- plotly::layout(
    plotly::plot_ly(x = colnames(cross_mat), y = rownames(cross_mat),
                    z = cross_mat,
                    colorscale = scales::col_numeric("Blues", domain = NULL)(unique(scales::rescale(c(volcano)))),
                    type = "heatmap"),
    title = "Cross-category effects",
    xaxis = list(title = "Influencial Category"),
    yaxis = list(title = "Influenced Category")
  )
  return(p)
}
