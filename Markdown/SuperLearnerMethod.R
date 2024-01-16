library(h2o)
library(ranger)
library(xgboost)
library(glmnet)

# h2o.init()
# Sys.setenv(https_proxy="") 
# Sys.setenv(http_proxy="") 
# Sys.setenv(http_proxy_user="") 
# Sys.setenv(https_proxy_user="")

#------------------------------------------------------------------------------#
# Running SuperLearner including training and CV
#------------------------------------------------------------------------------#
sl = function(input.data, target.feature, sl.library, cv.fold){
  # For reproducibility
  set.seed(20231207) 
  
  # split into X and Y variables
  traing.X = input.data[, setdiff(names(input.data), target.feature)]
  traing.Y = input.data[[target.feature]]
  
  # make sure feature name legit
  names(traing.X) = make.names(names(traing.X), unique = TRUE)
  
  # run SL 
  set.seed(20231207) 
  sl.model = SuperLearner(traing.Y, traing.X, family = gaussian(),
                          SL.library = sl.library, method = "method.NNLS")
  
  # Cross validation
  set.seed(20231207)
  sl.cv = CV.SuperLearner(X = traing.X,
                          Y = traing.Y,  
                          family = gaussian(),
                          cvControl = list(V = cv.fold), 
                          innerCvControl = list(list(V = cv.fold)),
                          SL.library = sl.library,
                          method = "method.NNLS")
  
  return(list(model = sl.model, cv = sl.cv, 
              traing.X = traing.X, traing.Y = traing.Y))
}

#------------------------------------------------------------------------------#
# Extract important features
#------------------------------------------------------------------------------#
extract_importance = function(sl){
  fit_sl = sl$model
  biggest_weight_idx = which.max(fit_sl$coef)
  fit_object = fit_sl$fitLibrary[[biggest_weight_idx]]$object
  
  imp_dt = data.frame()  # Initialize an empty dataframe
  
  if ("earth" %in% class(fit_object)) {
    earth_imp = earth::evimp(fit_object)
    imp_dt = data.frame(algo = "earth", Feature = rownames(earth_imp),
                        rank = earth_imp[,"rss"],
                        Importance = earth_imp[,"gcv"])
  }
  
  if ("ranger" %in% class(fit_object)) {
    # get importance of best ranger
    ranger_imp = ranger::importance(fit_object)
    ranger_imp_ranks = rank(-ranger_imp)
    imp_dt = data.frame(algo = "rf", Feature = names(ranger_imp),
                        rank = ranger_imp_ranks,
                        Importance = ranger_imp)
  }
  
  if (("xgboost" %in% class(fit_object)) | ("xgb.Booster" %in% class(fit_object))) {
    # get importance of best xgboost
    xgboost_imp = xgb.importance(model = fit_object)
    xgboost_imp$rank = seq_len(nrow(xgboost_imp))
    imp_dt = rbind(imp_dt, data.frame(algo = "xgboost", 
                                      Feature = xgboost_imp$Feature,
                                      rank = xgboost_imp$rank,
                                      Importance = xgboost_imp$Gain))
  }
  
  if (any(grepl("H2O", class(fit_object)))) {
    # get importance of best h2oboost
    h2oboost_imp_dt_init = fit_object@model$variable_importances
    h2oboost_imp_dt_init$rank = seq_len(nrow(h2oboost_imp_dt_init))
    imp_dt = data.frame(algo = "h2oboost", 
                        Feature = h2oboost_imp_dt_init$variable, 
                        rank = h2oboost_imp_dt_init$rank, 
                        Importance = h2oboost_imp_dt_init$relative_importance)
  }
  
  if ("cv.glmnet" %in% class(fit_object)) {
    # get importance of best glmnet
    glmnet_coef = fit_object$glmnet.fit$beta[, which(fit_object$lambda == fit_object$lambda.min)]
    glmnet_imp_rank = rank(-abs(glmnet_coef))
    imp_dt = data.frame(algo = "lasso", 
                        Feature = names(glmnet_coef),
                        rank = glmnet_imp_rank,
                        Importance = glmnet_coef)
  }
  
  if ("numeric" %in% class(fit_object)) {
    imp_dt = data.frame(algo = "mean", 
                        Feature = fit_sl$varNames, 
                        rank = mean(1:length(fit_sl$varNames)), 
                        Importance = 0)
  }
  
  row.names(imp_dt) = NULL
  return(imp_dt[order(imp_dt$rank), ])
}

#------------------------------------------------------------------------------#
# H2O algorithm
#------------------------------------------------------------------------------#
SL.h2oboost <- function(Y, X, newX, family, obsWeights = rep(1, length(Y)), ...)
  
{
  SuperLearner:::.SL.require("h2o")
  
  # Set GBM parameters to test
  hyper_parameters <- list(ntrees = list(1000),
                           max_depth = list(2,4,5,6),
                           learn_rate = list(0.05, 0.1, 0.2),
                           col_sample_rate = list(0.1, 0.2, 0.3))
  
  # Bind vector of outcome and covariate matrix together;
  # if Y is binary, make it a factor first
  Y <- switch((family$family == "binomial") + 1, Y, as.factor(Y))
  dat <- cbind(Y, X)
  
  # Convert dat to an h2o object
  dat.hex <- as.h2o(dat)
  
  # set up GBM hyperparameters
  if (family$family == "binomial") {
    h2o_dist <- "bernoulli"
    h2o_metric <- "AUC"
  } else if (family$family == "gaussian") {
    h2o_dist <- "gaussian"
    h2o_metric <- "MSE"
  } else {
    stop("The entered family isn't currently supported. Please enter one of 'binomial' or 'gaussian'.")
  }
  
  # search over the grid
  gbm.model <- h2o::h2o.grid("gbm",
                             hyper_params = hyper_parameters,
                             training_frame = dat.hex,
                             y = "Y",
                             distribution = h2o_dist,
                             nfolds = 5,
                             balance_classes = (h2o_dist == "bernoulli"),
                             max_after_balance_size = 5,
                             fold_assignment = ifelse(h2o_dist == "bernoulli",
                                                      "Stratified", "AUTO"),
                             stopping_metric = toupper(h2o_metric),
                             stopping_rounds = 3,
                             stopping_tolerance = 0.001,
                             max_runtime_secs = 60,
                             parallelism = 0)
  
  # get the models from the grid and sort by metric
  grid <- h2o::h2o.getGrid(gbm.model@grid_id,
                           sort_by = tolower(h2o_metric),
                           decreasing = (h2o_dist == "bernoulli"))
  
  # Save best parameters
  best.max_depth <- as.numeric(grid@summary_table[1, ]$max_depth)
  best.learn_rate <- as.numeric(grid@summary_table[1, ]$learn_rate)
  best.col_sample_rate <- as.numeric(grid@summary_table[1, ]$col_sample_rate)
  
  # Remove all models in grid to save memory
  h2o.removeAll(retained_elements = c(dat.hex))
  rm(gbm.model, grid)
  
  # Call garbage collection
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  h2o:::.h2o.garbageCollect()
  
  # Train the model with best hyperparameters
  gbm.final.model <- h2o::h2o.gbm(training_frame = dat.hex,
                                  y = "Y",
                                  distribution = h2o_dist,
                                  stopping_metric = toupper(h2o_metric),
                                  stopping_rounds = 3,
                                  stopping_tolerance = 0.001,
                                  ntrees = 1000,
                                  max_depth = best.max_depth,
                                  learn_rate = best.learn_rate,
                                  col_sample_rate = best.col_sample_rate)
  
  # Convert newdata to h2o object
  newX.hex <- as.h2o(newX)
  # Get predictions
  pred.raw <- h2o::h2o.predict(object = gbm.final.model,
                               newdata = newX.hex)
  
  # Extract predicted probabilities
  if(family$family == "gaussian"){
    pred <- as.numeric(as.vector(pred.raw))
  }
  else if(family$family == "binomial"){
    pred <- as.numeric(as.vector(pred.raw$p1))
  }
  # Make fit an object with the filepath we need to reload the h2o object
  fit <- list(object = gbm.final.model)
  class(fit) <- c("SL.h2oboost")
  out = list(pred = pred, fit = fit)
  
  return(out)
}

# Predict method for h2oboost
predict.SL.h2oboost <- function(object, newdata, ...)
{
  SuperLearner:::.SL.require("h2o")
  L <- list(...)
  # convert data to h2o object
  newdata.hex <- h2o::as.h2o(newdata)
  # Get predictions
  pred.raw <- h2o::h2o.predict(object = object$object,
                               newdata = newdata.hex)
  # Extract predicted probabilites
  if (L$family$family == "gaussian"){
    pred <- as.numeric(as.vector(pred.raw))
  }
  else if (L$family$family == "binomial"){
    pred <- as.numeric(as.vector(pred.raw$p1))
  }
  pred
}
descr_SL.h2oboost <- "boosted regression trees with (maximum depth, learning rate, column sampling rate) selected by 5-fold CV over the grid $(2, 4, 5, 6)\\times(.05, .1, .2)\\times(.1, .2, .3)$"

#------------------------------------------------------------------------------#
# xgboost algorithm
#------------------------------------------------------------------------------#
SL.xgboost.corrected <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, ntrees = 1000,
                                  max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(),
                                  nthread = 1, verbose = 0, save_period = NULL, ...)
{
  SuperLearner:::.SL.require("xgboost")
  if (packageVersion("xgboost") < 0.6)
    stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
  if (!is.matrix(X)) {
    X = model.matrix(~. - 1, X)
  }
  xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
  if (family$family == "gaussian") {
    model = xgboost::xgboost(data = xgmat, objective = "reg:squarederror",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  if (family$family == "binomial") {
    model = xgboost::xgboost(data = xgmat, objective = "binary:logistic",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, eval_metric = "logloss",
                             eta = shrinkage, verbose = verbose, nthread = nthread,
                             params = params, save_period = save_period)
  }
  if (family$family == "multinomial") {
    model = xgboost::xgboost(data = xgmat, objective = "multi:softmax",
                             nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode,
                             eta = shrinkage, verbose = verbose, num_class = length(unique(Y)),
                             nthread = nthread, params = params, save_period = save_period)
  }
  if (!is.matrix(newX)) {
    newX = model.matrix(~. - 1, newX)
  }
  pred = predict(model, newdata = newX)
  fit = list(object = model)
  class(fit) = c("SL.xgboost")
  out = list(pred = pred, fit = fit)
  return(out)
}

SL.xgboost.2 <- function(..., max_depth = 2){
  SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.4 <- function(..., max_depth = 4){
  SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.6 <- function(..., max_depth = 6){
  SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.8 <- function(..., max_depth = 8){
  SL.xgboost.corrected(..., max_depth = max_depth)
}
SL.xgboost.12 <- function(..., max_depth = 12){
  SL.xgboost.corrected(..., max_depth = max_depth)
}
descr_SL.xgboost <- "boosted regression trees with maximum depth of "
descr_SL.xgboost.2 <- paste0(descr_SL.xgboost, 2)
descr_SL.xgboost.4 <- paste0(descr_SL.xgboost, 4)
descr_SL.xgboost.6 <- paste0(descr_SL.xgboost, 6)
descr_SL.xgboost.8 <- paste0(descr_SL.xgboost, 8)
descr_SL.xgboost.12 <- paste0(descr_SL.xgboost, 12)

#------------------------------------------------------------------------------#
# random forests algorithm
#------------------------------------------------------------------------------#
SL.ranger.imp <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), num.trees = 500, mtry = floor(sqrt(ncol(X))),
                           write.forest = TRUE, probability = family$family == "binomial",
                           min.node.size = ifelse(family$family == "gaussian", 5, 1),
                           replace = TRUE, sample.fraction = ifelse(replace, 1, 0.632),
                           num.threads = 1, verbose = TRUE, ...) {
  SuperLearner:::.SL.require("ranger")
  if (family$family == "binomial") {
    Y = as.factor(Y)
  }
  if (is.matrix(X)) {
    X = data.frame(X)
  }
  fit <- ranger::ranger(`_Y` ~ ., data = cbind(`_Y` = Y, X),
                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                        replace = replace, sample.fraction = sample.fraction,
                        case.weights = obsWeights, write.forest = write.forest,
                        probability = probability, num.threads = num.threads,
                        verbose = verbose, importance = "impurity")
  pred <- predict(fit, data = newX)$predictions
  if (family$family == "binomial") {
    pred = pred[, "1"]
  }
  fit <- list(object = fit, verbose = verbose)
  class(fit) <- c("SL.ranger")
  out <- list(pred = pred, fit = fit)
  return(out)
}
SL.ranger.reg <- function(..., X, mtry = floor(sqrt(ncol(X)))){
  SL.ranger.imp(..., X = X, mtry = mtry)
}

SL.ranger.small <- function(..., X, mtry = floor(sqrt(ncol(X)) * 1/2)){
  SL.ranger.imp(..., X  = X, mtry = mtry)
}

SL.ranger.large <- function(..., X, mtry = floor(sqrt(ncol(X)) * 2)){
  SL.ranger.imp(..., X = X, mtry = mtry)
}
descr_SL.ranger.imp <- "random forest with `mtry` equal to "
descr_SL.ranger.reg <- paste0(descr_SL.ranger.imp, "square root of number of predictors")
descr_SL.ranger.small <- paste0(descr_SL.ranger.imp, "one-half times square root of number of predictors")
descr_SL.ranger.large <- paste0(descr_SL.ranger.imp, "two times square root of number of predictors")

#------------------------------------------------------------------------------#
# GLMnet algorithm
#------------------------------------------------------------------------------#
SL.glmnet.0 <- function (Y, X, newX, family, obsWeights = rep(1, length(Y)), id, alpha = 0, nfolds = 5,
                         nlambda = 100, useMin = TRUE, loss = "deviance", ...) {
  SuperLearner:::.SL.require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + ., X)
    newX <- model.matrix(~-1 + ., newX)
  }
  
  if(family$family == "binomial"){
    fold_id <- get_fold_id(Y)
    nfolds <- max(fold_id)
    if(nfolds != 0){
      fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                                 lambda = NULL, type.measure = loss, nfolds = nfolds,
                                 foldid = fold_id, family = family$family, alpha = alpha, nlambda = nlambda,
                                 ...)
      pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                        "lambda.min", "lambda.1se"))
      fit <- list(object = fitCV, useMin = useMin)
      class(fit) <- "SL.glmnet"
    }else{
      # if fewer than 3 cases, just use mean
      meanY <- weighted.mean(Y, w = obsWeights)
      pred <- rep.int(meanY, times = nrow(newX))
      fit <- list(object = meanY)
      out <- list(pred = pred, fit = fit)
      class(fit) <- c("SL.mean")
    }
  }else{
    fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                               lambda = NULL, type.measure = loss, nfolds = nfolds,
                               family = family$family, alpha = alpha, nlambda = nlambda,
                               ...)
    pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                      "lambda.min", "lambda.1se"))
    fit <- list(object = fitCV, useMin = useMin)
    class(fit) <- "SL.glmnet"
  }
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.glmnet.ridge <- function(..., alpha = 0){
  SL.glmnet.0(..., alpha = alpha)
}

SL.glmnet.25 <- function(..., alpha = 0.25){
  SL.glmnet.0(..., alpha = alpha)
}

SL.glmnet.50 <- function(..., alpha = 0.5){
  SL.glmnet.0(..., alpha = alpha)
}

SL.glmnet.75 <- function(..., alpha = 0.75){
  SL.glmnet.0(..., alpha = alpha)
}

SL.glmnet.lasso <- function(..., alpha = 1){
  SL.glmnet.0(..., alpha = alpha)
}

descr_SL.glmnet <- "elastic net with $\\lambda$ selected by 5-fold CV and $\\alpha$ equal to "
descr_SL.glmnet.50 <- paste0(descr_SL.glmnet, "0.5")
descr_SL.glmnet.25 <- paste0(descr_SL.glmnet, "0.25")
descr_SL.glmnet.75 <- paste0(descr_SL.glmnet, "0.75")
descr_SL.glmnet.0 <- paste0(descr_SL.glmnet, "1 (i.e., the lasso)")

descr_SL.mean <- "intercept only regression"
descr_SL.glm <- "main terms generalized linear model"

