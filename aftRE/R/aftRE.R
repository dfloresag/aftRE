TMB::compile("./src/aft_n_ub.cpp")
dyn.load(TMB::dynlib("./src/aft_n_ub"))

#' @title report_estimates:
#'
#' @param obj vector of observed durations
#' @param group_1_labels matrix containing the values of the covariates
#' @param group_2_labels vector of censoring indicators
#' @examples
#' #train the multiDA classifier using the SRBCT dataset, and find the resubstitution error rate
#'
#' y   <- SRBCT$y
#' X   <- SRBCT$X
#' res  <- multiDA(X, y, equal.var=TRUE, set.options="exhaustive", penalty="EBIC")
#' vals <- predict(res, newdata=X)$y.pred          #y.pred returns class labels
#' rser <- sum(vals!=y)/length(y)

#' @rdname aftRE
#' @export
#' @importFrom TMB sdreport
#' @importFrom TMB summary.sdreport
report_estimates <- function(obj, covariate_names){
  library(TMB)

  fixed <- try(unname(
    TMB::summary.sdreport(
      TMB::sdreport(
        obj = obj
        ),
      select = "fixed")
    ),
    silent = TRUE)

  report <- try((
    TMB::summary.sdreport(
      TMB::sdreport(
        obj = obj
        ),
      select = "report")),
    silent = TRUE)

  col_names  <- c(
    covariate_names,
    paste("log_sigma",0:2, sep="_"),
    paste("sigma",0:2, sep="_"),
    paste("sigma2",0:2, sep="_")
    )

  data.frame(
    Parameter = col_names,
    Estimates = c(fixed[,1],report[,1]),
    Std_Errors = c(fixed[,2],report[,2])
  )
}

#' @title report_predictions:
#'
#' @param obj vector of observed durations
#' @param group_1_labels matrix containing the values of the covariates
#' @param group_2_labels vector of censoring indicators
#' @examples
#'
#'
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom TMB sdreport
#' @importFrom TMB summary.sdreport
report_predictions <- function(obj,
                               group_1_labels,
                               group_2_labels){

  random <- try(unname(
    TMB::summary.sdreport(
      TMB::sdreport(
        obj = obj
        ),
      select = "random")),
    silent = TRUE)

  col_names  <- c(
    group_1_labels,
    group_2_labels
    )

  data.frame(
    Random_Effect = col_names,
    Prediction = random[,1],
    Conditonal_Variance = random[,2]
  )
}

#' @title report_predictions: Accelerated Failure Time Models with Random Effects
#'
#' @param obj vector of observed durations
#' @param group_1_labels matrix containing the values of the covariates
#' @param group_2_labels vector of censoring indicators
#' @examples
#'
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom TMB sdreport
#' @importFrom TMB summary.sdreport
select_column <- function(x, which_col){
  x[,which_col]
}

#' @title report_predictions: Accelerated Failure Time Models with Random Effects
#'
#' @param durations_t,
#' @param covariates_X,
#' @param censoring_ind,
#' @param grouping_factors, # list or matrix
#' @param conditional_distribution = c("log-normal", "log-logistic", "weibull"),
#' @param weights = NULL,
#' @param optimizer = "nlminb",
#' @param ...
#' @examples
#' #train the multiDA classifier using the SRBCT dataset, and find the resubstitution error rate
#'
#' y   <- SRBCT$y
#' X   <- SRBCT$X
#' res  <- multiDA(X, y, equal.var=TRUE, set.options="exhaustive", penalty="EBIC")
#' vals <- predict(res, newdata=X)$y.pred          #y.pred returns class labels
#' rser <- sum(vals!=y)/length(y)

#' @rdname multiDA
#' @export
#' @importFrom TMB compile
#' @importFrom TMB dynlib
#' @importFrom TMB MakeADFun
estimate_aft_nested <- function(
  durations_t,
  covariates_X,
  censoring_ind,
  grouping_factors, # list or matrix
  conditional_distribution = c(
    "log-normal",
    "log-logistic",
    "weibull"),
  weights = NULL,
  optimizer = "nlminb",
  ...
){

  # TODO : Verifying the inputs, coercing if necessary

  group_1_expanded <- grouping_factors[[1]] %>% as.factor()
  group_2_expanded <- grouping_factors[[2]] %>% as.factor()

  group_1_labels <- as.character(unique(group_1_expanded))
  group_2_labels <- as.character(unique(group_2_expanded))

  covariate_names  <- colnames(covariates_X)

  group_12  <- group_2_expanded %>%
    unique() %>% as.character() %>%
    strsplit(split = "\\.") %>%
    unlist() %>%
    matrix(ncol = ncol(grouping_factors), byrow = TRUE)%>%
    select_column(which_col=1) %>%
    as.factor()

  if(is.null(weights)){
    w = rep(1, nlevels(group_1_expanded))
  } else {
    w = weights
  }

  conditional_distribution <- match.arg(conditional_distribution)

  if(conditional_distribution == "log-normal"){
    cond_resp <- 1
  } else if(conditional_distribution == "log-logistic"){
    cond_resp <- 2
  } else if(conditional_distribution == "weibull") {
    cond_resp <- 3
  }

  dta_n   <- list(
    group_1  = group_1_expanded ,
    group_2  = group_2_expanded ,
    group_12 = group_12,
    X = covariates_X,
    w = w,
    t = durations_t,
    d = censoring_ind,
    cond_resp = cond_resp
  )

  par_n <- list(
    u1    = rep(0, nlevels(group_1_expanded)),
    u2    = rep(0, nlevels(group_2_expanded)),
    betas = rep(1, ncol(covariates_X)),
    log_sigma_0 = log(1),
    log_sigma_1 = log(1),
    log_sigma_2 = log(1)
  )

  obj  <- TMB::MakeADFun(
      data=dta_n,
      parameters=par_n,
      random=c("u1","u2"),
      type = c("ADFun", "ADGrad" , "Fun"),
      DLL= "aft_n_ub", silent = TRUE,
      hessian = FALSE, method = "CG",
      inner.control = list(maxit = 1000),
      control = list(maxit =1000)
    )

  if(optimizer == "nlminb"){
    opt <- try(
      nlminb(obj$par, obj$fn, obj$gr, control = obj$control, ...),
      silent = TRUE
    )
  } else if(optimizer == "optim") {
    opt <- try(
      do.call(optim, obj),
      silent = TRUE
    )
  }

  # TODO : Error Control

  if(inherits(opt, "try-error")){
    warning("estimation error: please review your initial values")
    out <- opt

  } else {
    if (opt$convergence == 1) message("issue: convergence")

    # TODO : Include elements of the frame as another argument

    out <- list(
      estimates      =
        report_estimates(
          obj,
          covariate_names),
      predictions    =
        report_predictions(
          obj,
          group_1_labels,
          group_2_labels
        ),
      distribution   = conditional_distribution,
      log_likelihood = -opt$objective,
      convergence    = opt$convergence, # TODO include option when the optimizer is `optim`
      frame = list(
        t = durations_t,
        X = covariates_X,
        d = censoring_ind,
        w = w,
        group_1 = group_1_expanded,
        group_2 = group_2_expanded
      )
    )
    class(out) = "aftRE"
  }
  invisible(out)
}

#' @title print.aftRE
#'
print.aftRE <- function(object, ...){
  # TODO: Call
  # cat("Call: \n")
  # print(x$call)
  cat("distribution: \n")
  print(object$distribution)
  # cat("method: \n")
  # print(x$method)
  cat("\n")
  cat("log-likelihood: ", object$log_likelihood, "\n")
  # if(!is.null(x$params$inv.phi)){ x$params$inv.phi <- NULL; }
  # crit <- inf.criteria(x)
  # df <- crit$k
  # cat("Residual degrees of freedom: ", length(x$y) - df, "\n")
  # cat("AIC: ", crit$AIC, "\n")
  # cat("AICc: ", crit$AICc, "\n")
  # cat("BIC: ", crit$BIC, "\n")
  # invisible(crit)
}

#' @title simulate_aft_nested: simulate a
#'
#' @param parameter_values
#' @param dimensions
#' @param censoring_scenario
#' @param conditional_distribution
#' @param max_trials
#' @param new_X = NULL,
#' @examples
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom TMB sdreport
#' @importFrom TMB summary.sdreport
#' @importFrom SpatialExtremes rgev
simulate_aft_nested <- function(
  parameter_values  = list(
    betas   = c(1, 1, -0.5, 1),
    sigma_0 = 1.21,
    sigma_1 = c(0.2),
    sigma_2 = c(0.2)
  ),
  dimensions= list(
    n_1 = 20,
    n_12 = 20,
    n_23 = 15
  ),
  censoring_scenario = 1,
  conditional_distribution =c(
    "log-normal",
    "log-logistic",
    "weibull"
  ),
  max_trials = 50,
  old_X = NULL,
  ...){

  # Beta
  betas   = parameter_values$betas
  sigma_0 = parameter_values$sigma_0
  sigma_1 = parameter_values$sigma_1
  sigma_2 = parameter_values$sigma_2

  n_1     = dimensions$n_1
  n_12    = dimensions$n_12
  n_23    = dimensions$n_23

  n = n_1*n_12*n_23

  group_1 <- gl(n_1, n_12*n_23)
  group_2 <- factor(paste(n_1,rep(gl(n_12, n_23), times = n_1), sep = "."))

  n_2 <- nlevels(group_2)

  m <- 1

  convergence_error <- TRUE
  estimation_error  <- TRUE

  while(m <= max_trials & convergence_error & estimation_error){
    if (is.null(old_X)){
      X  <-  cbind(
        intercept = rep(1, times = n),
        X_1_cont  = rnorm(n = n , mean = 0 , sd = 1),
        X_2_cont  = rnorm(n = n , mean = 0 , sd = 1),
        X_3_dich  = rbinom(n = n, size = 1 , prob = 0.5)
      ) %>%  as.matrix()
    } else {
      X <- old_X
    }

    u1  <- rnorm(n_1, mean = 0, sd = sigma_1)
    u2  <- rnorm(n_2, mean = 0, sd = sigma_2)

    e   <- as.vector(X%*%betas + u1[group_1] + u2[group_2])

    if(conditional_distribution == "log-normal"){
      t   <- exp(e + sigma_0*rnorm(n = n))
    } else if (conditional_distribution == "log-logistic"){
      t   <- exp(e + sigma_0*rlogis(n = n))
    } else if (conditional_distribution == "weibull"){
      t   <- exp(e - sigma_0*evd::rgumbel( n = n))
    }

    # TODO : A more precise censoring scheme

    if(censoring_scenario==0){
      cn  <- rep(1, n)
      dl  <- rep(0, n)
    } else {
      if (censoring_scenario==1){
        dl <- rexp(n = n, rate = 1/quantile(t, .975))
      } else if (censoring_scenario==2){
        dl <- rexp(n = n, rate = 1/quantile(t, .945))
      } else if(censoring_scenario==3){
        dl <- rexp(n = n, rate = 1/quantile(t, .905))
      } else if(censoring_scenario==4){
        dl <- rexp(n = n, rate = 1/quantile(t, .65))
      }
      cn   <- 1*(t<dl)
    }

    t <- (t*cn)+ (dl*(1-cn))

    out <- estimate_aft_nested(
      durations_t  = t,
      covariates_X = X,
      censoring_ind = cn,
      grouping_factors = data.frame(
        group_1 = group_1,
        group_2 = group_2
      ),
      conditional_distribution= conditional_distribution,
      optimizer = "nlminb",
      weights=NULL,
      ...
    )
    if(inherits(out, "try-error")) {
      m <- m+1
    } else if (out$convergence==1){
      m <- m+1
    } else {
      convergence_error <- FALSE
      estimation_error  <- FALSE
      out$trials <- m
      out$frame <- list(
        values= parameter_values,
        dimensions=dimensions,
        scenario=censoring_scenario,
        distribution=conditional_distribution,
        X = ifelse(is.null(old_X),
                   "not provided",
                   old_X)
      )
    }
  }
  # TODO : include the elements in the simulation frame, i.e. dimensions, betas, and X (if needed)
  out
}

#' @title gof.aftRE
#'
crit.aftRE <- function(object, ...){
  logLik  <- object$log_likelihood
  k       <- ncol(object$frame$X)
  AIC <-  -2*logLik+2*k
}

#' @title gof.aftRE
#'
gof.aftRE <- function(object, ...){
  logLik  <- object$log_likelihood
  k       <- ncol(object$frame$X)
  AIC <-  -2*logLik+2*k
}

#' @title summary.aftRE: simulate a
#'
#' @param obj
#' @param
#' @examples
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom
#' @importFrom
summary.aftRE <- function(object, ...){

  # TODO : Retrieve the call
  n <- NROW(object$y)
  p <- NCOL(object$y)
  nX <- dim(object$X)[2]
  nTR <- dim(object$TR)[2]
  num.lv <- object$num.lv
  family <- object$family

  M <- cbind(object$params$beta0, object$params$theta)
  sumry <- list()
  sumry$'log-likelihood' <- object$logL
  crit <- inf.criteria(object)
  sumry$df <- crit$k
  sumry$AIC <- crit$AIC
  sumry$AICc <- crit$AICc
  sumry$BIC <- crit$BIC

  crit <-
    newnams <- c("Intercept")

  if (num.lv > 0)
    newnams <- c(newnams, paste("theta.LV", 1:num.lv, sep = ""))
  colnames(M) <- newnams
  rownames(M) <- colnames(object$y)
  sumry$Call <- object$call
  sumry$family <- object$family
  sumry$Coefficients <- M

  if (!is.null(object$TR)) {
    if (!is.null(object$X)) {
      sumry$'Covariate coefficients' <- object$params$B
    }
  } else {
    if (!is.null(object$X)) {
      sumry$'Environmental coefficients' <- object$params$Xcoef
    }
  }
  if (!is.null(object$params$row.params)) {
    sumry$'Row intercepts' <- object$params$row.params
  }

  if (object$row.eff == "random") {
    object$params$sigma2 = object$params$sigma ^ 2
    names(object$params$sigma2) = "sigma^2"
    sumry$'Variance of random row intercepts' <- object$params$sigma2
  }

  if (object$family == "negative.binomial") {
    sumry$'Dispersion parameters' <- object$params$phi
  }
  if (object$family == "tweedie") {
    sumry$'Dispersion parameters' <- object$params$phi
  }
  if (object$family == "ZIP") {
    sumry$'Zero inflation p' <- object$params$phi
  }
  if(object$family == "gaussian"){
    sumry$'Standard deviations' <- object$params$phi
  }
  class(sumry) <- "summary.gllvm"
  return(sumry)
}

show_progress <- function(x, B) {
  if (x%%(B/20)==0) {
    message("Progress: ", paste(rep("=", times=(x/(B/20)))), paste(rep(" ", times=(20-x/(B/20)))), x/B*100,"%")
  }
}


#' @title ranef.aftRE: simulate a
#'
#' @param obj
#' @param
#' @examples
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom
#' @importFrom
ranef.aftRE <- function(object, ...){
  out <- object$predictions
  out
}

#' @title extract.aftRE: simulate a
#'
#' @param obj
#' @param
#' @examples
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom
#' @importFrom
extract.aftRE <- function(object,
                          what= c("coefficients",
                                  "fixef",
                                  "ranef",
                                  "fitted",
                                  "residuals"),
                          type,  ...){
  what = match.arg(what)
  pred <-  object$predictions
  if(what=="coefficients"){
    object$estimates
  } else if(what == "fixef"){
    object$estimates
  }

  if(what == "residuals"){
    if(type == "EAM"){

    } else if (type == "ITM"){

    }
  }
}

#' @title plot.aftRE: simulate a
#'
#' @param obj
#' @param
#' @examples
#'
#'
#' @rdname multiDA
#' @export
#' @importFrom
#' @importFrom
plot.aftRE <- function(object, what, type, ...){

  if(what == "residuals"){
   if(type == "EAM"){

   } else if (type == "ITM"){

   }
  }
}


# All the methods we can create for the class "aftRE"
# anova :
# coef  => extract (what = "fitted")
# confint
# df.residual
# extractAIC
# family
# fitted  => extract(what = "fitted")
# fixef   => extract(what = "fixed effects")
# formula
# getME
# isLMM
# logLik
# model.frame
# model.matrix
# nobs
# predict
# print
# profile
# ranef
# refit
# residuals => extract(what = "residuals")
# sigma
# simulate
# summary
# terms
# VarCorr  => extract(what = "Variance Components")
# vcov
# weights
# see '?methods' for accessing help and source code

