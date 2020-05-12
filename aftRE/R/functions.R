create_survival_data <- function(
  data,
  durations,
  covariates,
  censoring_ind,
  grouping_factors,
  slope_1=NULL,
  slope_2=NULL,
  rand_strc_1 = c("intercept", "intercept and slope"),
  rand_strc_2 = c("none", "intercept", "intercept and slope"),
  conditional_family = c("log-normal", "log-logistic", "weibull")
){

  group_1_expanded <- grouping_factors[[1]] %>% as.factor()


  if(rand_strc_1=="intercept"){
    if(rand_strc_1=="intercept"){

    }


  } else if(rand_strc_1=="intercept"){

  }
  if(rand_strc_)


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

  conditional_distribution <- match.arg(conditional_distribution)

  if(conditional_distribution == "log-normal"){
    cond_resp <- 1
  } else if(conditional_distribution == "log-logistic"){
    cond_resp <- 2
  } else if(conditional_distribution == "weibull") {
    cond_resp <- 3
  }
}

estimate_aft <- function(
  object,
  weights = NULL,
  optimizer = "nlminb",
  ...
){



  if(is.null(weights)){
    w = rep(1, nlevels(group_1_expanded))
  } else {
    w = weights
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
