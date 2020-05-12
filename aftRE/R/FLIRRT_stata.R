rm(list=ls())

# Log-Normal Results
library(dplyr)
library(tidyr)

parameter_names = c(
  "(Intercept)", "trt_PP", "sex_male", "age", "apache",
  "log_sigma_0", "log_sigma_1", "sigma_2",
  "sigma_0", "sigma_1", "sigma_2",
  "sigma2_0", "sigma2_1","sigma2_2"
)

tab_ln   = data.frame(

  Parameter = parameter_names,
  Estimates = c(
    3.620672, 0.6001863,  -.6774046, -.0266254, 0.0174304,
    0.1605357, NA, NA,
    NA, NA, NA,
    NA, 0.036606, 0.2532477
  ),
  Std_Errors = c(
    0.8347563, 0.2704065, 0.2871157, 0.0141525, 0.0046777,
    0.0685798, NA, NA,
    NA, NA, NA,
    NA, 0.178874, 0.2085091
  ),
  Distribution = "log-normal",
  Implementation  = "stata"
)

# Log-Logistic

tab_ll   = data.frame(
  Parameter = parameter_names,
  Estimates = c(
    4.619146, 0.7488391,  -0.5888637, -.0369474, 0.0120979,
    -0.5489914, NA, NA,
    NA, NA, NA,
    NA, 8.75e-33, .1407601
  ),
  Std_Errors =  c(
    0.670402, 0.2186811, 0.2249491, 0.0112046, 0.0037222,
    0.0759289, NA, NA,
    NA, NA, NA,
    NA, 3.22e-17, 0.0746872
  ),
  Distribution = "log-logistic",
  Implementation  = "stata"
)

# gamma
tab_gg   = data.frame(
  Parameter = parameter_names,
  Estimates = c(
    5.062561, 0.7376002,  -.6115868, -.0393207, 0.0116213,
    -0.1963199, NA, NA,
    NA, NA, NA,
    NA, 1.03e-34, .1329063
  ),
  Std_Errors = c(
    .6631211, .2022414, .2191173, .0110357, 0.0036367,
    .0546629, NA, NA,
    NA, NA, NA,
    NA, 1.37e-18, 0.0739107
  ),
  Distribution = "gamma",
  Implementation  = "stata"
)

# Weibull

tab_wb = data.frame(
  Parameter = parameter_names,
  Estimates = c(
    5.085003, 0.7305873,  -.5892463 , -.0385324, 0.0111244,
    -0.3246589, NA, NA,
    NA, NA, NA,
    NA, 2.21e-35, .1351321
  ),
  Std_Errors = c(
    0.6268485, 0.1991695 , 0.2105434 , 0.0104993, 0.0034611,
    0.075043, NA, NA,
    NA, NA, NA,
    NA, 6.86e-19, 0.0694978
  ),
  Distribution = "weibull",
  Implementation  = "stata"
)

stata_est <- bind_rows(tab_ln,tab_ll, tab_wb)

write.csv(stata_est, file = "./data/FLIRRT_estimates_stata.csv", row.names = FALSE)
read.csv("./data/FLIRRT_estimates_stata.csv")
