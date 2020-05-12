# FLIRRT Example

rm(list=ls())

# Log-Normal Results
library(dplyr)
library(tidyr)

# ML estimation of all parameters ####

# TODO : Store the FLIRRT Datasete in a format that is appropriate for analysis with this tool.

dta <- read.csv(
  "./data/FLIRRT.csv", na.strings = ""
) %>%
  na.omit() %>%
  arrange(flirrtid) %>%
  mutate(
    id       = factor(flirrtid),
    id.site  = factor(paste(flirrtid,site, sep=".")),
    id.side  = factor(paste(flirrtid,side, sep=".")),
    sex_male = ifelse(sex=="M", 1, 0),
    trt_PP   = ifelse(group_PP=="C1", 1, 0)
  )

t   <- dta$life   %>% as.vector()
cn  <- dta$censor %>% as.vector()

X   <- dta %>%
  mutate(
    '(Intercept)' = 1
  ) %>%
  select(
    '(Intercept)', trt_PP, sex_male, age, apache
  ) %>%
  as.matrix()

group_1_id      <- dta$id
group_2_id.site <- dta$id.site

source("./R/aftRE.R")

ln_aft <- estimate_aft_nested(
  durations_t  = t,
  covariates_X = X,
  censoring_ind = cn,
  grouping_factors = data.frame(
    group_1 = group_1_id,
    group_2 = group_2_id.site
  ),
  conditional_distribution ="log-normal",
  optimizer = "nlminb",
  weights=NULL
)

ll_aft <- estimate_aft_nested(
  durations_t  = t,
  covariates_X = X,
  censoring_ind = cn,
  grouping_factors = data.frame(
    group_1 = group_1_id,
    group_2 = group_2_id.site
  ),
  conditional_distribution ="log-logistic",
  weights = NULL,
  compile = TRUE
)

wb_aft <- estimate_aft_nested(
  durations_t  = t,
  covariates_X = X,
  censoring_ind = cn,
  grouping_factors = data.frame(
    group_1 = group_1_id,
    group_2 = group_2_id.site
  ),
  conditional_distribution ="weibull",
  weights = NULL,
  compile = TRUE
)

res_df <- dplyr::bind_rows(
  list(
    data.frame(ln_aft$estimates,
               Distribution = ln_aft$distribution,
               Implementation = "TMB"),
    data.frame(ll_aft$estimates,
               Distribution = ll_aft$distribution,
               Implementation = "TMB"),
    data.frame(wb_aft$estimates,
               Distribution = wb_aft$distribution,
               Implementation = "TMB"),
    read.csv("data/FLIRRT_estimates_stata.csv")
  )
)

# write.csv(res_df,"data/FLIRRT_results.csv" )
# write.csv(res_df,"/home/danitsp/Dropbox/[Uni] MONASH/[AFT-Documents]/[Reports]/Draft/datasets/FLIRRT_results.csv" )
# write.csv(res_df,"/home/danitsp/Dropbox/[Uni] MONASH/[AFT-Documents]/[Reports]/Datasets/FLIRRT_results.csv" )

library(ggplot2)

res_df %>%
  filter(Parameter %in% c("(Intercept)", "trt_PP", "sex_male", "age", "apache")) %>%
  ggplot()+
  geom_pointrange(aes(x = Parameter,
                      y = Estimates,
                      ymin = Estimates+qnorm(0.025)*Std_Errors,
                      ymax = Estimates+qnorm(0.975)*Std_Errors,
                      color= Implementation),
                  position = position_dodge(width=0.5))+
  facet_wrap(.~Distribution) +
  ggtitle("Comparison: Estimates + CI")

res_df %>%
  filter(Parameter %in% c("log_sigma_0", "sigma2_1","sigma2_2")) %>%
  ggplot()+
  geom_point(aes(Parameter, Estimates, color=Implementation), size = 4) +
  facet_wrap(.~Distribution) +
  ggtitle("Comparison: Standard Errors")



res_df %>%
  filter(Parameter %in% c("(Intercept)","trt_PP","sex_male","age","apache","log_sigma_0","sigma2_1","sigma2_2"))%>%
  select(Parameter, Distribution, Implementation, Estimates, Std_Errors) %>%
  pivot_wider(
    id_cols = c(Parameter,Distribution),
    names_from = Implementation,
    values_from = c(Estimates, Std_Errors)) %>%
  mutate(
    Diff_Estimates = Estimates_stata - Estimates_TMB,
    Diff_Std_Errors = Std_Errors_stata - Std_Errors_TMB) %>%
  ggplot() +
  geom_col(aes(x = Parameter, Diff_Estimates))+
  facet_wrap(.~Distribution)

res_df %>%
  filter(Parameter %in% c("(Intercept)","trt_PP","sex_male","age","apache","log_sigma_0","sigma2_1","sigma2_2"))%>%
  select(Parameter, Distribution, Implementation, Estimates, Std_Errors) %>%
  pivot_wider(
    id_cols = c(Parameter,Distribution),
    names_from = Implementation,
    values_from = c(Estimates, Std_Errors)) %>%
  mutate(
    Diff_Estimates = Estimates_stata - Estimates_TMB,
    Diff_Std_Errors = Std_Errors_stata - Std_Errors_TMB) %>%
  ggplot() +
  geom_col(aes(x = Parameter, Diff_Std_Errors))+
  facet_wrap(.~Distribution)




