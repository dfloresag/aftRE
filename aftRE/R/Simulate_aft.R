rm(list=ls())

source("./R/aftRE.R")
library(dplyr)

samples = 1000
simulations <- list()
start_time <-  Sys.time()

# dims =list(n_1 = 30, n_12 = 2, n_23 = 2)
# dims =list(n_1 = 30, n_12 = 2, n_23 = 50)
# dims =list(n_1 = 30, n_12 = 50, n_23 = 2)
# dims =list(n_1 = 30, n_12 = 50, n_23 = 50)

things <- lapply(
  X = as.list(1:samples),
  FUN = function(x){
    show_progress(x, B = samples)
    tmp <- simulate_aft_nested(
      dimensions = dims,
      conditional_distribution = "weibull",
      censoring_scenario = 4)
    if(tmp$convergence == 0){
    res <- data.frame(
      Sample= x,
      n_1  = tmp$frame$dimensions$n_1,
      n_12 = tmp$frame$dimensions$n_12,
      n_23 = tmp$frame$dimensions$n_23,
      tmp$estimates)
    }
    res
  }
  )

end_time <-  Sys.time() - start_time
things = bind_rows(things)

# summary(things)

filename <- paste(paste("sim",
                  "n_1" , dims$n_1,
                  "n_12", dims$n_12,
                  "n_23", dims$n_23, sep = "_"), ".csv", sep="")

write.csv(things, paste0("data/", filename))


dfs <- list()

# dims =list(n_1 = 30, n_12 = 2, n_23 = 2)
# dims =list(n_1 = 30, n_12 = 2, n_23 = 50)
# dims =list(n_1 = 30, n_12 = 50, n_23 = 2)
# dims =list(n_1 = 30, n_12 = 50, n_23 = 50)
# dfs[[4]] <-read.csv(file = paste(paste("data/sim",
#                      "n_1" , 30,
#                      "n_12", 50,
#                      "n_23", 50, sep = "_"), ".csv", sep=""), header =TRUE)
# tmp <- bind_rows(dfs)
# head(tmp)

# tmp %>% select(-X) %>% write.csv("data/sim_wb.csv", row.names = FALSE)


par(mfrow=c(2,2))
things %>% filter(Parameter=="intercept")%>%
  select(Estimates) %>% boxplot()
abline(h=1)
things %>% filter(Parameter=="X_1_cont")%>%
  select(Estimates) %>% boxplot()
abline(h=1)
things %>% filter(Parameter=="X_2_cont")%>%
  select(Estimates) %>% boxplot()
abline(h=-.5)
things %>% filter(Parameter=="X_3_dich")%>%
  select(Estimates) %>% boxplot()
abline(h=1)

par(mfrow=c(1,3))
things %>% filter(Parameter=="sigma_0")%>%
  select(Estimates) %>% boxplot()
abline(h=1.21)
things %>% filter(Parameter=="sigma_1")%>%
  select(Estimates) %>% boxplot()
abline(h=0.2)
things %>% filter(Parameter=="sigma_2")%>%
  select(Estimates) %>% boxplot()
abline(h=0.2)



# library(frailtyHL)
# search()
#
# data(cgd, package="frailtyHL")
# data_surv <- cgd
# mlmc1 <- jointmodeling(Model="mean",RespDist="AFT",Link="log",
#                        LinPred=Surv(tstop-tstart,status)~treat+(1|center)+(1|id),
#                        RandDist="gaussian")
# res_cgd <- mlmfit(jm1 = mlmc1, Maxiter=300)

# bootstrapping

B = 1000

start_time <-  Sys.time()

bootstrap.aftRE <- function(object, n_replicates=1000 , method = "rwlb",...){
  method = match.arg(method)
  if (method == "rwlb"){

    X <- object$frame$X
    n1 <- nlevels(object$frame$group_1)

    dfs <- lapply(
      X = as.list(1:n_replicates),
      FUN = function(x){
        tmp <- simulate_aft_nested(
          conditional_distribution = "weibull",
          # TODO : Modify this
          censoring_scenario = 4,
          weights = rexp(n=n1, rate=1), old_X = X)

        if(tmp$convergence == 0){
          res <- data.frame(Replicate = x,
                            n_1= tmp$frame$dimensions$n_1,
                            n_12= tmp$frame$dimensions$n_12,
                            n_23= tmp$frame$dimensions$n_23,
                            tmp$estimates)
        }
        res
      }
    )
  }
}

things <- lapply(
  X = as.list(1:B),
  FUN = function(x){
    tmp <- simulate_aft_nested(
      conditional_distribution = "weibull",
      censoring_scenario = 4,
      weights = rexp(n = 20, ), new_X = )

    if(tmp$convergence == 0){
      res <- data.frame(Replicate = x,
                        n_1= tmp$frame$dimensions$n_1,
                        n_12= tmp$frame$dimensions$n_12,
                        n_23= tmp$frame$dimensions$n_23,
                        tmp$estimates)
    }
    res
  }
)
