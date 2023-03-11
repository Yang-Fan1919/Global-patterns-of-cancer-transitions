library(foreach)
library(tidyverse)
library(haven)
library(poLCA)
library(lcmm)
library(ggsci)
library(ggpmisc)
library(nortest)
library(car)
library(readxl)
library(doParallel)


setwd("D:/paper/global cancer trend/data")
set.seed(2020)

###Load GBD data(including age-standardized incidence and mortality)
GBD1<-read_csv("IHME-GBD_2019_DATA-b0ace1f1-1.csv")
GBD2<-read_csv("IHME-GBD_2019_DATA-b0ace1f1-2.csv")
GBD3<-read_csv("IHME-GBD_2019_DATA-b0ace1f1-3.csv")

###Select ten most common cancer type
cancer_selet <- read_xlsx("IHME-GBD_2019_DATA-87064d1a-1.xlsx")
cancer_selet <- cancer_selet[,1] %>% unlist()
cancer_both <- cancer_selet[-c(3,4,7)] 
cancer_male <- cancer_selet[-c(3,7)] 
cancer_female <- cancer_selet[-4]
com_con <- c(paste("Incidence", cancer_both,"Both", sep = "_"),paste("Deaths", cancer_both,"Both", sep = "_"),
             paste("Incidence", cancer_male,"Male", sep = "_"),paste("Deaths", cancer_male,"Male", sep = "_"),
             paste("Incidence", cancer_female,"Female", sep = "_"),paste("Deaths", cancer_female,"Female", sep = "_"))
GBD<-bind_rows(GBD1,GBD2,GBD3)%>%
  mutate(label=group_indices(.,measure,location,sex),
         year = year-1990)

###Using parallel computing to speed up modelling
detectCores()
cl <- makeCluster(7)
registerDoParallel(cl)
eColi <- foreach(i=1:length(com_con),.multicombine = T,.packages = c("tidyverse","car","lcmm"),.errorhandling = "pass") %dopar% {
  set.seed(2022)
  inf <- strsplit(com_con[i], "_", fixed = TRUE) %>% unlist()
  data <- subset(GBD, (measure==inf[1])&(cause==inf[2])&sex==inf[3])
  pow <- summary(powerTransform(data$val))$result[2]
  if (abs(pow)<0.1){
    data[,ncol(data)+1] <- log(data$val)
  }else{
    data[,ncol(data)+1] <- (data$val)^pow}
  colnames(data) <- c(colnames(data)[-ncol(data)],"trsval")
  fit <- vector("list", 3)
  names(fit) <- c("fit_model", "fit_result","power")
  fit_model <- vector("list", 5)
  fit_result <- vector("list", 2)
  names(fit_result) <- c("fit statistics", "posterior class probabilities")
  fit_result$`posterior class probabilities` <- vector("list", 5)
  
  
  ###GMMs for one to five classes using random intercept (including a quadratic term to capture non-linearity trend)
  fit_model[[1]] <- hlme(trsval ~ year+I(year^2), subject = "label", random=~1, ng = 1, data = data)
  fit_model[[2]] <- gridsearch(rep = 50, maxiter = 10, minit = fit_model[[1]],
                               hlme(trsval ~ year+I(year^2), random=~1 , subject = "label", ng = 2, data = data, mixture = ~ year,nwg=T))
  fit_model[[3]] <- gridsearch(rep = 50, maxiter = 10, minit = fit_model[[1]],
                               hlme(trsval ~ year+I(year^2), random=~1 , subject = "label", ng = 3, data = data, mixture = ~ year,nwg=T))
  fit_model[[4]] <- gridsearch(rep = 50, maxiter = 10, minit = fit_model[[1]],
                               hlme(trsval ~ year+I(year^2), random=~1 , subject = "label", ng = 4, data = data, mixture = ~ year,nwg=T))
  fit_model[[5]] <- gridsearch(rep = 50, maxiter = 10, minit = fit_model[[1]],
                               hlme(trsval ~ year+I(year^2), random=~1 , subject = "label", ng = 5, data = data, mixture = ~ year,nwg=T))
  fit_result[[1]] <- summarytable(fit_model[[1]],fit_model[[2]], fit_model[[3]],fit_model[[4]],fit_model[[5]],
                                  which = c("G", "loglik",  "conv", "npm","AIC", "BIC", "SABIC", "entropy", "%class"))
  fit_result[[2]][[1]] <- postprob(fit_model[[1]])
  fit_result[[2]][[2]] <- postprob(fit_model[[2]])
  fit_result[[2]][[3]] <- postprob(fit_model[[3]])
  fit_result[[2]][[4]] <- postprob(fit_model[[4]])
  fit_result[[2]][[5]] <- postprob(fit_model[[5]])
  fit$fit_model <- fit_model
  fit$fit_result <- fit_result
  fit$power <- pow
  save(fit,file = paste(com_con[i],".Rdata",sep = ""))
}
stopImplicitCluster()
