library(poLCA)
library(tidyverse)
library(lcmm)
library(readxl)
library(patchwork)
library(cowplot)


setwd("D:/paper/global cancer trend/data/result")
set.seed(2022)


LoadToEnvironment <- function(RData, env = new.env()){
        load(RData, env)
        return(env) 
}


###Load GBD data
GBD1<-read_csv("D:/paper/global cancer trend/data/IHME-GBD_2019_DATA-b0ace1f1-1.csv")
GBD2<-read_csv("D:/paper/global cancer trend/data/IHME-GBD_2019_DATA-b0ace1f1-2.csv")
GBD3<-read_csv("D:/paper/global cancer trend/data/IHME-GBD_2019_DATA-b0ace1f1-3.csv")


###Final solutions for all GMMs ((see selection criteria in the Method section of manuscript))
n_class<-read_xlsx("D:/paper/global cancer trend/data/Noofclass.xlsx")
GBD<-bind_rows(GBD1,GBD2,GBD3)%>%
        mutate(label=group_indices(.,measure,location,sex),
               year = year-1990)
rm("GBD1","GBD2","GBD3")
SDI_raw <- read_xlsx("D:/paper/global cancer trend/data/IHME_GBD_2019_SDI_1990_2019_Y2020M10D15.xlsx",sheet = 2)
SDI_quintile <- read_xlsx("D:/paper/global cancer trend/data/IHME_GBD_2019_SDI_1950_2019_QUINTILES_Y2021M03D21.xlsx")


###Select ten most common cancer type
cancer_selet <- read_xlsx("D:/paper/global cancer trend/data/IHME-GBD_2019_DATA-87064d1a-1.xlsx")
cancer_selet <- cancer_selet[,1] %>% unlist()
cancer_both <- cancer_selet[-c(3,4,7)] 
cancer_male <- cancer_selet[-c(3,7)] 
cancer_female <- cancer_selet[-4]
both_name <- paste(c(paste("Incidence", cancer_both,"Both", sep = "_"),
                       paste("Deaths", cancer_both,"Both", sep = "_")),".Rdata",sep = "")
male_name <- paste(c(paste("Incidence", cancer_male,"Male", sep = "_"),
                       paste("Deaths", cancer_male,"Male", sep = "_")),".Rdata",sep = "")
female_name <- paste(c(paste("Incidence", cancer_female,"Female", sep = "_"),
                         paste("Deaths", cancer_female,"Female", sep = "_")),".Rdata",sep = "")

both_result <- vector("list",length = length(both_name))
names(both_result) <- both_name
for (i in 1:length(both_name)) {
        both_result[[i]] <- LoadToEnvironment(both_name[i])$fit
}

male_result <- vector("list",length = length(male_name))
names(male_result) <- male_name
for (i in 1:length(male_name)) {
        male_result[[i]] <- LoadToEnvironment(male_name[i])$fit
}

female_result <- vector("list",length = length(female_name))
names(female_result) <- female_name
for (i in 1:length(female_name)) {
        female_result[[i]] <- LoadToEnvironment(female_name[i])$fit
}


###Male
#Incidence
#####
male_in_class <- subset(n_class,measure == "Incidence"&type%in%cancer_male&sex=="Male")
male_in_class$type <- factor(male_in_class$type,levels = cancer_male)
male_in_class <- male_in_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
male_incidence <- male_result[[1]]$fit_model[[male_in_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(male_result[[2]]$fit_model[[male_in_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[3]]$fit_model[[male_in_class[3]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[4]]$fit_model[[male_in_class[4]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[5]]$fit_model[[male_in_class[5]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[6]]$fit_model[[male_in_class[6]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[7]]$fit_model[[male_in_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[8]]$fit_model[[male_in_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))
colnames(male_incidence) <- c("label","Lung","Colorectum","Prostate","Stomach","Leukemia","Esophageal",
                              "Liver","Pancreatic","location","SDI")
male_incidence <- male_incidence %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Prostate","Stomach","Leukemia","Esophageal",
                    "Liver","Pancreatic","SDI_c"),as.factor)
f_male <- cbind(Lung, Colorectum, Prostate, Stomach, Leukemia, Esophageal, Liver, Pancreatic) ~ 1
male_incidence_lca <- vector("list",7)
###Latent class modelling
for (i in 1:7) {
        male_incidence_lca[[i]] <- poLCA(f_male,male_incidence,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


#Mortality
#####
male_de_class <- subset(n_class,measure == "Death"&type%in%cancer_male&sex=="Male")
male_de_class$type <- factor(male_de_class$type,levels = cancer_male)
male_de_class <- male_de_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
male_death <- male_result[[9]]$fit_model[[male_de_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(male_result[[10]]$fit_model[[male_de_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[11]]$fit_model[[male_de_class[3]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[12]]$fit_model[[male_de_class[4]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[13]]$fit_model[[male_de_class[5]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[14]]$fit_model[[male_de_class[6]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[15]]$fit_model[[male_de_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[16]]$fit_model[[male_de_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))
colnames(male_death) <- c("label","Lung","Colorectum","Prostate","Stomach","Leukemia","Esophageal",
                              "Liver","Pancreatic","location","SDI")
male_death <- male_death %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Prostate","Stomach","Leukemia","Esophageal",
                    "Liver","Pancreatic","SDI_c"),as.factor)
male_death_lca <- vector("list",7)
###Latent class modelling
for (i in 1:7) {
        male_death_lca[[i]] <- poLCA(f_male,male_death,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


###Female
#Incidence
#####
female_in_class <- subset(n_class,measure == "Incidence"&type%in%cancer_female&sex=="Female")
female_in_class$type <- factor(female_in_class$type,levels = cancer_female)
female_in_class <- female_in_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
female_incidence <- female_result[[1]]$fit_model[[female_in_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(female_result[[2]]$fit_model[[female_in_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[3]]$fit_model[[female_in_class[3]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[4]]$fit_model[[female_in_class[4]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[5]]$fit_model[[female_in_class[5]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[6]]$fit_model[[female_in_class[6]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[7]]$fit_model[[female_in_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[8]]$fit_model[[female_in_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[9]]$fit_model[[female_in_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))
colnames(female_incidence) <- c("label","Lung","Colorectum","Breast","Stomach","Leukemia","Cervical","Esophageal",
                              "Liver","Pancreatic","location","SDI")
female_incidence <- female_incidence %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Breast","Stomach","Leukemia","Cervical","Esophageal",
                    "Liver","Pancreatic","SDI_c"),as.factor)
f_female <- cbind(Lung, Colorectum, Breast, Stomach, Leukemia, Cervical, Esophageal, Liver, Pancreatic) ~ 1
female_incidence_lca <- vector("list",7)
###Latent class modelling
for (i in 1:7) {
        female_incidence_lca[[i]] <- poLCA(f_female,female_incidence,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


#Mortality
#####
female_de_class <- subset(n_class,measure == "Death"&type%in%cancer_female&sex=="Female")
female_de_class$type <- factor(female_de_class$type,levels = cancer_female)
female_de_class <- female_de_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
female_death <- female_result[[10]]$fit_model[[female_de_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(female_result[[11]]$fit_model[[female_de_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[12]]$fit_model[[female_de_class[3]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[13]]$fit_model[[female_de_class[4]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[14]]$fit_model[[female_de_class[5]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[15]]$fit_model[[female_de_class[6]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[16]]$fit_model[[female_de_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[17]]$fit_model[[female_de_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[18]]$fit_model[[female_de_class[8]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))
colnames(female_death) <- c("label","Lung","Colorectum","Breast","Stomach","Leukemia","Cervical","Esophageal",
                                "Liver","Pancreatic","location","SDI")
female_death <- female_death %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Breast","Stomach","Leukemia","Cervical","Esophageal",
                    "Liver","Pancreatic","SDI_c"),as.factor)
female_death_lca <- vector("list",7)
###Latent class modelling
for (i in 1:7) {
        female_death_lca[[i]] <- poLCA(f_female,female_death,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


#####Total population
###Incidence
#####
b_index <- subset(n_class,measure == "Incidence"&type == "Breast cancer"&sex=="Female")$class_last
female_result$`Incidence_Breast cancer_Female.Rdata`$fit_model[[b_index ]]$pprob$label <- 
        female_result$`Incidence_Breast cancer_Female.Rdata`$fit_model[[b_index ]]$pprob$label-1
b_dedex <- subset(n_class,measure == "Death"&type == "Breast cancer"&sex=="Female")$class_last
female_result$`Deaths_Breast cancer_Female.Rdata`$fit_model[[b_dedex ]]$pprob$label <- 
        female_result$`Deaths_Breast cancer_Female.Rdata`$fit_model[[b_dedex ]]$pprob$label-1

c_index <- subset(n_class,measure == "Incidence"&type == "Cervical cancer"&sex=="Female")$class_last
female_result$`Incidence_Cervical cancer_Female.Rdata`$fit_model[[c_index]]$pprob$label <- 
        female_result$`Incidence_Cervical cancer_Female.Rdata`$fit_model[[c_index]]$pprob$label-1
c_dedex <- subset(n_class,measure == "Death"&type == "Cervical cancer"&sex=="Female")$class_last
female_result$`Deaths_Cervical cancer_Female.Rdata`$fit_model[[c_dedex]]$pprob$label <- 
        female_result$`Deaths_Cervical cancer_Female.Rdata`$fit_model[[c_dedex]]$pprob$label-1

p_index <- subset(n_class,measure == "Incidence"&type == "Prostate cancer"&sex=="Male")$class_last
male_result$`Incidence_Prostate cancer_Male.Rdata`$fit_model[[p_index]]$pprob$label <- 
        male_result$`Incidence_Prostate cancer_Male.Rdata`$fit_model[[p_index]]$pprob$label-2
p_dedex <- subset(n_class,measure == "Death"&type == "Prostate cancer"&sex=="Male")$class_last
male_result$`Deaths_Prostate cancer_Male.Rdata`$fit_model[[p_dedex]]$pprob$label <- 
        male_result$`Deaths_Prostate cancer_Male.Rdata`$fit_model[[p_dedex]]$pprob$label-2
        
both_in_class <- subset(n_class,measure == "Incidence"&type%in%cancer_both&sex=="Both")
both_in_class$type <- factor(both_in_class$type,levels = cancer_both)
both_in_class <- both_in_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
both_incidence <- both_result[[1]]$fit_model[[both_in_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(both_result[[2]]$fit_model[[both_in_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(both_result[[3]]$fit_model[[both_in_class[3]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[4]]$fit_model[[both_in_class[4]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[5]]$fit_model[[both_in_class[5]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[6]]$fit_model[[both_in_class[6]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[7]]$fit_model[[both_in_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[3]]$fit_model[[b_index]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[3]]$fit_model[[p_index]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[6]]$fit_model[[c_index]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))

colnames(both_incidence) <- c("label","Lung","Colorectum","Stomach","Leukemia","Esophageal","Liver",
                              "Pancreatic","Breast","Prostate","Cervical","location","SDI")
both_incidence <- both_incidence %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Stomach","Leukemia","Esophageal","Liver","Pancreatic","Breast",
                    "Prostate","Cervical","SDI_c"),as.factor)
f_both <- cbind(Lung, Colorectum, Stomach, Leukemia, Esophageal, Liver, Pancreatic,Breast,Prostate,Cervical) ~ 1
both_incidence_lca <- vector("list",5)
###Latent class modelling
for (i in 1:5) {
        both_incidence_lca[[i]] <- poLCA(f_both,both_incidence,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


###Mortality
#####
both_de_class <- subset(n_class,measure == "Death"&type%in%cancer_both&sex=="Both")
both_de_class$type <- factor(both_de_class$type,levels = cancer_both)
both_de_class <- both_de_class %>% 
        arrange(type) %>% 
        select(class_last) %>% 
        unlist
both_death <- both_result[[8]]$fit_model[[both_de_class[1]]]$pprob[,c(1,2)] %>% 
        left_join(both_result[[9]]$fit_model[[both_de_class[2]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(both_result[[10]]$fit_model[[both_de_class[3]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[11]]$fit_model[[both_de_class[4]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[12]]$fit_model[[both_de_class[5]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[13]]$fit_model[[both_de_class[6]]]$pprob[,c(1,2)],by="label") %>%
        left_join(both_result[[14]]$fit_model[[both_de_class[7]]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[12]]$fit_model[[b_dedex]]$pprob[,c(1,2)],by="label") %>% 
        left_join(male_result[[11]]$fit_model[[p_dedex]]$pprob[,c(1,2)],by="label") %>% 
        left_join(female_result[[15]]$fit_model[[c_dedex]]$pprob[,c(1,2)],by="label") %>% 
        left_join(GBD[,c(2,11)],by="label") %>% distinct() %>% 
        left_join(SDI_raw,by = c("location"="Location"))

colnames(both_death) <- c("label","Lung","Colorectum","Stomach","Leukemia","Esophageal","Liver","Pancreatic","Breast",
                          "Prostate","Cervical","location","SDI")
both_death <- both_death %>% 
        mutate(SDI_c = case_when(SDI<as.numeric(SDI_quintile[1,3])~1,
                                 SDI>=as.numeric(SDI_quintile[2,2])&SDI<as.numeric(SDI_quintile[2,3])~2,
                                 SDI>=as.numeric(SDI_quintile[3,2])&SDI<as.numeric(SDI_quintile[3,3])~3,
                                 SDI>=as.numeric(SDI_quintile[4,2])&SDI<as.numeric(SDI_quintile[4,3])~4,
                                 SDI>=as.numeric(SDI_quintile[5,2])~5)) %>% 
        mutate_at(c("Lung","Colorectum","Stomach","Leukemia","Esophageal","Liver","Pancreatic","Breast",
                    "Prostate","Cervical","SDI_c"),as.factor)
both_death_lca <- vector("list",6)
###Latent class modelling
for (i in 1:6) {
        both_death_lca[[i]] <- poLCA(f_both,both_death,nclass=i,nrep = 500,verbose=F,maxiter = 2000)
}


###Final solution for latent class modelling (see selection criteria in the Method section of manuscript)
lca_result <- vector("list",6)
lca_result[[1]] <- male_incidence_lca[[2]]
lca_result[[2]] <- male_death_lca[[2]]
lca_result[[3]] <- female_incidence_lca[[2]]
lca_result[[4]] <- female_death_lca[[2]]
lca_result[[5]] <- both_incidence_lca[[2]]
lca_result[[6]] <- both_death_lca[[2]]
names(lca_result) <- c("male_incidence","male_death","female_incidence","female_death","both_incidence","both_death")
