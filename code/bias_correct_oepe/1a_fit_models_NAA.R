# Brian Stock
# May 26, 2020
# 2D AR1 survival devs (Haikun's paper)
# Fit models

# use 2019 assessment ASAP data file
#  data til 2018
#  6 selectivity blocks for fleet
#  change from age-specific to logistic sel

# library(here); source(here("code","bias_correct_oepe","1a_fit_models_NAA.R"))

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)
# devtools::load_all("/home/bstock/Documents/wham")

# 2019 assessment from https://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php
# https://fish.nefsc.noaa.gov/saw/reviews_report_options.php
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))
res_dir <- here("results","revision2","NAA")
dir.create(res_dir, showWarnings=FALSE)

df.mods <- data.frame(NAA_cor = c('iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('rec','rec','rec+1','rec+1','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("NAA-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# run models
for(m in 1:n.mods){
  
  # blocks 1, 3, 9 have issues, try age-specific
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  input <- prepare_wham_input(asap3, recruit_model = 2, # no SR relationship
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"]),
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6)),
                              age_comp = "logistic-normal-pool0") 

  # Fit model with projections:
  #  - 3 years
  #  - F = 0
  #  - default: use WAA, MAA, maturity averaged over last 5 years
  mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=T, proj.opts=list(proj.F=rep(0.001, 3))) 

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=file.path(res_dir,paste0(df.mods$Model[m],".rds")))
}

# collect fit models into a list
mod.list <- file.path(res_dir,paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# calc results table
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=T)$tab)
df.aic$AIC[df.mods$pdHess==FALSE] <- NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

write.csv(df.mods, file=here("tables","revision2","NAA.csv"))

