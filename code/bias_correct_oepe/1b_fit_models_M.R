# Brian Stock
# May 26, 2020
# 2D AR1 survival devs (Haikun's paper)
# Fit models

# use 2019 assessment ASAP data file
#  data til 2018
#  6 selectivity blocks for fleet
#  change from age-specific to logistic sel

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/bias_correct_oepe/1b_fit_models_M.R")

# remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
# remotes::install_github("timjmiller/wham", ref="om_mode", dependencies=TRUE)
# remotes::install_github("timjmiller/wham", dependencies=TRUE)
library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)
# devtools::load_all("/home/bstock/Documents/wham")

# set up wham
# 2019 assessment from https://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php
# https://fish.nefsc.noaa.gov/saw/reviews_report_options.php
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))
res_dir <- here("results","revision2","M")
dir.create(res_dir, showWarnings=FALSE)

df.mods <- data.frame(M_re = c('none','iid','ar1_a','ar1_y','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("M-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods

# Base is same as in NAA, copy "Base.rds" to M results folder
file.copy(from=here("results","revision2","NAA","Base.rds"), to=here("results","revision2","M","Base.rds"))

# run models
for(m in 2:n.mods){
  # blocks 1, 3, 9 have issues, try age-specific
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  input <- prepare_wham_input(asap3, recruit_model = 2, # no SR relationship
                              model_name = df.mods$Model[m],                         
                              NAA_re = list(cor="iid", sigma="rec"),
                              M = list(re=df.mods$M_re[m]),
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6)),
                              age_comp = "logistic-normal-pool0") 

  # Fit model with projections:
  #  - 3 years
  #  - F = 0
  #  - WAA, MAA, maturity fixed at terminal year (2011)
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

write.csv(df.mods, file=here("tables","revision2","M.csv"))

