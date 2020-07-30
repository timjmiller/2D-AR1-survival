# Brian Stock
# July 30, 2020
# 2D AR1 smoother
# Fit final set of models with
#   CPI data (so can compare AIC), and
#   F=0 in projection years (so can compare SSB projections)

# use 2019 assessment ASAP data file, from https://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php
# https://fish.nefsc.noaa.gov/saw/reviews_report_options.php
#  data til 2018
#  removed larval indices
#  6 selectivity blocks for fleet
#  change from age-specific to logistic sel

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/bias_correct_oe/1i_fit_models_best_zeroFproj.R")

# remotes::install_github("timjmiller/wham", ref="om_mode", dependencies=TRUE)
library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

thedir = here("results","dat_2019","bias_correct_oe","best_zeroFproj")
dir.create(thedir, showWarnings=FALSE)

# assessment data
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

# CPI data
cpi <- read.csv(here("data","CPI_v3.csv"), header=TRUE) 
colnames(cpi) <- c("Year","CPI","CPI_sigma")
cpi$use <- 1
cpi$use[is.nan(cpi$CPI)] = 0 # don't use 2017 (NaN, fall survey missing)

# CPI model = AR1, limiting (2)
# don't estimate mean M
df.mods <- data.frame(NAA_cor = c('iid','iid','2dar1','iid','iid','iid','2dar1','iid'),
                      NAA_sigma = c('rec','rec+1','rec+1','rec','rec','rec+1','rec+1','rec+1'),
                      M_re = c('none','none','none','iid','2dar1','2dar1','iid','2dar1'),
                      CPI_how = c(rep(0,7),2), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods$lab <- c("NAA-1","NAA-3","NAA-6","M-2","M-5","NAA-M-2","NAA-M-5","NAA-M-CPI-2")
df.mods <- df.mods %>% select(Model, lab, everything()) # moves Model to first col
df.mods

# run models
for(m in 1:n.mods){
  # CPI-recruitment specification
  ecov <- list(
    label = "CPI",
    mean = as.matrix(cpi$CPI),
    logsigma = matrix(log(cpi$CPI_sigma), ncol=1, nrow=dim(cpi)[1]),
    year = cpi$Year,
    use_obs = matrix(cpi$use, ncol=1, nrow=dim(cpi)[1]), # use all obs except 2017 (missing)
    lag = 1, # CPI in year t affects Rec in year t + 1
    process_model = "ar1",
    where = "recruit", # CPI affects recruitment
    how = df.mods$CPI_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = "linear")  

  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = df.mods$lab[m],                         
                              NAA_re = list(cor=df.mods$NAA_cor[m], sigma=df.mods$NAA_sigma[m]),
                              M = list(re=df.mods$M_re[m]),
                              ecov = ecov,
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6))) 

  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # bias correct predicted CAA and IAA
  input$data$bias_correct_oe = 1

  # # start M-2 at optimal pars (convergence issues)
  # tmpmod <- fit_wham(input, do.fit=FALSE)
  # m2 = readRDS("/home/bstock/Documents/ms/2D-AR1-survival/results/dat_2019/bias_correct_oe/M/m2.rds")
  # tmppar <- tmpmod$env$parList(par=m2$env$last.par.best)
  # m5 = readRDS("/home/bstock/Documents/ms/2D-AR1-survival/results/dat_2019/bias_correct_oe/best_zeroFproj/m5.rds")
  # tmppar$Ecov_process_pars = m5$rep$Ecov_process_pars
  # tmppar$Ecov_re = m5$rep$Ecov_re
  # input$par <- tmppar

  # Fit model with projections:
  #  - 3 years, 2019-2021 (default)
  #  - WAA, maturity fixed at terminal year, 2018 (default)  
  #  - F = 0
  #  - continue M re process (default)
  mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=T, proj.opts=list(proj.F=rep(0,3), cont.Mre=TRUE))  

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=file.path(thedir, paste0(df.mods$Model[m],".rds")))
}

# didn't do projections... add after fit
for(m in c(1:3,5:8)){
  tmp <- readRDS(file.path(thedir, paste0("m",m,".rds")))
  mod <- project_wham(tmp, proj.opts=list(n.yrs=3, proj.F=rep(0,3), cont.ecov=TRUE))
  saveRDS(mod, file=file.path(thedir, paste0("m",m,".rds")))
}

