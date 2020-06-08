# Brian Stock
# May 11, 2020
# 2D AR1 survival devs (Haikun's paper)
# Fit models

# use 2019 assessment ASAP data file
#  data til 2018
#  6 selectivity blocks for fleet

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/1d_fit_models_M_GSI_v3.R")
# v4: fix M_re bug (no code here changed)

# remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
# remotes::install_github("timjmiller/wham", dependencies=TRUE)
library(here)
library(wham) # https://timjmiller.github.io/wham
library(ecodata) # EDAB package, hosts GSI and other env data, https://github.com/NOAA-EDAB/ecodata
library(tidyverse)

# # get SNEMAYT data from wham
# wham.dir <- find.package("wham")
# file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=here("data","ex2_SNEMAYT.dat"), overwrite=FALSE)

# set up wham
# 2019 assessment from https://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php
# https://fish.nefsc.noaa.gov/saw/reviews_report_options.php
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))

# GSI: rw + est obs error + 2
gsi <- read.csv(here("data","GSI_v4.csv"), header=TRUE) 

# NAA_re 4 options: 'iid','ar1_a','ar1_y','2dar1'
# M_re 5 options: 'none','iid','ar1_a','ar1_y','2dar1'
# M 2 options: fixed, estimate mean M
df.mods <- data.frame(NAA_re = rep(c('iid','ar1_a','ar1_y','2dar1'),10),
                      est_M = rep(c(FALSE,TRUE),each=20),
                      M_re = rep(rep(c('none','iid','ar1_a','ar1_y','2dar1'),each=4),2), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

# run models
for(m in 1:n.mods){
  NAA_list <- list(cor=df.mods[m,"NAA_re"], sigma='rec+1')

  if(df.mods[m,"est_M"]){
    M_list <- list(model="constant", est_ages=1, re=df.mods$M_re[m])
  } else {
    M_list <- list(re=df.mods$M_re[m])
  }

  ecov <- list(
    label = "GSI",
    mean = as.matrix(gsi$GSI.mean),
    logsigma = "est_1",
    year = gsi$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(gsi)[1]),
    lag = 1, 
    process_model = "rw",
    where = "recruit", 
    how = 2,
    link_model = "linear")

  # 2. use logistic sel
  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = df.mods$Model[m],                         
                              NAA_re = NAA_list,
                              M = M_list,
                              ecov=ecov,
                              selectivity=list(model=rep("logistic",9),
                                 initial_pars=rep(list(c(3,3)),9),
                                 fix_pars=rep(list(NULL),9)))  
  
  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  if(m %% 4 == 2) input$par$logit_selpars[1,7:8] = c(-1.5,-5)

# input$data$age_comp_model_indices = rep(7, input$data$n_indices)
# input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
# input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
# input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
# input$par$index_paa_pars = rep(0, input$data$n_indices)
# input$par$catch_paa_pars = rep(0, input$data$n_fleets)
# input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

  # Fit model with projections:
  #  - 3 years
  #  - F = 0
  #  - WAA, MAA, maturity fixed at terminal year 
  mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F)  
  # mod <- fit_wham(input, do.retro=T, do.osa=T, proj.opts = list(proj.F=rep(0,3), avg.yrs=2011))  

  if(exists("err")) rm("err") # need to clean this up
  opt_conv = as.logical(1-mod$opt$convergence)
  if(mod$na_sdrep==FALSE & !is.na(mod$na_sdrep)) ok_sdrep=TRUE else ok_sdrep=FALSE

  if(opt_conv){
    n.try <- 1
    new_input <- input
    while(!ok_sdrep & n.try < 11){
      n.try = n.try + 1
      # jitter pars and try again
      new_input <- input
      new_input$par <- mod$env$parList(mod$opt$par + rnorm(length(mod$opt$par), 0, 0.2)) 
      mod <- fit_wham(new_input, do.retro=F, do.osa=F, do.proj=F) 
      if(exists("err")) rm("err") # need to clean this up
      if(mod$na_sdrep==FALSE & !is.na(mod$na_sdrep)) ok_sdrep=TRUE else ok_sdrep=FALSE
    }

    if(ok_sdrep){ # do retro + proj if pdHess = TRUE
      # new_input <- input
      # new_input$par <- mod$env$parList(mod$opt$par)      
      mod <- fit_wham(new_input, do.retro=T, do.osa=F, do.proj=T)
      if(exists("err")) rm("err") # need to clean this up
    }
  }
  saveRDS(mod, file=here("results","dat_2019","NAA_M_GSI_v3",paste0(df.mods$Model[m],".rds"))) 
}

# # collect fit models into a list
# # df.mods <- df.mods[1:8,]
# mod.list <- here("results",paste0(df.mods$Model,".rds"))
# # mod.list <- here("results","index4_age_fix36_index1_age",paste0(df.mods$Model,".rds"))
# mods <- lapply(mod.list, readRDS)

# # t(sapply(mods, function(x) x$parList$logit_selpars[4,1:6]))
# # TMBhelper::Check_Identifiable(mods[[3]])
# # mods[[3]]$parList$logit_selpars

# opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
# ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
# df.mods$conv <- as.logical(opt_conv)
# df.mods$pdHess <- as.logical(ok_sdrep)

# # df.mods$na_sdrep <- sapply(mods, function(x) x$na_sdrep)
# df.mods$runtime <- sapply(mods, function(x) x$runtime)
# df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
# df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)
# df.aic$AIC[df.mods$pdHess==FALSE] <- NA
# minAIC <- min(df.aic$AIC, na.rm=T)
# df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
# df.mods <- cbind(df.mods, df.aic)
# rownames(df.mods) <- NULL

# # look at results table
# df.mods

# write.csv(df.mods, file=here("results","index4_age_fix36_index1_age.csv"))

# logit_selpars for models that didn't converge
# notconv <- which(!df.mods$pdHess)
# lapply(mods[notconv], function(x) x$parList$logit_selpars[,1:6])

# # notconv
# # 1  2  3  4  8 11 12 16 19 23 24 28 31
# # bad pars
# # 1-4, 8: high gradients
# # 11,12,16,19: index 1 logit_selpars
# # 23,24,28: mean_rec_pars
# # 31: index 1 logit_selpars AND mean_rec_pars

# TMBhelper::Check_Identifiable(mods[[11]])
# mods[[11]]$parList$logit_selpars

# TMBhelper::Check_Identifiable(mods[[3]])
# TMBhelper::Check_Identifiable(mods[[4]])


# # all of the ar1_y and 2dar1 models aren't converging
# # looks like selectivity of index 3 is the issue
# mods[[3]]$parList$logit_selpars
# TMBhelper::Check_Identifiable(mods[[3]])

# mods[[4]]$parList$logit_selpars
# TMBhelper::Check_Identifiable(mods[[4]])

# # plot output for all models
# for(m in 1:n.mods){
#   plot_wham_output(mod=mods[[m]], dir.main=file.path(getwd(),paste0("m",m)), out.type='html')
# }

# round 4
# age-specific for index 3 with ages 3 & 6 fixed, age-specific for index 1 no ages fixed
# mod.list <- here("results","index4_age_fix36_index1_age",paste0(df.mods$Model,".rds"))
# notconv <- which(!df.mods$pdHess)
# notconv
# # 1  2  3  4  7  8 15 20 24 28 32

# TMBhelper::Check_Identifiable(mods[[1]])
# # bad pars seem like estimating recruitment is the issue, not selectivity
# # 1-4,7-8: high gradients
# # 15,20,28,32: mean_rec_pars
# # 24: mean_rec_pars, logit_q[1], all logit_selpars for index 1

