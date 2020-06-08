# Brian Stock
# May 26, 2020
# 2D AR1 survival devs (Haikun's paper)
# Fit models

# use 2019 assessment ASAP data file
#  data til 2018
#  6 selectivity blocks for fleet
#  change from age-specific to logistic sel

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/1a_fit_models_NAA_GSI.R")

# remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
# remotes::install_github("timjmiller/wham", dependencies=TRUE)
library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)
devtools::load_all("/home/bstock/Documents/wham")

# # get SNEMAYT data from wham
# wham.dir <- find.package("wham")
# file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=here("data","ex2_SNEMAYT.dat"), overwrite=FALSE)

# created in 0_explore_data.R
gsi <- read.csv(here("data","GSI_v4.csv"), header=TRUE) 

# set up wham
# 2019 assessment from https://www.nefsc.noaa.gov/saw/sasi/sasi_report_options.php
# https://fish.nefsc.noaa.gov/saw/reviews_report_options.php
asap3 <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval.dat"))
# asap3_rmblocks <- read_asap3_dat(here("data","SNEMAYT_2019_rmlarval_rmblocks.dat"))

df.mods <- data.frame(NAA_cor = rep(c('iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),4),
                      NAA_sigma = rep(c('rec','rec','rec+1','rec+1','rec+1','rec+1'),4),
                      GSI_mod = c(rep("rw",12), rep("ar1",12)),
                      GSI_obserr = "est",
                      GSI_how = rep(c(rep(0,6),rep(2,6)),2), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

# run models
for(m in 1:n.mods){
# for(m in 17:n.mods){
# for(m in 1:8){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])

  if(df.mods$GSI_obserr[m] == "est"){
    gsi_obserr_val <- "est_1"
  } else {
    gsi_obserr_val <- as.matrix(log(gsi$GSI.se))
  }
  
  ecov <- list(
    label = "GSI",
    mean = as.matrix(gsi$GSI.mean),
    logsigma = gsi_obserr_val, # set above
    year = gsi$Year,
    use_obs = matrix(1, ncol=1, nrow=dim(gsi)[1]), # use all obs (=1)
    lag = 1, # GSI in year t affects Rec in year t + 1
    process_model = df.mods$GSI_mod[m], # "rw" or "ar1"
    where = "recruit", # GSI affects recruitment
    how = df.mods$GSI_how[m], # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
    link_model = "linear")
  

  # 1. selectivity as in ASAP file
  #      original: 6 blocks on fleet, age-specific fixed at 1 for ages 4-6
  #      rmblocks: 1 block on fleet, age-specific fixed at 1 for ages 4-5, 4, 2-4, none
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)

  # # prepare_wham_input doesn't read asap3 sel options correctly... need to fix. Specify here
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov,
  #                             selectivity=list(model=rep("age-specific",9),
  #                                initial_pars=list(c(.5,1,1,1,1,1), c(.01,.25,.75,1,1,1), c(.01,1,1,1,1,1), c(.01,.25,.75,1,1,1), c(.01,.25,.75,1,1,1), c(.01,.25,.75,1,1,1),
  #                                             c(.5,1,1,1,1,1), c(.5,1,1,1,1,1), c(.5,1,1,1,1,1)),
  #                                fix_pars=list(2:6, 4:6, 2:6, 4:6, 4:6, 4:6, 3:6, 4:6, 2:6)))

  # 2. use logistic sel
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov,
  #                             selectivity=list(model=rep("logistic",9),
  #                                initial_pars=rep(list(c(3,3)),9),
  #                                fix_pars=rep(list(NULL),9)))  
  
  # 2. blocks 1, 3, 9 have issues, try age-specific
  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = df.mods$Model[m],                         
                              NAA_re = NAA_list,
                              ecov=ecov,
                              selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
                                 initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
                                 fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6))) 

  # 3. like ICES papers, remove fleet selectivity blocks and use age-specific
  #      asap3_rmblocks: 1 block on fleet, age-specific fixed at 1 for ages 4-5, 4, 2-4, none
  # input <- prepare_wham_input(asap3_rmblocks, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)
  # input <- prepare_wham_input(asap3_rmblocks, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov,
  #                             selectivity=list(model=c("age-specific","age-specific","age-specific","age-specific"), re=c("none","none","none","none"),
  #                                             initial_pars=list(c(0.1,0.5,0.5,1,1,0.5), c(0.1,0.5,1,1,0.5,0.5), c(0.5,1,1,1,1,1), c(0.5,1,1,1,0.5,0.5)),
  #                                             fix_pars=list(4:5, 3:4, 2:5, 2:4)))


# sel_inits_age <- c(plogis(-4.0331660),0.5,0.5,1,1,0.5) # fix age 1 sel at mod 6 estimate (time constant, converged)
# fix_pars_age <- c(1,4,5)
# sel_inits_logistic <- c(2.1,0.15) # mods 2-4 converge with logit values ~ -0.63, -3.74
# fix_pars_logistic <- NULL  
# selectivity=list(model=c(df$sel_mod[m],"age-specific","age-specific"), re=c(df$sel_re[m],"none","none"),
#                     initial_pars=list(sel_inits, c(0.5,0.5,0.5,1,0.5,0.5), c(0.5,1,1,1,0.5,0.5)),
#                     fix_pars=list(fix_pars, 4, 2:4)))

#   input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
#                               model_name = df.mods$Model[m],                         
#                               NAA_re = NAA_list,
#                               ecov=ecov,
#                               selectivity=list(model=rep("age-specific",9),
#                                  initial_pars=list(c(.5,1,1,1,1,1), c(.01,.25,.75,1,1,1), c(.01,1,1,1,1,1), c(.01,.25,.75,1,1,1), c(.01,.25,.75,1,1,1), c(.01,.25,.75,1,1,1),
#                                               c(.5,1,1,1,1,1), c(.5,1,1,1,1,1), c(.5,1,1,1,1,1)),
#                                  fix_pars=list(2:6, 4:6, 2:6, 4:6, 4:6, 4:6, 3:6, 4:6, 2:6)))

  # # fit without ecov
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = "base",                         
  #                             NAA_re = list(cor="iid", sigma="rec"),
  #                             selectivity=list(model=c("age-specific","logistic","age-specific","logistic","logistic","logistic","logistic","logistic","age-specific"),
  #                                initial_pars=list(c(.01,1,1,1,1,1), c(3,3), c(.01,1,1,1,1,1), c(3,3), c(3,3), c(3,3), c(3,3), c(3,3), c(.01,.25,1,1,1,1)),
  #                                fix_pars=list(2:6, NULL, 2:6, NULL, NULL, NULL, NULL, NULL, 3:6))) 


  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # if(m == 10) input$par$logit_selpars[1,7:8] = c(-1.5,-5)
  # if(m == 19) input$par$Ecov_obs_logsigma = matrix(c(rep(-0.7,64),-1.3), ncol=1)
  if(m == 19){
    # input$par$Ecov_obs_logsigma = matrix(c(rep(-0.7,64),-1.3), ncol=1)
    input$par$log_NAA_sigma = -0.7
  }

# input$data$age_comp_model_indices = rep(7, input$data$n_indices)
# input$data$age_comp_model_fleets = rep(7, input$data$n_fleets)
# input$data$n_age_comp_pars_indices = rep(1, input$data$n_indices)
# input$data$n_age_comp_pars_fleets = rep(1, input$data$n_fleets)
# input$par$index_paa_pars = rep(0, input$data$n_indices)
# input$par$catch_paa_pars = rep(0, input$data$n_fleets)
# input$map = input$map[!(names(input$map) %in% c("index_paa_pars", "catch_paa_pars"))]

  # mod <- TMB::MakeADFun(input$data, input$par, DLL = "wham",map = input$map, random=input$random)
  # therep = mod$report()
  # sapply(grep("nll",names(therep),value=T), function(x) sum(therep[[x]]))

  # Fit model with projections:
  #  - 3 years
  #  - F = 0
  #  - WAA, MAA, maturity fixed at terminal year (2011)
  # mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F)  
  mod <- fit_wham(input, do.retro=T, do.osa=F, do.proj=T)  

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=here("results","dat_2019","NAA_GSI",paste0(df.mods$Model[m],".rds")))
}

# if(exists("err")) rm("err") # need to clean this up
# mod$parList$logit_selpars
# mod$sdrep
# TMBhelper::Check_Identifiable(mod)


# # # collect fit models into a list
# # # df.mods <- df.mods[1:8,]
# # mod.list <- here("results",paste0(df.mods$Model,".rds"))
# # # mod.list <- here("results","index4_age_fix36_index1_age",paste0(df.mods$Model,".rds"))
# # mods <- lapply(mod.list, readRDS)

# # # t(sapply(mods, function(x) x$parList$logit_selpars[4,1:6]))

# # 1  2  3  4  5  6  7 10 12 13 19 22
# # 13: last selpars
# # 19: rec_pars, sigma_R
# # 22: first selpars, rec_pars
# m=19
# mods[[m]]$parList$logit_selpars
# mods[[m]]$sdrep
# TMBhelper::Check_Identifiable(mods[[m]])

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

# # round 4
# # age-specific for index 3 with ages 3 & 6 fixed, age-specific for index 1 no ages fixed
# mod.list <- here("results","index4_age_fix36_index1_age",paste0(df.mods$Model,".rds"))
# notconv <- which(!df.mods$pdHess)
# notconv
# # 1  2  3  4  7  8 15 20 24 28 32

# TMBhelper::Check_Identifiable(mods[[1]])
# # bad pars seem like estimating recruitment is the issue, not selectivity
# # 1-4,7-8: high gradients
# # 15,20,28,32: mean_rec_pars
# # 24: mean_rec_pars, logit_q[1], all logit_selpars for index 1

