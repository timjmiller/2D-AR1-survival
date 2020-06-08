# Brian Stock
# May 11, 2020
# 2D AR1 survival devs (Haikun's paper)
# Fit models

# | Code | Paper | NAA devs | Recruitment |
# | --- | --- | --- | --- |
# | m1 | S0(base) | IID | Beverton-Holt |
# | m2 | S0(age) | AR1_a | Beverton-Holt |
# | m3 | S0(year) | AR1_y | Beverton-Holt |
# | m4 | S0(age,year) | 2D AR1 | Beverton-Holt |
# | m5 | S(base) | IID | Beverton-Holt + GSI (rw, limiting) |
# | m6 | S(age) | AR1_a | Beverton-Holt + GSI (rw, limiting) |
# | m7 | S(year) | AR1_y | Beverton-Holt + GSI (rw, limiting) |
# | m8 | S(age,year) | 2D AR1 | Beverton-Holt + GSI (rw, limiting) |

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/1_fit_models.R")

# remotes::install_github("noaa-edab/ecodata",build_vignettes=TRUE)
# remotes::install_github("timjmiller/wham", dependencies=TRUE)
library(here)
library(wham) # https://timjmiller.github.io/wham
library(ecodata) # EDAB package, hosts GSI and other env data, https://github.com/NOAA-EDAB/ecodata
library(tidyverse)

# get SNEMAYT data from wham
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=here("data","ex2_SNEMAYT.dat"), overwrite=FALSE)

# created in 0_explore_data.R
gsi <- read.csv(here("data","GSI_v4.csv"), header=TRUE) 

# edit -- don't use wham GSI bc no SE
# use wham GSI, fit obs_err
# gsi <- read.csv(here("data","GSI_wham.csv"), header=TRUE) 

# edit -- don't use ecodata GSI bc it only starts in 1993
# # get GSI data from ecodata
# # for now take annual mean... but should confirm
# gsi.orig <- ecodata::gsi
# gsi.orig$Year <- floor(gsi.orig$Time)
# gsi <- gsi.orig %>% group_by(Year) %>% summarize(GSI = mean(Value)) %>% as.data.frame
# write.csv(gsi, file=here("data","GSI.csv"), row.names=FALSE)

# edit -- don't use Haikun's GSI bc it doesn't have 1972 and has smaller obs_err than likely
# # get GSI from Haikun
# load(here("data","sneyt.Rdata"))
# years <- 1:x$n_years+1972
# gsi <- data.frame(Year=years, GSI=x$Ecov_obs, GSI.se=x$Ecov_obs_sigma)

# set up wham
asap3 <- read_asap3_dat(here("data","ex2_SNEMAYT.dat"))
df.mods <- data.frame(NAA_cor = rep(c('iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),4),
                      NAA_sigma = "rec+1",
                      GSI_mod = rep(c(rep("rw",8), rep("ar1",8)),2),
                      GSI_obserr = c(rep("est",16), rep("data",16)),
                      GSI_how = rep(c(rep(0,4),rep(2,4)),4), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# look at model table
df.mods

thedir <- "index4_age_fix36_index1_age"

# run models
for(m in 1:n.mods){
# for(m in 17:n.mods){
# for(m in 1:8){
  NAA_list <- list(cor=df.mods[m,"NAA_cor"], sigma=df.mods[m,"NAA_sigma"])
  if(NAA_list$sigma == '---') NAA_list = NULL

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
  

  # 1. use logistic selectivity for fleet and indices (as in ASAP file)
  #    mirror index 4 and 5 
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)
  
  # initial model fits did not converge for ar1_y or 2dar1 with logistic selectivity
  # logit_selpars for index 3 (row 4) were often the bad parameters
  # try age-specific selectivity for index 3

# prob should fix sel=1 for ages 3 and 6. start with 6 only.
# t(sapply(mods, function(x) x$parList$logit_selpars[4,1:6]))
#            [,1]       [,2]      [,3]      [,4]        [,5]      [,6]
# [1,] -0.9807071  0.8418645 297.01363 1.3462644  0.26049915 285.93756
# [2,] -0.8367982  1.0656324 231.51645 1.3447626  0.32959172 232.19362
# [3,] -0.5186264 27.1718403  94.94095 0.7773856 -0.02723632  62.69913
# [4,] -0.3801475 59.0730432 164.27797 0.7759044  0.02847899 131.34939
# [5,] -0.9609725  0.8603891 333.79518 1.3575478  0.26979792 415.31448
# [6,] -0.8213652  1.0644483 335.27989 1.3442526  0.33854667 345.23597
# [7,] -0.5214967 25.0453369 342.66289 0.8100341 -0.01660198 187.60379
# [8,] -0.3786363 10.7382261  17.13562 0.8132420  0.04774336  14.33564

  # # 2. age-specific for index 4, no fix
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             selectivity=list(model=c(rep("logistic",3),"age-specific","logistic","logistic"),
  #                                initial_pars=list(c(3,3), c(3,3), c(3,3), c(0.5,0.5,0.5,0.5,0.5,0.5), c(1.5,0.1), c(1.5,0.1)),
  #                                fix_pars=list(NULL, NULL, NULL, NULL, 1:2, 1:2)),
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)

  # # 3. fix age 6 sel = 1 for index 4
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             selectivity=list(model=c(rep("logistic",3),"age-specific","logistic","logistic"),
  #                                initial_pars=list(c(3,3), c(3,3), c(3,3), c(0.5,0.5,0.5,0.5,0.5,1), c(1.5,0.1), c(1.5,0.1)),
  #                                fix_pars=list(NULL, NULL, NULL, 6, 1:2, 1:2)),
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)

  # 4. also fix age 3
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             selectivity=list(model=c(rep("logistic",3),"age-specific","logistic","logistic"),
  #                                initial_pars=list(c(3,3), c(3,3), c(3,3), c(0.5,0.5,1,0.5,0.5,1), c(1.5,0.1), c(1.5,0.1)),
  #                                fix_pars=list(NULL, NULL, NULL, c(3,6), 1:2, 1:2)),
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)

  # 5. age-specific sel for index 1
  input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
                              model_name = df.mods$Model[m],                         
                              selectivity=list(model=c("logistic","age-specific","logistic","age-specific","logistic","logistic"),
               								   initial_pars=list(c(3,3), c(0.5,0.5,0.5,0.5,0.5,0.5), c(3,3), c(0.5,0.5,1,0.5,0.5,1), c(1.5,0.1), c(1.5,0.1)),
               								   fix_pars=list(NULL, NULL, NULL, c(3,6), 1:2, 1:2)),
                              NAA_re = NAA_list,
                              ecov=ecov)

  # cbind(df.mods$pdHess, t(sapply(mods, function(x) x$parList$logit_selpars[4,1:6])))
  # cbind(df.mods$pdHess, t(sapply(mods, function(x) x$parList$logit_selpars[2,1:6])))
  # apply(t(sapply(mods, function(x) x$parList$logit_selpars[2,1:6])), 2, mean)

  # 6. fix age 1 sel = 0 for index 1
  # input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment
  #                             model_name = df.mods$Model[m],                         
  #                             selectivity=list(model=c("logistic","age-specific","logistic","age-specific","logistic","logistic"),
  #                                initial_pars=list(c(3,3), c(0,0.5,0.5,0.5,0.5,0.5), c(3,3), c(0.5,0.5,1,0.5,0.5,1), c(1.5,0.1), c(1.5,0.1)),
  #                                fix_pars=list(NULL, 1, NULL, c(3,6), 1:2, 1:2)),
  #                             NAA_re = NAA_list,
  #                             ecov=ecov)  

  # age comp logistic normal pool obs (not multinomial, the default)
  input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
  input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
  input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
  n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
  n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
  input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
  input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

  # as in WHAM ex2, use selectivity from ASAP file (logistic)
  #   2 pars per block instead of n.ages
  #   sel pars of indices 4/5 fixed at 1.5, 0.1 (neg phase in .dat file)
  input$par$logit_selpars[1:4,7:8] <- 0 # original code started selpars at 0 (last 2 rows are fixed)
  # input$par$mean_rec_pars <- c(25,13) # help m15 converge

  # Fit model with projections:
  #  - 3 years (2012-2014)
  #  - F = 0
  #  - WAA, MAA, maturity fixed at terminal year (2011)
  mod <- fit_wham(input, do.retro=F, do.osa=F, do.proj=F)  
  # mod <- fit_wham(input, do.retro=T, do.osa=T, proj.opts = list(proj.F=rep(0,3), avg.yrs=2011))  

  # Save model
  if(exists("err")) rm("err") # need to clean this up
  saveRDS(mod, file=here("results",thedir,paste0(df.mods$Model[m],".rds")))
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
mod.list <- here("results","index4_age_fix36_index1_age",paste0(df.mods$Model,".rds"))
notconv <- which(!df.mods$pdHess)
notconv
# 1  2  3  4  7  8 15 20 24 28 32

TMBhelper::Check_Identifiable(mods[[1]])
# # bad pars seem like estimating recruitment is the issue, not selectivity
# # 1-4,7-8: high gradients
# # 15,20,28,32: mean_rec_pars
# # 24: mean_rec_pars, logit_q[1], all logit_selpars for index 1

