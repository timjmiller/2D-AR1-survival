# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/bias_correct_oepe/2_results_NAA.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

df.mods <- data.frame(NAA_cor = c('iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('rec','rec','rec+1','rec+1','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("NAA-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","bias_correct_oepe_rev","NAA",paste0(df.mods$Model,".rds"))
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

write.csv(df.mods, file=here("tables","dat_2019","bias_correct_oepe_rev","NAA.csv"))
