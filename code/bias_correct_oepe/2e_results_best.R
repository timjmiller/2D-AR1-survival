# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/bias_correct_oepe/2e_results_best.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

df.mods <- data.frame(NAA_cor = c('iid','iid','2dar1','iid','iid','iid','2dar1'),
                      NAA_sigma = c('rec','rec+1','rec+1','rec','rec','rec+1','rec+1'),
                      M_re = c('none','none','none','iid','2dar1','2dar1','iid'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods$lab <- c("Base","NAA-2","NAA-5","M-1","M-4","NAA-M-2","NAA-M-3")
df.mods <- df.mods %>% select(Model, lab, everything()) # moves Model to first col
df.mods

# load models
# mod.list <- here("results","dat_2019","bias_correct_oepe_rev","best",paste0(df.mods$Model,".rds"))
mod.list <- list("/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/NAA/Base.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/NAA/NAA-2.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/NAA/NAA-5.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/M/M-1.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/M/M-4.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/NAA_M/NAA-M-2.rds",
	"/home/bstock/Documents/ms/2D-AR1-survival/results/revision2/NAA_M/NAA-M-3.rds")
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

write.csv(df.mods, file=here("tables","revision2","best.csv"))
