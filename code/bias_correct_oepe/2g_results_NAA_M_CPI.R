# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/bias_correct_oepe/2g_results_NAA_M_CPI.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

df.mods <- data.frame(M_re = c('iid','2dar1','iid','2dar1'),
                      NAA_cor = c('iid','iid','2dar1','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("NAA-M-CPI-",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","bias_correct_oepe_rev","NAA_M_CPI",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# calc results table
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))

	have_peels <- sapply(mods, function(x) !is.null(x$peels))
	df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)
	df.aic$dAIC <- df.aic$rho_Fbar <- df.aic$rho_SSB <- df.aic$rho_R <- NA
	if(any(have_peels)){
		tmp <- as.data.frame(compare_wham_models(mods[have_peels], sort=FALSE, calc.rho=T)$tab)
		df.aic[have_peels==TRUE, 3:5] <- tmp[,3:5]
		df.aic[df.mods$pdHess==FALSE,] <- NA
	}

# df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=T)$tab)
# df.aic$AIC[df.mods$pdHess==FALSE] <- NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

write.csv(df.mods, file=here("tables","dat_2019","bias_correct_oepe_rev","NAA_M_CPI.csv"))
