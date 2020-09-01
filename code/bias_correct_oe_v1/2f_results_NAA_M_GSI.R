# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/2f_results_NAA_M_GSI.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

df.mods <- data.frame(NAA_re = c(rep('2dar1',4),rep('iid',4)),
                      M_re = c(rep('iid',4),rep('2dar1',4)),
                      est_M = c(FALSE,FALSE,TRUE,TRUE,FALSE,FALSE,TRUE,TRUE),
                      GSI_how = rep(c(0,2),4), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","NAA_M_GSI_ar1",paste0(df.mods$Model,".rds"))
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

write.csv(df.mods, file=here("tables","dat_2019","NAA_M_GSI_ar1.csv"))

# # dAIC boxplot, group by GSI-Rec model
# # point: 2D AR1 does best no matter how you model Recruitment/GSI
# df <- read.csv(here("tables","dat_2019","NAA.csv"), stringsAsFactors=F)[,-1]
# df$NAA_mod <- factor(paste(df$NAA_sigma,df$NAA_cor), levels=c("rec iid","rec ar1_y","rec+1 iid","rec+1 ar1_a","rec+1 ar1_y","rec+1 2dar1"))
# levels(df$NAA_mod) = c("Rec IID","Rec AR1 (year)","NAA IID", "NAA AR1 (age)", "NAA AR1 (year)", "NAA 2D AR1")
# df <- df %>% group_by(GSI_mod, GSI_how) %>% mutate(minAIC_group=min(AIC, na.rm=T)) %>% mutate(dAIC_group=AIC-minAIC_group) %>% as.data.frame

# png(here("plots","NAA_GSI_dAIC_boxplot.png"), units='in', height=5, width=7, res=300)
# ggplot(df, aes(x=NAA_mod,y=dAIC_group)) +
# 	geom_boxplot(fill='grey') +
# 	# facet_wrap(~sel_mod, nrow=1, scales="free_x") +
# 	xlab("NAA random effects") +
# 	ylab("dAIC") +
# 	theme_bw()
# dev.off()
