# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/2a_results_tables_2019.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

# create results table for each set with 2019 data
dir.names <- list.files(here("results","dat_2019"))
for(i in 1:length(dir.names)){
	if(grepl("M", dir.names[i], fixed=TRUE)){ # yes M
		df.mods <- data.frame(NAA_re = rep(c('iid','ar1_a','ar1_y','2dar1'),10),
		                      est_M = rep(c(FALSE,TRUE),each=20),
		                      M_re = rep(rep(c('none','iid','ar1_a','ar1_y','2dar1'),each=4),2), stringsAsFactors=FALSE)
	} else {
		if(grepl("estobserr", dir.names[i], fixed=TRUE)){
			df.mods <- data.frame(NAA_cor = rep(c('iid','ar1_a','ar1_y','2dar1'),8),
			                      NAA_sigma = "rec+1",
			                      GSI_mod = c(rep("rw",16), rep("ar1",16)),
			                      GSI_obserr = "est",
			                      GSI_how = rep(c(rep(0,4),rep(1,4),rep(2,4),rep(4,4)),2), stringsAsFactors=FALSE)

		} else {# no M, yes GSI
			df.mods <- data.frame(NAA_cor = rep(c('iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),4),
			                      NAA_sigma = "rec+1",
			                      GSI_mod = rep(c(rep("rw",8), rep("ar1",8)),2),
			                      GSI_obserr = c(rep("est",16), rep("data",16)),
			                      GSI_how = rep(c(rep(0,4),rep(2,4)),4), stringsAsFactors=FALSE)
		}
	}
	n.mods <- dim(df.mods)[1]
	df.mods$Model <- paste0("m",1:n.mods)
	df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

	# load models
	mod.list <- here("results","dat_2019",dir.names[i],paste0(df.mods$Model,".rds"))
	mods <- lapply(mod.list, readRDS)

	# calc results table
	opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
	ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
	df.mods$conv <- as.logical(opt_conv)
	df.mods$pdHess <- as.logical(ok_sdrep)
	df.mods$runtime <- sapply(mods, function(x) x$runtime)
	df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
	df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)
	df.aic$AIC[df.mods$pdHess==FALSE] <- NA
	minAIC <- min(df.aic$AIC, na.rm=T)
	df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
	df.mods <- cbind(df.mods, df.aic)
	rownames(df.mods) <- NULL

	write.csv(df.mods, file=here("results","tables","dat_2019",paste0(dir.names[i],".csv")))
}

# dAIC boxplot, group by GSI-Rec model
# point: 2D AR1 does best no matter how you model Recruitment/GSI
df <- read.csv(here("tables","dat_2019","NAA_GSI_estobserr.csv"))[,-1]
df$NAA_cor <- factor(as.character(df$NAA_cor), levels=c("iid","ar1_a","ar1_y","2dar1"))
levels(df$NAA_cor) = c("IID", "AR1 (age)", "AR1 (year)", "2D AR1")
df <- df %>% group_by(GSI_mod, GSI_how) %>% mutate(minAIC_group=min(AIC, na.rm=T)) %>% mutate(dAIC_group=AIC-minAIC_group) %>% as.data.frame

png(here("plots","NAA_GSI_dAIC_boxplot.png"), units='in', height=5, width=7, res=300)
ggplot(df, aes(x=NAA_cor,y=dAIC_group)) +
	geom_boxplot(fill='grey') +
	# facet_wrap(~sel_mod, nrow=1, scales="free_x") +
	xlab("NAA random effects") +
	ylab("dAIC") +
	theme_bw()
dev.off()
