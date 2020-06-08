# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/2b_results_2019.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(tidyverse)

# create results table for each set with 2019 data
# get all M + GSI runs
dir.names <- grep("M_GSI", list.files(here("results","dat_2019")), fixed=TRUE, value=TRUE)
df <- read.csv(text="Model,NAA_re,est_M,M_re,conv,pdHess,runtime,NLL,dAIC,AIC,rho_Fbar,rho_SSB,rho_R,ver")
for(i in 1:length(dir.names)){
	df.mods <- data.frame(NAA_re = rep(c('iid','ar1_a','ar1_y','2dar1'),10),
	                      est_M = rep(c(FALSE,TRUE),each=20),
	                      M_re = rep(rep(c('none','iid','ar1_a','ar1_y','2dar1'),each=4),2), stringsAsFactors=FALSE)
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
	# df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)

	have_peels <- sapply(mods, function(x) !is.null(x$peels))
	df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)
	df.aic$dAIC <- df.aic$rho_R <- df.aic$rho_SSB <- df.aic$rho_Fbar <- NA
	if(any(have_peels)){
		tmp <- as.data.frame(compare_wham_models(mods[have_peels], sort=FALSE, calc.rho=T)$tab)
		df.aic[have_peels==TRUE, 3:5] <- tmp[,3:5]
	}
	# minAIC <- min(df.aic$AIC, na.rm=T)
	# df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
	df.aic$AIC[df.mods$pdHess==FALSE] <- NA	
	df.mods <- cbind(df.mods, df.aic)
	rownames(df.mods) <- NULL
	df.mods$ver <- i

	df <- rbind(df, df.mods)
	# write.csv(df.mods, file=here("results","tables","dat_2019",paste0(dir.names[i],".csv")))
}

df$Model <- factor(as.character(df$Model))
df$NAA_re <- factor(as.character(df$NAA_re), levels=c("iid","ar1_a","ar1_y","2dar1"))
levels(df$NAA_re) = c("IID", "AR1 (age)", "AR1 (year)", "2D AR1")
df$est_M <- factor(df$est_M, levels=c(TRUE,FALSE))
levels(df$est_M) = c("Estimate mean M", "Fix age-specific M")
df$M_re <- factor(df$M_re, levels=c("none","iid","ar1_a","ar1_y","2dar1"))
levels(df$M_re) = c("None","IID", "AR1 (age)", "AR1 (year)", "2D AR1")

df.best <- df[1:n.mods,]
for(m in 1:n.mods){
	tmp <- df[df$Model == paste0("m",m),]
	ret <- tmp[5,]
	if(is.na(ret$rho_Fbar) & any(!is.na(tmp$rho_Fbar))){
		notNA <- tail(which(!is.na(tmp$rho_Fbar)),1)
		ret[,11:13] <- tmp[notNA,11:13]
	}
	if(is.na(ret$AIC) & any(!is.na(tmp$AIC))){
		lowAIC <- max(tmp$ver[which(tmp$AIC == min(tmp$AIC, na.rm=T))])
		ret[,c("conv","pdHess","runtime","NLL","AIC","ver")] <- tmp[lowAIC,c("conv","pdHess","runtime","NLL","AIC","ver")]
	}
	ret$NLL <- min(tmp$NLL[tmp$conv & tmp$pdHess], na.rm=T)
	df.best[m,] <- ret
}
minAIC <- min(df.best$AIC, na.rm=T)
df.best$dAIC <- round(df.best$AIC - minAIC,1)

png(here("plots","M_GSI_dAIC.png"), units='in', height=4, width=7, res=300)
ggplot(df.best, aes(x=NAA_re,y=M_re,fill=dAIC)) +
	geom_tile() +
	# scale_fill_gradient2(na.value = "grey80") +
	scale_fill_viridis(na.value = "grey80", direction=-1) +
	facet_wrap(~est_M, nrow=1) + # scales="free_x"
	xlab("NAA_re") +
	ylab("M_re") +
	theme_bw()
dev.off()

png(here("plots","M_GSI_rho_Fbar.png"), units='in', height=4, width=7, res=300)
ggplot(df.best, aes(x=NAA_re,y=M_re,fill=rho_Fbar)) +
	geom_tile() +
	scale_fill_gradient2(midpoint=0, na.value = "grey80", low = scales::muted("blue"),high = scales::muted("red")) +
	facet_wrap(~est_M, nrow=1) + 
	xlab("NAA_re") +
	ylab("M_re") +
	theme_bw()
dev.off()

png(here("plots","M_GSI_rho_SSB.png"), units='in', height=4, width=7, res=300)
ggplot(df.best, aes(x=NAA_re,y=M_re,fill=rho_SSB)) +
	geom_tile() +
	scale_fill_gradient2(midpoint=0, na.value = "grey80", low = scales::muted("blue"),high = scales::muted("red")) +
	facet_wrap(~est_M, nrow=1) + 
	xlab("NAA_re") +
	ylab("M_re") +
	theme_bw()
dev.off()

png(here("plots","M_GSI_rho_R.png"), units='in', height=4, width=7, res=300)
ggplot(df.best, aes(x=NAA_re,y=M_re,fill=rho_R)) +
	geom_tile() +
	scale_fill_gradient2(midpoint=0, na.value = "grey80", low = scales::muted("blue"),high = scales::muted("red")) +
	facet_wrap(~est_M, nrow=1) + 
	xlab("NAA_re") +
	ylab("M_re") +
	theme_bw()
dev.off()

write.csv(df.best, file=here("tables","dat_2019","NAA_M_GSI_best.csv"))
