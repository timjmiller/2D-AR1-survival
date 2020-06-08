# Brian Stock
# May 18, 2020
# 2D AR1 survival devs (Haikun's paper)
# Model results table - which converged?

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/2_results_table.R")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(ecodata) # EDAB package, hosts GSI and other env data, https://github.com/NOAA-EDAB/ecodata
library(tidyverse)

# create results table for each selectivity run
dir.names <- head(list.files(here("results")),6)
sel_ind3 <- c("Age-specific (3,6)", "Age-specific (3,6)", "Age-specific (3,6)", "Age-specific (6)", "Age-specific", "Logistic")
sel_ind1 <- c("Logistic", "Age-specific", "Age-specific (1)", "Logistic", "Logistic", "Logistic", "Logistic")
n.conv <- rep(NA, 6)

# df.mods <- read.csv(text="Model,NAA_cor,NAA_sigma,GSI_mod,GSI_obserr,GSI_how,sel_ind1,sel_ind3,conv,pdHess,runtime,NLL,dAIC,AIC")
df.mods <- read.csv(text="Model,NAA_cor,NAA_sigma,GSI_mod,GSI_obserr,GSI_how,sel_ind1,sel_ind3,pdHess,runtime,NLL,dAIC,AIC")
for(i in 1:length(dir.names)){
	# recreate df.mods
	tmp <- data.frame(NAA_cor = rep(c('iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),4),
	                      NAA_sigma = "rec+1",
	                      GSI_mod = rep(c(rep("rw",8), rep("ar1",8)),2),
	                      GSI_obserr = c(rep("est",16), rep("data",16)),
	                      GSI_how = rep(c(rep(0,4),rep(2,4)),4), stringsAsFactors=FALSE)
	n.mods <- dim(tmp)[1]
	tmp$Model <- paste0("m",1:n.mods)
	tmp <- tmp %>% select(Model, everything()) # moves Model to first col
	tmp$sel_ind1 <- sel_ind1[i]
	tmp$sel_ind3 <- sel_ind3[i]

	# load models
	mod.list <- here("results",dir.names[i],paste0(tmp$Model,".rds"))
	mods <- lapply(mod.list, readRDS)

	# calc results table
	opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
	ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
	# tmp$conv <- as.logical(opt_conv)
	tmp$pdHess <- as.logical(ok_sdrep)
	tmp$runtime <- sapply(mods, function(x) x$runtime)
	tmp$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
	df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F)$tab)
	df.aic$AIC[tmp$pdHess==FALSE] <- NA
	minAIC <- min(df.aic$AIC, na.rm=T)
	df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
	tmp <- cbind(tmp, df.aic)
	rownames(tmp) <- NULL

	n.conv[i] <- sum(tmp$pdHess)
	df.mods <- rbind(df.mods, tmp)
	# write results table
	# write.csv(tmp, file=here("results",paste0(dir.names[i],".csv")))
}

# "index4_age_fix36_index1_age" has the most converged models (21/32)
#   index 1 = age-specific
#   index 3 = age-specific (3,6)

# data.frame(file=dir.names, n.conv=n.conv)
#                               file n.conv
# 1                 index4_age_fix36     19
# 2      index4_age_fix36_index1_age     21
# 3 index4_age_fix36_index1_age_fix1     14
# 4                  index4_age_fix6     14
# 5                 index4_age_nofix     15
# 6                  index4_logistic     16

# minAIC approach = for each model, use whichever selectivity that converges and minimizes AIC
df.mods <- df.mods %>% group_by(Model) %>% mutate(minAIC = min(AIC, na.rm=T)) %>% as.data.frame 
df.mods$minAIC[df.mods$minAIC > -1500] <- NA
df <- df.mods[which(df.mods$AIC == df.mods$minAIC),]

notconv <- paste0("m",1:n.mods)[!paste0("m",1:n.mods) %in% df$Model]
df.res <- data.frame(Model = paste0("m",1:n.mods),
					  NAA_cor = rep(c('iid','ar1_a','ar1_y','2dar1','iid','ar1_a','ar1_y','2dar1'),4),
                      NAA_sigma = "rec+1",
                      GSI_mod = rep(c(rep("rw",8), rep("ar1",8)),2),
                      GSI_obserr = c(rep("est",16), rep("data",16)),
                      GSI_how = rep(c(rep(0,4),rep(2,4)),4),
                      sel_ind1 = NA,
                      sel_ind3 = NA, 
                      pdHess = FALSE,
                      runtime = NA, NLL = NA, dAIC=NA, AIC=NA, minAIC=NA, stringsAsFactors=FALSE)
for(m in 1:n.mods){
	if(paste0("m",m) %in% df$Model){
		ind <- which(df$Model == paste0("m",m))
		df.res[m,] <- df[ind,]
	}
}
df.res <- df.res[,-which(colnames(df.res)=="AIC")]
lowAIC <- min(df.res$minAIC, na.rm=T)
df.res$dAIC <- round(df.res$minAIC - lowAIC,1)

write.csv(df.res, file=here("results","tables","minAIC.csv"), row.names=FALSE)
