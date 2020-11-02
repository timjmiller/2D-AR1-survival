# Figure 5
# Adding CPI-Recruitment effect reduces $\rho_R$ and AIC. $\rho_{SSB}$ and $\rho_F$ unchanged.

# fig.height = 3, fig.width = 6

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig5_cpi_mohns_rho.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

df.mods <- data.frame(NAA_cor = rep(c('iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),4),
                      NAA_sigma = rep(c('rec','rec','rec+1','rec+1','rec+1','rec+1'),4),
                      CPI_mod = c(rep("rw",12), rep("ar1",12)),
                      CPI_how = rep(c(rep(0,6),rep(2,6)),2), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","bias_correct_oe_v1","NAA_CPI",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# calc results table
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))
df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=T, do.print=F)$tab)
df.aic$AIC[df.mods$pdHess==FALSE] <- NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

diff <- df.mods %>% filter(NAA_sigma=="rec+1",CPI_mod=='ar1' & CPI_how==2) %>% select(rho_R,rho_SSB,rho_Fbar) %>% abs - df.mods %>% filter(NAA_sigma=="rec+1",CPI_mod=='ar1' & CPI_how==0) %>% select(rho_R,rho_SSB,rho_Fbar) %>% abs
diff$NAA_re <- c("iid","ar1_a","ar1_y","2dar1")
df1 <- diff %>% pivot_longer(-NAA_re, names_to = "var", values_to = "val")
df1$M_re <- "none"

df.mods <- data.frame(NAA_re = rep(c(rep('iid',4),rep('2dar1',4),rep('iid',4)),2),
                      M_re = rep(c(rep('iid',4),rep('iid',4),rep('2dar1',4)),2),
                      est_M = rep(c(FALSE,FALSE,TRUE,TRUE),6),
                      CPI_how = rep(rep(c(0,2),6),2),
                      CPI_mod = c(rep("rw",12),rep("ar1",12)), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- paste0("m",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","bias_correct_oe_v1","NAA_M_CPI",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# calc results table
opt_conv = 1-sapply(mods, function(x) x$opt$convergence)
ok_sdrep = sapply(mods, function(x) if(x$na_sdrep==FALSE & !is.na(x$na_sdrep)) 1 else 0)
df.mods$conv <- as.logical(opt_conv)
df.mods$pdHess <- as.logical(ok_sdrep)
df.mods$runtime <- sapply(mods, function(x) x$runtime)
df.mods$NLL <- sapply(mods, function(x) round(x$opt$objective,3))

	have_peels <- sapply(mods, function(x) !is.null(x$peels))
	df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=F, do.print=F)$tab)
	df.aic$dAIC <- df.aic$rho_Fbar <- df.aic$rho_SSB <- df.aic$rho_R <- NA
	if(any(have_peels)){
		tmp <- as.data.frame(compare_wham_models(mods[have_peels], sort=FALSE, calc.rho=T, do.print=F)$tab)
		df.aic[have_peels==TRUE, 3:5] <- tmp[,3:5]
		df.aic[df.mods$pdHess==FALSE,] <- NA
	}

# df.aic <- as.data.frame(compare_wham_models(mods, sort=FALSE, calc.rho=T)$tab)
# df.aic$AIC[df.mods$pdHess==FALSE] <- NA
minAIC <- min(df.aic$AIC, na.rm=T)
df.aic$dAIC <- round(df.aic$AIC - minAIC,1)
df.mods <- cbind(df.mods, df.aic)
rownames(df.mods) <- NULL

diff <- df.mods %>% filter(NAA_re == "iid", CPI_how==2, CPI_mod=='ar1', pdHess==TRUE) %>% select(rho_R,rho_SSB,rho_Fbar) %>% abs - df.mods %>% filter(NAA_re == "iid", CPI_how==0, CPI_mod=='ar1', pdHess==TRUE) %>% select(rho_R,rho_SSB,rho_Fbar) %>% abs
diff$NAA_re <- "iid"
df2 <- diff %>% pivot_longer(-NAA_re, names_to = "var", values_to = "val")
df2$M_re <- c("iid / 2dar1")

df <- rbind(df1,df2)
df$var <- factor(df$var, levels=c("rho_R","rho_SSB","rho_Fbar"), labels=c("`Mohn's`~rho[R]","`Mohn's`~rho[SSB]","`Mohn's`~rho[F]"))
df$NAA_re <- factor(df$NAA_re, levels=c("iid","ar1_a","ar1_y","2dar1"), labels=c("IID","AR1[age]","AR1[year]","'2D'~AR1"))
df$M_re <- factor(df$M_re, levels=c("none","iid / 2dar1"), labels=c("None","IID / 2D AR1"))

ggplot(df, aes(x=NAA_re, y=val, fill=M_re)) +
  geom_point(pch = 21) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(vars(var), nrow=1, strip.position = "top", labeller = label_parsed) +
  ylab(expression(Reduction~"in"~"|"*"Mohn's"~rho*"|")) +
  xlab("NAA random effect structure") +
  scale_x_discrete(labels = function(l) parse(text=l)) +
  scale_fill_grey(name="M random effects") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))