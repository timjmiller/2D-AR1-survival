# figure focusing on specific effect of 2D AR1 structure on F and SSB, both for NAA and M

library(here)
library(wham)
library(tidyverse)
library(ggsci)

mod.list <- c(here("results","dat_2019","bias_correct_oe","NAA","m3.rds"),
              here("results","dat_2019","bias_correct_oe","NAA","m6.rds"),
              here("results","dat_2019","bias_correct_oe","NAA_M","m2.rds"),
              here("results","dat_2019","bias_correct_oe","M","m2.rds"),
              here("results","dat_2019","bias_correct_oe","M","m5.rds"),
              here("results","dat_2019","bias_correct_oe","NAA_M_iid","m5.rds"),
              here("results","dat_2019","bias_correct_oe","M","m7.rds"),
              here("results","dat_2019","bias_correct_oe","M","m10.rds"),
              here("results","dat_2019","bias_correct_oe","NAA_M_iid","m10.rds"))
mods <- lapply(mod.list, readRDS)
mod.cols <- c("IID","2D AR1", "NAA + M", "IID","2D AR1", "NAA + M", "IID","2D AR1", "NAA + M")
mod.labs <- c("NAA","NAA","NAA","M","M","M","M + mu_M", "M + mu_M", "M + mu_M")

years <- mods[[1]]$years_full # include projections
ny <- length(years)
ny.plot <- ny
ny.cut <- ny-ny.plot
na <- mods[[1]]$env$data$n_ages
alpha <- .05 # 95% CI
df <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(df) <- c("Year","var","val","lo","hi","Model","Col")
for(i in 1:length(mods)){
  std = summary(mods[[i]]$sdrep)
  ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
  ssb = exp(log.ssb)/1000
	ssb.cv <- std[ssb.ind,2]
  log.ssb.ci <- log.ssb + cbind(qnorm(1-alpha/2)*ssb.cv, -qnorm(1-alpha/2)*ssb.cv)
  ssb.ci = exp(log.ssb.ci)/1000
  df <- rbind(df, data.frame(Year=years, var="SSB", val=ssb, lo=ssb.ci[,1], hi=ssb.ci[,2], Model=mod.labs[i], Col=mod.cols[i]))

  n_ages = mods[[i]]$env$data$n_ages
	faa.ind <- which(rownames(std) == "log_FAA_tot")
	log.faa <- matrix(std[faa.ind,1], length(years), n_ages)
	faa.cv <- matrix(std[faa.ind,2], length(years), n_ages)
	age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
  full.f.ind = cbind(1:length(years), age.full.f)
  log.full.f <- log.faa[full.f.ind]
  full.f.cv <- faa.cv[full.f.ind]
  log.f.ci <- log.full.f + cbind(qnorm(1-alpha/2)*full.f.cv, -qnorm(1-alpha/2)*full.f.cv)
  full.f = exp(log.full.f)
  df <- rbind(df, data.frame(Year=years, var="F", val=full.f, lo=exp(log.f.ci[,1]), hi=exp(log.f.ci[,2]), Model=mod.labs[i], Col=mod.cols[i]))
}
df$Model <- factor(df$Model, levels=c("NAA","M","M + mu_M"), labels=c("NAA","M","M~+~mu[M]"))
df$Year <- as.integer(df$Year)
df$Col <- factor(df$Col, levels=c("IID","2D AR1", "NAA + M"), labels=c("IID","2D AR1", "2D AR1 + NAA/M"))
# dat <- data.table(df)
# dat[,y_min := 0, by = var]

# divide by IID within Model
dat <- df %>% group_by(Model) %>% mutate(val.rel=val/val[Col=="IID"]-1) %>% as.data.frame

png(here("plots","fig_2dar1_effect.png"),units='in',width=8.5,height=4.5,res=300)
ggplot(dat, aes(x=Year, y=val.rel, color=Col, group=Col)) +
  # geom_ribbon(aes(ymin=lo, ymax=hi, fill=Model), color=NA, alpha=.15) +
  geom_line(size=.8) +
  geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
  facet_grid(cols=vars(Model), rows=vars(var), labeller=label_parsed) +
  ylab("Relative difference") +
  # geom_blank(aes(y = y_min)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_color_jco() +
  # scale_fill_jco() +
  theme_bw() +
  theme(strip.placement = "outside", legend.title=element_blank())
dev.off()


