# Figure 7
# SSB and F (base, NAA, M, NAA + M, NAA + M + CPI)
# shows impact of model choice on stock status

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig7_stock_status.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

# fig.height = 5, fig.width = 7}
mod.list <- c(here("results","revision2","NAA","Base.rds"),
              here("results","revision2","NAA","NAA-2.rds"),
              here("results","revision2","NAA","NAA-5.rds"),
              here("results","revision2","M","M-4.rds"),
              here("results","revision2","NAA_M","NAA-M-2.rds"),
              here("results","revision2","NAA_M","NAA-M-3.rds"))
mods <- lapply(mod.list, readRDS)
mod.labs <- c("Base","NAA-2", "NAA-5", "M-4", "NAA-M-2","NAA-M-3")

years = mods[[1]]$years_full # include projections
alpha = .05 # 95% CI
df <- data.frame(matrix(NA, nrow=0, ncol=6))
colnames(df) <- c("Year","var","val","lo","hi","Model")
for(i in 1:length(mods)){
  std = summary(mods[[i]]$sdrep)
  ssb.ind <- which(rownames(std) == "log_SSB")
	log.ssb <- std[ssb.ind,1]
  ssb = exp(log.ssb)/1000
	ssb.cv <- std[ssb.ind,2]
  log.ssb.ci <- log.ssb + cbind(qnorm(1-alpha/2)*ssb.cv, -qnorm(1-alpha/2)*ssb.cv)
  ssb.ci = exp(log.ssb.ci)/1000
  df <- rbind(df, data.frame(Year=years, var="SSB", val=ssb, lo=ssb.ci[,1], hi=ssb.ci[,2], Model=mod.labs[i]))

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
  df <- rbind(df, data.frame(Year=years, var="F", val=full.f, lo=exp(log.f.ci[,1]), hi=exp(log.f.ci[,2]), Model=mod.labs[i]))
}
df$Model <- factor(df$Model, levels=mod.labs, labels=mod.labs)
df$Year <- as.integer(df$Year)
# trim ribbon bc y-axis limit can't be changed
ymax = 8
df$hi[df$hi > ymax] = ymax
df$val[df$val > ymax] = ymax
df$lo[df$lo > ymax] = ymax
df$var[df$var =="SSB"] = "SSB (x 1000 mt)"
dat <- data.table(df)
dat[,y_min := 0, by = var]

# exclude years before 1993 in dat bc ggplot doesn't recalculate y-axis limits
dat <- subset(dat, Year > 1992)

dev.new(width=7, height=5)
ggplot(dat, aes(x=Year, y=val, color=Model, group=Model)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill=Model), color=NA, alpha=.15) +
  geom_line(size=.8) +
  geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
  facet_wrap(vars(var), scales="free_y", ncol=1, strip.position = "left") +
  # facet_wrap(vars(var), scales="free", ncol=1, strip.position = "left") +
  ylab(NULL) +
  geom_blank(aes(y = y_min)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  coord_cartesian(xlim=c(1995,2021)) +
  scale_color_jco() +
  scale_fill_jco() +
  theme_bw() +
  theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(here("plots","revision2","fig6_stock_status.pdf"), device='pdf', width=7, height=5, units="in", dpi = 300)
dev.off()

