# figure focusing on specific effect of 2D AR1 structure on F and SSB, both for NAA and M

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/4_fig_2dar1_effect_v2.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)

mod.list <- c(here("results","dat_2019","bias_correct_oe_v2","NAA","NAA-2.rds"),
              here("results","dat_2019","bias_correct_oe_v2","NAA","NAA-5.rds"),
              here("results","dat_2019","bias_correct_oe_v2","NAA_M","NAA-M-3.rds"),
              here("results","dat_2019","bias_correct_oe_v2","M","M-1.rds"),
              here("results","dat_2019","bias_correct_oe_v2","M","M-4.rds"),
              here("results","dat_2019","bias_correct_oe_v2","NAA_M","NAA-M-2.rds"))
mods <- lapply(mod.list, readRDS)
mod.cols <- c("IID","2D AR1", "NAA + M", "IID","2D AR1", "NAA + M")
mod.labs <- c("NAA","NAA","NAA","M","M","M")

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
# df$Model <- factor(df$Model, levels=c("NAA","M","M + mu_M"), labels=c("NAA","M","M~+~mu[M]"))
df$Model <- factor(df$Model, levels=c("NAA","M"), labels=c("NAA","M"))
df$Year <- as.integer(df$Year)
df$Col <- factor(df$Col, levels=c("IID","2D AR1", "NAA + M"), labels=c("IID","2D AR1", "2D AR1 + NAA/M"))
df$lab <- "A"
df$lab[df$Model == "NAA" & df$var=="SSB"] = "B"
df$lab[df$Model == "M" & df$var=="F"] = "C"
df$lab[df$Model == "M" & df$var=="SSB"] = "D"

# dat <- data.table(df)
# dat[,y_min := 0, by = var]
df.labs <- data.frame(label=LETTERS[1:4], Model=factor(c("NAA","NAA","M","M"), levels=c("NAA","M")), var=factor(c("F","SSB","F","SSB"), levels=c("F","SSB")))

# divide by IID within Model
dat <- df %>% group_by(Model) %>% mutate(val.rel=val/val[Col=="IID"]-1) %>% as.data.frame
png(here("plots","fig_2dar1_effect_v2.png"),units='in',width=7,height=4.5,res=300)
print(ggplot(dat, aes(x=Year, y=val.rel, color=Col, group=Col)) +
  # geom_ribbon(aes(ymin=lo, ymax=hi, fill=Model), color=NA, alpha=.15) +
  geom_line(size=.8) +
  geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
  facet_grid(cols=vars(Model), rows=vars(var), labeller=label_parsed) +
  ylab("Relative difference") +
  # geom_blank(aes(y = y_min)) +
  scale_y_continuous(expand=c(0.01,0.01)) +
  coord_cartesian(ylim=c(-0.6,0.6)) +
  scale_color_jco() +
  geom_text(data=df.labs, aes(label=label), x=1974, y=0.53, color='black', size=6, inherit.aes=F) +  # scale_fill_jco() +
  theme_bw() +
  theme(strip.placement = "outside", legend.title=element_blank()))
dev.off()

# average abs(reldiff)
dat$proj <- 0
dat$proj[dat$Year > 2018] <- 1
dat$proj <- factor(dat$proj)
dat %>% group_by(Model, Col, var, proj) %>% summarize(meandiff = round(mean(abs(val.rel)), 2)) %>% as.data.frame
#   Model            Col var proj meandiff
# 5    NAA         2D AR1   F    0     0.12
# 7    NAA         2D AR1 SSB    0     0.09
# 8    NAA         2D AR1 SSB    1     0.22
# 9    NAA 2D AR1 + NAA/M   F    0     0.17
# 11   NAA 2D AR1 + NAA/M SSB    0     0.15
# 12   NAA 2D AR1 + NAA/M SSB    1     0.29
# 17     M         2D AR1   F    0     0.10
# 19     M         2D AR1 SSB    0     0.07
# 20     M         2D AR1 SSB    1     0.22
# 21     M 2D AR1 + NAA/M   F    0     0.14
# 23     M 2D AR1 + NAA/M SSB    0     0.09
# 24     M 2D AR1 + NAA/M SSB    1     0.49

# SSB in model years: .09, .15, .07, .09 / range 7-15% / mean 10%
# SSB in proj years: .22, .29, .22, .49 / range 22-49% / mean 31%

# F in model years: .12, .17, .10, .14 / range 10-17% / mean 13%

# mean(c(.09,.15,.07,.09))
# [1] 0.1
# mean(c(.22,.29,.22,.49))
# [1] 0.305
# mean(c(.12,.17,.1,.14))
# [1] 0.1325
