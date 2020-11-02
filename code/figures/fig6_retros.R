# Figure 6
# Retros for SSB and F (base, NAA, M, NAA + M, NAA + M + CPI)

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig6_retros.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

mod.list <- c(here("results","dat_2019","bias_correct_oepe_rev","NAA","Base.rds"),
              here("results","dat_2019","bias_correct_oepe_rev","NAA","NAA-2.rds"),
              here("results","dat_2019","bias_correct_oepe_rev","NAA","NAA-5.rds"),
              here("results","dat_2019","bias_correct_oepe_rev","M","M-4.rds"),
              here("results","dat_2019","bias_correct_oepe_rev","NAA_M","NAA-M-2.rds"),
              here("results","dat_2019","bias_correct_oepe_rev","NAA_M","NAA-M-3.rds"))
mods <- lapply(mod.list, readRDS)
mod.labs <- c("Base","NAA-2", "NAA-5", "M-4", "NAA-M-2","NAA-M-3")
years = mods[[1]]$years
nyears <- length(years) # don't use projections
npeels = length(mods[[1]]$peels)
nyears.plot <- 20
nyears.cut <- nyears-nyears.plot

dfrho <- data.frame(matrix(NA, nrow=0, ncol=3))
colnames(dfrho) <- c("val","var","Model")
df <- data.frame(matrix(NA, nrow=0, ncol=5))
colnames(df) <- c("Year","peel","val","var","Model")
for(i in 1:length(mods)){
  ssb = list(head(mods[[i]]$rep[["SSB"]],nyears))
  ssb[2:(npeels+1)] = lapply(mods[[i]]$peels, function(x) x$rep[["SSB"]])
  rel.ssb = lapply(1:length(ssb), function(x) ssb[[x]]/ssb[[1]][1:(nyears - x + 1)] - 1)
  tmp <- as.data.frame(t(plyr::ldply(rel.ssb, rbind))[(nyears.cut+1):(nyears-1),])[,-1]
  colnames(tmp) <- paste0("peel",1:npeels)
  tmp$Year <- tail(years, nyears.plot-1)
  tmp <- tmp %>% tidyr::pivot_longer(-Year, names_to = "peel", names_prefix = "peel",
                              names_transform = list(peel = as.integer), values_to = "val")
  tmp$var = "SSB"
  tmp$Model = mod.labs[i]
  df <- rbind(df, tmp)
  
  Fbar = list(head(mods[[i]]$rep[["Fbar"]],nyears))
  Fbar[2:(npeels+1)] = lapply(mods[[i]]$peels, function(x) x$rep[["Fbar"]])
  rel.Fbar = lapply(1:length(Fbar), function(x) Fbar[[x]]/Fbar[[1]][1:(nyears - x + 1)] - 1)  
  tmp <- as.data.frame(t(plyr::ldply(rel.Fbar, rbind))[(nyears.cut+1):(nyears-1),])[,-1]
  colnames(tmp) <- paste0("peel",1:npeels)
  tmp$Year <- tail(years, nyears.plot-1)
  tmp <- tmp %>% tidyr::pivot_longer(-Year, names_to = "peel", names_prefix = "peel",
                              names_transform = list(peel = as.integer), values_to = "val")
  tmp$var = "F"
  tmp$Model = mod.labs[i]
  df <- rbind(df, tmp)
  
  Rec = list(head(mods[[i]]$rep[["NAA"]][,1],nyears))
  Rec[2:(npeels+1)] = lapply(mods[[i]]$peels, function(x) x$rep[["NAA"]][,1])
  rel.Rec = lapply(1:length(Rec), function(x) Rec[[x]]/Rec[[1]][1:(nyears - x + 1)] - 1)
  tmp <- as.data.frame(t(plyr::ldply(rel.Rec, rbind))[(nyears.cut+1):(nyears-1),])[,-1]
  colnames(tmp) <- paste0("peel",1:npeels)
  tmp$Year <- tail(years, nyears.plot-1)
  tmp <- tmp %>% tidyr::pivot_longer(-Year, names_to = "peel", names_prefix = "peel",
                              names_transform = list(peel = as.integer), values_to = "val")
  tmp$var = "Recruitment"
  tmp$Model = mod.labs[i]
  df <- rbind(df, tmp)

  
  dfrho <- rbind(dfrho, data.frame(val=wham::mohns_rho(mods[[i]])[1:3], var=c("SSB","F","Recruitment"), Model=mod.labs[i]))
}

dfrho$Year = as.integer(years[nyears.cut]+1)
dfrho$yval = 1.75
dfrho$val = paste0("rho == ",sprintf("'%.2f'",round(dfrho$val,2)))

df$Model <- factor(df$Model, levels=mod.labs, labels=mod.labs)
dfrho$Model <- factor(dfrho$Model, levels=mod.labs, labels=mod.labs)

df$var <- factor(df$var, levels=c("Recruitment","SSB","F"))
dfpoints <- subset(df, peel == 1:7 & Year == 2018:2012)
df$peel <- as.character(df$peel)
df$Year <- as.integer(df$Year)
dfpoints$peel <- as.character(dfpoints$peel)
dfpoints$Year <- as.integer(dfpoints$Year)

dev.new(width=8, height=5)
ggplot(df, aes(x=Year,y=val, color=peel, group=peel)) +
  geom_line() +
  geom_point(data=dfpoints, aes(x=Year,y=val, color=peel)) +
  geom_label(data=dfrho, aes(x=Year, y=yval, label=val), inherit.aes=F, parse=T, hjust=0, label.r=unit(0, "lines"), label.size=NA, size=3.6) +
  ylab(bquote(paste("Mohn's ", rho))) +
  facet_grid(cols=vars(Model), rows=vars(var)) +
  coord_cartesian(ylim=c(-1,2)) +
  theme_bw() +
  theme(legend.position = "none", axis.text = element_text(size=8))
ggsave(here("plots","fig5_retros.pdf"), device='pdf', width=8, height=5, units="in", dpi = 300)
dev.off()
