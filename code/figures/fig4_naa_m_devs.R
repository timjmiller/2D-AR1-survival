# Figure 4
# NAA and M devs in final model (NAA + M + CPI random effects)
# Only showing NAA + M + CPI model with lowest rho (< 0.1 for all three), even though it has higher AIC. Could show all 3 in a multipanel

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig4_naa_m_devs.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

mod <- readRDS(here("results","dat_2019","bias_correct_oepe_rev","NAA_M","NAA-M-2.rds"))
n_ages <- mod$env$data$n_ages
df <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
colnames(df) <- c(paste0("Age_",1:n_ages),"Year","type")
tmp = as.data.frame(mod$rep$M_re)[1:length(mod$years_full),]
tmp$Year <- mod$years_full
colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
tmp$type = "log(M) deviations"
df <- rbind(df, tmp)

# add NAA devs
tmp = as.data.frame(mod$rep$NAA_devs)[1:(length(mod$years_full)-1),]
tmp$Year <- mod$years_full[-1]
colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
tmp$type = "log(NAA) deviations"
df <- rbind(df, tmp)

df.plot <- df %>% tidyr::pivot_longer(-c(Year,type),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "devs")

dev.new(width=7, height=7)
ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
      geom_tile(aes(fill=devs)) +
      geom_vline(xintercept = tail(mod$years,1), linetype=2, size=.4) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_wrap(vars(type), nrow=2, ncol=1) +
      scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))
ggsave(here("plots","fig4_naa_m_devs.pdf"), device='pdf', width=7, height=7, units="in", dpi = 300)
dev.off()

