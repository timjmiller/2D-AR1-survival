# Figure 4
# NAA and M devs in final model (NAA + M + CPI random effects)
# Only showing NAA + M + CPI model with lowest rho (< 0.1 for all three), even though it has higher AIC. Could show all 3 in a multipanel

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig4_naa_m_devs_2mods.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

# mod <- readRDS(here("results","dat_2019","bias_correct_oepe_rev","NAA_M","NAA-M-2.rds"))
mod.list <- here("results","dat_2019","bias_correct_oepe_rev","NAA_M",paste0("NAA-M-",2:3,".rds"))
mods <- lapply(mod.list, readRDS)

n_ages <- mods[[1]]$env$data$n_ages
df <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
colnames(df) <- c(paste0("Age_",1:n_ages),"Year","type","mod")
for(i in 1:length(mods)){
  tmp = as.data.frame(mods[[i]]$rep$NAA_devs)[1:(length(mods[[i]]$years_full)-1),]
  tmp$Year <- mods[[i]]$years_full[-1]
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$type = "log(NAA) deviations"
  tmp$mod = c("NAA-M-2","NAA-M-3")[i]
  df <- rbind(df, tmp)

  tmp = as.data.frame(mods[[i]]$rep$M_re)[1:length(mods[[i]]$years_full),]
  tmp$Year <- mods[[i]]$years_full
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$type = "log(M) deviations"
  tmp$mod = c("NAA-M-2","NAA-M-3")[i]
  df <- rbind(df, tmp)
}

df.plot <- df %>% tidyr::pivot_longer(-c(Year,type,mod),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "devs")

dev.new(width=7, height=4)
ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
      geom_tile(aes(fill=devs)) +
      geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_grid(rows=vars(type), cols=vars(mod)) +
      scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))
ggsave(here("plots","fig4_naa_m_devs_2mods.pdf"), device='pdf', width=7, height=4, units="in", dpi = 300)
dev.off()

