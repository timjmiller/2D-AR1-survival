# Figure 3
# M devs (in models where only M is random effect)

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig3_m_devs.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

df.mods <- data.frame(M_re = c('none','iid','ar1_a','ar1_y','2dar1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("M-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","revision2","M",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# df.mods$Model <- paste0("M-",1:n.mods)
df.mods$M_re <- factor(df.mods$M_re, levels=c('none','iid','ar1_a','ar1_y','2dar1'), labels=c("None","Indep.","AR(1)[age]","AR(1)[year]","'2D'~AR(1)"))
# df.mods$est_M <- factor(df.mods$est_M, levels=c(FALSE,TRUE), labels=c("mean M fixed","mean M estimated"))
# df.mods$est_M <- factor(df.mods$est_M, levels=c(FALSE,TRUE), labels=c("Fixed~mu[M]","Estimated~mu[M]"))
n_ages <- mods[[1]]$env$data$n_ages
df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+2))
# df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
# colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","M_re","est_M")
colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","M_re")
for(i in 1:length(mods)){
  tmp = as.data.frame(mods[[i]]$rep$M_re)[1:length(mods[[i]]$years_full),]
  tmp$Year <- mods[[i]]$years_full
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$M_re = df.mods$M_re[i]
  # tmp$est_M = df.mods$est_M[i]
  if(i %in% c(1,3,6,8)){
    for(y in (mods[[i]]$env$data$n_years_model+1):length(mods[[i]]$years_full)){
      tmp[y,1:6] <- tmp[length(mods[[i]]$years), 1:6]
    }
  }
  df.NAA <- rbind(df.NAA, tmp)
}
# df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,M_re,est_M),
df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,M_re),
          names_to = "Age",
          names_prefix = "Age_",
          names_transform = list(Age = as.integer),
          values_to = "M_devs")
# df.plot$M_re[df.plot$est_M == "RE: Recruit" & df.plot$Age > 1] = 0

dev.new(width=4, height=6)
ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
      geom_tile(aes(fill=M_devs)) +
      geom_vline(xintercept = tail(mods[[1]]$years,1), linetype=2, size=.4) +
      geom_label(aes(x=Year, y=Age, label=lab), size=4, alpha=1, label.r=unit(0, "lines"), label.size=NA,
          data=data.frame(Year=1976.5, Age=5.7, lab=df.mods$Model, M_re=df.mods$M_re)) +
      scale_x_continuous(expand=c(0,0)) +
      scale_y_continuous(expand=c(0,0)) +
      theme_bw() +
      facet_wrap(vars(M_re), ncol=1, labeller = label_parsed, strip.position="right") +
      scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))
ggsave(here("plots","revision2","fig3_m_devs.pdf"), device='pdf', width=4, height=6, units="in", dpi = 300)
dev.off()


