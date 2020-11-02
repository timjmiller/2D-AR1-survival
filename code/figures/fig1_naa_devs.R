# Figure 1
# NAA deviations
# Oct 30 2020 (revision)

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/figures/fig1_naa_devs.R")

library(here)
library(wham)
library(tidyverse)
library(ggsci)
library(data.table)

df.mods <- data.frame(NAA_cor = c('iid','ar1_y','iid','ar1_a','ar1_y','2dar1'),
                      NAA_sigma = c('rec','rec','rec+1','rec+1','rec+1','rec+1'), stringsAsFactors=FALSE)
n.mods <- dim(df.mods)[1]
df.mods$Model <- c("Base",paste0("NAA-",1:(n.mods-1)))
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col

# load models
mod.list <- here("results","dat_2019","bias_correct_oepe_rev","NAA",paste0(df.mods$Model,".rds"))
mods <- lapply(mod.list, readRDS)

# df.mods$Model <- paste0("NAA-",1:n.mods)
df.mods$NAA_cor <- factor(df.mods$NAA_cor, levels=c('iid','ar1_a','ar1_y','2dar1'), labels=c("Indep.","AR(1)[age]","AR(1)[year]","'2D'~AR(1)"))
df.mods$NAA_sigma <- factor(df.mods$NAA_sigma, levels=c('rec','rec+1'), labels=c("Random~effects:~Recruitment","Random~effects:~all~NAA"))
n_ages <- mods[[1]]$env$data$n_ages
df.NAA <- data.frame(matrix(NA, nrow=0, ncol=n_ages+3))
colnames(df.NAA) <- c(paste0("Age_",1:n_ages),"Year","NAA_cor","NAA_sigma")
for(i in 1:length(mods)){
  tmp = as.data.frame(mods[[i]]$rep$NAA_devs)[1:(length(mods[[i]]$years_full)-1),]
  tmp$Year <- mods[[i]]$years_full[-1]
  colnames(tmp) <- c(paste0("Age_",1:n_ages),"Year")
  tmp$NAA_cor = df.mods$NAA_cor[i]
  tmp$NAA_sigma = df.mods$NAA_sigma[i]
  df.NAA <- rbind(df.NAA, tmp)
}
df.plot <- df.NAA %>% tidyr::pivot_longer(-c(Year,NAA_cor,NAA_sigma),
                                          names_to = "Age",
                                          names_prefix = "Age_",
                                          names_transform = list(Age = as.integer),
                                          values_to = "NAA_re")
df.plot$NAA_re[df.plot$NAA_sigma == "Random~effects:~Recruitment" & df.plot$Age > 1] = 0
df.mods$xint = tail(mods[[1]]$years,1)

dev.new(width=7, height=6)
ggplot(df.plot, ggplot2::aes(x=Year, y=Age)) +
        geom_tile(aes(fill=NAA_re)) +
        geom_label(aes(x=Year, y=Age, label=lab), size=4, alpha=1, label.r=unit(0, "lines"), label.size=NA,
                   data=data.frame(Year=1978, Age=5.7, lab=df.mods$Model, NAA_cor=df.mods$NAA_cor, NAA_sigma=df.mods$NAA_sigma)) +
        geom_vline(data=df.mods, mapping=aes(xintercept = xint), linetype=2, size=.6) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        theme_bw() +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        facet_grid(rows=vars(NAA_cor), cols=vars(NAA_sigma), drop=F, labeller = label_parsed) +
        scale_fill_gradient2(name = "", low = scales::muted("blue"), mid = "white", high = scales::muted("red"))
ggsave(here("plots","fig1_naa_devs.pdf"), device='pdf', width=7, height=6, units="in", dpi = 300)
dev.off()
