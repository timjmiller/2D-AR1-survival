# Brian Stock
# Simulation test WHAM

# source("/home/bstock/Documents/ms/2D-AR1-survival/code/supp_all.R")

setwd("/home/bstock/Documents/ms/2D-AR1-survival")
res_dir <- "/home/bstock/Documents/ms/wham-sim/results"
simdata_dir <- "/home/bstock/Documents/ms/wham-sim/data/simdata"
plots_dir <- "/home/bstock/Documents/ms/2D-AR1-survival/plots/supp"
ids = c("SNEMAYT","SNEMAYT")
re = c("NAA","M")
bc.type = 2 # _oepe

# install.packages("ggplotFL", repos="http://flr-project.org/R")
# devtools::install_github("timjmiller/wham", dependencies=TRUE)
library(wham)
library(tidyverse)
library(ggplotFL)
library(ggsci)
library(cowplot)
library(data.table)
inv.rho.trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
rho_trans <- function(x) return(2/(1 + exp(-2*x)) - 1)
source("code/figures/supp_get_results.R")
source("code/figures/supp_plot_rel_err.R")
source("code/figures/supp_plot_rel_err_pars.R")

# ---------------------------------------------------------
# Relative error plots
for(j in 1:length(ids)){
  results <- get_results(stock.id=ids[j], re=re[j], bc.type=bc.type, res_dir=res_dir, simdata_dir=simdata_dir)
  plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=2, multipanel=TRUE, plots_dir=plots_dir)
  plot_rel_err(results, stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=2, multipanel=FALSE, plot.eps=FALSE, plots_dir=plots_dir)
  plot_rel_err_pars(stock.id=ids[j], re=re[j], bc.type=bc.type, sim.types=2, plot.eps=FALSE, res_dir=res_dir, simdata_dir=simdata_dir, plots_dir=plots_dir)
}

