# Get results for one stock

get_results <- function(stock.id="SNEMAYT", re="NAA", bc.type=2, 
                        res_dir=file.path(getwd(),"results"), 
                        simdata_dir=file.path(getwd(),"data","simdata")){ 
  id <- paste0(stock.id,"_",re)
  if(bc.type == 1){
    bc <- "bias_correct_oe"
    # plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oe")    
  }
  if(bc.type == 2){
    bc <- "bias_correct_oepe"
    # plots_dir <- file.path(plots_dir, bc, id)
    id <- paste0(id,"_oepe") 
  }
  # dir.create(plots_dir, showWarnings = FALSE)
  res_dir <- file.path(res_dir, bc, id)
  res.files <- list.files(path=res_dir, pattern = "results", full.names = TRUE)
  res.list <- lapply(res.files, readRDS)
  mod.files <- list.files(path=res_dir, pattern = "^m..rds", full.names = TRUE)
  n.mods <- length(mod.files)
  n.sim <- length(res.list[[1]][[1]])
  flatten.nested.list <- function(X) if(is.list(X)) Reduce(c, lapply(X, flatten.nested.list)) else list(X)
  results <- do.call(rbind, flatten.nested.list(res.list)) %>% as.data.frame
  results <- sapply(results, as.numeric)
  results <- as.data.frame(results)
  types <- c("OE","OEPE") # simulation type, not bias correction type
  tylabs = c("Simulated data: Obs error", "Simulated data: Obs + Process error (new NAA)")
  results$type <- factor(results$type, levels=1:2, labels=tylabs)
  if(re == "NAA" & n.mods == 4){
    mlabs = c("Base","NAA-1","NAA-2","NAA-5")
  }
  if(n.mods == 3){ # M or selectivity
    mlabs = c("Base","M-1","M-4")
  }
  results$om <- factor(results$om, levels=1:n.mods, labels=mlabs)
  results$em <- factor(results$em, levels=1:n.mods, labels=mlabs)
  # results$em.x <- fct_recode(results$em, !!!mlabs_short)
  # results$om2 <- factor(results$om, labels=mlabs_expr)
  # results$em2 <- factor(results$em, labels=mlabs_expr)	
  
  # calculate relative error
  results$SSB.rel = results$SSB_fit / results$SSB_sim
  results$SSB.rel.bc = results$SSB_fit_bc / results$SSB_sim
  results$F.rel = results$F_fit / results$F_sim
  results$F.rel.bc = results$F_fit_bc / results$F_sim
  results$relB.rel = results$relB_fit / results$relB_sim
  results$relB.rel.bc = results$relB_fit_bc / results$relB_sim
  results$relF.rel = results$relF_fit / results$relF_sim
  results$relF.rel.bc = results$relF_fit_bc / results$relF_sim
  results$catch.rel = results$catch_fit / results$catch_sim
  results$catch.rel.bc = results$catch_fit_bc / results$catch_sim
  
  simdata <- lapply(1:n.mods, function(x) readRDS(file.path(simdata_dir, bc, id, paste0("simdata_om",x,".rds"))))
  results$R.sim = NA
  for(om in 1:n.mods){
    for(em in 1:n.mods){
      for(i in 1:n.sim){
        for(ty in 1:2){
          res.ind <- which(results$om == mlabs[om] & results$em == mlabs[em] & results$sim == i & results$ty == tylabs[ty])
          results$R.sim[res.ind] <- simdata[[om]][[i]][[ty]]$NAA[,1]
        }
      }
    }
  }
  results$R.rel <- results$NAA1 / results$R.sim
  results$R.rel.bc <- results$NAA1_bc / results$R.sim

  results <- results[results$type == "Simulated data: Obs + Process error (new NAA)",]
  
  return(results)
}