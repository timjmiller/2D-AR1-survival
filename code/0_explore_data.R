# Brian Stock
# May 12, 2020
# 2D AR1 survival devs (Haikun's paper)

# 0) Look at data, compare Haikun's .RData with read_asap3_dat("ex2_SNEMAYT.dat")

library(here)
library(wham) # https://timjmiller.github.io/wham
library(ecodata) # EDAB package, hosts GSI and other env data, https://github.com/NOAA-EDAB/ecodata
library(tidyverse)

load(here("data","sneyt.Rdata"))

# get SNEMAYT data from wham
wham.dir <- find.package("wham")
file.copy(from=file.path(wham.dir,"extdata","ex2_SNEMAYT.dat"), to=here("data","ex2_SNEMAYT.dat"), overwrite=FALSE)
gsi.wham <- read.csv(here("data","GSI_wham.csv"), header=TRUE) 

# get GSI data from ecodata
# for now take annual mean... but should confirm
gsi.ecodata <- ecodata::gsi
gsi.ecodata$Year <- floor(gsi.ecodata$Time)
gsi.ecodata <- gsi.ecodata %>% group_by(Year) %>% summarize(Annual.mean = mean(Value)) %>% as.data.frame
write.csv(gsi.ecodata, file=here("data","GSI_ecodata.csv"), row.names=FALSE)

# set up wham
asap3 <- read_asap3_dat(here("data","ex2_SNEMAYT.dat"))
NAA_list <- list(cor="iid", sigma="rec+1")
ecov <- list(
	label = "GSI",
	mean = as.matrix(gsi$Annual.mean),
	logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
	year = gsi$Year,
	use_obs = matrix(1, ncol=1, nrow=dim(gsi)[1]), # use all obs (=1)
	lag = 1, # GSI in year t affects Rec in year t + 1
	process_model = 'rw', # "rw" or "ar1"
	where = "recruit", # GSI affects recruitment
	how = 2, # 0 = no effect (but still fit Ecov to compare AIC), 2 = limiting
	link_model = "linear")

input <- prepare_wham_input(asap3, recruit_model = 3, # Bev Holt recruitment                       
                          NAA_re = NAA_list,
                          ecov=ecov)

# age comp logistic normal pool obs (not multinomial, the default)
input$data$age_comp_model_fleets = rep(5, input$data$n_fleets) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
input$data$n_age_comp_pars_fleets = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets]
input$data$age_comp_model_indices = rep(5, input$data$n_indices) # 1 = multinomial (default), 5 = logistic normal (pool zero obs)
input$data$n_age_comp_pars_indices = c(0,1,1,3,1,2)[input$data$age_comp_model_indices]
n_catch_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_fleets[which(apply(input$data$use_catch_paa,2,sum)>0)]]
n_index_acomp_pars = c(0,1,1,3,1,2)[input$data$age_comp_model_indices[which(apply(input$data$use_index_paa,2,sum)>0)]]
input$par$catch_paa_pars = rep(0, sum(n_catch_acomp_pars))
input$par$index_paa_pars = rep(0, sum(n_index_acomp_pars))

# as in WHAM ex2, use selectivity from ASAP file (logistic)
#   2 pars per block instead of n.ages
#   sel pars of indices 4/5 fixed at 1.5, 0.1 (neg phase in .dat file)
input$par$logit_selpars[1:4,7:8] <- 0 # original code started selpars at 0 (last 2 rows are fixed)

# ----------------------------------------------------
# GSI
#  v1. wham. 1954-2017. no SE. Think I took the annual average.
#  v2. ecodata. 1993-2018. no SE. Need to adjust year to be numeric (have 1995.01 as Jan 1995, should be 1995.0416)
#  v3. Haikun. 1973-2011. has SE - how? How was it collapsed in time to get one value per year? Average whole year? Subset of months?
#        cites Joyce and Zhang 2010.
#  v4. from Vince. 1954-2017 quarterly. calculate annual average and SE.
#  Nye et al. 2011 uses "spring value" of GSI
#  O'Leary et al. 2019 also uses spring GSI but detrends, cites Xu 2018; Joyce and Zhang 2010.
years <- 1:x$n_years+1972
gsi.haikun <- data.frame(Year=years, GSI=x$Ecov_obs, GSI.se=x$Ecov_obs_sigma)
gsi.wham <- read.csv(here("data","GSI_wham.csv"), header=TRUE) 
gsi.v4 <- read.csv("/home/bstock/Documents/NRC/mapp/GSI.csv")
# se <- function(x) sqrt(var(x)/length(x))
gsi.v4$Year <- floor(gsi.v4$Year)
gsi.v4 <- gsi.v4 %>% group_by(Year) %>% summarize(GSI.mean = mean(GSI), GSI.se = sqrt(var(GSI)/n())) %>% as.data.frame
write.csv(gsi.v4, file=here("data","GSI_v4.csv"), row.names=FALSE)

# is Haikun's already lagged?
png(here("plots","GSI.png"), width=9, heigh=6, units='in', res=200)
plot(gsi.wham$Year, gsi.wham$GSI, lwd=2, type='b', xlab="Year", ylab="GSI")
lines(gsi.haikun$Year, gsi.haikun$GSI, lwd=2, col='red')
lines(gsi.ecodata$Year, gsi.ecodata$Annual.mean, lwd=2, col='blue')
legend("topleft", legend = c("wham", "Xu et al. (2018)", "ecodata"), col=c("black","red","blue"), lty=c(1,1,1), lwd=c(2,2,2), bty = "n", cex=1)
dev.off()

# --------------------------------------------------------------------
# Stock data 
shared.names <- names(x)[names(x) %in% names(input$data)]
not.match <- shared.names[!mapply(function(x,y) identical(x,y), x[shared.names], input$data[shared.names])]
for(i in 1:length(not.match)){
	print(x[not.match[i]])
	print(input$data[not.match[i]])
}

# not.match
#  [1] "n_selblocks"              "selblock_models"         
#  [3] "selblock_pointer_indices" "age_comp_model_fleets"   
#  [5] "age_comp_model_indices"   "catch_paa"               
#  [7] "index_paa"                "Ecov_obs"                
#  [9] "selpars_lower"            "selpars_upper"           
# [11] "n_NAA_sigma"              "NAA_sigma_pointers"      
# [13] "recruit_model"
print.both <- function(obj){
	print(x[[not.match[i]]])
	print(input$data[[not.match[i]]])
}
i=1
print.both(not.match[i])

# difference in dimensions
identical(x$catch_paa, input$data$catch_paa[1,,])
identical(x$index_paa[0*39+1:39,], input$data$index_paa[1,,])
identical(x$index_paa[1*39+1:39,], input$data$index_paa[2,,])
identical(x$index_paa[2*39+1:39,], input$data$index_paa[3,,])
identical(x$index_paa[3*39+1:39,], input$data$index_paa[4,,])
identical(x$index_paa[4*39+1:39,], input$data$index_paa[5,,])

# differences between Haikun's data/settings
#   index 5 selectivity set equal to index 4 selectivity (but mapped to NA, start at same values, so shouldn't be a difference)
#   GSI

# ----------------------------------------------------------------
# Cold Pool
#  CPI 			original file from Miller et al. 2016, 1973-2011, SE
#  CPI_ecodata 	2020 state of the ecosystem (ecodata), 1977-2018, no SE
#  CPI_v2 		stand-in for v3, add 2012-2018 + SE (use 0.25 for missing years)
#  CPI_v3 		from Chris Melrose, Miller et al. 2016 update, 1972-2018, SE

CPI <- read.csv(here("data","CPI.csv"), header=TRUE)
CPI_ecodata <- read.csv(here("data","CPI_ecodata.csv"), header=TRUE)
df <- data.frame(CPI=CPI[CPI$Year %in% 1977:2011, "CPI"],
				 CPI2=CPI_ecodata[CPI_ecodata$Year %in% 1977:2011, "VAR"])
plot(df$CPI, df$CPI2)
mod <- lm(CPI ~ CPI2, data=df)
df$predCPI = predict(mod)
# plot(df$predCPI, df$predCPI-df$CPI)
# abline(h=0, lty=2)

# 2017 missing, add NA
CPI_recent <- predict(mod, newdata=data.frame(CPI2=CPI_ecodata[CPI_ecodata$Year %in% 2012:2018, "VAR"]))
CPI_v2 <- rbind(CPI, data.frame(Year=c(2012:2016,2018,2017),
								CPI=c(CPI_recent,NA),
								CPI_sigma=median(CPI$CPI_sigma)))
CPI_v2 <- CPI_v2[order(CPI_v2$Year),]
CPI_v2$use <- as.numeric(!is.na(CPI_v2$CPI))

write.csv(CPI_v2, file=here("data","CPI_v2.csv"), row.names=FALSE)
