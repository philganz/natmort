#================================================================================================
#=================Get working directories
#================================================================================================

# set working directory to source file location
## Rstudio:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
## sourced in R:
# setwd(dirname(sys.frame(1)$ofile))
## run from terminal
# this_dir <- function(directory)
# setwd(file.path(getwd(),directory))

pathR <- getwd()
path   <- substr(pathR,1,nchar(pathR)-48)
pathOM <- paste(path,"/om/sablefish",sep="")
pathEM <- paste(path,"/em/sablefish",sep="")

#================================================================================================
#=================Load simulated data
#================================================================================================
if(grepl("high/one",pathR)){load(paste(pathOM,"/high_1.RData",sep="")); pathE <- paste(pathEM,"/high/one",sep="")
} else if(grepl("high/two",pathR)){load(paste(pathOM,"/high_2.RData",sep="")); pathE <- paste(pathEM,"/high/two",sep="")
} else if(grepl("low/one",pathR)){load(paste(pathOM,"/low_1.RData",sep="")); pathE <- paste(pathEM,"/low/one",sep="")
} else if(grepl("low/two",pathR)){load(paste(pathOM,"/low_2.RData",sep="")); pathE <- paste(pathEM,"/low/two",sep="")
} else {"data not found"}

# or load data with no observation error for calibration
# load(paste(pathOM,"/calibration.RData",sep="")); pathE <- paste(pathEM,"/calibration",sep=""); pathR <- pathE 

#================================================================================================
#=================Write .ctl file
#================================================================================================
M_case  <- 1
M_start <- 1

CTL <- c(
paste(1, "# Log recruitment (ph_logR)", sep=" "),
paste(1, "# Recruitment deviations phase (ph_Rdevs)", sep=" "),
paste(1, "# Initial abundance deviations phase (ph_Idevs)", sep=" "),
paste(1, "# Mean mortality phase (ph_M_0)", sep=" "),
paste(-1, "# Initial mortality phase (ph_M_1)", sep=" "),
paste(-1, "# Drift term phase (ph_a)", sep=" "),
paste(-1, "# Beta phase (ph_B)", sep=" "),
paste(1, "# Avg F phase (ph_avgF)", sep=" "),
paste(1, "# F deviations phase (ph_Fdevs)", sep=" "),
paste(1, "# Catchability phase (ph_q)", sep=" "),
paste(1, "# Survey Selectivity phase (ph_Ssel)", sep=" "),
paste(1, "# Fishery Selectivity phase (ph_Fsel)",sep=" "),
paste(-1, "# Correlation term phase (ph_phi)", sep=" "),
paste(-2, "# Mortality deviations phase (ph_Mdevs)", sep=" "),
paste(-2, "# Random effects sigma phase (ph_sig)", sep=" "),
paste(M_case, "# Natural mortality estimation case (M_case)", sep=" "),
paste(M_start, "# Start year for M_devs (M_start)", sep=" "))

write.table(CTL,file=paste(pathE,"/tem.ctl",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

#================================================================================================
#=================Write starting values to .pin file
#================================================================================================

PIN<-c(
"# logR:",
paste(0.4),
"# rec_devs:",
paste(as.vector(rep(0,nyears)), collapse=" "),
"# init_devs:",
paste(as.vector(rep(0,nages-1)), collapse=" "),
"# log_avg_F:",
paste(-2),
"# F_devs:",
paste(as.vector(rep(0,nyears)), collapse=" "),
"# log_M_0:",
paste(-2),
"# log_M_1:",
paste(-2),
"# log_phi:",
paste(-5000),
"# alpha:",
paste(0),
"# Beta:",
paste(0),
"# log_q_srv:",
paste(2),
"# a50_srv:",
paste(3),
"# delta_srv:",
paste(2),
"# a50_fish:",
paste(4),
"# delta_fish:",
paste(2),
"# log_sigma_M:",
paste(-4),
"# M_devs:",
paste(as.vector(rep(0,nyears)), collapse=" "))

write.table(PIN,file=paste(pathE,"/tem.pin",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

#================================================================================================
#=================Move master .tpl to estimation folder and run
#================================================================================================
file.copy(from=paste(pathEM,"/tem.tpl",sep=""), to=paste(pathE,"/tem.tpl",sep=""), overwrite = TRUE)
source(paste(path,"/em/estimate.R",sep=""))
