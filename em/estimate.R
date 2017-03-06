#================================================================================================
#===MCMC specifications
#================================================================================================
mcmc_N    <- 12000 #500000
mcmc_save <- 100 #500
burn_in   <- 100

#================================================================================================
#===For estimation methods that don't use the covariate
#================================================================================================
if(grepl("covariate",pathR)==0){cov_CV <- 0}

#================================================================================================
#===Combine weights for sablefish
#================================================================================================
if(grepl("sablefish",pathR)){
wt <- c(paste(as.vector(wt_F),collapse=" "),
        paste(as.vector(wt_M),collapse=" "))
} else if(grepl("pollock",pathR)){wt <- paste(as.vector(wt),collapse=" ")}

#================================================================================================
#===Set up results
#================================================================================================
Iter_base         <- read.delim(paste(pathE,"/iteration_base.rep",sep=""),sep="")
Results           <- array(NA,dim=c(R,length(Iter_base),m,length(cov_CV)))
colnames(Results) <- names(Iter_base)
STD               <- read.delim(paste(pathE,"/tem.std",sep=""),sep="")
mcmc_results      <- array(NA,dim=c(mcmc_N/mcmc_save,3,R,m,length(cov_CV)))
DIC               <- array(NA,dim=c(R,m,length(cov_CV)))

#================================================================================================
#===True parameters that change with simulated error
#================================================================================================
logR       <- array(NA,dim=c(R))
rec_devs   <- array(NA,dim=c(nyears,R))
init_devs  <- array(NA,dim=c(nages-1,R))
log_avg_F  <- array(NA,dim=c(R))
F_devs     <- array(NA,dim=c(nyears,R))
alpha      <- array(NA, dim=m)
M_devs     <- array(NA,dim=c(nyears,m))

#================================================================================================
#===Loop through M scenarios, covariate error, replicates
#================================================================================================
setwd(pathE)

#Compile estimation model
system("admb -r tem")

T_start <- Sys.time()
	
for (r in 2:5){
for (k in 1:1){
for (c in 1:1){

#================================================================================================
#=================Send true parameter values to .pin file (for calibration only)
#================================================================================================
logR[r]       <- mean(log(N[,1,r,]))
rec_devs[,r]  <- log(N[,1,r,k])-logR[r]
init_devs[,r] <- log(N[1,2:nages,r,k])
log_avg_F[r]  <- mean(log(F_year[,r]))
F_devs[,r]    <- log(F_year[,r]) - log_avg_F[r]
for(i in 2:nyears){alpha[k]      <- mean(M[i,k] - M[i-1,k])}
if(M_case<=2){sigma_M <- exp(log_sigma_M); M_devs[,k] <- (M[,k] - M_0)/sigma_M
} else{sigma_M <- 0; M_devs[,k] <- 0}

PIN<-c(
"# logR:",
as.character(logR[r]),
"# rec_devs:",
paste(as.vector(rec_devs[,r]), collapse=" "),
"# init_devs:",
paste(as.vector(init_devs[,r]), collapse=" "),
"# log_avg_F:",
as.character(log_avg_F[r]),
"# F_devs:",
paste(as.vector(F_devs[,r]), collapse=" "),
"# log_M_0:",
as.character(log_M_0),
"# log_M_1:",
paste(-1.8971200),
"# phi:",
paste(1),
"# alpha:",
paste(alpha[k]),
"# Beta:",
# paste(0),
as.character(Beta),
"# log_q_srv:",
as.character(log_q_srv),
"# a50_srv:",
as.character(a50_srv),
"# delta_srv:",
as.character(delta_srv),
"# a50_fish:",
as.character(a50_fish),
"# delta_fish:",
as.character(delta_fish),
"# sigma_M:",
as.character(sigma_M),
"# M_devs:",
paste(as.vector(M_devs[,k]), collapse=" "))

# write.table(PIN,file=paste(pathE,"/tem.pin",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

#================================================================================================
#===Send simulation data to ADMB
#================================================================================================

#Survey biomass
S_SB_sim<-paste(as.vector(obs_srv_biom[,r,k]),collapse=" ")

#Catch biomass
C_SB_sim<-paste(as.vector(obs_fish_biom[,r,k]),collapse=" ")

#Survey age comp
oac_srv_sim<-vector(length=nyears)
for(i in 1:nyears){
oac_srv_sim[i]<-paste(as.vector(obs_ac_srv[i,,r,k]),collapse=" ")}

#Fishery age comp
oac_fish_sim<-vector(length=nyears)
for(i in 1:nyears){
oac_fish_sim[i]<-paste(as.vector(obs_ac_fish[i,,r,k]),collapse=" ")}

#M covariate data
cov_sim<-paste(as.vector(obs_cov[,r,k,c]),collapse=" ")

DATs<-c(
"#==========================================================================================================================",
"# Sablefish .dat file",
"# Simulated data based on 2014 age-structured population model",
"#==========================================================================================================================",	
"",
"",
"#==========================================================================================================================",
"# Model input parameters/vectors",
"#==========================================================================================================================",
"# nyrs",
as.character(nyears),
"# nages",
as.character(nages),
"# spawn month",
as.character(spawn_month),
"# recage",
as.character(recage),
"",
"",
"#==========================================================================================================================",
"# Weight-at-age",
"#==========================================================================================================================",
wt,
"",
"",
"#==========================================================================================================================",
"# Maturity",
"#==========================================================================================================================",
paste(as.vector(mat),collapse=" "),
"",
"",
"#==========================================================================================================================",
"# M covariate",
"#==========================================================================================================================",
"# Observed covariate",
cov_sim,
"",
"",
"#==========================================================================================================================",
"# Domestic longline survey RPW",
"#==========================================================================================================================",
"# CV",
obs_srv_biom_CV,
"# Observed domestic RPW:",
S_SB_sim,
"",
"",
"#==========================================================================================================================",
"# Fixed gear fishery catch (mt):",
"#==========================================================================================================================",
"# CV",
obs_fish_biom_CV,
"# Observed catch (mt):",
C_SB_sim,
"",
"",
"#==========================================================================================================================",
"# Domestic LL survey Age Composition",
"#==========================================================================================================================",
"# Number of samples: nsamples_srv_age",
nsamples_srv_age,
"# Observed age compositions (proportions at age):",
oac_srv_sim,
"",
"",
"#==========================================================================================================================",
"# Fixed Gear Fishery Age Composition",
"#==========================================================================================================================",
"# Number of samples: nsamples_fish_age",
nsamples_fish_age,
"# Observed fishery age compositions (proportions at age):",
oac_fish_sim,
"",
"",
"#==========================================================================================================================",
"# EOF marker",
"#=======================================================================================================================",
42,
"#!")

write.table(DATs,file=paste(pathE,"/tem.dat",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)

#Run Model
system("./tem")

#Record iteration results
Iter_base <- read.delim(paste(pathE,"/iteration_base.rep",sep=""),sep="")
Results[r,,k,c] <- c(as.numeric(read.delim(paste(pathE,"/iteration_base.rep",sep=""),sep="")))

#Test for ADMB convergence using .std 
Results_base<-read.delim(paste(pathE,"/iteration_base.rep",sep=""),sep="")
Results[r,,k,c] <- c(as.numeric(read.delim(paste(pathE,"/iteration_base.rep",sep=""),sep="")))
if(length(scan(paste(pathE,"/tem.std",sep=""),what=character(0)))>0){
  STDi<-read.delim(paste(pathE,"/tem.std",sep=""),sep="")
  if(STD$value[length(STD$std)]!=STDi$value[length(STDi$std)]){
    Results[r,1,k,c]<-1
    STD<-STDi}
  else{Results[r,1,k,c]<-0}}else{STD<-STDi;Results[r,1,k,c]<-0}

# MCMC
# shell(paste("tem -mcmc ", mcmc_N, " -mcsave ", mcmc_save, sep=""))
# shell("tem -mceval")
# mcmc_results[,,r,k,c] <- as.matrix(read.csv(paste(pathE,"/mcmc_results.csv",sep="")))

# DIC
D_bar <- mean(2 * mcmc_results[-c(1:burn_in),,r,k,c])
D     <- 2 * Results[r,"obj_fun",k,c]
pD    <- D_bar - D
DIC[r,m,c] <- D_bar + pD

#End covariate error loop
}

#End M trend loop
}

#End replicate loop
}

T_end<-Sys.time()

runtime <- T_end-T_start
runtime

# save results
save(runtime,Results,mcmc_results,file=paste(pathR,"/Results_",R - sum(is.na(Results[,1,1,1])),".RData",sep=""))
