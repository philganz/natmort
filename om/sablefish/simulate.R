#================================================================================================
#=================Define universal working directories
#================================================================================================

# set working directory to source file location #

#================================================================================================
#=================Set seed for random number generator and load required packages
#================================================================================================
set.seed(269)

library(R2admb)

#================================================================================================
#===Define simulation cases
#================================================================================================

#Number of replicates
R <- 1000

#Error in simulated data
obs_srv_biom_CV   <- 0.05
nsamples_srv_age  <- 160
obs_fish_biom_CV  <- 0.01
nsamples_fish_age <- 200

#Number of natural mortality / covariate trend cases
m <- 4

#Ratio of fishing mortality to natural mortality
f2m <- 2

#Error in natural mortality covariate
cov_CV <- c(0.0, 0.2, 0.4, 0.6)

#================================================================================================
#===Define universal values
#================================================================================================

ages <- c(2:31)
years <- c(1985:2014)

nages <- length(ages)
nyears <-length(years)

recage <- 2
spawn_month <- 1

#================================================================================================
#=================Parameters from 2016 age-structured model, Alaska Fisheries Science Center
#================================================================================================
pars <- read_pars("tem")$coeflist

#Estimated within 2016 model
log_q_srv    <- pars$log_q_srv1
log_mean_rec <- pars$log_mean_rec
log_rec_dev  <- pars$log_rec_dev
a50_fish     <- exp((pars$log_a50_fish1_f + pars$log_a50_fish1_m)/2)
delta_fish   <- 2
a50_srv      <- exp((pars$log_a50_srv1_f + pars$log_a50_srv1_m)/2)
delta_srv    <- 2
log_avg_F    <- pars$log_avg_F_fish1
log_F_devs   <- pars$log_F_devs_fish1

#================================================================================================
#=================Natural mortality parameters
#================================================================================================
M_0         <- 0.1
log_M_0     <- log(M_0)
Beta        <- 0.030351642 
log_sigma_M <- log(Beta)

#================================================================================================
#=================Deterministic components estimated outside model, found in SAFE
#================================================================================================
#Maturity
mat <- 1/(1+exp(-0.84*(ages-6.6)))

#Weight-at-age
wt_F <- exp((log(5.47)+3.02*log(1-exp(-0.238*(ages+1.39)))))
wt_M <- exp(log(3.16)+2.96*log(1-exp(-0.356*(ages+1.13))))
wt <- (wt_F+wt_M)/2

#================================================================================================
#=================Create case arrays 
#================================================================================================

#Dependent values
cov        <- array(NA,dim=c(nyears,m))
M          <- array(NA,dim=c(nyears,m))
F_year     <- array(NA,dim=c(nyears,R))
F          <- array(NA,dim=c(nyears,nages,R))
Z          <- array(NA,dim=c(nyears,nages,R,m))
S          <- array(NA,dim=c(nyears,nages,R,m))
N          <- array(NA,dim=c(nyears,nages,R,m))
biom       <- array(NA,dim=c(nyears,R,m))
spawn_biom <- array(NA,dim=c(nyears,R,m))
rec        <- array(NA,dim=c(nyears,R,m))
srv        <- array(NA,dim=c(nyears,nages,R,m))
srv_biom   <- array(NA,dim=c(nyears,R,m))
ac_srv     <- array(NA,dim=c(nyears,nages,R,m))
fish       <- array(NA,dim=c(nyears,nages,R,m))
fish_biom  <- array(NA,dim=c(nyears,R,m))
ac_fish    <- array(NA,dim=c(nyears,nages,R,m))

#================================================================================================
#=================Organize recruitment and fishing mortality deviations
#================================================================================================
log_rec_dev <- data.frame(log_rec_dev=log_rec_dev)
rownames(log_rec_dev) <- as.character(1932:2015)
rec_sigma <- sd(log_rec_dev[as.character(1979:2015),"log_rec_dev"])

log_F_devs <- data.frame(log_F_devs=log_F_devs)
rownames(log_F_devs) <- as.character(1960:2016)
F_sigma <- sd(log_F_devs[,"log_F_devs"])

#================================================================================================
#=================Simulate the covariate
#================================================================================================
for (r in 1:R){
cov[,1] <- scale(1:nyears)[,1]
cov[,2] <- scale(nyears:1)[,1]
cov[,3] <- scale(sin(0.2*1:nyears))[,1]
cov[,4] <- 0}

#================================================================================================
#=================Get_Selectivity
#================================================================================================
srv_sel <- 1/(1+exp(-delta_srv*(ages-a50_srv)))
fish_sel <- 1/(1+exp(-delta_fish*(ages-a50_fish)))

#================================================================================================
#=================Get_Mortality_Rates
#================================================================================================
for (i in 1:nyears){
M[i,] <- exp(log_M_0)+Beta*cov[i,]}

for (r in 1:R){
F_year[,r] <-  exp(log_avg_F+log_F_devs[as.character(years),"log_F_devs"])*exp(rnorm(nyears,0,F_sigma)-0.5*F_sigma^2)
#Scale fishing intensity to 1 or 2 times natural mortality (with f2m)
F_year[,r] <- M_0/mean(F_year[,r])*F_year[,r]*f2m}

for (i in 1:nyears){
for (j in 1:nages){
for (r in 1:R){	
for (k in 1:m){
F[i,j,r] <- F_year[i,r]*fish_sel[j]
Z[i,j,r,k] <- M[i,k]+F[i,j,r]
S[i,j,r,k] <- exp(-(Z[i,j,r,k]))}}}}

#================================================================================================
#=================Get_Numbers_At_age
#================================================================================================
#Intitial abundance
for (r in 1:R){
for (k in 1:m){
for (j in 2:(nages-1)){
N[1,j,r,k] <- exp(rnorm(1, mean(log_mean_rec+log_rec_dev[as.character(1979:2013),"log_rec_dev"]), rec_sigma)-0.5*rec_sigma^2)*exp(-(j-1)*0.1)}
N[1,nages,r,k] <- exp(mean(log_mean_rec+log_rec_dev[as.character(1979:2013),"log_rec_dev"]))*exp(-(nages-1)*0.1)*1/(1-exp(-0.1))}

#Recruitment
N[,1,r,] <- exp(rnorm(nyears, mean(log_mean_rec+log_rec_dev[as.character(1979:2013),"log_rec_dev"]), rec_sigma)-0.5*rec_sigma^2)}

#Abundance after first age and year
for (r in 1:R){
for (k in 1:m){
for (i in 2:nyears){
for (j in 2:(nages-1)){
N[i,j,r,k] <- N[i-1,j-1,r,k]*S[i-1,j-1,r,k]}
N[i,nages,r,k] <- N[i-1,nages-1,r,k]*S[i-1,nages-1,r,k]+N[i-1,nages,r,k]*S[i-1,nages,r,k]}}}

#Total/spawning biomass
for (i in 1:nyears){
for (r in 1:R){
for (k in 1:m){
biom[i,r,k]       <- crossprod(N[i,,r,k],wt)
spawn_biom[i,r,k] <- sum(N[i,,r,k]*wt_F*mat)/2
rec[,r,k]         <- N[,1,r,k]}}}

#================================================================================================
#=================Get_Predicted_Values
#================================================================================================

#Survey
##Numbers-at-age
for (r in 1:R){
for (k in 1:m){
for (i in 1:nyears){
srv[i,,r,k] <- exp(log_q_srv) * srv_sel * N[i,,r,k]}}}

##Biomass
for (i in 1:nyears){
for (r in 1:R){
for (k in 1:m){
srv_biom[i,r,k] <- crossprod(srv[i,,r,k],wt)}}}

#Fishery
##Numbers-at-age
for (r in 1:R){
for (k in 1:m){
fish[,,r,k] <- N[,,r,k]*F[,,r]*(1-S[,,r,k])/Z[,,r,k]}}

##Biomass
for (i in 1:nyears){
for (r in 1:R){
for (k in 1:m){
fish_biom[i,r,k] <- crossprod(fish[i,,r,k],wt)}}}

#================================================================================================
#=================Get_Age_Comp
#================================================================================================

#Age composition
for (r in 1:R){
for (k in 1:m){
for (i in 1:nyears){
for (j in 1:nages){	
##Survey
ac_srv[i,j,r,k] <- srv[i,j,r,k]/sum(srv[i,,r,k])
##Fishery
ac_fish[i,j,r,k] <- fish[i,j,r,k]/sum(fish[i,,r,k])}}}}

#================================================================================================
#===Simulate data
#================================================================================================
obs_srv_biom  <- array(NA,dim=c(nyears,R,m)) 
obs_fish_biom <- array(NA,dim=c(nyears,R,m)) 
obs_ac_srv    <- array(NA,dim=c(nyears,nages,R,m)) 
obs_ac_fish   <- array(NA,dim=c(nyears,nages,R,m))
obs_cov       <- array(NA,dim=c(nyears,R,m,length(cov_CV)))

#Simulated observed survey biomass
for (i in 1:nyears){
for (r in 1:R){
for (k in 1:m){
#Observe data without error
obs_srv_biom[i,r,k] <- srv_biom[i,r,k]
obs_ac_srv[i,,r,k]<-ac_srv[i,,r,k]
obs_fish_biom[i,r,k] <- fish_biom[i,r,k]
obs_ac_fish[i,,r,k]<-ac_fish[i,,r,k]
}}}
      
#Observe data with error
# obs_srv_biom[i,r,k] <- srv_biom[i,r,k]*exp(rnorm(1,0,obs_srv_biom_CV)-0.5*obs_srv_biom_CV^2)
# obs_ac_srv[i,,r,k] <- rmultinom(1,nsamples_srv_age,ac_srv[i,,r,k])/nsamples_srv_age
# obs_fish_biom[i,r,k] <- fish_biom[i,r,k]*exp(rnorm(1,0,obs_fish_biom_CV)-0.5*obs_fish_biom_CV^2)
# obs_ac_fish[i,,r,k] <- rmultinom(1,nsamples_fish_age,ac_fish[i,,r,k])/nsamples_fish_age
# }}}

#Simulated observed covariate 
for (r in 1:R){
for (k in 1:m){
for (c in 1:length(cov_CV)){
#Observe data without error
obs_cov[,r,k,c] <- cov[,k]
}}}
      
#Observe data with error
# obs_cov[,r,k,c] <- cov[,k]+rnorm(nyears,0,cov_CV[c])}}}

#================================================================================================
#===Save values
#================================================================================================
# if(obs_srv_biom_CV==0.15  & f2m==1){save.image("high_1.RData")
# } else if(obs_srv_biom_CV==0.15  & f2m==2){save.image("high_2.RData")
# } else if(obs_srv_biom_CV==0.05  & f2m==1){save.image("low_1.RData")
# } else if(obs_srv_biom_CV==0.05  & f2m==2){save.image("low_2.RData")
# } else {"none of these conditions exist"}

# manually save data as calibration case
save.image("calibration.RData")
