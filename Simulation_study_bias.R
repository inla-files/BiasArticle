#------------------------------------
# Simulation code for Bias paper
#------------------------------------
#------------------------
# Set home directory:
setwd(".")

#------------------------------------------------

# Main program

NX <- 2000 # number of test takers test X
NY <- 2000 # number of test takers test Y
n  <- 80   # number of binary items
na <- 40   # number of anchor items
R  <- 500  # number of replicates


set.seed(6)
# Test form X
bX <- rnorm(n, 0.4, 1)   # binary item difficulty test X; For more difficult test rnorm(n, 0.9, 1)
aX <- rlnorm(n,0.3,0.4)  # discrimination
cX <- rbeta(n,1.6,6)     # guessing

# test form Y
bY <- rnorm(n, 0.4, 1)    
aY <- rlnorm(n,0.3,0.4)     
cY <- rbeta(n,1.6,6)    

# anchor items, X form
abX <- rnorm(na, 0.4, 1)   
aaX <- rlnorm(na,0.3,0.4)  
acX <- rbeta(na,1.6,6)    

# anchor items, Y form
abY <- rnorm(na, 0.4, 1)   
aaY <- rlnorm(na,0.3,0.4)  
acY <- rbeta(na,1.6,6)     

library(ltm)
library(mirt)
library(kequate)
library(equate)
library(MASS)
library(difR)

# Allocate space
FE_eq <- matrix(NA,n+1,R)
FESE_eq <- matrix(NA,n+1,R)
FEBIAS_eq <- matrix(NA,n+1,R)
FERMSE_eq <- matrix(NA,n+1,R)

FE_id <- matrix(NA,n+1,R)
FESE_id <- matrix(NA,n+1,R)
FEBIAS_id <- matrix(NA,n+1,R)
FERMSE_id <- matrix(NA,n+1,R)

FE_li <- matrix(NA,n+1,R)
FESE_li <- matrix(NA,n+1,R)
FEBIAS_li <- matrix(NA,n+1,R)
FERMSE_li <- matrix(NA,n+1,R)

FE_fe <- matrix(NA,n+1,R)
FESE_fe <- matrix(NA,n+1,R)
FEBIAS_fe <- matrix(NA,n+1,R)
FERMSE_fe <- matrix(NA,n+1,R)

FE_ce <- matrix(NA,n+1,R)
FESE_ce <- matrix(NA,n+1,R)
FEBIAS_ce <- matrix(NA,n+1,R)
FERMSE_ce <- matrix(NA,n+1,R)

CE_eq <- matrix(NA,n+1,R)
CESE_eq <- matrix(NA,n+1,R)
CEBIAS_eq <- matrix(NA,n+1,R)
CERMSE_eq <- matrix(NA,n+1,R)

CE_id <- matrix(NA,n+1,R)
CESE_id <- matrix(NA,n+1,R)
CEBIAS_id <- matrix(NA,n+1,R)
CERMSE_id <- matrix(NA,n+1,R)

CE_li <- matrix(NA,n+1,R)
CESE_li <- matrix(NA,n+1,R)
CEBIAS_li <- matrix(NA,n+1,R)
CERMSE_li <- matrix(NA,n+1,R)

CE_fe <- matrix(NA,n+1,R)
CESE_fe <- matrix(NA,n+1,R)
CEBIAS_fe <- matrix(NA,n+1,R)
CERMSE_fe <- matrix(NA,n+1,R)

CE_ce <-  matrix(NA,n+1,R)
CESE_ce <- matrix(NA,n+1,R)
CEBIAS_ce <-  matrix(NA,n+1,R)
CERMSE_ce <- matrix(NA,n+1,R) 

XSSm  <- matrix(NA,NX,R)
XSSam  <- matrix(NA,NX,R)
YSSm  <- matrix(NA,NY,R)
YSSam  <- matrix(NA,NY,R)

for (r in 1:R) {
  print(r)
  set.seed(r)
  # z.vals - abilities of test takers; Sigma can be used to fix correlations between the reg. test and anchor
  zdataX<- mvrnorm(n=NX, mu=c(0,0), Sigma=matrix(c(1,0.95,0.95,1),nrow=2),empirical=TRUE)
  zdataY<- mvrnorm(n=NY, mu=c(0,0), Sigma=matrix(c(1,0.95,0.95,1),nrow=2),empirical=TRUE)
  
  # Regular test scores
  XscoreB <-as.data.frame(rmvlogis(NX, cbind(bX,aX,cX),z.vals = zdataX[,1]))    # 3PL model (X test answers)  
  YscoreB <-as.data.frame(rmvlogis(NY, cbind(bY,aY,cY),z.vals = zdataY[,1]))  # Y test answers  
  
  #Anchor test scores
  XscoreBa <-rmvlogis(NX, cbind(abX,aaX,acX),z.vals = zdataX[,2])   # anchor test P
  YscoreBa <-rmvlogis(NY, cbind(abY,aaY,acY),z.vals = zdataY[,2])   # anchor test Q
  
  
  # Obtain sum scores for the four tests
  XSS  <- apply(XscoreB[,1:n],1,sum)  # Sum score for test X
  XSSa <- apply(XscoreBa[,1:na],1,sum) # Sum score for anchor test
  YSS  <- apply(YscoreB[,1:n],1,sum)  # Sum score for test X
  YSSa <- apply(YscoreBa[,1:na],1,sum) # Sum score for anchor test
  
  # Saving of scores
  XSSm[,r]  <- XSS
  XSSam[,r]  <- XSSa
  YSSm[,r]  <- YSS
  YSSam[,r]  <- YSSa
  
  #---------------------------------
  # Equating
  #--------------------------------
  # collect reg and anchor scores
  neat1k<-cbind(XSS,XSSa)
  neat2k<-cbind(YSS,YSSa)
  
  nex<-freqtab(x = neat1k, scales = list(0:n, 0:na))
  ney<-freqtab(x = neat2k, scales = list(0:n, 0:na))
  
  # Perform equating
  
  # Equiperc. criteria f-ion
  neq <- equate(nex,ney,type = "equipercentile", method = "none", bootse = TRUE)
  critF <- neq$conc$yx  
  
  # FREQ
  FEeq <- equate(nex,ney,type="equip", method ="frequency estimation")
  bootfe <- bootstrap(FEeq,crit = critF, reps = 100)
  # CHAIN
  CEeq <- equate(nex,ney,type="equip", method ="chain")
  bootce <- bootstrap(CEeq,crit = critF, reps = 100)
  #--------------------------------------------
  
  # Identity criteria f-ion
  neq <- equate(nex,ney,type = "identity", method = "none", bootse = TRUE)
  critF <- neq$conc$yx  
  
  # FREQ
  FEid <- equate(nex,ney,type="equip", method ="frequency estimation")
  bootfi <- bootstrap(FEid,crit = critF, reps = 100)
  # CHAIN
  CEid <- equate(nex,ney,type="equip", method ="chain")
  bootci <- bootstrap(CEid,crit = critF, reps = 100)
  #--------------------------------------------
  
  # Linear criteria f-ion
  neq <- equate(nex,ney,type = "linear", method = "none", bootse = TRUE)
  critF <- neq$conc$yx  
  
  # FREQ
  FEli <- equate(nex,ney,type="equip", method ="frequency estimation")
  bootfl <- bootstrap(FEli,crit = critF, reps = 100)
  # CHAIN
  CEli <- equate(nex,ney,type="equip", method ="chain")
  bootcl <- bootstrap(CEli,crit = critF, reps = 100)
  #--------------------------------------------
  
  # Frequency criteria f-ion
  # FREQ
  FEfe <- equate(nex,ney,type="equip", method ="frequency estimation") 
  critF <- FEfe$conc$yx 
  bootff <- bootstrap(FEfe,crit = critF, reps = 100)
  
  # CHAIN
  CEfe <- equate(nex,ney,type="equip", method ="chain")
  bootcf <- bootstrap(CEfe,crit = critF, reps = 100)
  #--------------------------------------------
  
  # Chain criteria f-ion
  # CHAIN
  CEce <- equate(nex,ney,type="equip", method ="chain")
  critF <- CEce$conc$yx 
  bootcc <- bootstrap(CEce,crit = critF, reps = 100)
  # FREQ
  FEce <- equate(nex,ney,type="equip", method ="frequency estimation")
  bootfc <- bootstrap(FEce,crit = critF, reps = 100)
  
  # Save the results to a matrix
  FE_eq[,r]  <- FEeq$conc$yx
  FESE_eq[,r] <- bootfe$se
  FEBIAS_eq[,r] <- bootfe$bias
  FERMSE_eq[,r] <- bootfe$rmse
  
  FE_id[,r]  <- FEid$conc$yx
  FESE_id[,r] <- bootfi$se
  FEBIAS_id[,r] <- bootfi$bias
  FERMSE_id[,r] <- bootfi$rmse
  
  FE_li[,r]  <- FEli$conc$yx
  FESE_li[,r] <- bootfl$se
  FEBIAS_li[,r] <- bootfl$bias
  FERMSE_li[,r] <- bootfl$rmse
  
  FE_fe[,r]  <- FEfe$conc$yx
  FESE_fe[,r] <- bootff$se
  FEBIAS_fe[,r] <- bootff$bias
  FERMSE_fe[,r] <- bootff$rmse
  
  FE_ce[,r]  <- FEce$conc$yx
  FESE_ce[,r] <- bootfc$se
  FEBIAS_ce[,r] <- bootfc$bias
  FERMSE_ce[,r] <- bootfc$rmse
  
  CE_eq[,r]  <- CEeq$conc$yx
  CESE_eq[,r] <- bootce$se
  CEBIAS_eq[,r] <- bootce$bias
  CERMSE_eq[,r] <- bootce$rmse
  
  CE_id[,r]  <- CEid$conc$yx
  CESE_id[,r] <- bootci$se
  CEBIAS_id[,r] <- bootci$bias
  CERMSE_id[,r] <- bootci$rmse
  
  CE_li[,r]  <- CEli$conc$yx
  CESE_li[,r] <- bootcl$se
  CEBIAS_li[,r] <- bootcl$bias
  CERMSE_li[,r] <- bootcl$rmse
  
  CE_fe[,r]  <- CEfe$conc$yx
  CESE_fe[,r] <- bootcf$se
  CEBIAS_fe[,r] <- bootcf$bias
  CERMSE_fe[,r] <- bootcf$rmse
  
  CE_ce[,r]  <- CEce$conc$yx
  CESE_ce[,r] <- bootcc$se
  CEBIAS_ce[,r] <- bootcc$bias
  CERMSE_ce[,r] <- bootcc$rmse  
} 

# scenario
t<-"S1"

# write results to .csv files
write.csv(FE_eq, file = paste("FE_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FESE_eq, file = paste("FESE_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FEBIAS_eq, file = paste("FEBIAS_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FERMSE_eq, file = paste("FERMSE_eq_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(FE_id, file = paste("FE_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FESE_id, file = paste("FESE_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FEBIAS_id, file = paste("FEBIAS_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FERMSE_id, file = paste("FERMSE_id_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(FE_li, file = paste("FE_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FESE_li, file = paste("FESE_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FEBIAS_li, file = paste("FEBIAS_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FERMSE_li, file = paste("FERMSE_li_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(FE_fe, file = paste("FE_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FESE_fe, file = paste("FESE_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FEBIAS_fe, file = paste("FEBIAS_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FERMSE_fe, file = paste("FERMSE_fe_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(FE_ce, file = paste("FE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FESE_ce, file = paste("FESE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FEBIAS_ce, file = paste("FEBIAS_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(FERMSE_ce, file = paste("FERMSE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
#-----------------------
write.csv(CE_eq, file = paste("CE_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CESE_eq, file = paste("CESE_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CEBIAS_eq, file = paste("CEBIAS_eq_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CERMSE_eq, file = paste("CERMSE_eq_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(CE_id, file = paste("CE_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CESE_id, file = paste("CESE_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CEBIAS_id, file = paste("CEBIAS_id_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CERMSE_id, file = paste("CERMSE_id_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(CE_li, file = paste("CE_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CESE_li, file = paste("CESE_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CEBIAS_li, file = paste("CEBIAS_li_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CERMSE_li, file = paste("CERMSE_li_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(CE_fe, file = paste("CE_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CESE_fe, file = paste("CESE_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CEBIAS_fe, file = paste("CEBIAS_fe_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CERMSE_fe, file = paste("CERMSE_fe_R500_",t,".csv",sep=""), row.names = FALSE)

write.csv(CE_ce, file = paste("CE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CESE_ce, file = paste("CESE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CEBIAS_ce, file = paste("CEBIAS_ce_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(CERMSE_ce, file = paste("CERMSE_ce_R500_",t,".csv",sep=""), row.names = FALSE)
#------------------------------
write.csv(XSSm, file = paste("XSS_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(XSSam, file = paste("XSSa_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(YSSm, file = paste("YSS_R500_",t,".csv",sep=""), row.names = FALSE)
write.csv(YSSam, file = paste("YSSa_R500_",t,".csv",sep=""), row.names = FALSE)
#===========================================================================================
