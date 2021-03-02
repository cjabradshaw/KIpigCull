## Rapid pig eradication assessment for Kangaroo Island
## PIRSA & Flinders University
## Corey Bradshaw
## March 2021

## remove everything
rm(list = ls())

library(gtools)

# stochastic beta sampler (single sample)
stoch.beta.func <- function(mu, var) {
  Sx <- rbeta(length(mu), (((1 - mu) / var - 1 / mu) * mu ^ 2), ((((1 - mu) / var - 1 / mu) * mu ^ 2)*(1 / mu - 1)))
  return(params=Sx)
}

sightBias <- 0.9
totA <- 1619.92 # km2
meanD <- (137+263)/totA * (1/sightBias)
SD.D <- 10 # %
habQ.class <- c(1,2,3,7,9,10)
habQ.A <- c(47.38, 787.14, 18.23, 372.47, 59.75, 334.96)
K.N <- 6000
maxD <- K.N/totA

# estimated total population (west end of KI only)
Kvec <- round(habQ.class/5 * maxD * habQ.A, 0)
Nvec1 <- round(habQ.class/5 * meanD * habQ.A, 0)
Nvec1/Kvec
totN <- round(sum(habQ.class/5 * meanD * habQ.A), 0)
totN

habQ <- data.frame(habQ.class, habQ.A, Kvec)
habQ

# pop dynamic model
rmax <- 0.34 # McMahon et al. 2010
theta <- 1.3 # McMahon et al. 2020
propMov <- 0.01 # proportion moving into another habitat class each time step

# kill distribution
# % pigs killed per hab quality class
killdist <- c(0,51,6,138,10,58)/2
killdistprop <- killdist/sum(killdist)

## cost info
mean.cost1 <- 418 # dollars per pig at start

N.start <- Nvec1 # - (263*killdistprop)
N.start <- ifelse(N.start < 0, 1, N.start)
N.start

# 1 step
N.start[2] * exp(rmax*(1-(N.start[2]/Kvec[2])^theta))


Nnew <- round(N.start * exp(rmax*(1-(N.start/Kvec)^theta)), 0)
Nmig <- round(Nnew * propMov, 0)
Nupd <- Nnew - Nmig + sample(Nmig, length(Nnew), replace=F)
Nupd <- ifelse(Nupd < 0, 0, Nupd)


###############################
## Scenario 1: No cull
## time to reach K
proj.int <- 20 # years

Nmat <- matrix(data = NA, nrow = (2*proj.int)+1, ncol=length(N.start))
Nmat[1,] <- N.start
Nmat

for (t in 1:(proj.int*2)) {
  Nmat[t+1,] <- round(Nmat[t,] * exp((rmax/2)*(1-(Nmat[t,]/Kvec)^theta)), 0)
  Nmig <- round(Nmat[t+1,] * propMov, 0)
  Nmat[t+1,] <- Nmat[t+1,] - Nmig + sample(Nmig, length(Nmat[t+1,]), replace=F)
  Nmat[t+1,] <- ifelse(Nmat[t+1,] < 0, 0, Nmat[t+1,])
}
Nmat

Npred <- rowSums(Nmat)
plot(1:dim(Nmat)[1], Npred, type="l",xlab="time step", ylab="N")




#########################################################################################
## stochastic, assuming 5% variation in following parameters: K, killdist, rmax, propMov
## Scenario 2: 90% winter kill rate
summer.cost1 <- 418 # dollars per pig at start
tc.hours <- 250 # total winter hours of thermal aerial culling
area.cph <- 500 # area culled per hour (ha) 
tot.area.tc <- (tc.hours * area.cph)/100 # total area covered by thermal culls (km2)
propcull.ta <- tot.area.tc/1619.92 # proportion of total area culled
costph <- 600000/250 # cost per hour in dollars
kill.efficiency <- 0.9
stochSD <- 0.05
iter <- 10000
itdiv <- iter/10

proj.int <- 3 # years
eradThresh <- 2 # eradication threshold
Npred.mat <- matrix(data=NA, nrow = iter, ncol = (2*proj.int)+1)
totCost <- erad <- rep(NA,iter)

for (i in 1:iter) {
  
  propMovStoch <- stoch.beta.func(propMov, stochSD*propMov)
  Kstoch <- round(rnorm(length(Kvec), mean=Kvec, sd=stochSD*Kvec),0)
  killdistStoch <- rpois(length(killdist), killdist)
  killtot <- round(sum(killdistStoch),0)
  rmaxStoch <- stoch.beta.func(rmax, stochSD*rmax)
    
  Nmat <- matrix(data = NA, nrow = (2*proj.int)+1, ncol=length(N.start))
  Nmat[1,] <- N.start
  Nmat
  cost.vec <- rep(NA,(2*proj.int)+1)
  
  for (t in 1:(proj.int*2)) {
    Nmat[t+1,] <- round(Nmat[t,] * exp((rmaxStoch/2)*(1-(Nmat[t,]/Kstoch)^theta)), 0)
    Nmig <- round(Nmat[t+1,] * propMovStoch, 0)
    Nmat[t+1,] <- Nmat[t+1,] - Nmig + sample(Nmig, length(Nmat[t+1,]), replace=F)
    Nmat[t+1,] <- ifelse(Nmat[t+1,] < 0, 0, Nmat[t+1,])
    Nmat[t+1,] <- ifelse(is.na(Nmat[t+1,]) == T, 0, Nmat[t+1,])
    
    ## summer culling
    if (t+1 > 2 & odd(t+1) == T) {
      Nkill1 <- round((killtot * Nmat[t+1,]/sum(Nmat[t+1,],na.rm=T)),0)
      Nkill1 <- ifelse(is.na(Nkill1) == T, 0, Nkill1)
      #Nkill2 <- round((sum.aerial.cull[t+1] * Nmat[t+1,]/sum(Nmat[t+1,],na.rm=T)),0)
      Nupd <- Nmat[t+1,] - Nkill1
      Nupd <- ifelse(Nupd < 0, 0, Nupd)
      Nmat[t+1,] <- Nupd
      
      cost.vec[t+1] <- (summer.cost1 * round(sum(Nkill1),0))
    }
    
    ## winter culling
    if (t+1 > 1 & even(t+1) == T) {
      Nupd <- Nmat[t+1,] - round((propcull.ta * kill.efficiency * Nmat[t+1,]), 0)
      Nupd <- ifelse(Nupd < 0, 0, Nupd)
      Nkilled <- Nmat[t+1,] - Nupd
      Nmat[t+1,] <- Nupd
      cost.vec[t+1] <- 600000
    }
    
  }
  
  Npred.mat[i, ] <- rowSums(Nmat, na.rm=T)
  totCost[i] <- sum(cost.vec,na.rm=T)
  erad[i] <- ifelse(Npred.mat[i,dim(Npred.mat)[2]] < eradThresh, 1, 0)  
  if (i %% itdiv==0) print(i) 
  
}

Pr.erad <- sum(erad, na.rm=T)/iter
Npred.mean <- apply(Npred.mat, MARGIN=2, mean, na.rm=T)
Npred.lo <- apply(Npred.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
Npred.up <- apply(Npred.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)

plot(1:dim(Npred.mat)[2], Npred.mean, type="l", xlab="time step", ylab="N", ylim=c(min(Npred.lo),max(Npred.up)))
lines(1:dim(Npred.mat)[2], Npred.lo, lty=2, col="red")
lines(1:dim(Npred.mat)[2], Npred.up, lty=2, col="red")
Pr.erad

stoch.out <- data.frame(1:dim(Npred.mat)[2], Npred.mean, Npred.up, Npred.lo)
colnames(stoch.out) <- c("t", "Nmn", "Nup", "Nlo")
write.csv(stoch.out, file="stoch1.csv")

mean(totCost)
quantile(totCost,probs=0.025,na.rm=T)
quantile(totCost,probs=0.975,na.rm=T)



#########################################################################################
## stochastic, assuming 5% variation in following parameters: K, killdist, rmax, propMov
## Scenario 3: additional, small aerial cull late summer
stochSD <- 0.05
iter <- 10000
itdiv <- iter/10

sum.aerial.cull <- c(0,0,130,0,65,0,33)

proj.int <- 3 # years
eradThresh <- 2 # eradication threshold
Npred.mat <- matrix(data=NA, nrow = iter, ncol = (2*proj.int)+1)
totCost <- erad <- rep(NA,iter)

for (i in 1:iter) {
  
  propMovStoch <- stoch.beta.func(propMov, stochSD*propMov)
  Kstoch <- round(rnorm(length(Kvec), mean=Kvec, sd=stochSD*Kvec),0)
  killdistStoch <- rpois(length(killdist), killdist)
  killtot <- round(sum(killdistStoch),0)
  rmaxStoch <- stoch.beta.func(rmax, stochSD*rmax)
  
  Nmat <- matrix(data = NA, nrow = (2*proj.int)+1, ncol=length(N.start))
  Nmat[1,] <- N.start
  Nmat
  cost.vec <- rep(NA,(2*proj.int)+1)
  
  for (t in 1:(proj.int*2)) {
    Nmat[t+1,] <- round(Nmat[t,] * exp((rmaxStoch/2)*(1-(Nmat[t,]/Kstoch)^theta)), 0)
    Nmig <- round(Nmat[t+1,] * propMovStoch, 0)
    Nmat[t+1,] <- Nmat[t+1,] - Nmig + sample(Nmig, length(Nmat[t+1,]), replace=F)
    Nmat[t+1,] <- ifelse(Nmat[t+1,] < 0, 0, Nmat[t+1,])
    Nmat[t+1,] <- ifelse(is.na(Nmat[t+1,]) == T, 0, Nmat[t+1,])
    
    ## summer culling
    if (t+1 > 2 & odd(t+1) == T) {
      Nkill1 <- round((killtot * Nmat[t+1,]/sum(Nmat[t+1,],na.rm=T)),0)
      Nkill1 <- ifelse(is.na(Nkill1) == T, 0, Nkill1)
      Nkill2 <- round((sum.aerial.cull[t+1] * Nmat[t+1,]/sum(Nmat[t+1,],na.rm=T)),0)
      Nkill2 <- ifelse(is.na(Nkill2) == T, 0, Nkill2)
      Nupd <- Nmat[t+1,] - Nkill1 - Nkill2 
      Nupd <- ifelse(Nupd < 0, 0, Nupd)
      Nmat[t+1,] <- Nupd
      
      cost.vec[t+1] <- (summer.cost1 * round(sum(Nkill1),0)) + 100000
    }
    
  }
  
  Npred.mat[i, ] <- rowSums(Nmat,na.rm=T)
  totCost[i] <- sum(cost.vec,na.rm=T)
  erad[i] <- ifelse(Npred.mat[i,dim(Npred.mat)[2]] < eradThresh, 1, 0)  
  if (i %% itdiv==0) print(i) 
  
}

Pr.erad <- sum(erad, na.rm=T)/iter
Npred.mean <- apply(Npred.mat, MARGIN=2, mean, na.rm=T)
Npred.lo <- apply(Npred.mat, MARGIN=2, quantile, probs=0.025, na.rm=T)
Npred.up <- apply(Npred.mat, MARGIN=2, quantile, probs=0.975, na.rm=T)

plot(1:dim(Npred.mat)[2], Npred.mean, type="l", xlab="time step", ylab="N", ylim=c(min(Npred.lo),max(Npred.up)))
lines(1:dim(Npred.mat)[2], Npred.lo, lty=2, col="red")
lines(1:dim(Npred.mat)[2], Npred.up, lty=2, col="red")
Pr.erad

stoch.out <- data.frame(1:dim(Npred.mat)[2], Npred.mean, Npred.up, Npred.lo)
colnames(stoch.out) <- c("t", "Nmn", "Nup", "Nlo")
write.csv(stoch.out, file="stoch2.csv")

mean(totCost)
quantile(totCost,probs=0.025,na.rm=T)
quantile(totCost,probs=0.975,na.rm=T)
