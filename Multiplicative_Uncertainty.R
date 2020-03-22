nTime <- 28
nEns <- 10^3

traj <- matrix(NA,nTime,nEns)
traj[1,] <- rep(1,nEns)

dub <- rnorm(nEns,2.5,0.15)

for(j in 1:nEns){
  for(i in 2:nTime){
    traj[i,j] <- traj[1,j]*(2^(i/dub[j]))
  }
}

#hist(traj[nTime,])

# Finding which trajectory has maximum growth rate 
# which(dub==min(dub)) maximum value 

# Finding trajectories near the mean 
# which(abs(dub - mean(dub)) < 0.001)
typical <- sample(which(abs(dub - mean(dub)) < 0.001),1)

par(mfcol=c(1,3))


# non-log scale
plot(seq(1,nTime,1),traj[,which(dub==min(dub))],type = 'l',lwd=0.1,
     ylab = "#Infected (as multiple of starting point)",xlab = "Day",
     main = "Uncertainty in doubling rate", 
     sub = "R_0 in [2.2,2.7]")

#lines(seq(1,nTime,1),traj[,which(dub==min(dub))])
for(i in 2:nEns){
  lines(seq(1,nTime,1),traj[,i],lwd=0.1)
}

lines(seq(1,nTime,1),traj[,typical],col=2)
abline(1000,0,col=3)




## Uncertainty in starting point 

nTime <- 28
nEns <- 10^3

traj <- matrix(NA,nTime,nEns)
traj[1,] <- rnorm(nEns,1,0.25)

dub <- 2.5

typical <- sample(which(abs(traj[1,] - mean(traj[1,])) < 0.005),1)


for(j in 1:nEns){
  for(i in 2:nTime){
    traj[i,j] <- traj[1,j]*(2^(i/dub))
  }
}

#hist(traj[nTime,])


# non-log scale

plot(seq(1,nTime,1),traj[,which(traj[1,]==max(traj[1,]))],type = 'l',lwd=0.1,
     ylab = "#Infected (as multiple of starting point)",xlab = "Day", main = 
       "Uncertainty in initial # of infected", 
     sub = "inital in [0.5,1.5]")

#lines(seq(1,nTime,1),traj[,which(dub==min(dub))])
for(i in 2:nEns){
  lines(seq(1,nTime,1),traj[,i],lwd=0.1)
}

lines(seq(1,nTime,1),traj[,typical],col=2)
abline(1000,0,col=3)


## Uncertainty in BOTH 

nTime <- 21
nEns <- 10^3

traj <- matrix(NA,nTime,nEns)
traj[1,] <- rnorm(nEns,750,125)

dub <- rnorm(nEns,2.5,0.15)


for(j in 1:nEns){
  for(i in 2:nTime){
    traj[i,j] <- traj[1,j]*(2^(i/dub[j]))
  }
}
typical <- sample(which(abs(traj[nTime,] - mean(traj[nTime,])) < 20),1)
max <- which(traj[nTime,]==max(traj[nTime,]))

plot(seq(1,nTime,1),traj[,max],type = 'l',lwd=0.1,
     ylab = "#Infected (as multiple of starting point)",xlab = "Day", main = 
       "Uncertainty in both # initial infections and rate of doubling", 
     sub = "R_0 in [2.2,2.7], initial in [0.5,1.5]")

#lines(seq(1,nTime,1),traj[,which(dub==min(dub))])
for(i in 2:nEns){
  lines(seq(1,nTime,1),traj[,i],lwd=0.1)
}

lines(seq(1,nTime,1),traj[,typical],col=2)
abline(100,0,col=3)

# FL statewide 

### Zooming in 
# non-log scale
plot(seq(30,nTime,1),traj[30:nTime,which(dub==min(dub))],type = 'l')

#lines(seq(1,nTime,1),traj[,which(dub==min(dub))])
for(i in 2:nEns){
  lines(seq(30,nTime,1),traj[30:nTime,i])
}

lines(seq(1,nTime,1),traj[,typical],col=2)
abline(10000,0,col=3)



# log-scale
plot(seq(1,nTime,1),log(traj[,which(dub==min(dub))]))
lines(seq(1,nTime,1),log(traj[,which(dub==min(dub))]))
for(i in 2:nEns){
  lines(seq(1,nTime,1),log(traj[,i]))
}


# Plot the discrepancy between the ensemble mean at each point in time and 
# the value of the "typical trajectory"

jensens <- rep(0,nTime)
sfree_jensens <- rep(0,nTime)
for(i in 1:nTime){
  jensens[i] = -(traj[i,typical] - mean(traj[i,]))
  sfree_jensens[i] <- jensens[i]/traj[i,typical]
}

plot(seq(1,nTime,1),sfree_jensens,type = 'l')

# Consider harm - below a certain threshold, harm from each case is linear, but above the 
# threshold, harm accelerates as it represents breakdown of health-care system, economy, etc. 

1*(2^(5/2.5))



