log(1-(0.02/20))
curve(exp(-0.001*x),0,200)
log(1-0.98/20)

### Basic SIR simulation 
# R_naught ~ 2 

days <- 360
init_pop <- 3*10^5
init_sick <- 500
k <- 5
p <- 0.02

alpha <- (k*p)/init_pop
CFR <- 0.02
sickDay <- 20 

c <- log(1-(CFR/sickDay))
v <- log(1-((1-CFR)/sickDay))

dt <- 1 

R_o <- -k*p/(c+v)


S <- rep(0,days/dt)
I <- rep(0,days/dt)
R <- rep(0,days/dt)
Dead <- rep(0,days/dt)
H <- rep(0,days/dt)

S[1] <- init_pop - init_sick
I[1] <- init_sick 
R[1] <- 0
Dead[1] <- 0 
H[1] <- 0 

for(i in 2:(days/dt)){
  S[i] <- S[i-1] - (alpha*I[i-1]*S[i-1])*dt
  I[i] <- I[i-1] + ((alpha*I[i-1]*S[i-1]) + (c+v)*I[i-1])*dt
  R[i] <- R[i-1] - (v*I[i-1])*dt
  Dead[i] <- Dead[i-1] - (c*I[i-1])*dt
 
}

plot(seq(1,days/dt,1),I)
I[1:10]
plot(seq(1,days/dt,1),Dead)
Dead[1:10]
plot(S,R)

####### Sensitivity analyses on isolation 
k_rep <- seq(0.1,10,0.1)
dead_fin <- rep(0,length(k_rep))

for(l in 1:length(k_rep)){
  
days <- 90 
init_pop <- 3*10^5
init_sick <- 500
k <- k_rep[l]
p <- 0.05

alpha <- (k*p)/init_pop
CFR <- 0.02
sickDay <- 20 

c <- log(1-(CFR/sickDay))
v <- log(1-((1-CFR)/sickDay))

dt <- 1 

S <- rep(0,days/dt)
I <- rep(0,days/dt)
R <- rep(0,days/dt)
Dead <- rep(0,days/dt)

S[1] <- init_pop - init_sick
I[1] <- init_sick 
R[1] <- 0
Dead[1] <- 0 

for(i in 2:(days/dt)){
  S[i] <- S[i-1] - (alpha*I[i-1]*S[i-1])*dt
  I[i] <- I[i-1] + ((alpha*I[i-1]*S[i-1]) + (c+v)*I[i-1])*dt
  R[i] <- R[i-1] - (v*I[i-1])*dt
  Dead[i] <- Dead[i-1] - (c*I[i-1])*dt
}

dead_fin[l] <- Dead[days/dt]
print(dead_fin[(days/dt)-10:(days/dt)])
}

plot(k_rep,dead_fin)




#### Adding a hospitalization intermediate 

days <- 360
init_pop <- 3*10^5
init_sick <- 500
k <- 5
p <- 0.04

alpha <- (k*p)/init_pop
hr <- 0.1 # Hospitalization as % of patients
CFR <- hr*0.2 # Overall CFR as product of % of hospitalized patients 
incubate <- 10 # Length of incubation before hospitalization 
stay <- 14 # Length of time from hospitalization to death 
sickDay <- 14 # average time to recovery/resistant status for non_hostpializaed

h <- -log(1-(hr/incubate)) # rate of hospitalization 
c <- -log(1-(0.2/stay)) # Mortality of hospitalized 
v_h <- -log(1-(0.8/stay)) # Recovery rate from hospitalized to well 
v_nh <- -log(1-((1-hr)/incubate)) # Recovery rate for non-hospitalized 

dt <- 1 

# R_o <- -k*p/(c+v)

S <- rep(0,days/dt)
I <- rep(0,days/dt)
H <- rep(0,days/dt)
R <- rep(0,days/dt)

Dead <- rep(0,days/dt)

S[1] <- init_pop - init_sick
I[1] <- init_sick 
R[1] <- 0
Dead[1] <- 0 
H[1] <- 0 

for(i in 2:(days/dt)){
  S[i] <- S[i-1] - (alpha*I[i-1]*S[i-1])*dt
  I[i] <- I[i-1] + ((alpha*I[i-1]*S[i-1]) - (v_nh+h)*I[i-1])*dt
  H[i] <- H[i-1] + (h*I[i-1] - (c+v_h)*H[i-1])*dt 
  R[i] <- R[i-1] + (v_nh*I[i-1] + v_h*H[i-1])*dt
  Dead[i] <- Dead[i-1] + (c*H[i-1])*dt
}

plot(seq(1,days/dt,1),I)
plot(seq(1,days/dt,1),H)


I[1:10]
plot(seq(1,days/dt,1),Dead)
Dead[1:10]
plot(S,R)
