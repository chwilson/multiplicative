##### The Wilson SIHR model :) 
##### Simulating SIR Model with a hospitalization intermediate to
##### predict demand on health services from SARS-COV2 

### Verity et al. 2020 for onset-to-death, onset-to-discharge

days <- 360
init_pop <- 243336
init_sick <- 1000
k <- 7 # #people exposed per infected per day 
p <- 0.02 # probability of transmission from infected to susceptible 

alpha <- (k*p)
hr <- 0.05 # Hospitalization as % of infected with SARS-COV2 
IFR <- hr*0.15 # Overall CFR as product of % of hospitalized patients 
incubate <- 21 # Length of incubation before hospitalization 
stay <- 15 # Length of time of hospitalization to death or discharge 
# population 

h <- -log(1-(hr/incubate)) # rate of hospitalization of infected 
c <- -log(1-(0.2/stay)) # mortality rate of hospitalized 
v_h <- -log(1-(0.8/(stay))) # Recovery rate from hospitalized to well 
v_nh <- -log(1-((1-hr)/incubate)) # Recovery rate for non-hospitalized 
# This is rather restrictive condition because it imposes rapid recovery 
# for those who will not be hospitalized, but is only way to ensure a 
# reasonable fraction end up hospitalized. 

dt <- 1 # Can toggle as needed if numerical instability crops up 

# Key ratios 
# Partitioning of infecteds to hospital: h/v_nh 
# Partitioning of hospitalized to dead/recovered c/v_h


# R_0 <- alpha/(h+v_nh) = 2.05 (I think!) 

S <- rep(0,days/dt)
I <- rep(0,days/dt)
H <- rep(0,days/dt)
R <- rep(0,days/dt)

Dead <- rep(0,days/dt)
sum <- rep(0,days/dt)
Hcount <- rep(0,days/dt)


S[1] <- init_pop - init_sick
I[1] <- init_sick 
R[1] <- 0
Dead[1] <- 0 
H[1] <- 2 
sum[1] <- S[1] + I[1] + R[1] + Dead[1] + H[1]
Hcount[1] <- H[1]

#alpha[i-1] = 0.15 - (i/360)*0.075
#(0.15 - (i/360)*0.15
for(i in 2:(days/dt)){
  S[i] <- S[i-1] - ((alpha)*I[i-1]*S[i-1]/init_pop)*dt
  I[i] <- I[i-1] + (((alpha)*I[i-1]*S[i-1]/init_pop) - (v_nh+h)*I[i-1])*dt
  H[i] <- H[i-1] + (h*I[i-1] - (c+v_h)*H[i-1])*dt 
  R[i] <- R[i-1] + (v_nh*I[i-1] + (v_h)*H[i-1])*dt
  Dead[i] <- Dead[i-1] + (c*H[i-1])*dt # Cumulative dead 
  sum[i] <- S[i] + I[i] + H[i] + R[i] + Dead[i] # Checking conservation 
  Hcount[i] <- Hcount[i-1] + h*I[i-1] # cumulative # hospitalized
}

plot(seq(1,days/dt,1),I,type = 'lwd')
plot(seq(1,days/dt,1),H,type = 'lwd')

round(H)
plot(seq(1,days/dt,1),Dead,type = 'lwd')
plot(seq(1,days/dt,1),Dead/I,type='lwd')

plot(seq(1,75,1),Dead[1:75]/I[1:75],type='lwd')


#plot(seq(1,days/dt,1),H,type = 'lwd')
Dead[360/dt]/init_pop # 0.7% overall mortality 
Hcount[360/dt]/init_pop # 4% hospitalized
Hcount[360/dt]/(Dead[360/dt]+Hcount[360/dt]) # 84% hospitalized recover

H[1:35]

##### Scenario planning for variations in R_0 

k <- 5 # #people exposed per infected per day 
p <- 0.03 # probability of transmission from infected to susceptible 
k*p

kp <- seq(0.01,1,0.05)
peak_hosp <- rep(0,length(kp))
when_peak <- rep(0,length(kp))
par(mfcol=c(1,3))
plot(seq(1,days,1),H,type='lwd',xlab="Day",ylab="#Hospitalized")
for(l in 1:length(kp)){

days <- 360
init_pop <- 3*10^5
init_sick <- 500

alpha <- kp[l]/init_pop
hr <- 0.025 # Hospitalization as % of infected with SARS-COV2 
CFR <- hr*0.2 # Overall CFR as product of % of hospitalized patients 
incubate <- 10 # Length of incubation before hospitalization 
stay <- 14 # Length of time of hospitalization to death or discharge 
sickDay <- 14 # average time to recovery/resistant status for non_hospitalized
# population 


h <- -log(1-(hr/incubate)) # rate of hospitalization of infected 
c <- -log(1-(0.2/stay)) # mortality rate of hospitalized 
v_h <- -log(1-(0.8/stay)) # Recovery rate from hospitalized to well 
v_nh <- -log(1-((1-hr)/sickDay)) # Recovery rate for non-hospitalized 

dt <- 1 # Can toggle as needed if numerical instability crops up 

# R_o <- k*p/(c+v)
# R_0 <- k*p/(CFR+v_nh) = 1.3 (I think!) 

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

peak_hosp[l] <- max(H)
when_peak[l] <- which(H==max(H))

lines(seq(1,days,1),H)

}



plot(kp,peak_hosp,ylab = "Peak # Hospitalized", xlab = "kp product")
lines(kp,peak_hosp)
plot(kp,when_peak, ylab = "Peak day", xlab = "kp product")
lines(kp,when_peak)

par(mfcol=c(1,1))
plot(peak_hosp,when_peak)



##### Global Sensitivity/Uncertainty Analysis 

k <- 5 # #people exposed per infected per day 
p <- 0.03 # probability of transmission from infected to susceptible 



peak_hosp <- rep(0,length(kp))
when_peak <- rep(0,length(kp))
par(mfcol=c(1,3))

plot(seq(1,days,1),H,type='lwd',xlab="Day",ylab="#Hospitalized")
for(l in 1:length(kp)){
  
  days <- 360
  init_pop <- 3*10^5
  init_sick <- 500
  
  alpha <- kp[l]/init_pop
  hr <- 0.025 # Hospitalization as % of infected with SARS-COV2 
  CFR <- hr*0.2 # Overall CFR as product of % of hospitalized patients 
  incubate <- 10 # Length of incubation before hospitalization 
  stay <- 14 # Length of time of hospitalization to death or discharge 
  sickDay <- 14 # average time to recovery/resistant status for non_hospitalized
  # population 
  
  
  h <- -log(1-(hr/incubate)) # rate of hospitalization of infected 
  c <- -log(1-(0.2/stay)) # mortality rate of hospitalized 
  v_h <- -log(1-(0.8/stay)) # Recovery rate from hospitalized to well 
  v_nh <- -log(1-((1-hr)/sickDay)) # Recovery rate for non-hospitalized 
  
  dt <- 1 # Can toggle as needed if numerical instability crops up 
  
  # R_0 <- k*p/(h+v_nh) = 1.3 (I think!) 
  
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
  
  peak_hosp[l] <- max(H)
  when_peak[l] <- which(H==max(H))
  
  lines(seq(1,days,1),H)
  
}



plot(kp,peak_hosp,ylab = "Peak # Hospitalized", xlab = "kp product")
lines(kp,peak_hosp)
plot(kp,when_peak, ylab = "Peak day", xlab = "kp product")
lines(kp,when_peak)

