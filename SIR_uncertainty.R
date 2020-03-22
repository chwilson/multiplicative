days <- 360
init_pop <- 243336
init_sick <- 6000
  

hr <- 0.05 # Hospitalization as % of infected with SARS-COV2 
IFR <- hr*0.15 # Overall CFR as product of % of hospitalized patients 
incubate <- 14 # Length of incubation before hospitalization 
stay <- 15 # Length of time of hospitalization to death or discharge 
# population 

sim <- 10^4


R_0 <- rnorm(sim,2.5,0.25)
h <- rnorm(sim,-log(1-(hr/incubate)),-0.25*log(1-(hr/incubate))) # rate of hospitalization of infected 
v_nh <- rnorm(sim,-log(1-((1-hr)/incubate)),-0.25*log(1-((1-hr)/incubate))) # Recovery rate for non-hospitalized 
c <- rnorm(sim,-log(1-(0.2/stay)),-0.25*log(1-(0.2/stay))) # mortality rate of hospitalized 
v_h <- rnorm(sim,-log(1-(0.8/(stay))),-0.25*log(1-(0.8/stay))) # Recovery rate from hospitalized to well 

hist(R_0)
hist(exp(-h*incubate))
hist(exp(-v_nh*incubate))


dt <- 1 # Can toggle as needed if numerical instability crops up 

# Key ratios 
# Partitioning of infecteds to hospital: h/v_nh 
# Partitioning of hospitalized to dead/recovered c/v_h


# R_0 <- alpha/(h+v_nh) = 2.05 (I think!) 

S <- matrix(0,sim,days/dt)
I <- matrix(0,sim,days/dt)
H <- matrix(0,sim,days/dt)
R <- matrix(0,sim,days/dt)

Dead <- matrix(0,sim,days/dt)
sum <- matrix(0,sim,days/dt)
Hcount <- matrix(0,sim,days/dt)


S[,1] <- init_pop - init_sick
I[,1] <- init_sick 
R[,1] <- 0
Dead[,1] <- 0 
H[,1] <- 2 
sum[,1] <- S[1] + I[1] + R[1] + Dead[1] + H[1]
Hcount[,1] <- H[1]

#alpha[i-1] = 0.15 - (i/360)*0.075
#(0.15 - (i/360)*0.15
for(i in 1:sim){
  for(j in 2:(days/dt)){
  S[i,j] <- S[i,j-1] - (R_0[i]*(v_nh[i]+h[i])*I[i,j-1]*S[i,j-1]/init_pop)*dt
  I[i,j] <- I[i,j-1] + ((R_0[i]*(v_nh[i]+h[i])*I[i,j-1]*S[i,j-1]/init_pop) - (v_nh[i]+h[i])*I[i,j-1])*dt
  H[i,j] <- H[i,j-1] + (h[i]*I[i,j-1] - (c[i]+v_h[i])*H[i,j-1])*dt 
  R[i,j] <- R[i,j-1] + (v_nh[i]*I[i,j-1] + (v_h[i])*H[i,j-1])*dt
  Dead[i,j] <- Dead[i,j-1] + (c[i]*H[i,j-1])*dt # Cumulative dead 
  sum[i,j] <- S[i,j] + I[i,j] + H[i,j] + R[i,j] + Dead[i,j] # Checking conservation 
  Hcount[i,j] <- Hcount[i,j-1] + h[i]*I[i,j-1] # cumulative # hospitalized
    }
}


plot(seq(1,days/dt,1),I[1,],type = 'lwd',ylim = c(0,150000))
for(i in 2:sim){
  lines(seq(1,days/dt,1),I[i,],size=0.00001)
}

# extracting peak hospitalized from each simulation 
peak_hosp <- apply(H,1,max)
peak_dead <- apply(Dead,1,max)

# function to standardize the inputs 

std_func <- function(x){
  x_std <- rep(0,length(x))
  for(i in 1:length(x)){
    x_std[i] = (x[i] - mean(x,na.rm=T))/sd(x,na.rm=T)
  }
  return(x_std)
}

library(dplyr)
peaks <- data.frame(hospital = peak_hosp, dead = peak_dead,
                    R_0 = std_func(R_0), v_nh = std_func(v_nh),
                    v_h = std_func(v_h), h = std_func(h),
                    c = std_func(c))

peaks <- peaks %>% filter(dead < init_pop, dead > 0, hospital < init_pop, hospital > 0,
                          h > -5, h < 5)

h_hist <- ggplot(peaks, aes(x=hospital)) + geom_histogram(stat="density") + xlab("Peak Hospitalized")

r_h <- ggplot(peaks, aes(x=R_0, y = hospital)) + geom_point() + xlab("R_0 (std)") +
  theme_bw() + ylab("Peak Hospitalized")

h_h <- ggplot(peaks, aes(x=h, y = hospital)) + geom_point() + xlab("Rate of hospitalization (std)") +
  theme_bw() + ylab("Peak Hospitalized")

vnh_h <- ggplot(peaks, aes(x=v_nh, y = hospital)) + geom_point() + xlab("Rate of recovery outside hospital (std)") +
  theme_bw() + ylab("Peak Hospitalized")

vh_h <- ggplot(peaks, aes(x=v_h, y = hospital)) + geom_point() + xlab("Rate of recovery inside hospital (std)") +
  theme_bw() + ylab("Peak Hospitalized")

c_h <- ggplot(peaks, aes(x=c, y = hospital)) + geom_point() + xlab("Mortality (in hospital,std)") +
  theme_bw() + ylab("Peak Hospitalized")


summary(lm(std_func(hospital)~R_0+v_nh+v_h + h + c,peaks))

library(cowplot)
plot_grid(r_h,h_h,vnh_h,vh_h,c_h,h_hist,nrow=2)

ggplot(peaks,aes(x=R_0,y=h)) + geom_point()
