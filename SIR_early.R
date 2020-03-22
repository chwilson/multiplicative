fev <- read.csv("alachua_fever_march.csv")
fev$exp <- seq(3.88,3.35,length=nrow(fev))
library(dplyr);library(ggplot2)
fev <- fev %>% mutate(anom_fev = Obs - exp) %>% mutate(anom = anom_fev - min(anom_fev)+0.001,
                                                       day = Day - 1)

ggplot(fev,aes(x=Day,y=anom)) + geom_point()

summary(lm(log(anom)~Day,fev))
summary(lm(anom~day,fev))


ggplot(data=fev,aes(x=Day,y=anom)) + geom_point() +
  stat_function(fun = function(x){(-0.45+ 0.24*x)})

seq(15,225,30)

28*25*205

library(lubridate)
nat_ts <- read.csv("covid_time_series.csv")
nat_ts$Date <- mdy(nat_ts$date)
nat_ts$day <- yday(nat_ts$Date) - min(yday(nat_ts$Date))

ggplot(nat_ts, aes(x=day,y=positive)) + geom_point()
summary(lm(log(positive)~day,nat_ts))

nat_ts$positive
curve(23000*exp(0.3*x),0,21)
23000*exp(0.3*21)

summary(lm(log(deaths)~day,nat_ts))
summary(lm(log(deaths/positive)~day,nat_ts))

ggplot(nat_ts, aes(x=positive,y=log(deaths))) + geom_point()

nat_ts$deaths
curve(292*exp(0.2*x),0,7)

292*exp(0.2*7)
292*exp(0.2*15)

str(nat_ts)
ggplot(nat_ts, aes(x=day,y=deaths)) + geom_point()
ggplot(nat_ts, aes(x=day,y=deaths/positive)) + geom_point()
ggplot(nat_ts, aes(x=day,y=positive/tested)) + geom_point()

# deaths/positive declining exponentially. 

# To do: simulate simple two-pool with early phase exponential growth of infecteds
# and a growing dead pool. My thought is that the accumulation of dead people is 
# the best indication of the "true" exponential growth rate of infection, assuming a 
# constant mortality rate. 

days <- 30
inf <- rep(0,days)
dead <- rep(0,days)
inf[1] <- 371 
dead[1] <- 19 
R_0 <- 4
v <- 0.12
m <- 0.003

# Need a prior on v, rate of recovery out of infectious
# and m/v, the ultimate IFR...
m/v


for(i in 2:days){
  inf[i] <- inf[i-1] + (v*(R_0-1)-m)*inf[i-1]
  dead[i] <- dead[i-1] + m*inf[i-1]
}
plot(inf[1:15])
plot(dead[1:15])

# Checking 
summary(lm(log(inf[1:15])~seq(1,15,1)))
summary(lm(log(inf)~seq(1,30,1)))
# Check! 

# Is the exponential growth rate in deaths homogeneous in time? 
summary(lm(log(dead[1:15])~seq(1,15,1)))
summary(lm(log(dead)~seq(1,30,1)))
# Answer is: NO! IN fact, it is increasing
# So...yea. All signs suggest, US death toll about to go up a LOT. 

