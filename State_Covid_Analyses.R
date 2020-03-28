# Florida Covid Analyses 

FLdata <- read.table("covid_ts_FL_032620.txt",header=T,sep=",")
str(FLdata)
FLdata$days <- round(FLdata$seconds_since_Epoch/(24*3600)) - min(round(FLdata$seconds_since_Epoch/(24*3600)))

ggplot(FLdata, aes(x=days,y=positive)) + geom_point()

### First testing elasticity 

summary(lm(log(positive)~log(tested),FLdata)) # 1.04 

# a 1% increase in tests is associated with a 1.04% increase in positive results 

summary(lm(log(positive)~ days + log(tested),FLdata)) #

# Most of the increase explained by temporal trend 
# day coef = 0.23 


summary(lm(log(positive)~ days,FLdata)) #
# day coef = 0.29 

# So, I would wager that the "true" increase in more like 0.23 
# which is a doubling time of 3 days 


# Now let's look at deaths
ggplot(FLdata, aes(x=days,y=deaths)) + geom_point()
# We truncate to day 4 since the early days are static
FLdataMort <- FLdata %>% filter(days>3)

summary(lm(log(deaths)~days,FLdataMort))
# Here we extracted an exponential coefficient of 0.18 
# doubling time: log(2)/0.18 = 3.85 days 

# Given lags in how long system would respond to even a drastic intervention, 
# I would feel comfortable extrapolating deaths out 4 weeks, or 28 days
curve(28*exp(0.18*x),1,28) 

# This gives us a MINIMUM of ~ 4K deaths statewide due to current Covid outbreak 
4000/(20*10^6) # Death rate of 0.0002 or 2 in 10K residents 

# 46 deaths in Alachua County. assuming 1% IFR, and 5% ICU rate - roughly 20% of ICU 
# patients have died, implying that 46*5 = 230 patients are in need of ICU in 
# alachua county within the month 


##### MI data 
MIdata <- read.table("covid_ts_MI_032620.txt",header=T,sep=",")
MIdata$days <- round(MIdata$seconds_since_Epoch/(24*3600)) - min(round(MIdata$seconds_since_Epoch/(24*3600)))
ggplot(MIdata, aes(x=days, y = positive)) + geom_point()

ggplot(MIdata, aes(x=days, y = deaths)) + geom_point()

summary(lm(log(deaths)~days,MIdata[-c(1:12),]))
exp(0.51*1.4)
curve(60*exp(0.51*x),1,7)

###### NY data 

NYdata <- read.table("covid_ts_NY_032620.txt",header=T,sep=",")
NYdata$days <- round(NYdata$seconds_since_Epoch/(24*3600)) - min(round(NYdata$seconds_since_Epoch/(24*3600)))
ggplot(NYdata, aes(x=days, y = positive)) + geom_point()

ggplot(NYdata, aes(x=days, y = deaths)) + geom_point()

summary(lm(log(deaths)~days,NYdata[-c(1:11),]))
# JFC: log(2)/0.44 - 1.6 days doubling rate 





