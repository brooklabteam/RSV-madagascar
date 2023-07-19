rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(mgcv)
library(reshape2)


homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)

# Figure 3
# calculate the lag of different climate predictors vs. the 
# case data. First, you will need to construct a time series, 
# with case data and the climate predictors for the same 
# timestep.

# load case data
dat <- read.csv(file=paste0(homewd,"/data/Tsiry_RSV_data_2011-2022_Update.csv"))

head(dat)

#clean and resort
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date")
dat <- dplyr::select(dat, -(X))
dat$DDN <- as.Date(dat$DDN)
dat$sampling_date <- as.Date(dat$sampling_date)

#calculate the epidemic week (week + year)
dat$epiweek <- as.Date(as.character(cut.Date(dat$sampling_date, breaks="week")))

#sum cases by month
dat.sum <- ddply(dat, .(epiweek), summarise, cases = sum(RSV), tested=length(RSV)) 

#calculate prevalence by month
dat.sum$prevalence <- dat.sum$cases/dat.sum$tested

#calculate year of sampling
dat.sum$year <- year(dat.sum$epiweek)

dat.sum <- subset(dat.sum, year>2010 & year<2022)

head(dat.sum)
length(dat.sum$epiweek[dat.sum$year==2011])

#load tsir data
tsir.dat <- read.csv(file = paste0(homewd, "/data/tsir_dat_beta.csv"), header = T, stringsAsFactors = F)

#now load the climate data and bring it in
#load the climate data and transform dates for epidemic year prior
clim.dat <- read.csv(file = paste0(homewd, "/data/POWER_Point_Daily_20110101_20220731_Tsiry.csv"))
head(clim.dat)

#make a date column
clim.dat$MO[clim.dat$MO<10] <- paste0("0", clim.dat$MO[clim.dat$MO<10])
unique(clim.dat$MO)
clim.dat$D[clim.dat$D<10] <- paste0("0", clim.dat$D[clim.dat$D<10])
clim.dat$date <- as.Date(paste0(clim.dat$YEAR, "-", clim.dat$MO, "-", clim.dat$DY))

head(clim.dat)

#get epiweek for the climate data
clim.dat$epiweek <- as.Date(cut.Date(clim.dat$date, breaks="week"))

#summarise by epiweek
clim.sum <- ddply(clim.dat, .(epiweek), summarise,
                                        meanTempMin = mean(T2M_MIN), 
                                        meanTempMax=mean(T2M_MAX), 
                                        sum_precip = sum(PRECTOTCORR),
                                        mean_H2M = mean(RH2M))

head(clim.sum)

#average the two temperatures
clim.sum$meanTemp <- rowMeans(cbind(clim.sum$meanTempMin, clim.sum$meanTempMax))

#drop the other min/max temp
clim.sum <- dplyr::select(clim.sum, -(meanTempMin), -(meanTempMax))

#make time column and interoplote
clim.sum$year <- year(clim.sum$epiweek)
clim.sum$doy <- yday(clim.sum$epiweek)
clim.sum$time <- clim.sum$year+(clim.sum$doy/365)

merge.dat <- tsir.dat
merge.dat$mean_H2M <- approx(x=clim.sum$time, y=clim.sum$mean_H2M, xout = dat.sum$time)$y
merge.dat$sum_precip <- approx(x=clim.sum$time, y=clim.sum$sum_precip, xout = dat.sum$time)$y
merge.dat$meanTemp <- approx(x=clim.sum$time, y=clim.sum$meanTemp, xout = dat.sum$time)$y


head(merge.dat)

# #melt your data to plot the climate with the case data
# merge.melt <- melt(merge.dat, id.vars = c("epiweek"))
# 
# head (merge.melt)
# #names(merge.melt)[names(merge.melt)=="variable"] <- "climate_variable"
# 
# # and plot these together - temp and cases
# case.dat1 = subset(merge.melt, variable=="cases")
# case.dat1$variable <- "meanTemp"
# case.dat2 = subset(merge.melt, variable=="cases")
# case.dat2$variable <- "sum_precip"
# case.dat3 = subset(merge.melt, variable=="cases")
# case.dat3$variable <- "mean_H2M"
# 
# 
# 
# case.dat <- rbind(case.dat1, case.dat2, case.dat3)
# head(case.dat)
# 
# case.dat$variable[case.dat$variable=="sum_precip"]<-  "sum precipitation (mm)"
# case.dat$variable[case.dat$variable=="meanTemp"] <- "mean temp (*C)"
# case.dat$variable[case.dat$variable=="mean_H2M"] <- "mean humidity (H2M)"
# 
# merge.melt$variable<- as.character(merge.melt$variable)
# merge.melt$variable[merge.melt$variable=="sum_precip"]<-  "sum precipitation (mm)"
# merge.melt$variable[merge.melt$variable=="meanTemp"] <- "mean temp (*C)"
# merge.melt$variable[merge.melt$variable=="mean_H2M"] <- "mean humidity (H2M)"
# 
# 
# unique(case.dat$variable)
# ptmp <- ggplot(data=subset(merge.melt, variable=="mean temp (*C)" |variable=="sum precipitation (mm)" |variable=="mean humidity (H2M)")) + 
#   geom_line(data=case.dat, aes(x=epiweek, y=value),  size=1, alpha=.2) +
#   geom_point(data=case.dat, aes(x=epiweek, y=value), size=3, alpha=.2) +
#   geom_point(aes(x=epiweek, y=value, color=variable),  size=3, show.legend = F) +
#   geom_line(aes(x=epiweek, y=value, color=variable),  size=1, show.legend = F) +
#   facet_grid(variable~., scales = "free") + ylim(c(0, NA)) +
#   theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
#                      legend.title = element_blank(),
#                      axis.title = element_blank(), 
#                      strip.text = element_text(size=14),
#                      strip.background = element_rect(fill="white"),
#                      legend.text = element_text(size=12),
#                      plot.margin = unit(c(.2,.1,1.3,1.1), "lines"),
#                      axis.text = element_text(size=14))
# 


#and look at climate through time

head(merge.dat)

#precip - increasing
gam1 <- gam(sum_precip~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam1) # slight pos trend
plot_model(gam1, type="pred", grid=T)

#humidity - also increasing
gam2 <- gam(mean_H2M~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam2) # slight pos tred
plot_model(gam2, type="pred", grid=T)

#temp - increasing
gam3 <- gam(meanTemp~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam3) # pos trend
plot_model(gam3, type="pred", grid=T)


#what about cases??? decreasing - 
#contradictory effects of increasing humidity (driving cases down) and increasing precip (driving cases up)
gam4 <- gam(cases~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam4) # neg trend
plot_model(gam4, type="pred", grid=T)

#plot each by year
clim.dat <- dplyr::select(merge.dat, -(births), -(pop), -(week), -(beta_low), -(beta_high))
clim.melt <- melt(clim.dat, id.vars = c("time", "year", "week_num"))
head(clim.melt)
clim.melt$year <- as.factor(clim.melt$year)
clim.melt$label <- as.character(clim.melt$variable)
clim.melt$label[clim.melt$label=="sum_precip"]<-  "sum precipitation (mm)"
clim.melt$label[clim.melt$label=="meanTemp"] <- "mean temp (*C)"
clim.melt$label[clim.melt$label=="mean_H2M"] <- "mean humidity (H2M)"

newFig2<- ggplot(data=subset(clim.melt, label!="beta")) + geom_line(aes(x=week_num, y=value, color=year)) + 
      facet_grid(label~., scales = "free_y") + 
      theme_bw() + theme(legend.position = "right", panel.grid = element_blank(),
      legend.title = element_blank(), axis.title = element_blank(), strip.text = element_text(size=14),
      strip.background = element_rect(fill="white"),legend.text = element_text(size=12),
      plot.margin = unit(c(.2,.1,1.3,1.1), "lines"), axis.text = element_text(size=14))
#and look for the cross correlations
#plot and get the printout

#beta lag - no lag really for the transmission rate -see below
#so instead we used the exact present timestep, consistent with Rachel
#regression of humidity and precip and temp on beta

head(merge.dat)
merge.dat$log_beta <- log(merge.dat$beta)

glm1 <- lm(log_beta~mean_H2M + sum_precip + meanTemp + year, data= merge.dat)
summary(glm1)
#tried year as a random effect and got nowhere.

#check best fit model
library(MuMIn)
glm1 <- lm(log_beta~mean_H2M + sum_precip + meanTemp + year, data= merge.dat, na.action = na.fail)
summary(glm1)
dredge(global.model = glm1) #best fit is just humidity and precip

glm2<- lm(log_beta~mean_H2M + sum_precip, data= merge.dat, na.action = na.fail)
summary(glm2)
AIC(glm1,glm2) #glm2 is better

library(sjPlot)

plot_model(glm2, type="est")

plot_model(glm2, type="pred", grid=T) #transmission decreases with increasing humidity and increases with increasing precip



#############
#climate lags below


#beta lag - no lag really for the transmission rate

print(ccf(merge.dat$sum_precip, merge.dat$beta))
# 95% CI at 0.09


#save as data
dat.lag <- cbind.data.frame(lag = print(ccf(merge.dat$sum_precip, merge.dat$beta))$lag, acf=print(ccf(merge.dat$sum_precip, merge.dat$beta))$acf)
dat.lag$variable <- "sum_precip"

#what is the optimal lag?
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)]
# 20 is maximized cross correlation: 
# so precip follows beta 20 epiweeks - no sig

# and what about temp ?
dat2 = cbind.data.frame(lag = print(ccf(merge.dat$meanTemp, merge.dat$beta))$lag, acf=print(ccf(merge.dat$meanTemp, merge.dat$beta))$acf)
dat2$variable <- "meanTemp"
dat2$lag[dat2$acf==max(dat2$acf)]
# -22 and -21 is maximized cross correlation
# cases follow temp by 22 weeks
# 95% CI at 0.09


#and humidity
dat3 = cbind.data.frame(lag = print(ccf(merge.dat$mean_H2M, merge.dat$beta))$lag, acf=print(ccf(merge.dat$mean_H2M, merge.dat$beta))$acf)
dat3$variable <- "mean_H2M"
dat3$lag[dat3$acf==max(dat3$acf)]
# 24 only
# transmission precedes  H2M by 24 week
# 95% CI at 0.09


#save together
dat.lag <- rbind(dat.lag, dat2, dat3)

#write.csv(dat.lag, file=paste0(homewd, "/data/lag_output.csv"), row.names = F)


#and plot acf
#include the optimal lag on plot
max.lag <- dlply(dat.lag, .(variable))
get.lag <- function(df){
  lag = min(abs(df$lag[df$acf==max(df$acf)]))
  df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
  return(df.out)
}
max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))



max.lag$label <- max.lag$variable
max.lag$label[max.lag$label=="sum_precip"]<-  "sum precipitation (mm)"
max.lag$label[max.lag$label=="meanTemp"] <- "mean temp (*C)"
max.lag$label[max.lag$label=="mean_H2M"] <- "mean humidity (H2M)"








#cases lag
print(ccf(merge.dat$sum_precip, merge.dat$cases))
# 95% CI at 0.09


#save as data
dat.lag <- cbind.data.frame(lag = print(ccf(merge.dat$sum_precip, merge.dat$cases))$lag, acf=print(ccf(merge.dat$sum_precip, merge.dat$cases))$acf)
dat.lag$variable <- "sum_precip"

#what is the optimal lag?
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)]
# -3 is maximized cross correlation: 
# so precip precedes cases by 3 epiweeks

# and what about temp ?
dat2 = cbind.data.frame(lag = print(ccf(merge.dat$meanTemp, merge.dat$cases))$lag, acf=print(ccf(merge.dat$meanTemp, merge.dat$cases))$acf)
dat2$variable <- "meanTemp"
dat2$lag[dat2$acf==max(dat2$acf)]
# -6 is maximized cross correlation
# cases follow temp by 6 weeks
# 95% CI at 0.09


#and humidity
dat3 = cbind.data.frame(lag = print(ccf(merge.dat$mean_H2M, merge.dat$cases))$lag, acf=print(ccf(merge.dat$mean_H2M, merge.dat$cases))$acf)
dat3$variable <- "mean_H2M"
dat3$lag[dat3$acf==max(dat3$acf)]
# 1 only
# cases follow mean H2M by 1 week
# 95% CI at 0.09


#save together
dat.lag <- rbind(dat.lag, dat2, dat3)

#write.csv(dat.lag, file=paste0(homewd, "/data/lag_output.csv"), row.names = F)


#and plot acf
#include the optimal lag on plot
max.lag <- dlply(dat.lag, .(variable))
get.lag <- function(df){
  lag = min(abs(df$lag[df$acf==max(df$acf)]))
  df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
  return(df.out)
}
max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))



max.lag$label <- max.lag$variable
max.lag$label[max.lag$label=="sum_precip"]<-  "sum precipitation (mm)"
max.lag$label[max.lag$label=="meanTemp"] <- "mean temp (*C)"
max.lag$label[max.lag$label=="mean_H2M"] <- "mean humidity (H2M)"

max.lag$x=47
max.lag$y=90
max.lag$y[max.lag$variable=="meanTemp"] <- 12
max.lag$y[max.lag$variable=="sum_precip"] <- 350
max.lag$text = paste0("optimal lag=", max.lag$lag, " wks")
max.lag$text[max.lag$variable=="mean_H2M"] <- "optimal lag concurrent"

newFig2<- ggplot(data=clim.melt) + geom_line(aes(x=week, y=value, color=year)) + 
  geom_label(data=max.lag, aes(x=x,y=y, label=text), label.size = 0) +
  facet_grid(label~., scales = "free_y") + 
  theme_bw() + theme(legend.position = "right", panel.grid = element_blank(),
                     legend.title = element_blank(), axis.title = element_blank(), strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"),legend.text = element_text(size=12),
                     plot.margin = unit(c(.2,.1,1.3,1.1), "lines"), axis.text = element_text(size=14))
#and loo

ggsave(file = paste0(homewd, "/figures/Fig2-new.png"),
       plot = newFig2,
       units="mm",  
       width=90, 
       height=80, 
       scale=3, 
       dpi=300)


head(merge.dat)

#make a new dataset with the lagged versions of these climate variables
# we can't start cases at timestep 1 because they are only predicted by all three variables starting
# at timestep 7 (because the lag for temp is 6).

#temp first
merge.shift <- merge.dat[7: length(merge.dat$epiweek),]
merge.shift$meantempLag <- merge.dat$meanTemp[1:(length(merge.dat$meanTemp)-6)]

#humid stays as is

#lag precip by 3
merge.shift$precip_lag <- merge.dat$sum_precip[4:(length(merge.dat$sum_precip)-3)]


head(merge.shift)

#now, make a linear model predicting cases
merge.shift$n_hospitals <- 1
merge.shift$n_hospitals[merge.shift$epiweek>="2020-01-01" & merge.shift$epiweek< "2022-08-01"] <- 2
merge.shift$n_hospitals[merge.shift$epiweek>= "2022-08-01"] <- 3

#save new lagged data
write.csv(merge.shift, file = paste0(homewd, "/data/lagged-dat.csv"), row.names = F)

#first try no predictor by year

m1 <- glm(cases~n_hospitals + mean_H2M + meantempLag + precip_lag, data = merge.shift, family = "poisson")
summary(m1)

library(sjPlot)

#then try with year as a random effect
merge.shift$year <- as.factor(merge.shift$year)
library(lme4)

m2 <- glmer(cases~n_hospitals + mean_H2M + meantempLag + precip_lag + (1|year),data = merge.shift, family = "poisson")
summary(m2)

AIC(m1, m2) #m2 is better

#ideally, you would add package sjPlot
#plot model output

plot_model(m2)
