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

#add the cases by hospital
# 1 before 2020
dat.sum$cases_by_hospital <- dat.sum$cases
# 2 from 2020 to end July 2022
dat.sum$cases_by_hospital[dat.sum$epiweek>"2020-01-01" & dat.sum$epiweek<"2022-07-31"] <- dat.sum$cases_by_hospital[dat.sum$epiweek>"2020-01-01" & dat.sum$epiweek<="2022-07-31"]/2
# 3 from August 2022 til now
dat.sum$cases_by_hospital[dat.sum$epiweek>"2022-07-31"] <- dat.sum$cases_by_hospital[dat.sum$epiweek>"2022-07-31"]/3

#look at your data
head(dat.sum)

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
                                        mean_precip = mean(PRECTOTCORR),
                                        mean_H2M = mean(RH2M))

head(clim.sum)

#now merge with your case data
merge.dat <- merge(dat.sum, clim.sum, by="epiweek")

head(merge.dat)

#melt your data to plot the climate with the case data
merge.melt <- melt(merge.dat, id.vars = c("epiweek"))

head (merge.melt)
#names(merge.melt)[names(merge.melt)=="variable"] <- "climate_variable"

# and plot these together - temp and cases
p1 <- ggplot(data=subset(merge.melt, variable=="meanTempMin" | variable=="cases_by_hospital")) + 
      geom_point(aes(x=epiweek, y=value, color=variable),  size=3) +
      geom_line(aes(x=epiweek, y=value, color=variable),  size=1) +
      scale_color_manual(values = c("tomato", "cornflowerblue")) +
      theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                         legend.title = element_blank(),
                         axis.title = element_blank(), 
                         legend.text = element_text(size=12),
                         plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                         axis.text = element_text(size=14))

#humidity and cases
p2 <- ggplot(data=subset(merge.melt, variable=="mean_H2M" | variable=="cases_by_hospital")) + 
  geom_point(aes(x=epiweek, y=value, color=variable),  size=3) +
  geom_line(aes(x=epiweek, y=value, color=variable),  size=1) +
  scale_color_manual(values = c("tomato", "purple")) + ylim(c(0,150)) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     legend.text = element_text(size=12),
                     axis.title = element_blank(), 
                     plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                     axis.text = element_text(size=14))

#and precip and cases
p3 <- ggplot(data=subset(merge.melt, variable=="mean_precip" | variable=="cases_by_hospital")) + 
  geom_point(aes(x=epiweek, y=value, color=variable),  size=3) +
  geom_line(aes(x=epiweek, y=value, color=variable),  size=1) +
  scale_color_manual(values = c("tomato", "forestgreen")) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     legend.text = element_text(size=12),
                     axis.title = element_blank(), 
                     plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                     axis.text = element_text(size=14))


#and all together
Fig3 <- cowplot::plot_grid(p1,p2,p3, ncol=1, nrow = 3, labels = c("A", "B", "C"))

#save 


ggsave(file = paste0(homewd, "/figures/Fig3.png"),
       plot = Fig3,
       units="mm",  
       width=60, 
       height=100, 
       scale=3, 
       dpi=300)

#and look for the cross correlations
#plot and get the printout
print(ccf(merge.dat$mean_precip, merge.dat$cases_by_hospital))

#what is the optimal lag
out.precip = print(ccf(merge.dat$mean_precip, merge.dat$cases_by_hospital))
out.precip$lag[which(out.precip$acf==max(out.precip$acf))]
# -3 is maximized cross correlation: 
# so cases follow precip by 3 epiweeks

# and what about temp ?
out.minTemp = print(ccf(merge.dat$meanTempMin, merge.dat$cases_by_hospital))
out.minTemp$lag[which(out.minTemp$acf==max(out.minTemp$acf))]
# -5 is maximized cross correlation
# cases follow min temp by 5 weeks

#max temp too?
out.maxTemp = print(ccf(merge.dat$meanTempMax, merge.dat$cases_by_hospital))
out.maxTemp$lag[which(out.maxTemp$acf==max(out.maxTemp$acf))]
# -14 is maximized cross correlation
# cases follow max temp by 14 weeks

#and humidity
out.H2M = print(ccf(merge.dat$mean_H2M, merge.dat$cases_by_hospital))
out.H2M$lag[which(out.H2M$acf==max(out.H2M$acf))]
# 1 only
# cases follow max temp by 1 weeks

#try shifting and plotting
merge.dat$humid_shift = c(merge.dat$mean_H2M[2:length(merge.dat$mean_H2M)], merge.dat$mean_H2M[1])
out.H2shift = print(ccf(merge.dat$humid_shift, merge.dat$cases_by_hospital))
out.H2shift$lag[which(out.H2shift$acf==max(out.H2shift$acf))] # 0! the lag is gone

#plot as regression with random effect of year
library(lme4)
merge.dat$year <- as.factor(merge.dat$year)
m1a <- lm(cases_by_hospital~mean_H2M, data=merge.dat)
summary(m1a)

m1 <- lmer(cases_by_hospital~mean_H2M + (1|year), data=merge.dat)
summary(m1)

m2a <- lm(cases_by_hospital~humid_shift, data=merge.dat)
summary(m2a)

m2 <- lmer(cases_by_hospital~humid_shift + (1|year), data=merge.dat)
summary(m2)

AIC(m1a, m1, m2a, m2) #m2 best: mixed effects with regression

ggplot(data=merge.dat) + 
  geom_point(aes(x=humid_shift, y=cases_by_hospital))
