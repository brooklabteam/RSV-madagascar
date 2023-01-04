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
                                        sum_precip = sum(PRECTOTCORR),
                                        mean_H2M = mean(RH2M))

head(clim.sum)

#average the two temperatures
clim.sum$meanTemp <- rowMeans(cbind(clim.sum$meanTempMin, clim.sum$meanTempMax))

#drop the other min/max temp
clim.sum <- dplyr::select(clim.sum, -(meanTempMin), -(meanTempMax))
#now merge with your case data
merge.dat <- merge(dat.sum, clim.sum, by="epiweek")

head(merge.dat)

#melt your data to plot the climate with the case data
merge.melt <- melt(merge.dat, id.vars = c("epiweek"))

head (merge.melt)
#names(merge.melt)[names(merge.melt)=="variable"] <- "climate_variable"

# and plot these together - temp and cases
case.dat1 = subset(merge.melt, variable=="cases_by_hospital")
case.dat1$variable <- "meanTemp"
case.dat2 = subset(merge.melt, variable=="cases_by_hospital")
case.dat2$variable <- "sum_precip"
case.dat3 = subset(merge.melt, variable=="cases_by_hospital")
case.dat3$variable <- "mean_H2M"

case.dat <- rbind(case.dat1, case.dat2, case.dat3)
head(case.dat)
unique(case.dat$variable)
Fig2left <- ggplot(data=subset(merge.melt, variable=="meanTemp" |variable=="sum_precip" |variable=="mean_H2M")) + 
  geom_line(data=case.dat, aes(x=epiweek, y=value),  size=1, alpha=.2) +
  geom_point(data=case.dat, aes(x=epiweek, y=value), size=3, alpha=.2) +
  geom_point(aes(x=epiweek, y=value, color=variable),  size=3, show.legend = F) +
  geom_line(aes(x=epiweek, y=value, color=variable),  size=1, show.legend = F) +
  facet_grid(variable~., scales = "free") + ylim(c(0, NA)) +
  theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                     legend.title = element_blank(),
                     axis.title = element_blank(), 
                     strip.text = element_text(size=14),
                     strip.background = element_rect(fill="white"),
                     legend.text = element_text(size=12),
                     plot.margin = unit(c(.2,.1,1.3,1.1), "lines"),
                     axis.text = element_text(size=14))


#and look for the cross correlations
#plot and get the printout
print(ccf(merge.dat$sum_precip, merge.dat$cases_by_hospital))
# 95% CI at 0.09


#save as data
dat.lag <- cbind.data.frame(lag = print(ccf(merge.dat$sum_precip, merge.dat$cases_by_hospital))$lag, acf=print(ccf(merge.dat$sum_precip, merge.dat$cases_by_hospital))$acf)
dat.lag$variable <- "sum_precip"

#what is the optimal lag?
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)]
# -3 is maximized cross correlation: 
# so precip precedes cases by 3 epiweeks

# and what about temp ?
dat2 = cbind.data.frame(lag = print(ccf(merge.dat$meanTemp, merge.dat$cases_by_hospital))$lag, acf=print(ccf(merge.dat$meanTemp, merge.dat$cases_by_hospital))$acf)
dat2$variable <- "meanTemp"
dat2$lag[dat2$acf==max(dat2$acf)]
# -6 is maximized cross correlation
# cases follow temp by 6 weeks
# 95% CI at 0.09


#and humidity
dat3 = cbind.data.frame(lag = print(ccf(merge.dat$mean_H2M, merge.dat$cases_by_hospital))$lag, acf=print(ccf(merge.dat$mean_H2M, merge.dat$cases_by_hospital))$acf)
dat3$variable <- "mean_H2M"
dat3$lag[dat3$acf==max(dat3$acf)]
# 1 only
# cases follow mean H2M by 1 week
# 95% CI at 0.09


#save together
dat.lag <- rbind(dat.lag, dat2, dat3)

write.csv(dat.lag, file=paste0(homewd, "/data/lag_output.csv"), row.names = F)


#and plot acf
#include the optimal lag on plot
max.lag <- dlply(dat.lag, .(variable))
get.lag <- function(df){
  lag = df$lag[df$acf==max(df$acf)]
  df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
  return(df.out)
}
max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))
max.lag$label = paste0("lag=", max.lag$lag, " epiwks")

Fig2right <- ggplot(dat.lag) + geom_label(data=max.lag, aes(x=18,y=.4, label=label), label.size = 0) +
             geom_bar(aes(x=lag, y=acf), stat = "identity") + ylim(c(NA,.45)) +
             geom_hline(aes(yintercept=0.09), color="blue", linetype=2) +
             geom_hline(aes(yintercept=-0.09), color="blue", linetype=2) +
             facet_grid(variable~.) + theme_bw() + theme(legend.position = c(.2,.87), panel.grid = element_blank(),
                                                         legend.title = element_blank(),
                                                         axis.title = element_text(size=16),
                                                         strip.background = element_rect(fill="white"),
                                                         strip.text = element_text(size=14),
                                                         legend.text = element_text(size=12),
                                                         plot.margin = unit(c(.2,.1,.1,1.1), "lines"),
                                                         axis.text = element_text(size=14))


Fig2 <- cowplot::plot_grid(Fig2left, Fig2right, rel_widths = c(1,1.1), nrow = 1, ncol = 2, labels = c("A", "B"), label_size = 22)


ggsave(file = paste0(homewd, "/figures/Fig2.png"),
       plot = Fig2,
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
