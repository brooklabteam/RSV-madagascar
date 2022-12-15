rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(mgcv)


homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)

# Figure 2
# Plot date of case onset on the y axis, per year and averaged
# climate parameters from the year before on the X-axis

# load data
dat <- read.csv(file=paste0(homewd,"/data/Tsiry_RSV_data_2011-2022_Update.csv"))

head(dat)

#clean and resort
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date")
dat <- dplyr::select(dat, -(X))
dat$DDN <- as.Date(dat$DDN)
dat$sampling_date <- as.Date(dat$sampling_date)

#calculate the day of year
dat$doy <- yday(dat$sampling_date)
#and year
dat$year <- year(dat$sampling_date)
sort(unique(dat$year)) #2010-2022

yday("2010-10-27") #300

# You observe that cases start to take off at the end of 
# one calendar year and then increase into early months 
# in the next year, beginning around the 300th day of year.
# So we set our epidemic year as beginning on Oct 28 of the
# year prior and ending on Oct 27 of the true year.

# So first, we restructure the data to assign a day of epidemic year.
# 2020 is a leap year
dat$doy[dat$doy>365] <- 365
dat$doy_new <- NA
# for all doy above 300, subtract 300, so these
# become the start of the following epidemic year
dat$doy_new[dat$doy>300] <- dat$doy[dat$doy>300]-300
# for those that are under day 300, add 65, so these become
# the end of the current epidemic year
dat$doy_new[is.na(dat$doy_new)] <- dat$doy[is.na(dat$doy_new)] +65
min(dat$doy_new)
max(dat$doy_new)

# now, epi year is the current year for values >65 and the current
# year +1 for values 65 and under
dat$epi_year <- NA
dat$epi_year[dat$doy_new>65] <- dat$year[dat$doy_new>65]
dat$epi_year[dat$doy_new<=65] <- dat$year[dat$doy_new<=65]+1

sort(unique(dat$epi_year)) #2010 to 2022
#but 2010 and 2022 should be incomplete
max(dat$doy_new[dat$epi_year==2022]) #246
min(dat$doy_new[dat$epi_year==2022]) #1

max(dat$doy_new[dat$epi_year==2010]) #120
min(dat$doy_new[dat$epi_year==2010]) #111

#so eliminate 2010 and 2022

dat.epi = subset(dat, epi_year>2010 & epi_year<2022)
sort(unique(dat.epi$epi_year)) #2011 - 2021

#now find date of case onset per year
dat.pos = subset(dat.epi, RSV==1)
dat.onset <- ddply(dat.pos, .(epi_year), summarise, case_onset=min(doy), sum_cases = sum(RSV))

dat.onset

dat.doy <- ddply(dat.epi, .(epi_year, doy_new), summarise, cases=sum(RSV))

# and plot
dat.onset$epi_year <- as.factor(dat.onset$epi_year)
ggplot(data=dat.onset) + geom_point(aes(x=epi_year, case_onset))

dat.doy$epi_year <- as.factor(dat.doy$epi_year)
ggplot(data=dat.doy) + geom_point(aes(x=doy_new, y=cases, color=epi_year)) + 
  geom_line(aes(x=doy_new, y=cases, color=epi_year))


#plot first date of onset against the climate data


#load the climate data and transform dates for epidemic year prior
clim.dat <- read.csv(file = paste0(homewd, "/data/POWER_Point_Daily_20110101_20220731_Tsiry.csv"))
head(clim.dat)

#make a date column
clim.dat$MO[clim.dat$MO<10] <- paste0("0", clim.dat$MO[clim.dat$MO<10])
unique(clim.dat$MO)
clim.dat$D[clim.dat$D<10] <- paste0("0", clim.dat$D[clim.dat$D<10])
clim.dat$date <- as.Date(paste0(clim.dat$YEAR, "-", clim.dat$MO, "-", clim.dat$DY))

# and transform to epidemic year data to get average from prior year
clim.dat$doy_new <- NA
clim.dat$doy <- yday(clim.dat$date)
clim.dat$doy_new[clim.dat$doy>300] <- clim.dat$doy[clim.dat$doy>300]-300
clim.dat$doy_new[is.na(clim.dat$doy_new)] <- clim.dat$doy[is.na(clim.dat$doy_new)] + 65
min(clim.dat$doy_new)
max(clim.dat$doy_new)

# and calculate epi_year - then also do epi_year predictive and epi_year minus one
clim.dat$epi_year <- NA
clim.dat$epi_year[clim.dat$doy_new>65] <- clim.dat$YEAR[clim.dat$doy_new>65]
clim.dat$epi_year[clim.dat$doy_new<=65] <- clim.dat$YEAR[clim.dat$doy_new<=65]+1

sort(unique(clim.dat$epi_year)) #20111 to 2022
#but 2011 and 2022 should be incomplete
max(clim.dat$doy_new[clim.dat$epi_year==2022]) #277
min(clim.dat$doy_new[clim.dat$epi_year==2022]) #1

max(clim.dat$doy_new[clim.dat$epi_year==2011]) #365
min(clim.dat$doy_new[clim.dat$epi_year==2011]) #66

#so eliminate 2011 and 2022
clim.epi = subset(clim.dat, epi_year>2011 & epi_year<2022)

sort(unique(clim.epi$epi_year)) #2012-2021

# and predict year (the climate in the current year predicts the 
# cases in the following year)
clim.epi$predict_year <- clim.epi$epi_year +1

#and summarise by year
clim.sum <- ddply(clim.epi, .(predict_year), summarise, 
                  min_temp = mean(T2M_MIN), max_temp = mean(T2M_MAX), 
                  preci_mean = sum(PRECTOTCORR), 
                  H2M_mean = sum(RH2M))

names(clim.sum)[names(clim.sum)=="predict_year"] <- "epi_year"

clim.sum <- ddply(clim.epi, .(epi_year), summarise, 
                  min_temp = mean(T2M_MIN), max_temp = mean(T2M_MAX), 
                  preci_mean = sum(PRECTOTCORR), 
                  H2M_mean = sum(RH2M))



#now merge with case onset data
merge.dat <- merge(dat.onset, clim.sum, by="epi_year")
head(merge.dat)

#melt
library(reshape2)
merge.melt <- melt(merge.dat, id.vars = c("epi_year", "sum_cases", "case_onset"))
#plot as regression
merge.melt$epi_year <- as.factor(merge.melt$epi_year)
ggplot(data=merge.melt) + geom_point(aes(x=value, y=case_onset, color=epi_year), size=5) +
  facet_wrap(~variable, scales = "free", ncol=2) + ylab("case onset day")


ggplot(data=subset(merge.melt, variable=="preci_sum")) + geom_point(aes(x=value, y=case_onset, color=epi_year), size=5) +
        facet_grid(~variable, scales = "free") + ylab("case onset day")

ggplot(data=subset(merge.melt, variable=="H2M_sum")) + geom_point(aes(x=value, y=case_onset, color=epi_year), size=5) +
  facet_grid(~variable, scales = "free") + ylab("case onset day")

ggplot(data=subset(merge.melt, variable=="min_temp")) + geom_point(aes(x=value, y=case_onset, color=epi_year), size=5) +
  facet_grid(~variable, scales = "free") + ylab("case onset day")

ggplot(data=subset(merge.melt, variable=="max_temp")) + geom_point(aes(x=value, y=case_onset, color=epi_year), size=5) +
  facet_grid(~variable, scales = "free") + ylab("case onset day")


# and look at total cases
ggplot(data=merge.melt) + geom_point(aes(x=value, y=sum_cases, color=epi_year), size=5) +
  facet_wrap(~variable, scales = "free", ncol=2) + ylab("total cases")


#and total cases by case onset?
ggplot(data=merge.melt) + geom_point(aes(x=case_onset, y=sum_cases, color=epi_year), size=5) +
   ylab("total cases") + xlab("day of case onset")
