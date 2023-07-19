rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tsiR)


homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)


#first, population data for madagascar
pop.dat <- read.csv(file=paste0(homewd, "/data/WorldBankMada.csv"), header = T, stringsAsFactors = F)
head(pop.dat)
pop.vec <- pop.dat[2,5:ncol(pop.dat)]
names(pop.vec) <- seq(2009,2022,1) #assume this is census at the end of each year,


#so one timestep before gives you the population at the beginning of the year
pop.vec <- c(unlist(pop.vec[which(names(pop.vec)=="2011"):which(names(pop.vec)=="2021")]))

#load pop.vector from the region
pop.vec <- read.csv(file=paste0(homewd, "/data/catchment_pop_by_year.csv"), header = T, stringsAsFactors = F)
#load pop.vector from the region
pop.vec <- pop.vec$tot_pop[2:length(pop.vec$year)]
names(pop.vec) <- 2011:2021


# this is at the national level - now scale down to just 1 hospital catchment in Antananarivo,
# assuming ~1% of the country pop and constant through the time series
#pop.vec <- pop.vec*.001

birth.vec <- pop.dat[1,5:ncol(pop.dat)]
names(birth.vec) <- seq(2009,2022,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="2011"):which(names(birth.vec)=="2021")]
#birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year


#now scale up by population size to get total births per year
birth.vec = c(unlist(birth.vec*(pop.vec/1000)))

#and load case data - first the age-structured
dat <- read.csv(file=paste0(homewd,"/data/rsv_data_update_2011_2022.csv"))

head(dat)
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date", 'hospital')
dat <- dplyr::select(dat, -(X))
unique(dat$DDN)
dat$DDN <- as.Date(dat$DDN,  format = "%m/%d/%y")
sort(unique(dat$DDN))
unique(dat$sampling_date)


#correct DDN for those with dates over 2022
dat$DDN[dat$DDN>as.Date("2021-12-31") & !is.na(dat$DDN)] <- as.Date(paste0((as.numeric(sapply(strsplit(as.character(dat$DDN[dat$DDN>as.Date("2021-12-31") & !is.na(dat$DDN)]), split="-"),'[',1))-100), "-", sapply(strsplit(as.character(dat$DDN[dat$DDN>as.Date("2021-12-31")& !is.na(dat$DDN)]), split="-"),'[',2), "-", sapply(strsplit(as.character(dat$DDN[dat$DDN>as.Date("2021-12-31")& !is.na(dat$DDN)]), split="-"),'[',3)))

dat$sampling_date <- as.Date(dat$sampling_date, format = "%m/%d/%y")

head(dat)
unique(dat$hospital)
dat$hospital[dat$hospital=="CENHOSOA "] <- "CENHOSOA"
dat$hospital[dat$hospital=="CSMI TSL"] <- "CSMI-TSL"
nrow(subset(dat, hospital== "")) #25 with no hospital ID
#remove from dataset
dat = subset(dat, hospital!="") #4152
#now plot cases by time and fit a gam describing seasonality

#calculate the day of year
dat$doy <- yday(dat$sampling_date)


dat <- dplyr::select(dat, num_viro, sampling_date, RSV)
names(dat) <- c("num_viro", "date", "RSV")

#get epieek and sum by that
dat$epiweek <-  cut.Date(dat$date, breaks = "week")
dat$year <- year(as.Date(dat$epiweek))
dat$doy <- yday(as.Date(dat$epiweek))

dat$time = dat$year + round(dat$doy/365,2)
head(dat) 
tail(dat)

dat <- arrange(dat, date)
dat$week <- week(as.Date(dat$epiweek))
head(dat) 
tail(dat) #stops end of June 2022



sub.tsir.dat <- ddply(dat, .(year, epiweek, week, time), summarise, cases = sum(RSV))
head(sub.tsir.dat)

sub.tsir.dat = subset(sub.tsir.dat, year>2010 & year<2022)

sub.tsir.dat <- arrange(sub.tsir.dat, time)

#need to fill in the 
sub.tsir.dat$cases[sub.tsir.dat$cases==0] <- 1

#there are gaps in the time series
ggplot(sub.tsir.dat) + geom_point(aes(x=week, y=cases, color=year, group=year)) + 
  geom_line(aes(x=week, y=cases, color=year, group=year)) + ylab("cases in 2 hospital catchment") +
  facet_wrap(~year)


tsir.year <- dlply(sub.tsir.dat, .(year))

get.tsir.data <- function(df, births, pop){
  intbirths <- rep(births[names(births)==unique(df$year)],52)/52
  intpop <- rep(pop[names(pop)==unique(df$year)]-births[names(births)==unique(df$year)], 52) 
  intpop <- intpop+ cumsum(intbirths)
  
  
  if(nrow(df)==52){
    df$week <- 1:52
    df$time <- trunc(df$time) + round(cumsum(rep(.019,52)),2)
    df <- dplyr::select(df, -(epiweek))
    df$births = intbirths
    df$pop = intpop
  }else{
    #interpolate
    intcase = round(approx(x=(df$time-trunc(df$time)), y=df$cases, xout = round(cumsum(rep(.019,52)),2))$y,0)
    df <- cbind.data.frame(year=rep(unique(df$year), 52), week=1:52, time= (unique(df$year)+round(cumsum(rep(.019,52)),2)), cases=intcase, births=intbirths, pop=intpop)
  }
  
  df <- dplyr::select(df, year, week, time, cases, births,pop)
  return(df)
}

tsir.dat <- data.table::rbindlist(lapply(tsir.year, get.tsir.data, births=birth.vec, pop=pop.vec))

#tsir.dat = tsiRdata(time=sub.tsir.dat.new$time, cases = sub.tsir.dat.new$cases, births=birth.vec, pop = pop.vec, IP=1)

head(tsir.dat)
#correct by reporting hospital - no need if all two hospitals (no 2022)


ggplot(tsir.dat) + geom_point(aes(x=time, y=cases)) + geom_line(aes(x=time, y=cases)) + ylab("cases")

tsir.dat$cases[is.na(tsir.dat$cases)]<- 1

tsir.dat$year <- as.factor(trunc(tsir.dat$time))
tsir.dat$week <- tsir.dat$time- as.numeric(as.character(tsir.dat$year))
head(tsir.dat)
ggplot(tsir.dat) + geom_point(aes(x=week, y=cases, color=year, group=year)) + 
  geom_line(aes(x=week, y=cases, color=year, group=year)) + ylab("weekly cases") #consider replotting figure 1b


ggplot(tsir.dat) + geom_point(aes(x=week, y=cases, color=year, group=year)) + 
  geom_line(aes(x=week, y=cases, color=year, group=year)) + ylab("weekly cases") + facet_wrap(~year)
  

#now fit tsir model to data
# # 
  fittedpars <- estpars(data=tsir.dat,
                            IP=1, alpha=.97, sbar=NULL, xreg = "cumcases",
                            regtype='lm',family='poisson',link='log')

# fittedpars <- estpars(data=tsir.dat, 
#                       IP=1, alpha=.97, sbar=NULL, xreg = "cumcases",
#                       regtype='lm',family='gaussian')

beta.df <- cbind.data.frame(biweek = rep(1:26, each=2), week=1:52, beta = fittedpars$contact$beta, beta_low = fittedpars$contact$betalow, beta_high = fittedpars$contact$betahigh)
#beta.df <- beta.df[!duplicated(beta.df),]

ggplot(data=beta.df) + geom_point(aes(x=biweek, y=beta)) +geom_line(aes(x=biweek, y=beta)) + geom_linerange(aes(x=biweek, ymin=beta_low, ymax=beta_high))
ggplot(data=beta.df) + geom_point(aes(x=week, y=beta)) +geom_line(aes(x=week, y=beta)) + geom_linerange(aes(x=week, ymin=beta_low, ymax=beta_high))
plot(fittedpars$contact$beta, type="b")

simfitted <- simulatetsir(data=tsir.dat,
                           IP = 1,
                           #epidemics = "break",
                           parms=fittedpars,
                           nsim=100)

plotres(simfitted)


#and do regression with climate variables
names(beta.df)[names(beta.df)=="week"] <- "week_num"
beta.df <- dplyr::select(beta.df, -(biweek))

tsir.dat$week_num <- rep(1:52, length(unique(tsir.dat$year)))

tsir.dat <- merge(tsir.dat, beta.df, by="week_num", all.x = T)
head(tsir.dat)

#now take to climate and look for predictors with lags - save here to load later with figure 2
write.csv(tsir.dat, file = paste0(homewd, "/data/tsir_dat_beta.csv"), row.names = F)

