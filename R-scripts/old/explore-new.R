rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)

homewd="/Users/carabrook/Developer/RSV-madagascar"
setwd(homewd)

#and load data - first the age-structured
dat <- read.csv(file=paste0(homewd,"/data/Tsiry_RSV_data_2011-2022_Update.csv"))

head(dat)
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date")
dat <- dplyr::select(dat, -(X))
dat$DDN <- as.Date(dat$DDN)
dat$sampling_date <- as.Date(dat$sampling_date)

#replace those commas
dat$age <- sub(",", ".", dat$age)
unique(dat$age)

dat$age<- as.numeric(dat$age)
unique(dat$age[!is.na(dat$DDN)])
unique(dat$DDN[!is.na(dat$age)])


#fill in the missing ages
dat$age[is.na(dat$age)] <- (dat$sampling_date[is.na(dat$age)] - dat$DDN[is.na(dat$age)])/365

#a bunch still don't have ages...
subset(dat, is.na(age))


# look at age structured prevalence by year
dat$age_round <- round(dat$age,0)
dat$year <- year(dat$sampling_date)
dat.pos <- subset(dat, RSV==1 & !is.na(age))
max(dat.pos$age_round) #75

unique(dat.pos$year)

dat.split <- dlply(dat.pos,.(year))

#and get cumulative incidence by age
get.cuminc <- function(df){
  df <- arrange(df, age)
  
  age.list = 0:75
  
  df.out <- cbind.data.frame(age=age.list)
  
  df.sum <- ddply(df,.(age_round), summarise, cases = sum(RSV))
  names(df.sum)[names(df.sum)=="age_round"] <- "age"
  
  df.out <- merge(df.out, df.sum, by="age", all.x = T)
  df.out$cases[is.na(df.out$cases)] <- 0
  
  Ntot=sum(df.out$cases)
  
  df.out$cumcases = cumsum(df.out$cases)
  
  df.out$cumprop <- df.out$cumcases/Ntot
  
  #and add
  df.out$year <- unique(df$year)
  
  return(df.out)
}

dat.split <- lapply(dat.split, get.cuminc)

dat.inc <- data.table::rbindlist(dat.split)

dat.inc$year <- as.factor(dat.inc$year)

p1 <- ggplot(data=dat.inc) + 
      geom_point(aes(x=age, y=cumprop, color=year)) +
      geom_line(aes(x=age, y=cumprop, color=year)) + facet_grid(year~.)
print(p1)


p2 <- ggplot(data=dat.inc) + theme_bw() +
  geom_point(aes(x=age, y=cumprop, color=year)) +
  geom_line(aes(x=age, y=cumprop, color=year)) +
  ylab("cumulative proportion of cases")


ggsave(file = paste0(homewd, "/figures/cum-prop-cases.png"),
       plot = p2,
       units="mm",  
       width=60, 
       height=50, 
       scale=3, 
       dpi=200)


#does the age distribution change with time? can you calculate 
#foi changes with time? 


head(dat.inc)

head(dat)

#get average age by year???
dat.age <- ddply(dat, .(year), summarise, mean_age = mean(age, na.rm=T))

ggplot(dat.age) + geom_point(aes(x=year, y=mean_age)) + geom_line(aes(x=year, y=mean_age))

# lots of "risk factor" analysis to do
# what is the odds ratio of infection by year, by sex, by age, by month?
# first, need to know purpose of sampling set
# and need to find those years where we are missing age data and DDN

#binomial
head(dat)
dat$month <- month(dat$sampling_date)
dat$RSV <- as.numeric(dat$RSV)
m1 <- glm(RSV ~year + month + age + sex, data=dat, family="binomial")
summary(m1)


#now estimate FOI


library(mgcv)
dat$year
dat$month <- as.numeric(dat$month)
dat$sex <- as.factor(dat$sex)
gam1 <- gam(RSV ~ year + 
                  s(month, bs="cc")  + 
                  #s(age, bs="tp")+ 
                  s(sex, bs="re"), 
                  data=dat, 
                  family="binomial")
summary(gam1)

plot(gam1)

df.pred <- cbind.data.frame(pred= predict.gam(gam1, type="response", exclude=c("s(month)", "s(age)", "s(sex)")))
df.pred$year <- dat$year

with(df.pred, plot(year, pred, type="l"))

# then, plot the climate data
clim.dat <- read.csv(file = paste0(homewd, "/data/POWER_Point_Daily_20110101_20220731_Tsiry.csv"))
head(clim.dat)

#make a date column
clim.dat$MO[clim.dat$MO<10] <- paste0("0", clim.dat$MO[clim.dat$MO<10])
unique(clim.dat$MO)
clim.dat$D[clim.dat$D<10] <- paste0("0", clim.dat$D[clim.dat$D<10])
clim.dat$date <- as.Date(paste0(clim.dat$YEAR, "-", clim.dat$MO, "-", clim.dat$DY))

#and get epiweek to match the data
clim.dat$epiweek = cut.Date(clim.dat$date, breaks = "week", start.on.monday = T) 

# and split by epiweek and take average temp and sum precip and humid
clim.split <- dlply(clim.dat, .(epiweek))

sum.dat <- function(df){
  avg_tmp = mean(c((df$T2M_MIN), (df$T2M_MAX)))
  sum_precip = sum(c(df$PRECTOTCORR))
  sum_H2M = sum(c(df$RH2M))
  df.out = cbind.data.frame(epiweek=unique(df$epiweek), avg_tmp=avg_tmp, sum_precip=sum_precip, sum_H2M=sum_H2M)
  return(df.out)
  
}

clim.plot <- data.table::rbindlist(lapply(clim.split, sum.dat))  
head(clim.plot)

#and melt
library(reshape2)
clim.long <- melt(clim.plot)
head(clim.long)
clim.long$epiweek <- as.Date(clim.long$epiweek)

#and plot

p3 <- ggplot(data=clim.long) + geom_line(aes(x=epiweek, y=value, color=variable), show.legend = F) +
      facet_grid(variable~., scales = "free_y") + theme_bw()


ggsave(file = paste0(homewd, "/figures/climate-data-plot.png"),
       plot = p3,
       units="mm",  
       width=60, 
       height=90, 
       scale=3, 
       dpi=200)
