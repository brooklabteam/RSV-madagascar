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
dat <- read.csv(file=paste0(homewd,"/data/rsv_data_update_2011_2022.csv"))

head(dat)


head(dat)
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date", "hospital")
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


#calculate the epidemic week
dat$epiweek <- as.Date(as.character(cut.Date(dat$sampling_date, breaks="week")))

#sum cases by week
dat.sum <- ddply(dat, .(epiweek), summarise, cases = sum(RSV), tested=length(RSV)) 

#calc how many cases tested at each hospital each year
dat$year <- year(dat$epiweek)

#calculate year of sampling
dat.sum$year <- year(dat.sum$epiweek)

#limit years 
dat.sum <- subset(dat.sum, year>2010 & year<2022)

#look at your data
head(dat.sum)

# Now also try to plot the cases by week of year to visualize seasonality
# These will instead go into Figure 2

dat.sum$week_of_year <- week(dat.sum$epiweek)
dat.sum$year <- as.factor(dat.sum$year)

# by week - scaling by test amount
Fig2Aa <- ggplot(data=dat.sum) + geom_point(aes(x=week_of_year, y=cases, color=year, size=tested)) +
  geom_line(aes(x=week_of_year, y=cases, color=year))+ ylab("cases in catchment") +
  xlab("month of year") + theme_bw() +
  scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.6), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))
print(Fig2Aa)

#no scaling - use this one
Fig2A <- ggplot(data=dat.sum) + #geom_point(aes(x=week_of_year, y=cases, color=year)) +
  geom_line(aes(x=week_of_year, y=cases, color=year), size=1)+ ylab("weekly cases") +
  xlab("month of year") + theme_bw() +
  scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.6), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))
print(Fig2A)

#now add a time column for merging later
dat.sum$time <- yday(dat.sum$epiweek)/365 + year(dat.sum$epiweek)

#Then, load beta for transmission by week of year (or biweek)
#load tsir data
tsir.dat <- read.csv(file = paste0(homewd, "/data/tsir_dat_beta.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)

beta.dat = subset(tsir.dat, year == 2011)

#Fig 2B - beta by biweek
tsir.dat$year <- as.factor(tsir.dat$year)
Fig2B <- ggplot(data=beta.dat) + 
  geom_ribbon(aes(x=week, ymin=beta_low, ymax=beta_high), alpha=.3) +
  geom_line(aes(x=week, y=beta), size=1) +
  ylab(bquote(beta~', transmission')) +
  xlab("week of year") + theme_bw() +
  scale_x_continuous(breaks=c(seq(1,52,4.5)/52), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))
print(Fig2B)


Fig2AB = cowplot::plot_grid(Fig2A, Fig2B, ncol=1, nrow=2, align = "hv", labels = c("A", "B"), label_size = 22)

 
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
head(merge.dat)
merge.dat$mean_H2M <- approx(x=clim.sum$time, y=clim.sum$mean_H2M, xout = merge.dat$time)$y
merge.dat$sum_precip <- approx(x=clim.sum$time, y=clim.sum$sum_precip, xout = merge.dat$time)$y
merge.dat$meanTemp <- approx(x=clim.sum$time, y=clim.sum$meanTemp, xout = merge.dat$time)$y
#and look at climate through time

head(merge.dat)

library(sjPlot)
merge.dat$year <- as.numeric(as.character(merge.dat$year))

#precip - increasing
gam1 <- gam(sum_precip~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam1) # slight pos trend
plot_model(gam1, type="pred", grid=T)

precip.df <- get_model_data(gam1, type="pred", grid=T)
precip.df$x[precip.df$group_col=="week"] <- precip.df$x[precip.df$group_col=="week"]*52
precip.df$group_col <- as.character(precip.df$group_col)
precip.df$group_col[precip.df$group_col=="week"] <- "week of year"
precip.df$group_col <- factor(precip.df$group_col, levels=c("year", "week of year"))

FigS3A <- ggplot(data=precip.df) + facet_grid(~group_col,scales = "free_x") +
          geom_line(aes(x=x, y=predicted), size=1) + ylab("predicted sum precipitation (mm)") +
          geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), alpha=.3) + theme_bw() +
          theme(panel.grid = element_blank(), strip.text = element_text(size=18),
          strip.background = element_rect(fill="white"),
          axis.title.y = element_text(size = 16),
          axis.title.x = element_blank(), axis.text = element_text(size = 13),
          plot.margin = unit(c(.1,.1,.1,.9), "cm")) 


#humidity - also increasing
gam2 <- gam(mean_H2M~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam2) # slight pos tred
plot_model(gam2, type="pred", grid=T)

humid.df <- get_model_data(gam2, type="pred", grid=T)
humid.df$x[humid.df$group_col=="week"] <- humid.df$x[humid.df$group_col=="week"]*52
humid.df$group_col <- as.character(humid.df$group_col)
humid.df$group_col[humid.df$group_col=="week"] <- "week of year"
humid.df$group_col <- factor(humid.df$group_col, levels=c("year", "week of year"))

FigS3B <- ggplot(data=humid.df) + facet_grid(~group_col,scales = "free_x") +
  geom_line(aes(x=x, y=predicted), size=1) + ylab("predicted mean humidity (%)") +
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), alpha=.3) + theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size=18),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(), axis.text = element_text(size = 13),
        plot.margin = unit(c(.1,.1,.1,.9), "cm")) 


#temp - increasing
gam3 <- gam(meanTemp~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam3) # pos trend
plot_model(gam3, type="pred", grid=T)


temp.df <- get_model_data(gam3, type="pred", grid=T)
temp.df$x[temp.df$group_col=="week"] <- temp.df$x[temp.df$group_col=="week"]*52
temp.df$group_col <- as.character(temp.df$group_col)
temp.df$group_col[temp.df$group_col=="week"] <- "week of year"
temp.df$group_col <- factor(temp.df$group_col, levels=c("year", "week of year"))

FigS3C <- ggplot(data=temp.df) + facet_grid(~group_col,scales = "free_x") +
  geom_line(aes(x=x, y=predicted), size=1) + 
  ylab(bquote('predicted mean temp ( '^0~'C)'))+
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), alpha=.3) + theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size=18),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(), axis.text = element_text(size = 13),
        plot.margin = unit(c(.1,.1,.1,.9), "cm")) 

#what about cases??? decreasing - 
#contradictory effects of increasing humidity (driving cases down) and increasing precip (driving cases up)
gam4 <- gam(cases~year +  s(week, k=7, bs="cc"), data=merge.dat)
summary(gam4) # neg trend
case.df <- get_model_data(gam4, type="pred", grid=T)

case.df$x[case.df$group_col=="week"] <- case.df$x[case.df$group_col=="week"]*52
case.df$group_col <- as.character(case.df$group_col)
case.df$group_col[case.df$group_col=="week"] <- "week of year"
case.df$group_col <- factor(case.df$group_col, levels=c("year", "week of year"))


FigS3D <- ggplot(data=case.df) + facet_grid(~group_col,scales = "free_x") +
  geom_line(aes(x=x, y=predicted), size=1) + 
  ylab("predicted cases")+
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), alpha=.3) + theme_bw() +
  theme(panel.grid = element_blank(), strip.text = element_text(size=18),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(), axis.text = element_text(size = 13),
        plot.margin = unit(c(.1,.1,.1,.9), "cm")) 


#and compile
FigS3 <- cowplot::plot_grid(FigS3A, FigS3B, FigS3C, FigS3D, ncol=1, nrow=4, labels=c("A", "B", "C", "D"), label_size = 22, align = "hv")


ggsave(file = paste0(homewd, "/figures/FigS3.png"),
       plot = FigS3,
       units="mm",  
       width=70, 
       height=115, 
       scale=3, 
       dpi=300)




#plot each by year
clim.dat <- dplyr::select(merge.dat, -(births), -(pop), -(week), -(beta_low), -(beta_high))
clim.melt <- melt(clim.dat, id.vars = c("time", "year", "week_num"))
head(clim.melt)
clim.melt$year <- as.factor(clim.melt$year)
clim.melt$label <- as.character(clim.melt$variable)
clim.melt$label[clim.melt$label=="sum_precip"]<-  "sum precipitation (mm)"
clim.melt$label[clim.melt$label=="meanTemp"] <- "mean temp (*C)"
clim.melt$label[clim.melt$label=="mean_H2M"] <- "mean humidity (H2M)"

Fig2C <-  ggplot(data=subset(clim.melt, label=="mean temp (*C)")) + 
                    geom_line(aes(x=week_num, y=value, color=year), size=1, show.legend = F) + 
                    theme_bw() + scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
                    ylab(bquote('mean temp ( '^0~'C)'))+
                    theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
                    axis.title.x = element_blank(), axis.text = element_text(size = 13), 
                    #legend.title = element_blank(),
                    #legend.direction = "vertical", legend.position = c(.8,.7), 
                    #legend.background  = element_rect(color="black"),
                    plot.margin = unit(c(.1,.1,1.1,.9), "cm")) 

Fig2D <-  ggplot(data=subset(clim.melt, label=="sum precipitation (mm)")) + 
  geom_line(aes(x=week_num, y=value, color=year), size=1, show.legend = F) + 
  theme_bw() + scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  ylab("sum precipitation (mm)")+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        #legend.title = element_blank(),
        #legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) 
  

Fig2E <-  ggplot(data=subset(clim.melt, label=="mean humidity (H2M)")) + 
  geom_line(aes(x=week_num, y=value, color=year), size=1, show.legend = F) + 
  theme_bw() + scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  ylab("mean relative humidity (%)")+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        #legend.title = element_blank(),
        #legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) 


Fig2 <- cowplot::plot_grid(Fig2A, Fig2B, Fig2C, Fig2D, Fig2E, ncol=1, nrow=5, align = "hv", labels = c("A", "B","C",  "D", "E"), label_size = 22)


ggsave(file = paste0(homewd, "/figures/Fig2.png"),
       plot = Fig2,
       units="mm",  
       width=90, 
       height=130, 
       scale=3, 
       dpi=300)


# Now, test whether lagged climate variables offer a better predictor of 
# beta


head(merge.dat)
merge.dat$log_beta <- log(merge.dat$beta)

#save this and plot it in Figure 3
write.csv(merge.dat, file = paste0(homewd, "/data/lagged_climate_transmission_RSV.csv"), row.names = F)

glm1 <- lm(log_beta~mean_H2M + sum_precip + meanTemp + year, data= merge.dat)
summary(glm1) #all three climate variables are significant
#tried year as a random effect and got nowhere.

#check best fit model
library(MuMIn)
glm1 <- lm(log_beta~mean_H2M + sum_precip + meanTemp + year, data= merge.dat, na.action = na.fail)
summary(glm1)
dredge(global.model = glm1) 
#best fit is humidity, temp, and precip (no year)

glm2<- lm(log_beta~mean_H2M + sum_precip+ meanTemp, data= merge.dat, na.action = na.fail)
summary(glm2)
AIC(glm1,glm2) #glm2 is better

library(sjPlot)

plot_model(glm2, type="est")

plot_model(glm2, type="pred", grid=T) #transmission decreases with increasing humidity and increases with increasing precip

# Figure 3 will be model comparison + some variation on above

#############
#Test whether climate lags are appropriate predictors


#beta lag - no lag really for the transmission rate

print(ccf(merge.dat$sum_precip, merge.dat$log_beta))
# 95% CI at 0.09


#save as data
dat.lag <- cbind.data.frame(lag = print(ccf(merge.dat$sum_precip, merge.dat$log_beta))$lag, acf=print(ccf(merge.dat$sum_precip, merge.dat$log_beta))$acf)
dat.lag$variable <- "sum_precip"

#what is the optimal lag?
dat.lag$lag[dat.lag$acf==max(dat.lag$acf)]
# 20 is maximized cross correlation: 
# so precip follows beta 20 epiweeks - not a sig lag

# and what about temp ?
dat2 = cbind.data.frame(lag = print(ccf(merge.dat$meanTemp, merge.dat$log_beta))$lag, acf=print(ccf(merge.dat$meanTemp, merge.dat$log_beta))$acf)
dat2$variable <- "meanTemp"
dat2$lag[dat2$acf==max(dat2$acf)]
#  -21 is maximized cross correlation
# cases follow temp by 21 weeks - pretty feeble signature
# 95% CI at 0.09


#and humidity
dat3 = cbind.data.frame(lag = print(ccf(merge.dat$mean_H2M, merge.dat$log_beta))$lag, acf=print(ccf(merge.dat$mean_H2M, merge.dat$log_beta))$acf)
dat3$variable <- "mean_H2M"
dat3$lag[dat3$acf==max(dat3$acf)]
# 24 only
# transmission precedes % humidity by 24 weeks - also no sig signature
# 95% CI at 0.09



#save together
dat.lag <- rbind(dat.lag, dat2, dat3)

write.csv(dat.lag, file=paste0(homewd, "/data/lag_output.csv"), row.names = F)


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


#We decided not to use these because they aren't really significant
LagPlot<- ggplot(dat.lag) + geom_label(data=max.lag, aes(x=18,y=.4, label=label), label.size = 0) +
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






