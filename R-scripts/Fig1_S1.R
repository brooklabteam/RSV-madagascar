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


# Figure 1
# panel A: plot case prevalence, organized by month OR
# plot raw cases, divided by number of reporting hospitals 
# (2 for all years but 3 for 2022)

#and load data - first the age-structured
dat <- read.csv(file=paste0(homewd,"/data/rsv_data_update_2011_2022.csv"))

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
dat$year <- year(dat$sampling_date)
sort(unique(dat$year))
subset(dat, year==2022)

dat.ts <- ddply(dat,.(year), summarise, mean_age = mean(age, na.rm=T) )
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
calc.hosp <- ddply(dat, .(year, hospital), summarise, cases=sum(RSV), tested=length(RSV))
N.year <- ddply(dat, .(year), summarise, yearly_cases=sum(RSV), yearly_tested=length(RSV))

calc.hosp <- merge(calc.hosp, N.year, by="year", all.x = T)
calc.hosp$percent_cases = calc.hosp$cases/calc.hosp$yearly_cases
calc.hosp$percent_cases = calc.hosp$tested/calc.hosp$yearly_tested

#most cases drom CENHOSOA

#what percent of testing belonged to each hospital each year

percent.hosp <- ddply(calc.hosp,.(year), summarise, cases=sum(cases), tested=sum(tested))

percent.hosp$prevalence = percent.hosp$cases/percent.hosp$tested



percent.hosp <- ddply(calc.hosp, .(hospital), summarise, cases=sum(cases), tested=sum(tested))
#only one data point from 2010 with no hospital identifier

calc.hosp = subset(calc.hosp, year>2010 & year<2022)
#most of these cases are from CENHOSOA but a few are not - this is figure S1
FigS1 <- ggplot(data=calc.hosp) + 
          geom_line(aes(x=year, y=tested, color=hospital), size=1) +
          geom_line(aes(x=year, y=cases, color=hospital), linetype=2, size=1) + 
          scale_x_continuous(breaks=seq(2011,2021, by=2)) +
          theme_bw() + ylab("febrile patients") +
          theme(panel.grid = element_blank(),
          axis.title.y = element_text(size=18),
          axis.title.x = element_blank(),
          axis.text = element_text(size=14),
          plot.margin = unit(c(.1,.1,.5,1), "cm"),
          legend.position = c(.25,.9),
          legend.direction = "horizontal") +
          guides(color=guide_legend(ncol=2))


ggsave(file = paste0(homewd, "/figures/FigS1.png"),
       plot = FigS1,
       units="mm",  
       width=65, 
       height=50, 
       scale=3, 
       dpi=300)


#catchment by hospital
#: BHK(28574), CSMI TSL(43222), MJR(8000), CENHOSOA(89000), CHUMET(43222)
pop.hospital <- cbind.data.frame(hospital= c("BHK", "CSMI-TSL", "MJR", "CENHOSOA", "CHUMET"), pop=c(28572, 43222, 8000, 89000, 43222))

#and join
calc.hosp <- merge(calc.hosp, pop.hospital, by="hospital", all.x = T)

# what was population served catchment by year?
pop.year <- ddply(calc.hosp, .(year), summarise, tot_pop= sum(pop))
pop.year <- pop.year[complete.cases(pop.year),]
ggplot(pop.year) + geom_line(aes(year, tot_pop))
#save this for the tsir
write.csv(pop.year, file = paste0(homewd, "/data/catchment_pop_by_year.csv"), row.names = F)

#calculate prevalence by week
dat.sum$prevalence <- dat.sum$cases/dat.sum$tested

#calculate year of sampling
dat.sum$year <- year(dat.sum$epiweek)
sort(unique(dat.sum$year))
#limit years 
dat.sum <- subset(dat.sum, year>2010 & year<2022)

#look at your data
head(dat.sum)


# plot prevalence by week
Fig1A <- ggplot(data=dat.sum) + geom_point(aes(x=epiweek, y=prevalence, size=tested)) +
         geom_line(aes(x=epiweek, y=prevalence)) + theme_bw() + 
         scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("prevalence RSV") +
         theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
               axis.title.x = element_blank(), axis.text = element_text(size = 14), 
               legend.direction = "horizontal", legend.position = c(.7,.9), 
               legend.background  = element_rect(color="black"), 
               plot.margin = unit(c(.5,.1,.1,.1), "cm"))
print(Fig1A)

# 
# # Also plot raw cases - these don't account for increased sampling effort in some years
# Fig1Ab<- ggplot(data=dat.sum) + geom_point(aes(x=epiweek, y=cases, size=tested)) +
#   geom_line(aes(x=epiweek, y=cases)) + theme_bw() + 
#   scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("cases in catchment") +
#   theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
#         axis.title.x = element_blank(), axis.text = element_text(size = 14), 
#         legend.direction = "horizontal", legend.position = c(.15,.9), 
#         legend.background  = element_rect(color="black"),
#         plot.margin = unit(c(.1,.1,1.1,.9), "cm"))
# 
# print(Fig1Ab)


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
        legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))
print(Fig2Aa)

#no scaling - use this one
Fig2A <- ggplot(data=dat.sum) + geom_point(aes(x=week_of_year, y=cases, color=year)) +
  geom_line(aes(x=week_of_year, y=cases, color=year))+ ylab("cases in catchment") +
  xlab("month of year") + theme_bw() +
  scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.2,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))
print(Fig2A)

## Now build a GAM to test for predictors of positive cases, especially seasonality
## Start by including all predictors

# Look at data
head(dat) 
#limit years 
dat$year <- year(dat$sampling_date)
dat <- subset(dat, year>2010 & year<2022)
sort(unique(dat$year))
dat$RSV


# Clean data to remove your incorrect variables
# Sex can only be M or F
dat$sex[dat$sex=="H"] <- "M"
dat$sex[ dat$sex=="I" | dat$sex=="ND"] <- NA


# Make sex a factor
dat$sex <- as.factor(dat$sex)

# Calculate year of sampling and keep as continuous
dat$year <- year(dat$sampling_date)


# Clean age 
# Can we fill in the blank here for the missing data?
unique(dat$age)

# where you have DDN but not age, you can fill in age
dat$age[is.na(dat$age) & !is.na(dat$DDN)] <- (dat$sampling_date[is.na(dat$age) & !is.na(dat$DDN)] - dat$DDN[is.na(dat$age) & !is.na(dat$DDN)])/365
unique(dat$age)
subset(dat, age<0) 
length(dat$age[is.na(dat$age)]) #173 missing
length(dat$sex[is.na(dat$sex)]) #90 missing
length(dat$sex[is.na(dat$sex) & is.na(dat$age)]) #73 missing both

dat$year <- as.factor(dat$year)


head(dat)
dat$hospital <- as.factor(dat$hospital)
dat$year <- as.numeric(as.character(dat$year))
unique(dat$hospital)

dat$year <- as.factor(dat$year)

# Run GAM with all predictors - here year as a random effect
gam1 <- gam(RSV~s(doy, bs="cc") + 
                #s(year, bs="tp") +
                s(year, bs="re") +
                s(age, bs="tp") +
                s(hospital, bs="re") +
                s(sex, bs="re"),
                data=dat,
                family = "binomial")

summary(gam1)

#n=3245 because a few are missing sex or age.
#total N=3432


# doy, year, age and hospital have significant effects here
# sex is not significant


# Now, load prior script to calculate and plot random effects:
source(paste0(homewd, "/R-scripts/mollentze-streicker-2020-functions.R"))

#load plotting functions
plot.partial <- function(df, var, response_var){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  p2 <- ggplot(data=df2, aes(var,  Residual)) +
       geom_boxplot(aes(var~Residual))
  
  p1 <- ggplot(data=df1, aes(var, y)) + 
    geom_crossbar(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                  alpha=.4, show.legend = F) +
    #geom_point(aes(x=var, y=y, color=var), size=5) +
    #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
    scale_fill_manual(values = fillz) +
    geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
    theme(panel.grid = element_blank(), axis.title.x = element_blank(),
          axis.text = element_text(size=14, angle = 90),
          axis.title.y = element_text(size=18, angle = 90),
          legend.text = element_text(size=10),
          plot.margin = unit(c(.1,.1,.5,1), "cm"))+
    ylab(paste0("partial effect on ", response_var)) 
  
  #print(p1)
  
  return(p1)
}
plot.partial.cont <- function(df, log, var, response_var, alt_var, legend.on){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  if(legend.on==TRUE){
    if(log==F){
      
      p1 <- ggplot(data=df1, aes(var, y)) + 
        geom_line(aes(color=IsSignificant), size=3)+
        geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                    alpha=.4, show.legend = F) +
        #geom_point(aes(x=var, y=y, color=var), size=5) +
        #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
        scale_fill_manual(values = fillz) +
        scale_color_manual(values = fillz) +
        scale_x_continuous(labels=scales::comma) +
        geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=18),
              axis.text = element_text(size=14),
              plot.margin = unit(c(.5,.1,.1,1), "cm"),
              legend.position = c(.15,.15))+
        ylab(paste0("partial effect on ", response_var)) + xlab(alt_var)
      
      
      
    }else{
      df1$var <- 10^(df1$var)
      
      p1 <- ggplot(data=df1, aes(var, y)) + 
        geom_line(aes(color=IsSignificant), size=3)+
        geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                    alpha=.4, show.legend = F) +
        #geom_point(aes(x=var, y=y, color=var), size=5) +
        #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
        scale_fill_manual(values = fillz) +
        scale_color_manual(values = fillz) +
        scale_x_continuous(labels=scales::comma) +
        geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
        theme(panel.grid = element_blank(), #axis.title.x = element_blank(),
              #axis.text.x = element_text(size=8, angle = 45),
              plot.margin = unit(c(.5,.1,.1,1), "cm"),
              legend.position = c(.15,.15))+
        ylab(paste0("partial effect on ", response_var)) + xlab(alt_var)
    }}else{
      if(log==F){
        
        p1 <- ggplot(data=df1, aes(var, y)) + 
          geom_line(aes(color=IsSignificant), size=3, show.legend = F)+
          geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                      alpha=.4, show.legend = F) +
          #geom_point(aes(x=var, y=y, color=var), size=5) +
          #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
          scale_fill_manual(values = fillz) +
          scale_color_manual(values = fillz) +
          scale_x_continuous(labels=scales::comma) +
          geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
          theme(panel.grid = element_blank(),
                axis.title = element_text(size=18),
                axis.text = element_text(size=14),
                plot.margin = unit(c(.5,.1,.1,1), "cm"),
                legend.position = c(.15,.15))+
          ylab(paste0("partial effect on ", response_var)) + xlab(alt_var)
        
        
        
      }else{
        df1$var <- 10^(df1$var)
        
        p1 <- ggplot(data=df1, aes(var, y)) + 
          geom_line(aes(color=IsSignificant), size=3, show.legend = F)+
          geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                      alpha=.4, show.legend = F) +
          #geom_point(aes(x=var, y=y, color=var), size=5) +
          #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
          scale_fill_manual(values = fillz) +
          scale_color_manual(values = fillz) +
          scale_x_continuous(labels=scales::comma) +
          geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
          theme(panel.grid = element_blank(), #axis.title.x = element_blank(),
                #axis.text.x = element_text(size=8, angle = 45),
                plot.margin = unit(c(.5,.1,.1,1), "cm"), 
                legend.position = c(.15,.15))+
          ylab(paste0("partial effect on ", response_var)) + xlab(alt_var)
      }
    }
  
  print(p1)
  
  return(p1)
}
plot.partial.year <- function(df, log, var, response_var, alt_var, legend.on){
  df1 = df$effects
  df2= df$partialResiduals
  #head(df2)
  
  #head(df1)
  names(df1)[names(df1)==var] <- "var"
  names(df2)[names(df2)==var] <- "var"
  
  fillz = c("No"="gray70", "Yes" = "skyblue3")
  
  
  #p2 <- ggplot(data=df2, aes(var,  Residual)) +
  #     geom_boxplot(aes(var~Residual))
  if(legend.on==TRUE){
    if(log==F){
      
      p1 <- ggplot(data=df1, aes(var, y)) + 
        geom_line(aes(color=IsSignificant), size=3)+
        geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                    alpha=.4, show.legend = F) +
        #geom_point(aes(x=var, y=y, color=var), size=5) +
        #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
        scale_fill_manual(values = fillz) +
        scale_color_manual(values = fillz) +
        scale_x_continuous(breaks=seq(2011,2021, by=2)) +
        geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              axis.text = element_text(size=14),
              plot.margin = unit(c(.2,.1,.5,1), "cm"),
              legend.position = c(.85,.85))+
        ylab(paste0("partial effect on ", response_var)) + xlab(alt_var) 
        
      
      
      
    }else{
      df1$var <- 10^(df1$var)
      
      p1 <- ggplot(data=df1, aes(var, y)) + 
        geom_line(aes(color=IsSignificant), size=3)+
        geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                    alpha=.4, show.legend = F) +
        #geom_point(aes(x=var, y=y, color=var), size=5) +
        #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
        scale_fill_manual(values = fillz) +
        scale_color_manual(values = fillz) +
        scale_x_continuous(breaks=seq(2011,2021, by=2)) +
        geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
        theme(panel.grid = element_blank(),
              axis.title.y = element_text(size=18),
              axis.title.x = element_blank(),
              axis.text = element_text(size=14),
              plot.margin = unit(c(.1,.1,.5,1), "cm"),
              legend.position = c(.15,.15))+
        ylab(paste0("partial effect on ", response_var)) + xlab(alt_var) 
      
    }}else{
      if(log==F){
        
        p1 <- ggplot(data=df1, aes(var, y)) + 
          geom_line(aes(color=IsSignificant), size=3, show.legend = F)+
          geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                      alpha=.4, show.legend = F) +
          #geom_point(aes(x=var, y=y, color=var), size=5) +
          #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
          scale_fill_manual(values = fillz) +
          scale_color_manual(values = fillz) +
          scale_x_continuous(breaks=seq(2011,2021, by=2)) +
          geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
          theme(panel.grid = element_blank(),
                axis.title.y = element_text(size=18),
                axis.title.x = element_blank(),
                axis.text = element_text(size=14),
                plot.margin = unit(c(.1,.1,.5,1), "cm"),
                legend.position = c(.15,.15))+
          ylab(paste0("partial effect on ", response_var)) + xlab(alt_var) 
        
        
        
        
      }else{
        df1$var <- 10^(df1$var)
        
        p1 <- ggplot(data=df1, aes(var, y)) + 
          geom_line(aes(color=IsSignificant), size=3, show.legend = F)+
          geom_ribbon(aes(ymin=ylower, ymax=yupper, fill=IsSignificant), 
                      alpha=.4, show.legend = F) +
          #geom_point(aes(x=var, y=y, color=var), size=5) +
          #geom_jitter(data=df2, aes(x=var, y=Residual), width=.1, alpha=.2, size=.3)+
          scale_fill_manual(values = fillz) +
          scale_color_manual(values = fillz) +
          scale_x_continuous(breaks=seq(2011,2021, by=2)) +
          geom_hline(aes(yintercept=0), linetype=2) + theme_bw() +
          theme(panel.grid = element_blank(),
                axis.title.y = element_text(size=18),
                axis.title.x = element_blank(),
                axis.text = element_text(size=14),
                plot.margin = unit(c(.1,.1,.5,1), "cm"),
                legend.position = c(.15,.15))+
          ylab(paste0("partial effect on ", response_var)) + xlab(alt_var) 
        
      }
    }
  
  print(p1)
  
  return(p1)
}

dat$year <- as.numeric(as.character(dat$year)) 
#now, fit the same gam with year as a numeric thinplate
gam2 <- gam(RSV~s(doy, bs="cc") + 
              s(year, bs="tp") +
              #s(year, bs="re") +
              s(age, bs="tp") +
              s(hospital, bs="re") +
              s(sex, bs="re"),
            data=dat,
            family = "binomial")

summary(gam2)

AIC(gam1, gam2) #model is better with year as a random effect, so get all the partial effects from this 



# Use functions from above script to calculate random effects for each
age.df <- get_partial_effects_continuous(gamFit = gam1, var="age")
doy.df <- get_partial_effects_continuous(gamFit = gam1, var="doy")
#look at deviations by specific year:
year.df <- get_partial_effects(fit=gam1, var = "year")
sex.df <- get_partial_effects(fit=gam1, var = "sex")
hosp.df <- get_partial_effects(fit=gam1, var = "hospital")


# Sex
plot.partial(df=sex.df, var="sex", response_var = "RSV positivity") 
# Sex here is not significant
# save as panel E
Fig1E <- plot.partial(df=sex.df, var="sex", response_var = "RSV positivity") 

#Hospital
#remove the one with o 
plot.partial(df=hosp.df, var="hospital", response_var = "RSV positivity") 
#no sig terms by partial effect which is good
# save as panel F
Fig1F <- plot.partial(df=hosp.df, var="hospital", response_var = "RSV positivity") 

#Year
#save as panel C
Fig1C <- plot.partial(df=year.df, var="year", response_var = "RSV positivity") 

#Age
plot.partial.cont(df=age.df, log=F, var="age", response_var = "RSV positivity", alt_var ="age", legend.on = TRUE) 
# Significant trend with age shows that higher ages have decreased association 
# with RSV positivity (disease is a disease of kids)

#can see specifically where it becomes neg/pos in the data
age.df$effects 
#relationship is only positive in ages 0-21. Much less certainty - and, in fact, no directionality in ages >85

#save as panel D
Fig1D <- plot.partial.cont(df=age.df, log=F, var="age", response_var = "RSV positivity", alt_var ="age", legend.on = FALSE) 

#DOY
plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern
# Here is the pattern of positivity by day of year - 
# we see that positive cases are found between day 1 and ~90 and negative the rest of the
# year.
#save as panel B
Fig1B <- plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern


# a few years with significant deviations - consistent with linear trends

dat$year <- as.factor(dat$year)

# Now, because a bunch of data were missing age and sex, so
# go ahead and redo the GAM with only day of year and year as predictors!
gam3 <- gam(RSV~s(doy, bs="cc") +
                s(year, bs="re"),
                data=dat,
                family = "binomial")

summary(gam3) # N= 3435 all of the data here

doy.df <- get_partial_effects_continuous(gamFit = gam3, var="doy")
year.df <- get_partial_effects(fit=gam3, var = "year")

plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = F) #seasonal pattern
plot.partial(df=year.df, var="year", response_var = "RSV positivity") #replicates 2012 & 2013. 2019 also sig here



#But, compare models
AIC(gam1, gam2, gam3)# gam1 is A LOT better, so let's use it


Fig1AB <- cowplot::plot_grid(Fig1A, Fig1B, nrow = 1, ncol = 2, labels=c("A", "B"), label_size = 22, align = "hv")
Fig1CD <- cowplot::plot_grid(Fig1C, Fig1D, nrow = 1, ncol = 2, labels=c("C", "D"), label_size = 22, align = "hv")
Fig1EF <- cowplot::plot_grid(Fig1E, Fig1F, nrow = 1, ncol = 2, labels=c("E", "F"), label_size = 22, align = "hv")

Fig1 <- cowplot::plot_grid(Fig1AB,Fig1CD, Fig1EF, nrow = 3, ncol = 1, label_size = 22, rel_heights = c(1,1.1,1.1))



ggsave(file = paste0(homewd, "/figures/Fig1.png"),
       plot = Fig1,
       units="mm",  
       width=100, 
       height=100, 
       scale=3, 
       dpi=300)
