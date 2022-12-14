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
dat <- read.csv(file=paste0(homewd,"/data/Tsiry_RSV_data_2011-2022_Update.csv"))

head(dat)
names(dat) <-c("num_viro", "DDN", "age", "sex", "RSV", "X", "sampling_date")
dat <- dplyr::select(dat, -(X))
dat$DDN <- as.Date(dat$DDN)
dat$sampling_date <- as.Date(dat$sampling_date)


#now plot cases by time and fit a gam describing seasonality

#calculate the day of year
dat$doy <- yday(dat$sampling_date)


#calculate the epidemic month (month + year)
dat$epimonth <- as.Date(as.character(cut.Date(dat$sampling_date, breaks="month")))

#sum cases by month
dat.sum <- ddply(dat, .(epimonth), summarise, cases = sum(RSV), tested=length(RSV)) 

#calculate prevalence by month
dat.sum$prevalence <- dat.sum$cases/dat.sum$tested

#calculate year of sampling
dat.sum$year <- year(dat.sum$epimonth)


#look at your data
head(dat.sum)


# plot prevalence by month
Fig1Aa<- ggplot(data=dat.sum) + geom_point(aes(x=epimonth, y=prevalence, size=tested)) +
         geom_line(aes(x=epimonth, y=prevalence)) + theme_bw() + 
         scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("prevalence RSV") +
         theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
               axis.title.x = element_blank(), axis.text = element_text(size = 14), 
               legend.direction = "horizontal", legend.position = c(.7,.9), 
               legend.background  = element_rect(color="black"), 
               plot.margin = unit(c(.1,.1,1.1,.1), "cm"))
print(Fig1Aa)


# Also calculate cases by reporting hospital
# There were 2 hospitals for most of the time series, 
# so this is just the number of cases divided by 2
dat.sum$cases_by_hospital <- dat.sum$cases/2

# Then, in 2022, you had 3 hospitals so you replace here
dat.sum$cases_by_hospital[dat.sum$year==2022] <- dat.sum$cases[dat.sum$year==2022]/3

# Now plot cases by hospital
Fig1Ab<- ggplot(data=dat.sum) + geom_point(aes(x=epimonth, y=cases_by_hospital, size=tested)) +
  geom_line(aes(x=epimonth, y=cases_by_hospital)) + theme_bw() + 
  scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("cases per reporting hospital") +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 14), 
        legend.direction = "horizontal", legend.position = c(.35,.9), 
        legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm"))

print(Fig1Ab)


# Now also try to plot the cases by week of year to visualize seasonality
# Add year, week, month
dat$year <- as.factor(year(dat$sampling_date))
dat$week <- week(dat$sampling_date)
dat$month <- month(dat$sampling_date)

# Now summarize by month or week of year
dat.wk <- ddply(dat, .(year, week), summarise, cases = sum(RSV), tested=length(RSV))
dat.wk$prevalence <- dat.wk$cases/dat.wk$tested
dat.wk$cases_by_hospital <- dat.wk$cases/2
dat.wk$cases_by_hospital[dat.wk$year==2022] <- dat.wk$cases[dat.wk$year==2022]/3

dat.month <- ddply(dat, .(year, month), summarise, cases = sum(RSV), tested=length(RSV))
dat.month$prevalence <- dat.month$cases/dat.month$tested
dat.month$cases_by_hospital <- dat.month$cases/2
dat.month$cases_by_hospital[dat.month$year==2022] <- dat.month$cases[dat.month$year==2022]/3

# by month (best)
Fig1Ba <- ggplot(data=dat.month) + geom_point(aes(x=month, y=cases_by_hospital, color=year), size=1.5) +
  geom_line(aes(x=month, y=cases_by_hospital, color=year)) + ylab("cases per reporting hospital") +
  xlab("month of year") + theme_bw() +
  scale_x_continuous(breaks=c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))


# by week
Fig1Bb <- ggplot(data=dat.wk) + geom_point(aes(x=week, y=cases_by_hospital, color=year)) +
  geom_line(aes(x=week, y=cases_by_hospital, color=year))


## Now build a GAM to test for predictors of positive cases, especially seasonality
## You should add data of reporting hospital when you find int because
## each hospital might have different seasonal patterns.
## Start by including all predictors

# Look at data
head(dat) 
dat$RSV

# Record the date as a date class
dat$sampling_date <- as.Date(dat$sampling_date)

# Clean data to remove your incorrect variables
# Sex can only be M or F
dat$sex[dat$sex=="H" | dat$sex=="I" | dat$sex=="ND"] <- NA

# Make sex a factor
dat$sex <- as.factor(dat$sex)

# Calculate year of sampling and make a factor
dat$year <- year(dat$sampling_date)
dat$year <- as.factor(dat$year)

# Clean age to get rid of the commas
# Can we fill in the blank here for the missing data?
dat$age <- gsub(x=dat$age, pattern= ",", replacement =  ".", fixed = T )
dat$age <- as.numeric(dat$age)

# Run GAM with all predictors
gam1 <- gam(RSV~s(doy, bs="tp") + 
                s(year, bs="re") +
                s(age, bs="tp") +
                s(sex, bs="re"),
                data=dat,
                family = "binomial")

summary(gam1)

#n=1245 because a bunch of daat are missing sex or age


# The summary shows a significant effect of  doy, year (as a random effect),
# age, and sex is minimally significant


# Now, load prior script to calculate and plot random effects:
source(paste0(homewd, "/R-scripts/mollentze-streicker-2020-functions.R"))

# Use functions from above script to calculate random effects for each
age.df <- get_partial_effects_continuous(gamFit = gam1, var="age")
doy.df <- get_partial_effects_continuous(gamFit = gam1, var="doy")
sex.df <- get_partial_effects(fit=gam1, var = "sex")
year.df <- get_partial_effects(fit=gam1, var = "year")

# And plot each partial effect
plot.partial(df=sex.df, var="sex", response_var = "RSV positivity") 
# Sex here is not significant

plot.partial(df=year.df, var="year", response_var = "RSV positivity")
# 2 years show significant deviations. We are missing later years because they are not accompanied by age data

plot.partial.cont(df=age.df, log=F, var="age", response_var = "RSV positivity", alt_var ="age", legend.on = TRUE) 
# Significant trend with age shows that higher ages have decreased association 
# with RSV positivity (disease is a disease of kids)

plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern
# Here is the pattern of positivity by day of year - 
# we see that positive cases are found between day 1 and ~90 and negative the rest of the
# year.
#


# Now, because a bunch of data were missing age and sex, so
# go ahead and redo the GAM with only day of year and year as predictors!
gam2 <- gam(RSV~s(doy, bs="tp") +
                s(year, bs="re"),
                data=dat,
                family = "binomial")

summary(gam2) 
# both are significant

# Calculate partial effects
doy.df <- get_partial_effects_continuous(gamFit = gam2, var="doy")
year.df <- get_partial_effects(fit=gam2, var = "year")

# Here is a better plot of seasonality within a single year:
plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern

#let's save this panel for the final plot
Fig1D <- plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern

# And here is a measure of the years that have significant impacts on RSV positivity:
plot.partial(df=year.df,  var="year", response_var = "RSV positivity") #a few high years, a few low years

#This looks great too! Are there climatic reasons for the years that are significantly 
# higher or lower????
Fig1C <- plot.partial(df=year.df,  var="year", response_var = "RSV positivity") 

# Now combine together to make Figure 1
# use cases per hospital version of panel A and month version of panel B

Fig1 <- cowplot::plot_grid(Fig1Ab, Fig1Ba, Fig1C, Fig1D, nrow = 2, ncol = 2, labels=c("A", "B", "C", "D"), label_size = 22)



ggsave(file = paste0(homewd, "/figures/Fig1.png"),
       plot = Fig1,
       units="mm",  
       width=90, 
       height=70, 
       scale=3, 
       dpi=300)
