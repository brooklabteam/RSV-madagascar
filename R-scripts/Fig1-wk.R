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
dat$epiweek <- as.Date(as.character(cut.Date(dat$sampling_date, breaks="week")))

#sum cases by week
dat.sum <- ddply(dat, .(epiweek), summarise, cases = sum(RSV), tested=length(RSV)) 

#calculate prevalence by biweek
dat.sum$prevalence <- dat.sum$cases/dat.sum$tested

#calculate year of sampling
dat.sum$year <- year(dat.sum$epiweek)

#limit years 
dat.sum <- subset(dat.sum, year>2010 & year<2022)

#look at your data
head(dat.sum)


# plot prevalence by week
Fig1Aa<- ggplot(data=dat.sum) + geom_point(aes(x=epiweek, y=prevalence, size=tested)) +
         geom_line(aes(x=epiweek, y=prevalence)) + theme_bw() + 
         scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("prevalence RSV") +
         theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
               axis.title.x = element_blank(), axis.text = element_text(size = 14), 
               legend.direction = "horizontal", legend.position = c(.7,.9), 
               legend.background  = element_rect(color="black"), 
               plot.margin = unit(c(.1,.1,1.1,.1), "cm"))
print(Fig1Aa)


# Also plot raw cases

# Now plot cases by hospital
Fig1Ab<- ggplot(data=dat.sum) + geom_point(aes(x=epiweek, y=cases, size=tested)) +
  geom_line(aes(x=epiweek, y=cases)) + theme_bw() + 
  scale_size_continuous(name="specimens\ntested", breaks = c(10,100,200)) + ylab("cases in catchment") +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 14), 
        legend.direction = "horizontal", legend.position = c(.15,.9), 
        legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm"))

print(Fig1Ab)


# Now also try to plot the cases by week of year to visualize seasonality
# Add year, week, month
head(dat)
dat$year <- year(dat$epiweek)
dat$week <- week(dat$epiweek)
dat$month <- month(dat$epiweek)

dat <- subset(dat, year>2010 & year<2022)

dat$year <- as.factor(dat$year)

# Now summarize by  week of year
dat.wk <- ddply(dat, .(year, week), summarise, cases = sum(RSV), tested=length(RSV))
dat.wk$prevalence <- dat.wk$cases/dat.wk$tested


# by week
Fig1B <- ggplot(data=dat.wk) + geom_point(aes(x=week, y=cases, color=year)) +
  geom_line(aes(x=week, y=cases, color=year))+ ylab("cases in catchment") +
  xlab("month of year") + theme_bw() +
  scale_x_continuous(breaks=seq(1,52,4.5), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13), 
        legend.title = element_blank(),
        legend.direction = "vertical", legend.position = c(.8,.7), 
        #legend.background  = element_rect(color="black"),
        plot.margin = unit(c(.1,.1,1.1,.9), "cm")) + 
  guides(color=guide_legend(ncol = 2))



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
dat$sex[dat$sex=="H"] <- "M"
dat$sex[ dat$sex=="I" | dat$sex=="ND"] <- NA


# Make sex a factor
dat$sex <- as.factor(dat$sex)

# Calculate year of sampling and make a factor
dat$year <- year(dat$sampling_date)
dat$year <- as.factor(dat$year)

# Clean age to get rid of the commas
# Can we fill in the blank here for the missing data?
dat$age <- gsub(x=dat$age, pattern= ",", replacement =  ".", fixed = T )
dat$age <- as.numeric(dat$age)

dat$DDN <- as.Date(dat$DDN)
# where you have DDN but not age, you can fill in age
dat$age[is.na(dat$age) & !is.na(dat$DDN)] <- (dat$sampling_date[is.na(dat$age) & !is.na(dat$DDN)] - dat$DDN[is.na(dat$age) & !is.na(dat$DDN)])/365

#dat$year <- as.numeric(as.character(dat$year))
dat$year <- as.factor(dat$year)

# Run GAM with all predictors
gam1 <- gam(RSV~s(doy, bs="cc") + 
                #s(year, bs="tp") +
                s(year, bs="re") +
                s(age, bs="tp") +
                s(sex, bs="re"),
                data=dat,
                family = "binomial")

summary(gam1)

#n=3246 because a few are missing sex or age

#then with year as a linear term
# Run GAM with all predictors
dat$year <- as.numeric(as.character(dat$year))
gam1b <- gam(RSV~s(doy, bs="cc") + 
               s(year, bs="tp") +
               #s(year, bs="re") +
               s(age, bs="tp") +
               s(sex, bs="re"),
             data=dat,
             family = "binomial")

summary(gam1b)





# The summary shows a significant effect of  doy, year (as a random effect),
# age, and sex is minimally significant


# Now, load prior script to calculate and plot random effects:
source(paste0(homewd, "/R-scripts/mollentze-streicker-2020-functions.R"))

# Use functions from above script to calculate random effects for each
age.df <- get_partial_effects_continuous(gamFit = gam1, var="age")
doy.df <- get_partial_effects_continuous(gamFit = gam1, var="doy")
sex.df <- get_partial_effects(fit=gam1, var = "sex")
#year.df <- get_partial_effects_continuous(gamFit=gam1, var = "year")
year.df <- get_partial_effects(fit=gam1, var = "year")

# And plot each partial effect
plot.partial(df=sex.df, var="sex", response_var = "RSV positivity") 
# Sex here is not significant

# save as panel E
Fig1E <- plot.partial(df=sex.df, var="sex", response_var = "RSV positivity") 

plot.partial(df=year.df, var="year", response_var = "RSV positivity")
#plot.partial.cont(df=year.df, log=F, var="year", response_var = "RSV positivity", alt_var ="year", legend.on = TRUE) 
# 3 years show significant deviations. 

#save as panel C
Fig1C <- plot.partial(df=year.df, var="year", response_var = "RSV positivity")

plot.partial.cont(df=age.df, log=F, var="age", response_var = "RSV positivity", alt_var ="age", legend.on = TRUE) 
# Significant trend with age shows that higher ages have decreased association 
# with RSV positivity (disease is a disease of kids)

#save as panel F
Fig1F <- plot.partial.cont(df=age.df, log=F, var="age", response_var = "RSV positivity", alt_var ="age", legend.on = TRUE) 

plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern
# Here is the pattern of positivity by day of year - 
# we see that positive cases are found between day 1 and ~90 and negative the rest of the
# year.
#

#save as panel D
Fig1D <- plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = F) #seasonal pattern


# Now, because a bunch of data were missing age and sex, so
# go ahead and redo the GAM with only day of year and year as predictors!
gam2 <- gam(RSV~s(doy, bs="cc") +
                s(year, bs="re"),
                data=dat,
                family = "binomial")

summary(gam2) 
# both are significant

#But, compare models
AIC(gam1, gam2, gam1b)# gam1 is A LOT better, so let's use it
# interesting no effect of covid really on cases in 2020/2021

#year when included as a numeric
year.df <- get_partial_effects_continuous(gamFit = gam1b, var="year")
extraPlot <- plot.partial.cont(df=year.df, log=F, var="year", response_var = "RSV positivity", alt_var ="year", legend.on = TRUE) 
#slight declining trend through time

# # Calculate partial effects
# doy.df <- get_partial_effects_continuous(gamFit = gam2, var="doy")
# year.df <- get_partial_effects(fit=gam2, var = "year")
# 
# # Here is a better plot of seasonality within a single year:
# plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern
# 
# #let's save this panel for the final plot
# Fig1D <- plot.partial.cont(df=doy.df, log=F, var="doy", response_var = "RSV positivity", alt_var ="day of year", legend.on = T) #seasonal pattern
# 
# # And here is a measure of the years that have significant impacts on RSV positivity:
# plot.partial(df=year.df,  var="year", response_var = "RSV positivity") #a few high years, a few low years
# 
# #This looks great too! Are there climatic reasons for the years that are significantly 
# # higher or lower????
# Fig1C <- plot.partial(df=year.df,  var="year", response_var = "RSV positivity") 

# Now combine together to make Figure 1
# use cases per hospital version of panel A and month version of panel B

Fig1 <- cowplot::plot_grid(Fig1Ab, Fig1B, Fig1C, Fig1D, Fig1E, Fig1F, nrow = 3, ncol = 2, labels=c("A", "B", "C", "D", "E", "F"), label_size = 22)



ggsave(file = paste0(homewd, "/figures/Fig1week.png"),
       plot = Fig1,
       units="mm",  
       width=90, 
       height=100, 
       scale=3, 
       dpi=300)
