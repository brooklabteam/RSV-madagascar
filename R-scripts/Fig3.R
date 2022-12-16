rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(mgcv)
library(reshape2)
library(sjPlot)
library(lme4)


homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)

#load data 
merge.shift <- read.csv(file = paste0(homewd, "/data/lagged-dat.csv"), header = T, stringsAsFactors = F)
head(merge.shift)

#first try no predictor by year

#lags + no random effects
m1 <- glm(cases~n_hospitals + mean_H2M + meantempLag + precip_lag, data = merge.shift, family = "poisson")
summary(m1)

library(sjPlot)

#then try with year as a random effect
merge.shift$year <- as.factor(merge.shift$year)
library(lme4)

#lags + random effects
m2 <- glmer(cases~n_hospitals + mean_H2M + meantempLag + precip_lag + (1|year),data = merge.shift, family = "poisson")
summary(m2)


#compare against the unlagged data
# no lags + no random effects
m3 <- glm(cases~n_hospitals + mean_H2M + meanTemp + sum_precip,data = merge.shift, family = "poisson")
summary(m3)


AIC(m1, m2, m3) #m3 is better


#ideally, you would add package sjPlot
#plot model output
plot_model(m2, type="est") # fixed effects
plot_model(m2, type="re") # random effects

