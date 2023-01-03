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

#load data - this is the data that has been time-shifted from the lags in Figure 2
merge.shift <- read.csv(file = paste0(homewd, "/data/lagged-dat.csv"), header = T, stringsAsFactors = F)
head(merge.shift)

#try some models to predict cases from the lagged climate variables,
#with the offset of number of hospitals

#first try no predictor by year

#lags + no random effects
m1 <- glm(cases~ mean_H2M + meantempLag + precip_lag + offset(n_hospitals), data = merge.shift, family = "poisson")
summary(m1)

#library(sjPlot)

#then try with year as a random effect
merge.shift$year <- as.factor(merge.shift$year)
library(lme4)

#lags + random effects
m2 <- glmer(cases ~ mean_H2M + meantempLag + precip_lag + offset(n_hospitals) + (1|year),data = merge.shift, family = "poisson")
summary(m2)

#compare against the unlagged climate data
# no lags + no random effects
m3 <- glm(cases~n_hospitals + mean_H2M + meanTemp + sum_precip + offset(n_hospitals ),data = merge.shift, family = "poisson")
summary(m3)


AIC(m1, m2, m3) #m2 is best


#below, tested GAM but you can ignore this. better to use glm
# library(mgcv)
# m4 <- gam(cases_by_hospital~#s(n_hospitals,  bs="re") + 
#                 s(mean_H2M, k=7, bs="tp") + 
#                 s(meantempLag, k=7, bs="tp") + 
#                 s(precip_lag, k=7, bs="tp") + 
#                 s(year, bs="re"),
#                 #family = "poisson",
#                 data = merge.shift)
# 
# 
# summary(m4)
# 
# plot.gam(m4)
# 

#save the output from best fit model 2

out.m2 <- as.data.frame(summary(m2)$coefficients)
head(out.m2)

#compute confidence intervals 
#(all effects are significant except)

out.m2$lci <- out.m2[,1] - 1.96*out.m2[,2]
out.m2$uci <- out.m2[,1] + 1.96*out.m2[,2]
out.m2$effect <- row.names(out.m2)

out.m2 <- subset(out.m2, effect!="(Intercept)")
#out.m2$effect[out.m2$effect=="n_hospitals"] <- "number\nhospitals surveyed"
out.m2$effect[out.m2$effect=="meantempLag"] <- "mean temperature,\nlagged 6 weeks"
out.m2$effect[out.m2$effect=="precip_lag"] <- "sum precipitation,\nlagged 3 weeks"
out.m2$effect[out.m2$effect=="mean_H2M"] <- "mean humidity,\nconcurrent"
names(out.m2)[names(out.m2)=="effect"] <- "predictor"
out.m2$predictor <- factor(out.m2$predictor, levels=c("sum precipitation,\nlagged 3 weeks","mean humidity,\nconcurrent", "mean temperature,\nlagged 6 weeks",  "number\nhospitals surveyed"))
out.m2$isSig <- "yes"
out.m2$isSig[out.m2$predictor=="number\nhospitals surveyed"] <- "no"



Fig3a <- ggplot(out.m2) + geom_point(aes(x=predictor, y=Estimate), size=3) +
        geom_errorbar(aes(x=predictor, ymin=lci, ymax=uci), linewidth=1) +  
        geom_label(aes(x=predictor, y=0.315, label="***"), label.size = 0, size=8) +
        geom_hline(aes(yintercept=0), linetype=2, color="cornflowerblue", linewidth=1) +
        theme_bw() +  coord_flip(ylim=c(-.03,.335)) +
        ylab("Slope") +
        theme(panel.grid = element_blank(),
              axis.title = element_text(size=16),
              plot.margin = unit(c(.5,.5,.1,.5), "lines"),
              axis.text = element_text(size=14))


# Now plot the random effects using sjPlot

#ideally, you would add package sjPlot
#plot model output
plot_model(m2, type="est") # fixed effects - you've already shown these above in Fig 3a
plot_model(m2, type="re") # random effects - this is what you want to use now

#save the random effects data 
random.dat = plot_model(m2, type="re")$data 

#look at data
head(random.dat)

#add a significance term
random.dat$sig <- "yes"

random.dat$sig[random.dat$group=="neg" & random.dat$conf.high>=1] <- "no"
random.dat$sig[random.dat$group=="pos" & random.dat$conf.low<=1] <- "no"

random.dat$group[random.dat$sig=="no"] <- "not_sig"

#and make your own custom plot
colz = c('pos' = "red", 'neg' = "navy", 'not_sig' = "gray60")

Fig3b <- ggplot(random.dat) + geom_point(aes(x=term, y=estimate, color=group), size=3, show.legend = F) +
  geom_errorbar(aes(x=term, ymin=conf.low, ymax=conf.high, color=group), linewidth=1, show.legend = F) +  
  #geom_label(aes(x=predictor, y=0.315, label="***"), label.size = 0, size=8) +
  geom_hline(aes(yintercept=1), linetype=2, color="cornflowerblue", linewidth=1) +
  theme_bw() + scale_color_manual(values=colz) + coord_flip() +
  ylab("Partial Effects") +
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size=16),
        axis.title.y =element_blank(),
        plot.margin = unit(c(.5,.5,.1,1.8), "lines"),
        axis.text = element_text(size=14))


# Finally, plot cases by the lagged climate predictors

# Include predicted cases from the glm... 
# choose predictions
predict.dat = cbind.data.frame(precip_lag = 0:200, mean_H2M = seq(70,90, length=201), meantempLag=seq(10,22, length=201))
predict.dat$n_hospitals = 1
predict.dat$cases = NA

predict.dat$cases = exp(predict(m1, newdata = predict.dat, type="response"))
predict.dat$cases_lci = exp(predict(m1, newdata = predict.dat, type="response")-1.96*predict(m1, newdata = predict.dat, type="response", se.fit = T)$se)
predict.dat$cases_uci = exp(predict(m1, newdata = predict.dat, type="response")+1.96*predict(m1, newdata = predict.dat, type="response", se.fit = T)$se)

# # excluding random effects (so predict from m1 instead of m2)
# merge.shift$predict_case <- exp(predict(m1, type="response"))
# merge.shift$predict_cases_by_hospital <- merge.shift$predict_case/merge.shift$n_hospitals


merge.melt <- dplyr::select(merge.shift, year, cases_by_hospital,  precip_lag,  meantempLag,mean_H2M )
merge.melt <- melt(merge.melt, id.vars = c("cases_by_hospital",  "year"))

#now rename predictors to match those in Fig3b
merge.melt$variable <- as.character(merge.melt$variable)
merge.melt$variable[merge.melt$variable=="meantempLag"] <- "mean temperature,\nlagged 6 weeks"
merge.melt$variable[merge.melt$variable=="precip_lag"] <- "sum precipitation,\nlagged 3 weeks"
merge.melt$variable[merge.melt$variable=="mean_H2M"] <- "mean humidity,\nconcurrent"

merge.melt$variable <- factor(merge.melt$variable, levels=c("sum precipitation,\nlagged 3 weeks","mean humidity,\nconcurrent", "mean temperature,\nlagged 6 weeks",  "number\nhospitals surveyed"))

predict.melt <-  melt(predict.dat, id.vars = c("cases", "cases_lci", "cases_uci", "n_hospitals"))

predict.melt$variable <- as.character(predict.melt$variable)
predict.melt$variable[predict.melt$variable=="meantempLag"] <- "mean temperature,\nlagged 6 weeks"
predict.melt$variable[predict.melt$variable=="precip_lag"] <- "sum precipitation,\nlagged 3 weeks"
predict.melt$variable[predict.melt$variable=="mean_H2M"] <- "mean humidity,\nconcurrent"

predict.melt$variable <- factor(predict.melt$variable, levels=c("sum precipitation,\nlagged 3 weeks","mean humidity,\nconcurrent", "mean temperature,\nlagged 6 weeks",  "number\nhospitals surveyed"))
predict.melt = subset(predict.melt, cases < (max(merge.melt$cases_by_hospital)+5) & cases_lci < (max(merge.melt$cases_by_hospital)+5) & cases_uci < (max(merge.melt$cases_by_hospital)+5))


merge.melt$year <- as.numeric(as.character(merge.melt$year))

Fig3c <- ggplot(merge.melt) + 
          geom_point(aes(x=value, y=cases_by_hospital, color=year)) +
          facet_grid(~variable, scales = "free", switch = "x") + theme_bw() + ylab("cases per hospital") +
          scale_color_continuous(high="darkblue", low="lightblue") +
          geom_ribbon(data=predict.melt, aes(x=value, ymin=cases_lci, ymax=cases_uci), alpha=.3) +
          geom_line(data=predict.melt, aes(x=value, y=cases), size=1) +
          #geom_point(aes(x=value, y=predict_cases_by_hospital)) +
          theme(panel.grid = element_blank(), axis.title.x = element_blank(),
                axis.title = element_text(size=16),
                legend.position = c(x=.74,y=.8),
                strip.placement = "outside",
                strip.text = element_text(size=14),
                plot.margin = unit(c(.5,.5,.5,.5), "lines"),
                axis.text = element_text(size=14),
                strip.background = element_blank())


#and plot together
Fig3top <- cowplot::plot_grid(Fig3a, Fig3b, ncol = 2, nrow = 1, labels = c("A", "B"), rel_widths = c(1.25, 1), label_size = 22, label_x = c(0,-.03))
Fig3 <- cowplot::plot_grid(Fig3top, Fig3c, ncol = 1, nrow = 2, labels = c("", "C"), label_size = 22, rel_heights = c(1,1.3))



ggsave(file = paste0(homewd, "/figures/Fig3.png"),
       plot = Fig3,
       units="mm",  
       width=100, 
       height=80, 
       scale=3, 
       dpi=300)
