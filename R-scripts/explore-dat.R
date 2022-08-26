rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#set working directory - Tsiry, you can write your own here
homewd = "/Users/carabrook/Developer/RSV-madagascar"

setwd(homewd)

#read in data and plot it
dat <- read.csv(file= paste0(homewd, "/data/RSV_2011_2022_Tsiry.csv"), header = T, stringsAsFactors = F)
head(dat)

#rename without the "X"
names(dat) <- c("week", 2011:2022)

#and make into long format
dat.long <- melt(dat, measure.vars = 2:ncol(dat), variable.name = "year")

head(dat.long)
names(dat.long)[names(dat.long)=="value"] <- "count"

dat.long$year <- as.factor(dat.long$year)

#and plot as a time series within a year to get at season
p1 <- ggplot(data=dat.long) + 
      geom_point(aes(x=week, y=count, color=year)) +
      geom_line(aes(x=week, y=count, color=year)) + theme_bw() +
      ylab("RSV cases") + xlab("week of year") +
      theme(plot.margin =unit(c(.2,.2,.2,1), "lines"),
            axis.text = element_text(size=12), 
            axis.title = element_text(size=16))


#and join together as one long time series
dat.long$date <- as.Date(paste(dat.long$week, dat.long$year, 'Mon'), '%U %Y %a')
head(dat.long)

subset(dat.long, is.na(date))

#and plot the long time series
p2 <- ggplot(data=dat.long) + geom_line(aes(x=date, y=count)) + theme_bw() +
      theme(axis.title.x = element_blank(), 
            axis.text = element_text(size=12), axis.title.y = element_text(size=16),
            plot.margin =unit(c(.2,5,.2,1), "lines")) + ylab("RSV cases")


#plot together and label
pAll <- cowplot::plot_grid(p1, p2, ncol=1, nrow=2, labels = c("A", "B"), label_size = 22)


ggsave(file = paste0(homewd, "/figures/RSV-time-series.png"),
       plot = pAll,
       units="mm",  
       width=90, 
       height=110, 
       scale=2.5, 
       dpi=200)
