rm(list=ls())

homewd= "/Users/carabrook/Developer/RSV-madagascar"


setwd(homewd)

library(readxl)
library(lubridate)


#read in data and plot it
dat_tsiry <- read.csv(paste0(homewd, "/data/tsiry_edit.csv"),header=T, stringsAsFactors = F)
head(dat_tsiry)
dat_tsiry$prelevement_date_prelevement<-as.Date(dat_tsiry$prelevement_date_prelevement,"%d/%M/%Y")
unique(dat_tsiry$prelevement_date_prelevement)
subset(dat_tsiry, is.na(prelevement_date_prelevement))

dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00091-13"] <- as.Date("2013-10-01")
dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00092-13"] <- as.Date("2013-10-01")
dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00093-13"] <- as.Date("2013-10-01")
dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00175-13"] <- as.Date("2013-10-01")
dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00091-13"] <- as.Date("2013-10-01")
dat_tsiry$prelevement_date_prelevement[dat_tsiry$Num..Viro=="00091-13"] <- as.Date("2013-10-01")




dat_tsiry$semaine_epid<-epiweek(dat_tsiry$prelevement_date_prelevement)
subset(dat_tsiry, is.na(semaine_epid))
dat_tsiry$annee_epid<-epiyear(dat_tsiry$prelevement_date_prelevement)
dat_tsirypos<-dat_tsiry%>%filter(RSV=="1")
head(dat_tsirypos)
dat_tsir_count<-dat_tsirypos%>%group_by(semaine_epid,annee_epid,RSV)%>%dplyr::summarise(count=n())
dat_tsir_count$annee_epid <- as.factor(dat_tsir_count$annee_epid)

subset(dat_tsir_count, is.na(annee_epid))
#and plot as a time series within a year to get at season
p1<-ggplot(data=dat_tsir_count) +
  geom_point(aes(x=semaine_epid, y=count, color=annee_epid)) +
  geom_line(aes(x=semaine_epid, y=count, color=annee_epid)) + theme_bw() +
  ylab("RSV cases") + xlab("week of year") +
  theme(plot.margin =unit(c(.2,.2,.2,1), "lines"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=16))
