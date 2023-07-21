rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tsiR)
library(reshape2)


homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)

# 
# #first, birth data for madagascar
 birth.dat <- read.csv(file=paste0(homewd, "/data/WorldBankMada.csv"), header = T, stringsAsFactors = F)
 head(birth.dat)

#load pop.vector from the region
pop.vec <- read.csv(file=paste0(homewd, "/data/catchment_pop_by_year.csv"), header = T, stringsAsFactors = F)
#load pop.vector from the region
pop.vec <- pop.vec$tot_pop
names(pop.vec) <- 2011:2021


#select only the birth row from the national data
birth.vec <- birth.dat[1,5:ncol(birth.dat)]
names(birth.vec) <- seq(2009,2022,1) #assume this is census at the end of each year,

#here, select the range of dates for your data
birth.vec <- birth.vec[which(names(birth.vec)=="2011"):which(names(birth.vec)=="2021")]


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

# need to fill in the gaps in the data -
# TSIR breaks when there are 'extinctions' in the dataset
# So we make an assumption of 'low' cases only
# it does not change the results
sub.tsir.dat$cases[sub.tsir.dat$cases==0] <- 1

# This fills in the gaps in the time series
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

head(tsir.dat)


ggplot(tsir.dat) + geom_point(aes(x=time, y=cases)) + geom_line(aes(x=time, y=cases)) + ylab("cases")

tsir.dat$cases[is.na(tsir.dat$cases)]<- 1

tsir.dat$year <- as.factor(trunc(tsir.dat$time))
tsir.dat$week <- tsir.dat$time- as.numeric(as.character(tsir.dat$year))
head(tsir.dat)
ggplot(tsir.dat) + geom_point(aes(x=week, y=cases, color=year, group=year)) + 
  geom_line(aes(x=week, y=cases, color=year, group=year)) + ylab("weekly cases") #consider replotting figure 1b


ggplot(tsir.dat) + geom_point(aes(x=week, y=cases, color=year, group=year)) + 
  geom_line(aes(x=week, y=cases, color=year, group=year)) + ylab("weekly cases") + facet_wrap(~year)
  

# Now fit tsir model to data usig the tSIR package
 
fittedpars <- estpars(data=tsir.dat,
                            IP=1, alpha=.97, sbar=NULL, xreg = "cumcases",
                            regtype='lm',family='poisson',link='log')



beta.df <- cbind.data.frame(biweek = rep(1:26, each=2), week=1:52, beta = fittedpars$contact$beta, beta_low = fittedpars$contact$betalow, beta_high = fittedpars$contact$betahigh)


ggplot(data=beta.df) + geom_point(aes(x=biweek, y=beta)) +geom_line(aes(x=biweek, y=beta)) + geom_linerange(aes(x=biweek, ymin=beta_low, ymax=beta_high))
#ggplot(data=beta.df) + geom_point(aes(x=week, y=beta)) +geom_line(aes(x=week, y=beta)) + geom_linerange(aes(x=week, ymin=beta_low, ymax=beta_high))
#plot(fittedpars$contact$beta, type="b")

simfitted <- simulatetsir(data=tsir.dat,
                           IP = 1,
                           #epidemics = "break",
                           parms=fittedpars,
                           nsim=1000)

#and have the TSIR output diagnostics as Fig S1
#edit plot res function to have labels

plotres <- function (dat, max.plot = 10) {
  if (is.null(dat$SIRS) == TRUE) {
    dat$SIRS <- FALSE
  }
  if (dat$SIRS == TRUE) {
    ll.melt <- dat$ll.melt
    p1 <- ggplot(ll.melt, aes_string(x = "X1", y = "X2", 
                                     z = "loglik")) + geom_tile(aes_string(fill = "loglik")) + 
      scale_fill_gradient(low = "white", high = "red") + 
      theme_bw() + geom_contour(col = "black") + geom_point(aes(x = dat$k, 
                                                                y = dat$m), col = "black") + xlab("strength of immunity") + 
      ylab("duration of immunity") + ggtitle(sprintf("duration = %g, strength = %g", 
                                                     round(dat$m), signif(dat$k, 2)))
    update.ll <- ll.melt[, which(!apply(ll.melt == 0, 2, 
                                        all))]
    if (nrow(update.ll) * ncol(update.ll) != nrow(ll.melt) * 
        ncol(ll.melt)) {
      lldf <- NULL
      lldf$loglik <- update.ll$loglik
      lldf$Var <- update.ll[names(update.ll)[which(names(update.ll) != 
                                                     "loglik")]]
      lldf <- as.data.frame(lldf)
      names(lldf) <- c("loglik", "Var")
      p1 <- ggplot(lldf, aes_string("Var", "loglik")) + 
        geom_line(size = 2) + theme_bw() + ggtitle(sprintf("duration = %g, strength = %g", 
                                                           round(dat$m), signif(dat$k, 2)))
    }
    Sdf <- NULL
    Sdf$time <- head(dat$res$time, -1)
    Sdf$S <- dat$S
    Sdf <- as.data.frame(Sdf)
    p2 <- ggplot(Sdf, aes_string("time", "S")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(rho) == .(signif(mean(1/dat$rho), 
                                                       2))))
    p3 <- ggplot(data = dat$contact, aes_string("time", "beta")) + 
      geom_line(size = 2) + geom_ribbon(aes_(ymin = ~betalow, 
                                             ymax = ~betahigh), alpha = 0.5, col = "dodgerblue", 
                                        fill = "dodgerblue") + ylim(c(min(dat$contact$betalow), 
                                                                      max(dat$contact$betahigh))) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                                                                                                  .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == .(signif(dat$alpha, 
                                                                                                                                                                                         3)))) + ylab(bquote(beta))
    if (sum(sum(is.na(dat$contact))) > 0) {
      p3 <- ggplot(betadf, aes_string("time", "beta")) + 
        geom_line(size = 2) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                            .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                            .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
    }
    p4 <- logcorr(dat) + geom_abline(slope = 1, colour = "dodgerblue")
    drops <- c("mean", "sd", "error", "cases", "time")
    sim.only <- dat$res[, !(names(dat$res) %in% drops)]
    n <- ncol(sim.only)
    error <- qt(0.975, df = n - 1) * dat$res$sd/sqrt(n)
    dat$res$error <- error
    eb <- aes(ymax = mean + error, ymin = mean - error)
    p6 <- ggplot(data = dat$res, aes_string("time")) + theme(legend.position = "none") + 
      geom_line(aes_string(y = "cases"), colour = "dodgerblue", 
                size = 1) + xlab("year") + ylab("cases") + geom_line(aes_string(y = "mean"), 
                                                                     colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                                                                    alpha = 0.3) + theme_bw()
    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases
    p7 <- ggplot(data = inversecases, aes_string("time")) + 
      theme(legend.position = "none") + geom_line(aes_string(y = "cases"), 
                                                  colour = "dodgerblue", size = 1) + xlab("time") + 
      ylab("cases") + geom_line(aes_string(y = "mean"), 
                                colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                               alpha = 0.3) + theme_bw()
    if (dat$nsim > max.plot) {
      sampledat <- sample(sim.only, max.plot)
      sampledat$time <- dat$res$time
    }
    else {
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }
    meltdf <- melt(sampledat, id = "time")
    p8 <- ggplot(meltdf, aes_string(x = "time", y = "value", 
                                    fill = "variable")) + geom_line(alpha = 0.6, colour = "orangered4") + 
      xlab("time") + ylab("cases") + geom_line(data = dat$res, 
                                               aes_string(x = "time", y = "cases", fill = NA), colour = "dodgerblue", 
                                               size = 1) + theme_bw()
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 2)))
    print(p1, vp = vplayout(1, 1))
    print(p2, vp = vplayout(1, 2))
    print(p3, vp = vplayout(2, 1))
    print(p4, vp = vplayout(2, 2))
    print(p8, vp = vplayout(3, 1:2))
    print(p7, vp = vplayout(4, 1:2))
  }
  else {
    regdf <- NULL
    regdf$X <- dat$X
    regdf$Y <- dat$Y
    regdf$Yhat <- dat$Yhat
    regdf$time <- dat$res$time
    regdf <- as.data.frame(regdf)
    meltdf <- melt(regdf, id.vars = "time")
    p1 <- ggplot(meltdf, aes_string("time", "value", col = "variable")) + 
      geom_line(size = 2) + theme_bw()
    rhodf <- NULL
    rhodf$time <- dat$res$time
    rhodf$rho <- 1/dat$rho
    rhodf <- as.data.frame(rhodf)
    p2 <- ggplot(rhodf, aes_string("time", "rho")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(rho) == .(signif(mean(1/dat$rho), 
                                                       2)))) + ylab(bquote(1/rho))
    resdf <- NULL
    resdf$time <- dat$res$time
    resdf$Z <- dat$Z
    resdf$S <- dat$Z + dat$sbar
    resdf <- as.data.frame(resdf)
    meltres <- melt(resdf, id.vars = "time")
    p3 <- ggplot(meltres, aes_string("time", "value", col = "variable")) + 
      geom_line(size = 2) + theme_bw()
    loglikdf <- NULL
    loglikdf$sbar <- dat$Smean
    loglikdf$loglik <- dat$loglik
    loglikdf <- as.data.frame(loglikdf)
    p9 <- ggplot(loglikdf, aes_string("sbar", "loglik")) + 
      geom_line() + geom_point() + theme_bw() + geom_vline(xintercept = dat$sbar, 
                                                           linetype = "longdash") + ggtitle(bquote(bar(S) == 
                                                                                                     .(signif(dat$sbar, 2)) ~ "," ~ .(signif(dat$sbar/mean(dat$pop) * 
                                                                                                                                               100, 2)) ~ "%")) + xlab(bquote(bar(S)))
    betadf <- NULL
    betadf <- dat$contact
    betadf <- as.data.frame(betadf)
    p4 <- ggplot(betadf, aes_string("time", "beta")) + geom_line(size = 2) + 
      theme_bw() + ggtitle(bquote(bar(beta) == .(signif(mean(dat$beta), 
                                                        2)) ~ "," ~ alpha == .(signif(dat$alpha, 3)))) + 
      ylab(bquote(beta))
    if ("contact" %in% names(dat)) {
      p4 <- ggplot(data = dat$contact, aes_string("time", 
                                                  "beta")) + geom_line(size = 2) + geom_ribbon(aes_(ymin = ~betalow, 
                                                                                                    ymax = ~betahigh), alpha = 0.5, col = "dodgerblue", 
                                                                                               fill = "dodgerblue") + ylim(c(min(dat$contact$betalow), 
                                                                                                                             max(dat$contact$betahigh))) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                                                                                                                                                         .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                                                                                                                                                         .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
      if (sum(sum(is.na(dat$contact))) > 0) {
        p4 <- ggplot(betadf, aes_string("time", "beta")) + 
          geom_line(size = 2) + theme_bw() + ggtitle(bquote(bar(beta) == 
                                                              .(signif(mean(dat$beta), 2)) ~ "," ~ alpha == 
                                                              .(signif(dat$alpha, 3)))) + ylab(bquote(beta))
      }
    }
    p4 <- p4 + xlab(sprintf("time mod %g", nrow(dat$contact)))
    if (dat$inits.fit == TRUE) {
      inits.grid <- dat$inits.grid
      p5 <- ggplot(inits.grid, aes_string(x = "S0", y = "I0", 
                                          z = "log10LS")) + geom_tile(aes_string(fill = "log10LS")) + 
        scale_fill_gradient(low = "white", high = "black") + 
        theme_bw() + geom_contour(col = "black") + geom_point(aes(x = dat$inits[1]/mean(dat$pop), 
                                                                  y = dat$inits[2]/mean(dat$pop)), col = "red") + 
        xlab("prop. init. sus.") + ylab("prop. init. inf.")
    }
    drops <- c("mean", "sd", "error", "cases", "time")
    sim.only <- dat$res[, !(names(dat$res) %in% drops)]
    n <- ncol(sim.only)
    error <- qt(0.975, df = n - 1) * dat$res$sd/sqrt(n)
    dat$res$error <- error
    eb <- aes(ymax = mean + error, ymin = mean - error)
    p6 <- ggplot(data = dat$res, aes_string("time")) + theme(legend.position = "none") + 
      geom_line(aes_string(y = "cases"), colour = "dodgerblue", 
                size = 1) + xlab("year") + ylab("cases") + geom_line(aes_string(y = "mean"), 
                                                                     colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                                                                    alpha = 0.3) + theme_bw()
    inversecases <- dat$res
    inversecases$cases <- -dat$res$cases
    p7 <- ggplot(data = inversecases, aes_string("time")) + 
      theme(legend.position = "none") + geom_line(aes_string(y = "cases"), 
                                                  colour = "dodgerblue", size = 1) + xlab("time") + 
      ylab("cases") + geom_line(aes_string(y = "mean"), 
                                colour = "orangered4", size = 1) + geom_ribbon(eb, 
                                                                               alpha = 0.3) + theme_bw()
    if (dat$nsim > max.plot) {
      sampledat <- sample(sim.only, max.plot)
      sampledat$time <- dat$res$time
    }
    else {
      sampledat <- sim.only
      sampledat$time <- dat$res$time
    }
    meltdf <- melt(sampledat, id = "time")
    p8 <- ggplot(meltdf, aes_string(x = "time", y = "value", 
                                    fill = "variable")) + geom_line(alpha = 0.6, colour = "orangered4") + 
      xlab("time") + ylab("cases") + geom_line(data = dat$res, 
                                               aes_string(x = "time", y = "cases", fill = NA), colour = "dodgerblue", 
                                               size = 1) + theme_bw()
    #grid.newpage()
    #pushViewport(viewport(layout = grid.layout(5, 2)))
    if (dat$inits.fit == TRUE) {
      if (all(is.na(dat$loglik)) == T) {
        print(p1, vp = vplayout(1, 1))
        print(p2, vp = vplayout(1, 2))
        print(p3, vp = vplayout(2, 1:2))
        print(p4, vp = vplayout(3, 1))
        print(p5, vp = vplayout(3, 2))
        print(p8, vp = vplayout(4, 1:2))
        print(p7, vp = vplayout(5, 1:2))
      }
      else {
        print(p1, vp = vplayout(1, 1))
        print(p2, vp = vplayout(1, 2))
        print(p3, vp = vplayout(2, 1))
        print(p9, vp = vplayout(2, 2))
        print(p4, vp = vplayout(3, 1))
        print(p5, vp = vplayout(3, 2))
        print(p8, vp = vplayout(4, 1:2))
        print(p7, vp = vplayout(5, 1:2))
      }
    }
    else {
      if (all(is.na(dat$loglik)) == T) {
        print(p1, vp = vplayout(1, 1))
        print(p2, vp = vplayout(1, 2))
        print(p3, vp = vplayout(2, 1:2))
        print(p4, vp = vplayout(3, 1:2))
        print(p8, vp = vplayout(4, 1:2))
        print(p7, vp = vplayout(5, 1:2))
      }
      else {
        
        return(list(p1,p2,p3,p9,p4,p8,p7))
        # print(p1, vp = vplayout(1, 1))
        # print(p2, vp = vplayout(1, 2))
        # print(p3, vp = vplayout(2, 1))
        # print(p9, vp = vplayout(2, 2))
        # print(p4, vp = vplayout(3, 1:2))
        # print(p8, vp = vplayout(4, 1:2))
        # print(p7, vp = vplayout(5, 1:2))
      }
    }
  }
}
out.plot <- plotres(dat=simfitted)

pAB <- cowplot::plot_grid(out.plot[[1]], out.plot[[2]], nrow=1, ncol=2, align = "vh", labels = c("A", "B"))
pCD <- cowplot::plot_grid(out.plot[[3]], out.plot[[4]],  nrow=1, ncol=2, align = "vh", labels = c("C", "D"))
pEFG <- cowplot::plot_grid(out.plot[[5]], out.plot[[6]], out.plot[[7]],  nrow=3, ncol=1, align = "vh", labels = c("E", "F", "G"))


#and all together
pS2top <- cowplot::plot_grid(pAB, pCD, nrow = 2, ncol=1,align = "vh")
FigS2 <- cowplot::plot_grid(pS2top, pEFG, nrow=2, ncol=1, align = "vh", rel_heights = c(1,1.3))

ggsave(file = paste0(homewd, "/figures/FigS2.png"),
       plot = FigS2,
       units="mm",  
       width=100, 
       height=110, 
       scale=3, 
       dpi=300)


#now save the beta output to plot with the climate variables
names(beta.df)[names(beta.df)=="week"] <- "week_num"
beta.df <- dplyr::select(beta.df, -(biweek))

tsir.dat$week_num <- rep(1:52, length(unique(tsir.dat$year)))

tsir.dat <- merge(tsir.dat, beta.df, by="week_num", all.x = T)
head(tsir.dat)

#now take to climate and look for predictors with lags - save here to load later with figure 2
write.csv(tsir.dat, file = paste0(homewd, "/data/tsir_dat_beta.csv"), row.names = F)

