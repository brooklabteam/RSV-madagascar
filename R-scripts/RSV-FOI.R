rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(matrixcalc)

homewd="/Users/carabrook/Developer/RSV-madagascar"
setwd(homewd)

#load functions
model.age.incidence.series <- function(par.dat, age_vect, year.start){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and 35 ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  #N_sero = rep(N.sero.guess, length=lts)
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  #pexposed= as.list(rep(NA, (length(age_vect)*lts))) 
  #pexposed= rep(NA, (max(age_vect)))
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  #pnaive= as.list(rep(NA, (length(age_vect)*lts))) 
  #pnaive= rep(NA, (max(age_vect)))
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker[1] <- year.start
  year_tracker = rep(list(year_tracker),lts)
  
  year.start= min(par.dat$year)
  
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ #for-loops over all possible ages in our data range across all the years
      
      
      
      #remake the age vector for each year
      #print(a)
      age = age_vect_tmp[a]
      age_trunc = floor(age)
      age_current = age-age_trunc
      age_current[age_current==0] <- 1 #if you are at the end of the year, you spent the whole year here
      #age_class = ceiling(age)
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero#[2:nrow(par.dat)]
      #lambda = rep(exp(as.numeric(log.lambda.guess)), length=(lts))
      lambda =  year.par$lambda#[2:nrow(par.dat)] #one per year
      
      #make a vector of durations for the lambdas across the time series
      if (age>1){
        dur = rep(1, length(lambda))  
      }else if(age<=1){
        dur = rep(0, length(lambda))  
      }
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>2){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      dur[length(dur)] <- age_current
      
      
      # #age should never exceed i under our new regime
      # dur1 <- age - i
      # if(dur1 <=0){
      #   dur[1] <- 0
      # }else{
      #   dur[1] <- dur[1] + dur1
      # }
      #dur1[dur1<0] <- 0
      
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      
      #if you are < 1 year in age, we still want it to be possible for you
      #to experience a multitypic infection (we know that cases are reported <1 year in
      #our natl data, and we make the assumption that all reported data signifies secondary
      #infections)
      
      #therefore, if you are <1 year, we have to account for this by breaking down the
      #foi across each quarter
      
      #if (age<1){
      #this just means, we split the current year into time before this timestep and now
      lambda <- c(lambda, lambda[i])
      dur <- c(dur, dur[i])
      
      #and we account for the current timestep accordingly
      dur[length(dur)-1] <-   dur[length(dur)]-.25
      dur[length(dur)] <- .25 #always just a quarter year in here
      
      if(sum(dur) != age){
        warning("mismatch in age and integration")
      }
      #}else if (age>=1){
      #still get the bonus of more time in this timstep
      # dur[i] <-  dur[i]
      #}
      
      #
      
      #now, we integrate the hazard of becoming infected over all the years in each
      
      
      inte_pre = sum(dur[1:(length(dur)-1)]*lambda[1:(length(lambda)-1)])
      inte_now = sum(dur[length(dur)]*lambda[length(dur)])  
      
      #now sum them
      inte_all <-  inte_pre + inte_now
      
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_all))
      
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- 1-(1-exp(-inte_all))
      
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive = c(unlist(pnaive)))
  
  #any_primary_inf=c(unlist(pprim)), 
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(p.out[,3:4])
  
  
  # #head(p.out)
  # p.add <- cbind.data.frame(year = seq((year.start), (year.start+lts-1), 1),
  #                           age=rep(0, length=lts), 
  #                           exposed=rep(0, length=lts), naive=rep(1, length=lts), all_prim=rep(0, length=lts), 
  #                           multi=rep(0, length=lts), sum_exp_naive=rep(1, length=lts), sum_naive_prim_multi=rep(1, length=lts))
  # p.out <- rbind(p.add, p.out)
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=exposed))  + 
  # geom_line(aes(x=age, y=exposed)) +facet_wrap(~year)
  # 
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$exposed
  #ggplot(data=p.out) + geom_point(aes(x=age, y=cum_prop_cases))  + geom_line(aes(x=age, y=cum_prop_cases)) +facet_wrap(~year)
  p.out <- p.out[complete.cases(p.out),]
  #and split by age category
  p.out$age_trunc = ceiling(p.out$age)
  
  #p.split <- dlply(p.out, .(year, age_trunc))
  
  #get the mean value in each
  p.sum <- ddply(p.out, .(year, age_trunc), summarise, exposed = mean(exposed), naive=mean(naive), sum_exp_naive= mean(sum_exp_naive), cum_prop_cases = mean(cum_prop_cases) )
  
  names(p.sum)[names(p.sum)=="age_trunc"] <- "age"
  
  p.add <- cbind.data.frame(year = unique(p.sum$year), age = rep(0, length(unique(p.sum$year))), exposed=0, naive=1, sum_exp_naive=1, cum_prop_cases =0)
  
  p.sum <- rbind(p.add, p.sum)
  p.sum <- arrange(p.sum, year, age)
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=cum_prop_cases))  + 
    geom_line(aes(x=age, y=cum_prop_cases)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.fit.all <- function(p, par.dat, dat, year.start){
  
  par.dat$lambda <- exp(p)
  
  age_vect=seq(0, max(dat$age), by=1/4)
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, year.start,
                                        age_vect=age_vect)  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  #now, select only the years for fitting for which there are data
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  #plot(out.mod$cum_prop_cases, type="b")
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  dat.merge <- dplyr::select(dat, age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  # ggplot(data=out.merge) + geom_line(aes(x=age,y= cum_prop_cases), color="tomato") + facet_wrap(~year) + 
  # geom_point(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age, y=cum_prop_cases_data))
  # 
  # # # # 
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],p=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
    
  }
  
  if (ll==-Inf){
    ll= -1000000
  }
  
  #par.dat$lambda_high <- 0
  #par.dat$lambda_high[par.dat$lambda>1] <- 1
  #tot_over = sum(par.dat$lambda_high)
  #if(tot_over>0){
    #add a penalty for any lambda value over 1
    #   print("correcting for high lambda")
    #   print("original:")
    #   print(-ll)
   # ll=ll-(tot_over*1000000)
    #   
    #   print("corrected:")
    #   print(-ll)
  #}#
  
  return(-ll)
}
sum.yr.all <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(df$age),year=rep(unique(df$year), max(df$age)))
  
  df.out <- merge(df.out, df.sum, by=c("age", "year"), all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
fit.all.yrs.seq.yr.nlm <- function(dat,  lambda.guess, N.sero.fix,  fit.CI){
  
  
  #get dist back here for the oldest individual in the first year of the dataset
  dist.back =   max(dat$age[dat$year==min(dat$year)]) 
  #dist.back =  min(dat$year) - min(dat$year_of_first_FOI)
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  #min.year <- min(as.numeric(as.character(dat$year)))
  #dist.back <- max(dat$age[dat$year==min(dat$year)])#22
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  #age_vect_year = floor(age_vect)[!duplicated(floor(age_vect))]
  
  year.dat.sum <- lapply(year.dat, sum.yr.all)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  
  # # #head(df.out)
  #       ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #          geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # # # # # # # # 
  # # #make your guess parameters
  #lambda is takes data from the previous year and creates infections in this year
  
  #now, make the appropriate number of serotypes and lambdas
  
  if( length(lambda.guess)==1){ #here, number of serotypes is fixed across the time series
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = rep(lambda.guess, length((min(dat$year)-dist.back +1):max(dat$year))))
    
    
  }else if (length(lambda.guess)>1){ # here you have pre-prepped only the lambda but are going to fix serotypes across the dataset
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = lambda.guess)
    
  }
  
  #now, for all lambda prior to 1970, assume a very, very, very low FOI
  #par.dat$lambda[par.dat$year<1970] <- 0.0000001
  #and fit it cumulatively across the entire time series
  
  
  #if this is the first year, you need
  log.lambda.guess <- log(par.dat$lambda)
  
  out.NS <- nlm(f=log.lik.fit.all, 
                p = log.lambda.guess, 
                par.dat=par.dat, 
                year.start = min(par.dat$year),#this is the year the model will start iterating with the first birth cohort
                #age_vect=age_vect, 
                dat=df.out, 
                iterlim=1000,
                hessian = T)
  
  
  par.dat$lambda <- exp(out.NS$estimate)
  par.dat$llik <- out.NS$minimum
  par.dat$convergence <- out.NS$code
  
  if(fit.CI==TRUE){
    
    if (is.positive.definite(out.NS$hessian)==TRUE){
      hess <- solve(out.NS$hessian)
      prop_sigma <-sqrt(diag(hess))
      upper<-exp(out.NS$par)+1.96*prop_sigma
      lower<-exp(out.NS$par)-1.96*prop_sigma
      CI <-data.frame(lower=lower, upper=upper)
      CI$lower[CI$lower<0] <- 0
      CI$upper[CI$upper<0] <- 0
      
      par.dat$lci <- CI$lower
      par.dat$uci <- CI$upper
      
    }else{
      par.dat$lci <- "not yet"
      par.dat$uci <- "not yet"
      
      # #now apply over all the parameters to get your likelihoods across this range of lambda values
      # index.list <- as.list(1:nrow(par.dat))
      # par.dat$lambda_min <- 0.000001
      # par.dat$lambda_min[par.dat$lambda_min>par.dat$lambda] <- par.dat$lambda[par.dat$lambda_min>par.dat$lambda]/10
      # par.dat$lambda_max <- 1
      # par.dat$lambda_max[par.dat$lambda_max<par.dat$lambda] <- par.dat$lambda[par.dat$lambda_max<par.dat$lambda]*10
      # 
      # out.list <- lapply(index.list, get.lliks, par.dat=par.dat, df= df.out, 
      #                    year.start = (min(dat$year)-dist.back),
      #                    n.iterations = 100)
      # 
      # #out.list returns the parameter estimates and CIs by year
      
      #bind and return back
      
      #par.dat <- data.table::rbindlist(out.list)
    }
    
  }
  
  #and return
  
  return(par.dat)
  
}


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
dat$age_new <- NA
dat$age_new[!is.na(dat$DDN)] <-(dat$sampling_date[!is.na(dat$DDN)] - dat$DDN[!is.na(dat$DDN)])/365

dat$age_new[is.na(dat$age_new)] <- dat$age[is.na(dat$age_new)] 

dat$age <- dat$age_new
#a bunch still don't have ages...
subset(dat, is.na(age))


# look at age structured prevalence by year

dat$year <- year(dat$sampling_date)
dat.pos <- subset(dat, RSV==1 & !is.na(age))
dat.pos$age <- ceiling(dat.pos$age)

#everyone gets 1 added
#0s become 1s
dat.pos$age[dat.pos$age==0] <- 1#dat.pos$age+1

max(dat.pos$age) #76
min(dat.pos$age) #1


unique(dat.pos$year)
dat.pos = subset(dat.pos, year>2010 & year<2022)

dat.pos <- arrange(dat.pos, year, sampling_date)

year.dat <- dlply(dat.pos, .(year))

year.dat.sum <- lapply(year.dat, sum.yr.all)
df.out <- data.table::rbindlist( year.dat.sum)



df.out$year <- as.factor(df.out$year)


p2 <- ggplot(data=df.out) + theme_bw() +
  geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year)) + coord_cartesian(xlim=c(0,15)) +
  ylab("cumulative proportion of cases") + facet_wrap(~year)

p2

#seems to be accelerating


p3 <-ggplot(data=df.out) + theme_bw() +
  geom_point(aes(x=age, y=cum_prop_cases, color=year)) +
  geom_line(aes(x=age, y=cum_prop_cases, color=year)) + coord_cartesian(xlim=c(0,15)) +
  ylab("cumulative proportion of cases") 


p3
#does the age distribution change with time? can you calculate 
#foi changes with time? 

#get average age by year???
dat.age <- ddply(dat.pos, .(year), summarise, mean_age = mean(age, na.rm=T))

#mea age is definitely decreasing
p4 <- ggplot(dat.age) + geom_point(aes(x=year, y=mean_age)) + geom_line(aes(x=year, y=mean_age))
p4


dat.pos <- arrange(dat.pos, sampling_date)

dat.pos$year <- as.numeric(as.character(dat.pos$year))


#plot age distribution of cases
head(dat.pos)
ggplot(dat=dat.pos) + geom_violin(aes(x=year, y=age, group=year),draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_point(aes(x=year, y=age, group=year, color=year), position = position_jitterdodge(jitter.width = .8), alpha=.2) + scale_y_log10()

#just estimate FOI below 10
#goes from 924 to 891
dat.pos = subset(dat.pos,  age<=10) #lose 

#now estimate FOI 
out <- fit.all.yrs.seq.yr.nlm(dat=dat.pos,
                              lambda.guess=0.1,
                              fit.CI = F)


