rm(list=ls())

library(ggplot2)
library(lubridate)
library(plyr)
library(dplyr)
library(mgcv)
library(matrixcalc)



homewd="/Users/carabrook/Developer/RSV-madagascar"
# Tsiry, here add your own directory in place of mine:
#homewd="path_to_Tsiry_directory"
setwd(homewd)


# first, prep the data
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

dat$sampling_date <- as.Date(dat$sampling_date, format = "%m/%d/%y")
dat$year <- year(dat$sampling_date)

#plot prevalence by age 

age.dat <- subset(dat, !is.na(age))
age.dat$age <- ceiling(age.dat$age)
age.dat$age[age.dat$age==0] <- 1 
# 
# #plot
# 
age.count <- ddply(age.dat, .(age), summarise, N=sum(RSV))

all.age <- cbind.data.frame(age=0:89)
all.age <- merge(all.age, age.count, by="age", all.x = T)
all.age$N[is.na(all.age$N)] <- 0

all.age$cum_RSV <- cumsum(all.age$N)
all.age$Ntot <- max(all.age$cum_RSV)
all.age$cum_prop_RSV <- all.age$cum_RSV/max(all.age$cum_RSV)
ggplot(all.age) + geom_point(aes(x=age, y=cum_prop_RSV)) + geom_line(aes(x=age, y=cum_prop_RSV)) + theme_bw()


#now get foi
model.age.prev <- function(lambda, data, cate){ 
  #takes in a value for lambda from each age bin, the data itself, and a vector that gives the lower bound of every age bin
  dur=c(diff(cate), 0) #calculates the length of each age bin
  p=as.list(rep(NA,max(data$age))) #make a list to store the prevalence by age predicted from the model
  for(a in 1:max(data$age)){ #for-loops over all possible ages in our data range
    dummy1=a>cate  #creates dummy variables to identify which age class this individual falls in
    dummy2 = a>cate & !c(a>cate[-1], FALSE) #adds a false at the end and takes the inverse to make true only those values above you
    dummy1=c(a>cate, FALSE)[-1]# and dummy 1 gives you a true for all age bins below you
    inte=sum(dur*lambda*dummy1) + lambda[dummy2]*(a-cate[dummy2]) 
    #integrates over the amount of time in previous age classes + how long you've been in this one as a susceptible
    #to give you the hazard of infection
    p[[a]]=1-exp(-inte) #returns the probability of infection for a given age.
  }
  return(p) #returns prevalence by age for the length of all possible ages in our dataset
}
log.lik.fit <- function(lambda.guess, data.in, lower_age_bin){
  
  out.prev = cbind.data.frame(age=unique(data.in$age), model_prop=c(0,unlist(model.age.prev(lambda = lambda.guess, data=data.in, cate=lower_age_bin))))
  
  #head(out.prev)
  #ggplot(data=out.prev) + geom_point(aes(x=age,y= model_prop))  + geom_line(aes(x=age,y= model_prop)) +theme_bw()
  
  #now merge with data
  
  merge.dat <- merge(data.in, out.prev, by="age", all.x = T)
  
  merge.dat <- arrange(  merge.dat, age)
  #ggplot(data=merge.dat) + geom_point(aes(x=age,y= cum_prop_RSV))  + geom_line(aes(x=age,y= cum_prop_RSV)) +geom_point(aes(x=age,y= model_prop), color="red")  + geom_line(aes(x=age,y= model_prop), color="red") + theme_bw()
  #head(merge.dat)
  
  merge.dat$Ntot = sum(merge.dat$N)
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(merge.dat$age)){
    ll=ll+dbinom(merge.dat$cum_RSV[i],merge.dat$Ntot[i],prob=merge.dat$model_prop[i],log=T)
    #print(ll)
  }
  
  
  if(ll==-Inf){
    ll <- -1000000
  }
  
  return(-ll)
}
fit.foi <- function(lambda.guess, data.in, lower_age_bin_vector, hyp, do.plot){
  
  
  data.in <- arrange(data.in, age)
  
  out.NS <- optim(par = lambda.guess, 
                  fn=log.lik.fit, 
                  method = "L-BFGS-B",
                  lower = 0.0000001,
                  upper=10,
                  data.in=data.in, 
                  control = list(maxit=500),
                  lower_age_bin = lower_age_bin_vector,
                  hessian = T)
  
  #get confidence intervals
  
  if (is.positive.definite(out.NS$hessian)==TRUE){
    hess <- solve(out.NS$hessian)
    prop_sigma <-sqrt(diag(hess))
    upper<-out.NS$par+1.96*prop_sigma
    lower<-out.NS$par-1.96*prop_sigma
    CI <-data.frame(lower=lower, upper=upper)
    CI$lower[CI$lower<0] <- 0
    CI$upper[CI$upper<0] <- 0
  
    #and return estimates and fit
    est.dat <- cbind.data.frame(age_bin_lower = lower_age_bin_vector, age_bin_upper = c(lower_age_bin_vector[-1], max(data.in$age)), foi = out.NS$par, foi_lci = CI$lower, foi_uci=CI$upper, neg_llik = rep(out.NS$value, length(lambda.guess)), convergence = rep(out.NS$convergence, length(lambda.guess)))
    
  }else{
    est.dat <- cbind.data.frame(age_bin_lower = lower_age_bin_vector, age_bin_upper = c(lower_age_bin_vector[-1], max(data.in$age)), foi = out.NS$par, foi_lci = rep(NA, length(lambda.guess)), foi_uci=rep(NA, length(lambda.guess)), neg_llik = rep(out.NS$value, length(lambda.guess)), convergence = rep(out.NS$convergence, length(lambda.guess)))
  }
    
  
  est.dat$npar <- rep(length(lambda.guess), length(lambda.guess))
  est.dat$hyp <- hyp
  
  est.dat$AIC <- (2*nrow(est.dat))+(2*est.dat$neg_llik)
  
  # also, plot with data
  if (do.plot==TRUE){
    
    out.prev = cbind.data.frame(age=unique(data.in$age), model_prop=c(0,unlist(model.age.prev(lambda = out.NS$par, data=data.in, cate=lower_age_bin_vector))))
    out.prev.lci = cbind.data.frame(age=unique(data.in$age), model_prop_lci=c(0,unlist(model.age.prev(lambda = est.dat$foi_lci, data=data.in, cate=lower_age_bin_vector))))
    out.prev.uci = cbind.data.frame(age=unique(data.in$age), model_prop_uci=c(0,unlist(model.age.prev(lambda = est.dat$foi_uci, data=data.in, cate=lower_age_bin_vector))))

    out.prev <- merge(out.prev, out.prev.lci, by="age")    
    out.prev <- merge(out.prev, out.prev.uci, by="age")
    #head(out.prev)
    #ggplot(data=out.prev) + geom_point(aes(x=age,y= model_prop))  + geom_line(aes(x=age,y= model_prop)) +theme_bw()
    
    #now merge with data
    
    merge.dat <- merge(data.in, out.prev, by="age", all.x = T)
    
    merge.dat <- arrange(  merge.dat, age) 
    p1 <- ggplot(data=merge.dat) + 
          geom_line(aes(x=age,y= cum_prop_RSV)) +
          geom_ribbon(aes(x=age,ymin= model_prop_lci, ymax=model_prop_uci), fill="red", alpha=.3)  + 
          geom_line(aes(x=age,y= model_prop), color="red") + theme_bw() +
          ylab("cumulative proportion of RSV cases")+
          geom_point(aes(x=age,y= cum_prop_RSV))  
    
    foi.dat = cbind.data.frame(age=c(rep(lower_age_bin_vector, each=2)[-1], max(data.in$age)), foi=rep(out.NS$par, each=2), foi_lower = rep(est.dat$foi_lci, each=2), foi_upper = rep(est.dat$foi_uci, each=2))

    p2 <- ggplot(foi.dat) + geom_ribbon(aes(x=age, ymin=foi_lower, ymax=foi_upper), fill="red", alpha=.3) +
          geom_line(aes(x=age, y=foi), color="red")  + theme_bw()
    
    pall <- cowplot::plot_grid(p1,p2, labels = c("A", "B"), label_size = 22)
    print(pall)
  }
  
  return(est.dat)
}
midpt <- function(vector){
  mid <- list()
  for (i in 1:(length(vector)-1)){
    mid[[i]] <-  (vector[i+1] -vector[i])/2 + vector[i]
  }
  mid <- c(unlist(mid))
  return(mid)
}

# and compare the following hypotheses
h0 <- fit.foi(lambda.guess = c(.1),
                data.in=all.age,
                hyp = 0,
                lower_age_bin_vector=c(0),
                do.plot = TRUE)


#then, try two classes, with a focus on the younger years

h1 <- fit.foi(lambda.guess = c(1,.1),
                     data.in=all.age,
                     hyp = 1,
                     lower_age_bin_vector=c(0,5),
                     do.plot = TRUE)

h2 <- fit.foi(lambda.guess = c(1,.1),
                  data.in=all.age,
                  hyp = 2,
                  lower_age_bin_vector=c(0,4),
                  do.plot = TRUE)

h3 <- fit.foi(lambda.guess = c(1,.1),
                  data.in=all.age,
                  hyp = 3,
                  lower_age_bin_vector=c(0,3),
                  do.plot = TRUE)

h4 <- fit.foi(lambda.guess = c(1,.1),
                  data.in=all.age,
                  hyp = 4,
                  lower_age_bin_vector=c(0,2),
                  do.plot = TRUE)

#not as good - needs the two
h5 <- fit.foi(lambda.guess = c(1,.1),
                  data.in=all.age,
                  hyp = 5,
                  lower_age_bin_vector=c(0,1),
                  do.plot = TRUE)

#move to three classes
h6 <- fit.foi(lambda.guess = c(1,.1,.1),
                  data.in=all.age,
                  hyp = 6,
                  lower_age_bin_vector=c(0,1,2),
                  do.plot = TRUE)


h7 <- fit.foi(lambda.guess = c(1,.1,.1,.1),
                   data.in=all.age,
                   hyp = 7,
                   lower_age_bin_vector=c(0,1,2,3),
                   do.plot = TRUE)

h8 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1),
                    data.in=all.age,
                    hyp = 8,
                    lower_age_bin_vector=c(0,1,2,3,4),
                    do.plot = TRUE)

#not as good - no need for class age 5
h9 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1),
                     data.in=all.age,
                     hyp = 9,
                     lower_age_bin_vector=c(0,1,2,3,4,5),
                     do.plot = TRUE)

#now, add to the rest of the ages from 4
#gets a lot better
h10 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
                      data.in=all.age,
                      hyp = 10,
                      lower_age_bin_vector=c(0,1,2,3,4,10,20,30, 40, 50, 60,70,80),
                      do.plot = TRUE)

#80 is too high - 70 is sufficient
h11 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
                            data.in=all.age,
                            hyp = 11,
                            lower_age_bin_vector=c(0,1,2,3,4,10,20,30, 40, 50,60, 70),
                            do.plot = TRUE)

#not enough if you stop at 50
h12 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
                            data.in=all.age,
                            hyp = 12,
                            lower_age_bin_vector=c(0,1,2,3,4,10,20,30, 40, 50, 60),
                            do.plot = TRUE)


#what about losing some classes in the middle?
#even better
h13 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 13,
                                lower_age_bin_vector=c(0,1,2,3,4,10,30, 40, 50, 60,70),
                                do.plot = TRUE)
#and keep at it
h14 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 14,
                                lower_age_bin_vector=c(0,1,2,3,4,10,30, 50, 60,70),
                                do.plot = TRUE)

#and again - same as h14 - but fewer pars
h15 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 15,
                                lower_age_bin_vector=c(0,1,2,3,4,10,30,60,70),
                                do.plot = TRUE)

#not as good
h16 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 16,
                                lower_age_bin_vector=c(0,1,2,3,4,10,30,70),
                                do.plot = TRUE)

#best yet - fwer pars (8)
h17 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 17,
                                lower_age_bin_vector=c(0,1,2,3,4,10,60,70),
                                do.plot = TRUE)

#not as good - but pretty close
h18 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 18,
                                lower_age_bin_vector=c(0,1,2,3,4,60,70),
                                do.plot = TRUE)

#not as good - but close
h19 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 19,
                                lower_age_bin_vector=c(0,1,2,3,4,5,60,70),
                                do.plot = TRUE)

#best yet
h20 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 20,
                                lower_age_bin_vector=c(0,1,2,3,4,20,60,70),
                                do.plot = TRUE)

#too many, no need to add
h21 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 21,
                                lower_age_bin_vector=c(0,1,2,3,4,10,20,60,70),
                                do.plot = TRUE)

#is 20 too high? 15 is worse actually
h22 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 22,
                                lower_age_bin_vector=c(0,1,2,3,4,15,60,70),
                                do.plot = TRUE)
#is 20 too low? - 25 the same as 20
h23 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 23,
                                lower_age_bin_vector=c(0,1,2,3,4,25,60,70),
                                do.plot = TRUE)

#30 too high
h24 <- fit.foi(lambda.guess = c(1,.1,.1,.1,.1,.1,.1,.1),
                                data.in=all.age,
                                hyp = 24,
                                lower_age_bin_vector=c(0,1,2,3,4,30,60,70),
                                do.plot = TRUE)

# and compare against all ages - did not converge - 
# but gives a sense of a similar shape
h25 <- fit.foi(lambda.guess = rep(.1,90),
                                data.in=all.age,
                                hyp = 25,
                                lower_age_bin_vector=0:89,
                                do.plot = TRUE)


# and combine the hypotheses into a table
dat.all <- rbind(h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18,h19,h20,h21,h22,h23,h24,h25)

#simplify 
dat.table <- ddply(dat.all,.(hyp), summarise, npar = unique(npar), neg_llik=unique(neg_llik), AIC=unique(AIC), convergence = unique(convergence))
dat.table$convergence[dat.table$convergence==0] <- "yes"
dat.table$convergence[dat.table$convergence==1] <- "no"
dat.table$convergence[dat.table$convergence==52] <- "no"
write.csv(dat.table, file = paste0(homewd, "/data/Table_SX_foi.csv"), row.names = F)

#and the final shape for the actual figure
fig.2 <- function(lambda, lambda_lci, lambda_uci, data.in, lower_age_bin_vector,do.save, filename){
  
  
    data.in <- arrange(data.in, age)
  
    out.prev = cbind.data.frame(age=unique(data.in$age), model_prop=c(0,unlist(model.age.prev(lambda = lambda, data=data.in, cate=lower_age_bin_vector))))
    out.prev.lci = cbind.data.frame(age=unique(data.in$age), model_prop_lci=c(0,unlist(model.age.prev(lambda = lambda_lci, data=data.in, cate=lower_age_bin_vector))))
    out.prev.uci = cbind.data.frame(age=unique(data.in$age), model_prop_uci=c(0,unlist(model.age.prev(lambda = lambda_uci, data=data.in, cate=lower_age_bin_vector))))
    
    out.prev <- merge(out.prev, out.prev.lci, by="age")    
    out.prev <- merge(out.prev, out.prev.uci, by="age")
    #head(out.prev)
    #ggplot(data=out.prev) + geom_point(aes(x=age,y= model_prop))  + geom_line(aes(x=age,y= model_prop)) +theme_bw()
    
    #now merge with data
    
    merge.dat <- merge(data.in, out.prev, by="age", all.x = T)
    
    merge.dat <- arrange(  merge.dat, age) 
    p1 <- ggplot(data=merge.dat) + 
      geom_line(aes(x=age,y= cum_prop_RSV)) +
      geom_ribbon(aes(x=age,ymin= model_prop_lci, ymax=model_prop_uci), fill="red", alpha=.3)  + 
      geom_line(aes(x=age,y= model_prop), color="red") + theme_bw() +
      ylab("cumulative proportion of RSV cases")+
      geom_point(aes(x=age,y= cum_prop_RSV))  +
      theme(axis.title = element_text(size=16), axis.text = element_text(size=12))
    
    foi.dat = cbind.data.frame(age=c(rep(lower_age_bin_vector, each=2)[-1], max(data.in$age)), foi=rep(lambda, each=2), foi_lower = rep(lambda_lci, each=2), foi_upper = rep(lambda_uci, each=2))
    
    p2 <- ggplot(foi.dat) + geom_ribbon(aes(x=age, ymin=foi_lower, ymax=foi_upper), fill="red", alpha=.3) +
      geom_line(aes(x=age, y=foi), color="red")  + theme_bw() + ylab(bquote(lambda~', force of infection')) +
      theme(axis.title = element_text(size=16), axis.text = element_text(size=12))
    
    pall <- cowplot::plot_grid(p1,p2, labels = c("A", "B"), label_size = 22)
  
    
    if(do.save==TRUE){
      ggsave(file = filename,
             plot = pall,
             units="mm",  
             width=90, 
             height=45, 
             scale=3, 
             dpi=300)
    }
  
  
}

fig.2(lambda = h20$foi,
      lambda_lci = h20$foi_lci,
      lambda_uci = h20$foi_uci,
      data.in = all.age,
      lower_age_bin_vector = c(0,1,2,3,4,20,60,70),
      do.save=TRUE,
      filename = paste0(homewd, "/figures/Fig2.png"))

