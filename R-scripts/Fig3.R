rm(list=ls())


library(ggplot2)
library(dplyr)
library(plyr)
library(MuMIn)
library(relaimpo)
library(sjPlot)


homewd="/Users/carabrook/Developer/RSV-madagascar"


dat <- read.csv(file = paste0(homewd, "/data/lagged_climate_transmission_RSV.csv"), header = TRUE, stringsAsFactors = F)


# First, make a global model of all possible predictors of transmission
# That is here
global.model <- lm(log_beta~mean_H2M + sum_precip + meanTemp + year, data= dat, na.action = na.fail)
summary(global.model)

# Now, use dredge to compare and rank all possible model combinations
out.comp <- dredge(global.model = global.model) 
print(out.comp)

# Top performing model has 3 predictors: humidity, temp, and precip (no year)

# For the first panel of Fig 3, we will compare all model predictors in a heatmap


# First, load functions
calc.one.glm <- function(mod, dat, delta_AICc, model_num){
  
  #run the model
  
  m.out <- glm(mod$call[2], data=dat, na.action = na.fail)
  
  #get the R2 of the whole model and the relative contributions of each predictor,
  #both forced to 100% and as a subset
  
  if(length(attr(mod$terms, "dataClasses"))>2){
    
    
    #out1 <- domin(formula_overall = formula(mod$call[[2]]), reg=glm, fitstat=list(pscl::pR2, "r2CU"), data=dat, family = poisson())
    out1 <-calc.relimp(m.out,rela=F)
    out2 <-calc.relimp(m.out,rela=T)
    
    
    
    # com.df <- cbind.data.frame(predictor= names(out1$General_Dominance),
    #                            modelRsq = rep(out1$Fit_Statistic_Overall, length(names(out1$General_Dominance))),
    #                            lmg_raw = out1$General_Dominance,
    #                            lmg_percent = out1$Standardized)
    # # 
    # 
    # 
     com.df <- cbind.data.frame(predictor= rownames(out1@ave.coeffs),
                                modelRsq = rep(out1@R2, length(names(out1@lmg))),
                                lmg_raw = out1@lmg,
                                lmg_percent = out2@lmg)
    
  }else{
    
    sum.df = summary(mod)
    
    com.df <- cbind.data.frame(predictor= names(attr(mod$terms, "dataClasses")[2]),
                               modelRsq = sum.df$r.squared,
                               lmg_raw = sum.df$r.squared,
                               lmg_percent = 1)
    
    # com.df <- cbind.data.frame(predictor= names(attr(mod$terms, "dataClasses")[2]),
    #                            modelRsq = pR2(mod)['McFadden'],
    #                            lmg_raw = pR2(mod)['McFadden'],
    #                            lmg_percent = 100)
    
  }
  
  com.df$AICc <- AICc(mod)
  com.df$delta_AICc <- delta_AICc
  com.df$model_num <- model_num
  #and mark which predictors are significant
  
  #sum.mod <- summary(mod)
  return(com.df)
}
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Subselect the top 10 models to compare
AIC.max <-sort(out.comp$AICc)[10] # could change if you wanted more or fewer 
out.sum = get.models(out.comp, subset=AICc<=AIC.max) # This function is from MuMIN package - this selects the top 10 models from dredge and formats as a list

# Then calculate the relative contribution of each predictor within each model
# To the total R-squared within each mode
# This uses functions I wrote above. make sure '1:10' after 'model_num' and 'delta' terms
# matches the number of models selected under AIC.max (ends in 10 here but could be 5 or 15, etc)
rel.comps <- mapply(calc.one.glm, mod=out.sum, model_num=as.list(1:10), 
                    delta_AICc = as.list(out.comp$delta[1:10]), 
                    MoreArgs = list(dat=dat), SIMPLIFY = F)

# Bind output into a data table
comp.df <- data.table::rbindlist(rel.comps)

# NA predictors are random intercepts, so rename these
comp.df$predictor[is.na(comp.df$predictor)] <- "random\nintercept"

# Rename your other predictors for plotting purposes
comp.df$predictor[comp.df$predictor=="mean_H2M"] <- "mean\nhumidity"
comp.df$predictor[comp.df$predictor=="meanTemp"] <- "mean\ntemperature"
comp.df$predictor[comp.df$predictor=="sum_precip"] <- "sum\nprecipitation"

# Plot the output of all the models as a heatmap. Here is the main body of the plot panel
p1a <- ggplot(data=comp.df) + geom_tile(aes(x=predictor, y=model_num, fill=lmg_percent), color="gray", size=1) +
  scale_fill_viridis_c(limits=c(0,1), direction = -1) + scale_y_reverse() + theme_bw() + 
  #facet_grid(variable~DENV.serotype) +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), #axis.title.y = element_text(size=14), 
        legend.position = "top",
        axis.text.x = element_text(size=10),
        #axis.text.y = element_text(size=12)
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + labs(fill="Relative contribution\nto standardized R-sq.") 


# Here, we rank by delta AIC
p1b <- ggplot(data=comp.df) + geom_label(aes(x=1, y=model_num, label=specify_decimal(delta_AICc,2))) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_blank(), axis.ticks=element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_blank(),
                     plot.margin = unit(c(1.8,0,2.1,0), "cm"),
                     plot.title = element_text(hjust=0.5, face = "bold", vjust=9),
                     panel.border = element_blank())+labs(title="delta\nAICc") + scale_y_reverse()

# Here we rank by R-squared
p1c <- ggplot(data=comp.df) + geom_label(aes(x=0, y=model_num, label=specify_decimal(modelRsq,2))) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_blank(), axis.ticks=element_blank(),
                     panel.background = element_blank(),
                     plot.background = element_blank(),
                     plot.margin = unit(c(2.3,0,2.1,0), "cm"),
                     plot.title = element_text(hjust=0.5, face = "bold", vjust=10),
                     panel.border = element_blank()) +labs(title="R-sq.") + scale_y_reverse()

# Here we combine in panel A of Fig 3
p1all <- cowplot::plot_grid(p1a, p1b, p1c, rel_widths = c(1,.1,.1), ncol = 3, nrow = 1)
p1all <- cowplot::plot_grid(p1b, p1c, p1a, rel_widths = c(.1,.1, 1), ncol = 3, nrow = 1)

Fig3A <- p1all + theme(panel.border = element_rect(color="white", fill=NA), panel.background = element_rect(fill="white", colour = "white"),  
                       plot.margin = unit(c(.2,.2,.4,.1), "cm")) 

print(Fig3A)

# Now, next, we plot the elements of the top model

# Here is the top model
top_model <- lm(log_beta~mean_H2M + sum_precip + meanTemp, data= dat, na.action = na.fail)
summary(top_model)

# First, let's plot the main predictors
est.df <- get_model_data(top_model, type="est")
est.df$term <- as.character(est.df$term)
est.df$term[est.df$term=="mean_H2M"] <- "mean\nhumidity"
est.df$term[est.df$term=="sum_precip"] <- "sum\nprecipitation"
est.df$term[est.df$term=="meanTemp"] <- "mean\ntemperature"

est.df$term <- factor(est.df$term, levels=rev(c("sum\nprecipitation", "mean\nhumidity", "mean\ntemperature")))
est.df$posneg <- "pos"
est.df$posneg[est.df$estimate<0 & est.df$conf.high<0] <- "neg"
est.df$posneg[est.df$estimate<0 & est.df$conf.high>0] <- "notsig"
est.df$posneg[est.df$estimate>0 & est.df$conf.low<0] <- "notsig"
est.df$p.stars[est.df$p.stars==""] <- "."

colz <- c('pos' = "red3", 'neg'="cornflowerblue", 'notsig'="gray50")


Fig3B <- ggplot(data=est.df) + coord_flip(ylim=c(-.007,0.016)) + scale_color_manual(values=colz) +
         geom_label(aes(x=term, y=0.015, label=p.stars), label.size = 0, size=8) +
         geom_point(aes(x=term, y = estimate, color=posneg), show.legend = F, size=3) +
         ylab("slope") + geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() +
         geom_linerange(aes(x=term, ymin=conf.low, ymax=conf.high, color=posneg), show.legend = F, size=1) +
         theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18),
         axis.title.y = element_blank(), axis.text = element_text(size = 13),
         plot.margin = unit(c(2.2,.2,.2,.2), "cm")) + scale_y_continuous(breaks=c(-0.005,0, 0.005, 0.01))

# Now, show predictions based on where eahc predictor drives transmission (up or down)

pred.df <- get_model_data(top_model, type="pred")

pred.df <- data.table::rbindlist(pred.df)
head(pred.df)
pred.df$posneg <- "pos"
pred.df$posneg[pred.df$group_col=="mean_H2M"] <- "neg"
pred.df$posneg[pred.df$group_col=="meanTemp"] <- "notsig"

pred.df$group_col <- as.character(pred.df$group_col)
pred.df$group_col[pred.df$group_col=="mean_H2M"] <- "atop('mean humidity','(%)')"
pred.df$group_col[pred.df$group_col=="sum_precip"] <- "atop('sum precipitation','(mm)')"
pred.df$group_col[pred.df$group_col=="meanTemp"] <- "atop('mean temperature', '('^0~'C)')"

pred.df$predicted <- exp(pred.df$predicted)
pred.df$conf.low <- exp(pred.df$conf.low)
pred.df$conf.high <- exp(pred.df$conf.high)

Fig3C <- ggplot(data=pred.df) + facet_grid(~group_col,scales = "free_x", labeller = label_parsed) +
  geom_line(aes(x=x, y=predicted, color=posneg), size=1, show.legend = F) + ylab(bquote('predicted transmission,'~beta)) +
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high, fill=posneg), alpha=.3, show.legend = F) + theme_bw() +
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  theme(panel.grid = element_blank(), strip.text = element_text(size=18),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_blank(), axis.text = element_text(size = 13),
        plot.margin = unit(c(.2,.2,.2,.9), "cm")) 
Fig3C 

# And combine into Figure 3

Fig3AB <- cowplot::plot_grid(Fig3A, Fig3B, ncol = 2, nrow = 1, labels=c("A", "B"), label_size = 22, rel_widths = c(1,.63))

Fig3 <- cowplot::plot_grid(Fig3AB, Fig3C, ncol = 1, nrow = 2, labels = c("", "C"), label_size = 22, rel_heights = c(1,.9)) + theme(plot.background = element_rect(fill="white"))


ggsave(file = paste0(homewd, "/figures/Fig3.png"),
       plot = Fig3,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)
