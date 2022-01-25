

##### ####

# processData.R
# should be called from SPIKEanalysis.Rmd

require(magrittr)
require(readxl)
library(tidyr)
library(dplyr)

library(ggplot2)

library(nlme)

setwd("O:\\HSL\\HSL_COVID-19\\Chunming Zhu\\SARS\\ELISA\\Chunming_result_28Sep21\\spike_sep21\\")


############
# cutpoint #
############

# read i n and process LLOQ data


cutpoint <- read.csv("spike_all_linear.csv", header=T) #, na = c('', 'Sample did not dilute down properly-can not use data')) # %>%

cutpoint$lg_Acon <- log(cutpoint$Acon)
####################################
### try ###

IgG <- cutpoint[cutpoint$Assay=="Spike IgG",]
IgM <- cutpoint[cutpoint$Assay=="Spike IgM",]





data1 <- IgG[lo <= IgG$OD & IgG$OD <= hi, ]
model0<- lme(lg_Acon ~ OD , random = ~1 |SID, data=data1, na.action = na.omit)
out <- cbind(coef=summary(model0)$coefficients$fixed ,  anova(model0)[4]) [2,] 
 out$low_point <-lo
 out$hi_point<-hi
 
#######################################
### screen cut off range IgG
#######################################
 
 d<-NULL

 for (i in 1:14) {
    lo <- i/20
     for (j in 12:38 ) {
     
      hi <- j/10
      data1 <- IgG[lo <= IgG$OD & IgG$OD <= hi, ]
      
      model<-  lme(lg_acon ~ lg_OD , random = ~1 |SID, data=data1, na.action = na.omit, method="ML")
      df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
      data.frame(df) 
      df <- cbind(low_point= lo,df) 
      df <-  cbind(hi_point= hi, df)
     df_ij<- data.frame(df) 
     d<-rbind(d,df_ij) 
  
    }
   
 }
 final_igg <- d[substr(rownames(d),2,5) != "Inter" & d$p.value >0.05, ]  #0.05--2.1
 
 final$hi_od <- 4*exp(final$hi_point)/(1 + exp(final$hi_point))
 final$lo_od <- 4*exp(final$low_point)/(1 + exp(final$low_point))
 
 write.csv(d,"linearity_random_IgG.csv")
 
 
 
 ####################################################################################
 
 #######################################
 ### screen cut off range IgM
 #######################################
 
 d<-NULL
 
 for (i in 1:14) {
   lo <- i/20
   for (j in 12:38 ) {
     
     hi <- j/10
     data1 <- IgM[lo <= IgM$OD & IgM$OD <= hi, ]
     
     model<-  lme(lg_acon ~ lg_OD , random = ~1 |SID, data=data1, na.action = na.omit, method="ML")
     df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
     data.frame(df) 
     df <- cbind(low_point= lo,df) 
     df <-  cbind(hi_point= hi, df)
     df_ij<- data.frame(df) 
     d<-rbind(d,df_ij) 
     
   }
   
 }
 
 final_igM <- d[substr(rownames(d),2,5) != "Inter" & d$p.value >0.05, ]  #0.05--2.1
 
 final$hi_od <- 4*exp(final$hi_point)/(1 + exp(final$hi_point))
 final$lo_od <- 4*exp(final$low_point)/(1 + exp(final$low_point))
 write.csv(d,"linearity_random_IgM.csv")
 
 
 ###################################################################################
 ############# plot the fitted model ##############################################
 ##################################################################################
 
 IgG  %>% ggplot(aes(x=lg_OD, y = Result)) + 
    geom_point(alpha=0.1) + geom_smooth(method='lm') +
    xlab('Log_OD') +
    scale_y_log10() +
    
    ylab('Titer(AU/mL)')   +  
    theme_bw()+
   scale_x_continuous(breaks = seq(-7, 5, 1.0)) +
    theme(legend.position  = 'bottom',legend.title = element_blank()) 

 
 
 IgM %>% filter(lg_OD < 5) %>% ggplot(aes(x=lg_OD, y = Result)) + 
    geom_point(alpha=0.1) + geom_smooth(method='lm') +
    xlab('Log_OD') +
    scale_y_log10() +
    ylab('Titer(AU/mL)')   +  
    scale_x_continuous(breaks = seq(-7, 5, 1.0)) +
    theme_bw()+
    theme(legend.position  = 'bottom',legend.title = element_blank()) 
 
 
 ###################################################################################
 ############# plot the stability data ###########################################
 ##################################################################################
   freeze <-read.csv("stability_freeze.csv", header=T)  
 freeze$Treatment =  factor(freeze$Treatment, levels= c("1X", "5X","10X"))
 
 
 freeze %>% ggplot(aes(x=Treatment, y = geomean,group=Sample_ID)) + 
    geom_line(aes(colour=Sample_ID), lwd=1) +
    xlab('Treatment : Freeze and Thraw times') +
    scale_y_log10() +
    facet_wrap(~Assay) +
    ylab('Titer (AU/mL)')   +  
       theme_bw()+
    theme(legend.position  = 'null',legend.title = element_blank()) 
 
 
  