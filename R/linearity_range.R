

##### ####

# processData.R
# should be called from SPIKEanalysis.Rmd

require(magrittr)
require(readxl)
library(tidyr)
library(dplyr)

library(ggplot2)

library(nlme)

setwd("O:\\HSL\\HSL_COVID-19\\Chunming Zhu\\SARS\\Chunming_result\\trying\\")

############
# cutpoint #
############

# read i n and process LLOQ data


cutpoint <- read.csv("COVID19 neucleocapsid_lin_all820.csv", header=T) #, na = c('', 'Sample did not dilute down properly-can not use data')) # %>%

cutpoint$lg_Acon <- log(cutpoint$Acon)
####################################
### try ###

IgG <- cutpoint[cutpoint$Assay=="Nucleocapsid IgG",]
IgM <- cutpoint[cutpoint$Assay=="Nucleocapsid IgM",]

data1 <- IgG[lo <= IgG$OD & IgG$OD <= hi, ]
model0<- lme(lg_acon ~ OD , random = ~1 |SID, data=data1, na.action = na.omit)
out <- cbind(coef=summary(model0)$coefficients$fixed ,  anova(model0)[4]) [2,] 
 out$low_point <-lo
 out$hi_point<-hi
 
#######################################
### screen cut off range IgG
#######################################
 
 d<-NULL

 for (i in 1:35) {
    lo <- i/50
     for (j in 60:180 ) {
     
      hi <- j/50
      data1 <- IgG[lo <= IgG$OD & IgG$OD <= hi, ]
      
      model<-  lme(lg_acon ~ OD , random = ~1 |sample_idd, data=data1, na.action = na.omit, method="ML")
      df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
      data.frame(df) 
      df <- cbind(low_point= lo,df) 
      df <-  cbind(hi_point= hi, df)
     df_ij<- data.frame(df) 
     d<-rbind(d,df_ij) 
  
    }
   
 }
 final <- d[substr(rownames(d),2,5) != "Inter" & d$p.value >0.05, ]
 write.csv( final,"linearity_plusOQ_random_neucleo_IgG.csv")
 
 ##############################################################
 
 
 d<-NULL
 
 for (i in 5:35) {
    lo <- i/50
    for (j in 60:180 ) {
       
       hi <- j/50
       data1 <- IgG[lo <= IgG$OD & IgG$OD <= hi, ]
       
       model<-  lme(lg_acon ~ OD , random = ~1 + Analyst|sample_idd, data=data1, na.action = na.omit, method="ML")
       df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
       data.frame(df) 
       df <- cbind(low_point= lo,df) 
       df <-  cbind(hi_point= hi, df)
       df_ij<- data.frame(df) 
       d<-rbind(d,df_ij) 
       
    }
    
 }
 final <- d[substr(rownames(d),2,5) != "Inter" & d$p.value >0.05, ]
 
 write.csv( final,"linearity_plusOQ_random_slope_IgG.csv")
 
 
 
 ####################################################################################
 
 #######################################
 ### screen cut off range IgM
 #######################################
 
 d<-NULL
 
 for (i in 1:35) {
    lo <- i/50
    for (j in 60:180 ) {
       
       hi <- j/50
     data1 <- IgM[lo <= IgM$OD & IgM$OD <= hi, ]
     
     model<-  lme(lg_acon ~ OD , random = ~1 |sample_idd, data=data1, na.action = na.omit, method="ML")
     df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
     data.frame(df) 
     df <- cbind(low_point= lo,df) 
     df <-  cbind(hi_point= hi, df)
     df_ij<- data.frame(df) 
     d<-rbind(d,df_ij) 
     
   }
   
 }
 
 final <- d[substr(rownames(d),2,5) != "Inter" & d$p.value >0.05, ]
 write.csv(d,"linearity_plusOQ_random__neucleoIgM.csv")
 
 
 ###################################################################################
 ############# plot the fitted model ##############################################
 ##################################################################################
 
 
 lo <- 0.5
 hi <- 2.2
 data1 <- cutpoint[lo <= cutpoint$bg_OD & cutpoint$bg_OD <= hi, ]
 
 model<- lme(AIC.5pf ~ bg_OD , random = ~ 1 |SID, data=data1)
 df<-cbind(coef=summary(model)$coefficients$fixed ,  anova(model)[4]) 
 data.frame(df) 
 df <- cbind(low_point= lo,df) 
 df <-  cbind(hi_point= hi, df)
 data.frame(df)  
 
 
 data1$fit <- predict(model) 
 
 summary(model)
 
 
##### plot the points
  ggplot(data=data1, aes(x=bg_OD, y=AIC.5pf , color=SID)) + geom_point( alpha = 0.9) +

### fitted regression lines for each SID    
 geom_line(aes(y=fit,   color=SID),lwd=1) +    
  
#####Combined fitted line          
 geom_abline(intercept = 3.540365, slope = 0.098693,  lwd = 1.5, color="black") 
 

   
 

 
  