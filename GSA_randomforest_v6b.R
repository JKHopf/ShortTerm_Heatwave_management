# Global Sensitivity Analysis of Cascade Sorte's competition model
# Using procedures described in Harper et al. 2011 Ecological Applications

#install.packages('rpart')
#install.packages('randomForest')

library('rpart'); library(tidyverse)
library('randomForest') #Liaw and Wiener 

#A. Liaw and M. Wiener (2002). Classification and Regression by #randomForest. R News 2(3), 18--22.

# load in data
D_readin <- read.csv('Implicitv6b_MngtScen_GSA_none_20250324b.csv',head=F)

#name columns
DF <- D_readin %>%
  rename(kelp_pers_t=V1,PBEsa=V2,kelp.RK=V3,kelp.RKstdv=V4,kelp.mu=V5,kelp.ddD=V6,kelp.muvar=V7,kelp.g=V8,
         kelp.rS=V9,kelp.c=V10,kelp.rD=V11,kelp.d=V12,kelp.ajae=V13,kelp.adc=V14,kelp.h=V15,
         urchin.RU=V16,urchin.RUstdv=V17,urchin.gJ=V18,urchin.MJ=V19,
         urchin.MH=V20,urchin.ME=V21,urchin.PE=V22,urchin.PH=V23, 
         urchin.PLD=V24,urchin.w2=V25,urchin.kmin=V26)


# Change in para values
  # set up table
  DF_change = DF 
  
  for (i in 2:26){
    # DF_change[,i] <- mean(D_readin[,i])-D_readin[,i]
    DF_change[,i] <- (mean(D_readin[,i])-D_readin[,i])/sd(D_readin[,i])
  }


# remove extremes of 0 & 100% persistence
# DF_change <- DF_change %>% filter(kelp_pers_t >100 & kelp_pers_t <9900)


# Do the random forest analysis
DFA.rf <- randomForest(as.factor(kelp_pers_t)~.,data=DF_change,importance=TRUE)
DFA.impt <- importance(DFA.rf,type=1) # type 1 is mean decrease in accuracy

# Table
DFA.impt


# plot forest
DF.fit <- rpart(as.factor(kelp_pers_t)~.,data=DF_change,method='class',
                control=rpart.control(cp=0.005))
plot(DF.fit)
text(DF.fit,digits=1,xpd=TRUE)



# find values


DF_change$id <- c(1:10000)

id <- DF_change %>% filter(urchin.PH>=0.6975) %>% 
      arrange(urchin.PH) %>% 
      pull(id)

DF$urchin.PH[id[1]]



DF_change$id <- c(1:10000)

id <- DF_change %>% filter(PBEsa>=-0.4067) %>% 
  arrange(PBEsa) %>% 
  pull(id)

DF$PBEsa[id[1]]














