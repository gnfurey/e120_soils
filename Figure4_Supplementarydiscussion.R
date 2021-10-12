################
library(tidyverse)
#custom functions
source("e120soils_functions.R")
#read in trait data 
newdat <- read.csv("SubmissionData/Dataset_S2.csv")
#
mymod0 <- lm(RootBiomass~1,data=newdat)
mymod1 <- lm(RootBiomass~Potassium,data=newdat)
mymod2 <- lm(RootBiomass~Nitrogen,data=newdat)
mymod3 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat)
summary(mymod3)
NROW(newdat)
#
AIC(mymod0,mymod1,mymod2,mymod3)
#
mymod.step <- step(mymod3)
summary(mymod.step)
#tests 
car::Anova(mymod3,type="III")
shapiro.test(residuals(mymod3))
lmtest::bptest(mymod3)
car::ncvTest(mymod3)
#for the SI
mymod4 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat[newdat$Species!="Poapr",])
summary(mymod4)
mymod5 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat[newdat$Species!="Monfi",])
summary(mymod5)
mymod6 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat[newdat$Species %ni% 
                                                          c("Monfi","Poapr"),])
summary(mymod6)
#################
#Monoculture only       
#################
library(tidyverse)
#custom functions
source("e120soils_functions.R")
#read in trait data 
dat <- read.csv("SubmissionData/Dataset_S7.csv")
###
#these two species are planted together
#we do not have good monoculture data for them
#so we are using a duplicate value for both for root biomass
tmp <- dat %>% filter(Species =="Monfi Solri")
tmp1 <- rbind(tmp,tmp)
tmp1$Species <- c("Monfi","Monfi","Solri","Solri")
tmp2 <- dat %>% filter(Species !="Monfi Solri")
tmp1$Species
#
dat1<- rbind(tmp1,tmp2)

#get values in long format 
dat.w <- dat1 %>% 
  filter(NumSp==1) %>% 
  select(Plot,Species,nut,val,func) %>% 
  pivot_wider(names_from="nut",values_from="val") %>% 
  filter(Plot!=308)#no biomass for asctu
#get functional groups
dat.w$func <- get_funcgroup(get_specidfrom5(dat.w$Species))
#create a new data set 
newdat <- dat.w
#run regression
residplots <- function(mod){
  foo <- data.frame(fitted=fitted(mod),resid=residuals(mod))
  p1 <- ggplot(foo,aes(x=fitted(mod),y=residuals(mod)))+
    geom_point(shape=21)+
    theme_bw()
  return(p1)
}
#SI, Supplemental Discussion 
mymod0 <- lm(RootBiomass~1,data=newdat)
mymod1 <- lm(RootBiomass~Potassium,data=newdat)
mymod2 <- lm(RootBiomass~Nitrogen,data=newdat)
mymod3 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat)
summary(mymod3)
mymod4 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat[newdat$Species %ni% 
                                                          c("Poapr","Monfi"),])
summary(mymod4)
#