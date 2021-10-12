#
library(tidyverse)
source("e120soils_functions.R")
bd <- read.csv("Data/data_e120soils_bulk.density_gcm3.csv")
base <- read.csv("Data/data_e120soils_baseset.csv")
dat <- left_join(bd,base) %>% 
  filter(is.na(NumSp)==FALSE)#we measured some bareground plots 
#
#OR
dat1 <- read.csv("SubmissionData/Dataset_S4.csv")
dat1 <- dat1 %>% select(Plot,NumSp,bd.pred,Bulk.Density.g.cm3_0.20cm_2017) %>% 
  filter(bd.pred=="measured") %>% 
  distinct()  
#
####################################
#regression 
mod1 <- lm(Bulk.Density.g.cm3_0.20cm_2017~
           log(NumSp),
           data=dat)
summary(mod1)
#
mod1 <- lm(Bulk.Density.g.cm3_0.20cm_2017~
             log(NumSp),
           data=dat1)
summary(mod1)
#
#get means 
datm <- dat %>% 
  group_by(NumSp) %>% 
  summarise(
    bd=mean(Bulk.Density.g.cm3_0.20cm_2017))
#get percent difference 
datmw <- datm %>% 
  pivot_wider(names_from = "NumSp",values_from = "bd")
pchange(new = datmw$`16`,old = datmw$`1`)
#range
range(datm$bd)
range(dat$Bulk.Density.g.cm3_0.20cm_2017)
#get mean SE 
datm1 <- dat %>% 
  group_by(NumSp) %>% 
  summarise(
    se=my.stand(Bulk.Density.g.cm3_0.20cm_2017),
    bd=mean(Bulk.Density.g.cm3_0.20cm_2017))
datm1
#show the data 
ggplot(datm1,aes(x=NumSp,y=bd))+
  geom_jitter(data=dat,aes(y=Bulk.Density.g.cm3_0.20cm_2017),
             shape=21,alpha=0.5,
             height=0,width=0.1)+
  geom_point()+
  geom_errorbar(aes(ymax=bd+se,ymin=bd-se))+
  theme_bw()+
  xlab("Number of Species")+
  scale_x_continuous(breaks=c(1,2,4,8,16))+
  ylab(expression(paste("Soil Bulk Density (0-20 cm) ","(","g"%.%"cm"^-3,")")))+
  theme(panel.grid.minor = element_blank())
#
