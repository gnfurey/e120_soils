#
source("e120soils_functions.R")
library(tidyverse)
dat <- read.csv("SubmissionData/Dataset_S3.csv")
#
base <- read.csv("Data/data_e120soils_baseset.csv")
#
dat <- left_join(dat,base)
#there is also nitrate data, but not pre-treatment nitrate data 
#nitrate also showed low values. 
#these soils are low in N 
dat %>% group_by(Year,Depth,nut) %>%
  filter(nut=="Nitrogen") %>% 
  mutate(val=val/10000) %>% 
  summarise(se=my.stand(val,na.rm = TRUE),
            val=mean(val,na.rm = TRUE)) %>% 
  filter(Depth=="0.20cm") %>% 
  filter(Year==1994)
#
dat %>% group_by(Year,Depth,nut) %>%
  mutate(val=val/10000) %>% 
  filter(nut=="Organic.Matter") %>% 
  summarise(se=my.stand(val,na.rm = TRUE),
            val=mean(val,na.rm = TRUE)) %>% 
  filter(Depth=="0.20cm") %>% 
  filter(Year==1994)
#
dat %>% group_by(Year,Depth,nut) %>%
  filter(nut=="Phosphorus") %>% 
  summarise(se=my.stand(val,na.rm = TRUE),
            val=mean(val,na.rm = TRUE)) %>% 
  filter(Depth=="0.20cm") %>% 
  filter(Year==1994)
#
dat %>% group_by(Year,Depth,nut) %>%
  filter(nut=="Potassium") %>% 
  summarise(se=my.stand(val,na.rm = TRUE),
            val=mean(val,na.rm = TRUE)) %>% 
  filter(Depth=="0.20cm") %>% 
  filter(Year==1994)

