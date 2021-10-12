#
library(tidyverse)
#read in data for conversions
dat <- read.csv("SubmissionData/Dataset_S4.csv")
##
AppendixS3col <- colnames(dat)
write.csv(x = AppendixS3col,file = "Data/S3_DatasetS4.columns.csv",row.names = FALSE)

#functions 
source("e120soils_functions.R")
#get means 
datm <- dat %>%
  group_by(NumSp) %>% 
  summarise(across(c(gm2_2017_u,gm2_2017_add,gm2_2017_0.20cm,Depth_add_cm,gm2_2017_add),
                   ~mean(.x))) %>% 
  mutate(diff=pchange(new = gm2_2017_0.20cm,old =gm2_2017_u))
datm
####when contraction occurred
test <- dat[1,]
test
#get soil mass 
initial <- test$Bulk.Density.g.cm3_0.20cm_1994*20*10000
final <- test$Bulk.Density.g.cm3_0.20cm_2017*20*10000
diff <- test$Bulk.Density.g.cm3_0.20cm_2017*test$Depth_add_cm*10000
#
(initial-final)
(diff)
#
round((initial-final),3)==round((diff),3)#should be equal
#
sub <- test$X0.20cm_2017*
  test$Bulk.Density.g.cm3_0.20cm_2017*
  test$Depth_add_cm*10000*(1/1000)*(1/1000)
sub#should equal amount added 
round(sub,3)==round(test$gm2_2017_add,3)
#when expansion occurred 
test1 <- dat[3,]
initial <- test1$Bulk.Density.g.cm3_0.20cm_1994*20*10000
final <- test1$Bulk.Density.g.cm3_0.20cm_2017*20*10000
diff <- test1$Bulk.Density.g.cm3_20.40cm_2017*test1$Depth_add_cm*10000
#
(initial-final)
(diff)
#
round((initial-final),3)==round((diff),3)
#
add <- test1$X20.40cm_2017*
  test1$Bulk.Density.g.cm3_20.40cm_2017*
  test1$Depth_add_cm*10000*(1/1000)*(1/1000)
add
round(add,3)==round(test1$gm2_2017_add,3)
#
#
#
#
#
#