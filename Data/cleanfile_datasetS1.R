library(tidyverse)
#make file for all data 
pos <- read.csv("Data/data_e120soils_baseset.csv")
#SI_Appendix2.Rmd
diffnuts <- read.csv("Data/data_e120soils_bulknutrients.csv")
unique(diffnuts$type)
####wide format
diffnuts.w <- diffnuts %>% spread(type,val)
unique(diffnuts.w$nut)
#join in data with base set 
dat1 <- left_join(pos,diffnuts.w)
colnames(dat1)
#spread to wide format 
dat.t1 <- dat1 %>% pivot_wider(names_from = c(nut,extract),
                               values_from = c(
                                 soil_rawdiffgm2_0.20cm,
                                 gm2_1994_0.20cm,
                                 gm2_2017_0.20cm))
#reorganize 
colnames(dat.t1)
dat.t1 <- dat.t1 %>% select(Plot:Species,Fgset,Fgset2,soil_rawdiffgm2_0.20cm_Calcium_aa:gm2_2017_0.20cm_Potassium_aa)
colnames(dat.t1)
colnames(dat.t1) <- str_replace_all(string = colnames(dat.t1),
                                    pattern = "raw",
                                    replacement = "")
colnames(dat.t1)
#get ANPP data 
anpp <- read.csv("Data/data_e120soils_ANPP_root.csv")
colnames(anpp)
#
shortdat <- left_join(dat.t1,anpp)
colnames(shortdat)
######
#get plant gm2 data 
plant <- read.csv("SubmissionData/Dataset_S5.csv")
colnames(plant)
plant <- plant %>% select(Plot,Nutrient,Type,val.g)
#wide format 
plant.w <- plant %>% 
  pivot_wider(names_from = c("Nutrient","Type"),
              values_from="val.g")
##rename
colnames(plant.w)
colnames(plant.w) <- str_replace_all(string = colnames(plant.w),
                                     pattern = "p.Abv",
                                     replacement = "gm2_Abv")
colnames(plant.w) <- str_replace_all(string = colnames(plant.w),
                                     pattern = "p.Root",
                                     replacement = "gm2_Root.")
colnames(plant.w)
shortdat <- left_join(shortdat,plant.w)
##############
#check with anpp data
biomass <- read.csv("SubmissionData/Dataset_S5.csv")
biomass <- biomass %>% select(Plot,Type,Biomass) %>% 
  distinct() %>% 
  pivot_wider(names_from = "Type",values_from="Biomass")
colnames(anpp)
colnames(biomass) <- str_replace_all(string = colnames(biomass),
                                     pattern = "p.Abv",
                                     replacement = "AbvBioAnnProd.2015.2017")
colnames(biomass) <- str_replace_all(string = colnames(biomass),
                                     pattern = "p.Root",
                                     replacement = "Root030cm.2015.2017")
#no need to upload redundant data. ANPP and Rootbiomass are available in Dataset_S5
all(biomass$AbvBioAnnProd.2015.2017==shortdat$AbvBioAnnProd.2015.2017)
all(biomass$Root030cm.2015.20170.30cm==shortdat$Root030cm.2015.2017)
##############
e120 <- shortdat
##
#colnames(e120)
#get plant nutrient totals 
plant.tot <- e120 %>% 
  select(Plot,NumSp,Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm) %>% 
  pivot_longer(c(Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm),
               names_to = "nut",values_to = "val.g") %>% 
  separate(nut,into=c("nut","unit","type"),sep="_") %>% 
  select(-unit) 
#check data 
table(plant.tot$nut)#154*2 for each nut
#get total amounts per plot of both roots and shoots
plant.tot <- plant.tot %>% 
  group_by(Plot,nut) %>% 
  summarise(tot=sum(val.g))
table(plant.tot$nut)
####
#test that nutrients total are right
test1 <- plant.tot %>% filter(Plot==2)
test1_dat <- e120 %>% filter(Plot==2)
(test1_dat$Nitrogen_gm2_Abv+test1_dat$Nitrogen_gm2_Root.0.30cm)==
test1[test1$nut=="Nitrogen","tot"]
####
#colnames(e120)
#get long format for soils 
soil <- e120 %>% select(Plot,NumSp,
                        Fgset,Fgset2,
                        soil_diffgm2_0.20cm_Calcium_aa:
                          soil_diffgm2_0.20cm_Potassium_aa) %>% 
  pivot_longer(soil_diffgm2_0.20cm_Calcium_aa:soil_diffgm2_0.20cm_Potassium_aa,
               names_to="nut",values_to="val") %>% 
  separate(nut,into = c("type","Year","Depth","nut","extract"),sep = "_")
#get nutrients for figures 
nuts <- c("Nitrogen","Calcium","Magnesium","Potassium","Phosphorus")
#filter only macronutrients 
soil <- soil %>% 
  filter(nut %in% c(nuts))
#simplify data 
soil <- soil %>% 
  select(Plot,nut,val)
#############
#join soils and plants 
soil1 <- left_join(soil,plant.tot)
colnames(soil1)
#rename 
colnames(soil1)[3] <- "soil"
colnames(soil1)[4] <- "biomass"
#get total pool change
soil1$total <- soil1$soil+soil1$biomass
soil1 <- soil1 %>% select(Plot,nut,total)
#############
#test that nutrients total are right
soiltest1 <- soil1 %>% filter(Plot==2)
soiltest1_dat <- e120 %>% filter(Plot==2)
(soiltest1_dat$Nitrogen_gm2_Abv+soiltest1_dat$Nitrogen_gm2_Root.0.30cm+
  soiltest1_dat$soil_diffgm2_0.20cm_Nitrogen_ea)==
soiltest1[soiltest1$nut=="Nitrogen","total"]
##############
soil1.w <- soil1 %>% spread(nut,total)
##
colnames(soil1.w)[2:6] <- paste(colnames(soil1.w)[2:6],"_total.soil.plant.g.m2",sep="")
colnames(soil1.w)
#
e120 <- left_join(e120,soil1.w)
##
write.csv(x = e120,file = "SubmissionData/Dataset_S1.csv",row.names = FALSE)
##test
tots <- e120$soil_diffgm2_0.20cm_Calcium_aa+e120$Calcium_gm2_Abv+e120$Calcium_gm2_Root.0.30cm
tots <- round(tots,2)
tots==round(e120$Calcium_total.soil.plant.g.m2,2)

diffs <- e120$gm2_2017_0.20cm_Calcium_aa-e120$gm2_1994_0.20cm_Calcium_aa
diffs <- round(diffs,2)
diffs==round(e120$soil_diffgm2_0.20cm_Calcium_aa,2)
#within rounding error
tmp <- data.frame(diffs=diffs,soildiff=e120$soil_diffgm2_0.20cm_Calcium_aa)
#
#
