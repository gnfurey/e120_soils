#
library(tidyverse)
e120 <- read.csv("SubmissionData/Dataset_S1.csv")
colnames(e120)
dat <- e120 %>%
  select(Plot,TotalProd.2015.2017,NumSp,gm2_2017_0.20cm_Calcium_aa:
           gm2_2017_0.20cm_Potassium_aa)
colnames(dat)
#rename
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "gm2_2017_0.20cm_",replacement = "")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "aa|ea|m3|LOI|bray|icp|1.1w",replacement = "")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "_",replacement = "")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Organic.Matter",replacement = "OM")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Calcium",replacement = "Ca")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Magnesium",replacement = "Mg")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Potassium",replacement = "K")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Phosphorus",replacement = "P")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Carbon",replacement = "C")
colnames(dat) <- str_replace_all(string = colnames(dat),pattern = "Nitrogen",replacement = "N")
#get logNumSp
dat$logNumSp <- log(dat$NumSp)
colnames(dat)
#Plot=273 # this plot is an outlier, but there is not a 
#non statistical reason to remove it. It is a real value and shows that biomass 
#has not reached its maximum on these soils
ggplot(dat,aes(x=logNumSp,y=TotalProd.2015.2017))+
  geom_point()+
  geom_smooth(method="lm")
#log Y does not look much better 
ggplot(dat,aes(x=logNumSp,y=log(TotalProd.2015.2017)))+
  geom_point()+
  geom_smooth(method="lm")
#resid
residplot <- function(mod,type1){
  fit <- fitted(mod)
  resid <- residuals(mod,type=type1)
  dat <- data.frame(fit=fit,resid=resid)
  p1 <- ggplot(dat,aes(x=fit,y=resid))+geom_point(shape=21)+
    theme_bw()+
    ylab("Residuals")+
    xlab("Fitted Values")
  return(p1)
}
##OLS reg 
mod1 <- lm(TotalProd.2015.2017~logNumSp,
           data=dat)
summary(mod1)
#examine 
residplot(mod1,"pearson")
#tests 
shapiro.test(residuals(mod1))
library(nlme)
#OLS
mod.lm <- lm(TotalProd.2015.2017~
              logNumSp+
              K+
              Ca+
              Mg+
              N+
              C+
              pH+
              P,
            data=dat)
summary(mod.lm)
anova(mod.lm)
tmpmod <- step((mod.lm))
summary(tmpmod)
car::Anova(mod.lm,type="III")
car::vif(mod.lm)#CN as expected are correlated 
#gls 
mod0 <- gls(TotalProd.2015.2017~
             logNumSp+
             K+
             Ca+
             Mg+
             N+
             C+
             pH+
             P,
            method="ML",
           data=dat)
#error weights 
mod1 <- update(mod0,
            weights=varPower(form=~fitted(.)),
            data=dat)
#
summary(mod1)
anova(mod0,mod1)#increased the fit 
tmpmod2 <- MASS::stepAIC(mod1)#backwards selection
summary(tmpmod2)
#tests
shapiro.test(residuals(mod0,type="normalized"))
shapiro.test(residuals(mod1,type="normalized"))
#examine 
gridExtra::grid.arrange(residplot(mod0,"normalized"),
                        residplot(mod1,"normalized"),
                        ncol=1)
#there is a slight non-linearity in the residuals due to the relationship with some soil variables
#the effect of these nutrients on plant productivity is not linear, but it is a reasonable approx.
#I chose to leave these terms as linear to provide the effect size in g of biomass per g of nutrient
####
library(MuMIn)
library(flextable)
library(officer)
options(na.action = "na.fail")
dd <- dredge(mod1,
             rank="BIC")#BIC as we want a parsimonious model
##
dd#see all models 
####model average
mod5 <- model.avg(dd, subset = delta < 4)
summary(mod5)
####
#general comment
#please email me if you find something interesting here that I missed
#this model lines up well with the soil test's agronomic report so it has external validity
#There is probably a paper here examining the correlations among these variables
#I am not interested in pursuing that route, but it is open for anyone to run
###############
tab_sub <- as.data.frame(summary(mod5)$coefmat.subset)
tab_sub$Nutrient <- rownames(tab_sub)
tab_sub <- tab_sub %>% arrange(`Pr(>|z|)`) %>% 
  select(Nutrient,everything()) %>% 
  mutate(across(Estimate:`Pr(>|z|)`,~round(.x,3)))
#
colnames(tab_sub)[2] <- "Coefficient"
tmp <- flextable(tab_sub) %>% fontsize(size=12)
set_table_properties(tmp, width = 1, layout = "autofit")
tmp <- font(tmp,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 6", style = "Normal") %>%
  body_add_flextable(value = tmp)
print(doc, target = "Tables/SupplementalTable6_anppreg_raw.docx")
############
