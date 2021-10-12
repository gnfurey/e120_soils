#############################
library(tidyverse)
e120 <- read.csv("SubmissionData/Dataset_S1.csv")
#rename the data 
colnames(e120)
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "soil_diffgm2_0.20cm_",replacement ="" )
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_ea",replacement ="" )
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_aa",replacement ="" )
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_LOI",replacement ="" )
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_1.1w",replacement ="" )
colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_bray",replacement ="" )
#not the best solution, but to keep the breaks 
colnames(e120)[6:14] <- paste(colnames(e120)[6:14],"soil_soil",sep="_")
colnames(e120)
#long format 
dat <- e120 %>% select(Plot,NumSp,
                       Fgset,Fgset2,
                       Calcium_soil_soil:
                         Potassium_soil_soil,
                       Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm) %>% 
  pivot_longer(Calcium_soil_soil:Calcium_gm2_Root.0.30cm,
               names_to="nut",values_to="val") %>% 
  separate(nut,into = c("nut","drop","pool"),sep = "_") %>% 
  mutate(drop=NULL) %>% 
  filter(nut %in% c("Nitrogen","Potassium","Calcium","Magnesium"))
###get soil C 
carbon <- e120 %>% select(Plot,NumSp,
                          Fgset,Fgset2,
                          Calcium_soil_soil:
                            Potassium_soil_soil,
                          Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm) %>% 
  pivot_longer(Calcium_soil_soil:Calcium_gm2_Root.0.30cm,
               names_to="nut",values_to="val") %>% 
  separate(nut,into = c("nut","drop","pool"),sep = "_") %>% 
  mutate(drop=NULL) %>% 
  filter(nut %in% c("Carbon"))
###get biomass 
biomass <- e120 %>% select(Plot,NumSp,Fgset,Fgset2,
                           AbvBioAnnProd.2015.2017,
                           Root030cm.2015.2017) %>% 
  pivot_longer(AbvBioAnnProd.2015.2017:Root030cm.2015.2017,
               names_to="nut",values_to="val") %>% 
  mutate(pool = ifelse(nut=="AbvBioAnnProd.2015.2017","Abv","Root.0.30cm")) %>% 
  mutate(nut="Carbon")
###biomass and carbon together 
carb.both <- rbind(biomass,carbon)
###
dat <- rbind(dat,carb.both)
###
dat$Fgset2 <- as.factor(dat$Fgset2)
dat$Fgset2 = factor(dat$Fgset2,
                    levels(dat$Fgset2)[c(3,1,7,4,2,6,5)])
#############################
allstuff <- function(x){
  type1=x
  tabs <- function(x,type1){
    # x="Nitrogen"
    # type1="soil"
    tmp <- dat %>% 
      filter(nut==x) %>% 
      filter(pool==type1)
    library(nlme)
    library(emmeans)
    form1 <- as.formula("val~Fgset2")
    if(x=="Carbon"&type1=="Root.0.30cm"){
      form1 <- as.formula("log(val) ~ Fgset2")
    }else{
      form1 <- as.formula("val~Fgset2")
    }
    # form1
    #the FL combination only has 5 points and it fails 
    #when using a separate weight. I gave it a fixed weight at the cost
    #of increasing its degrees of freedom with the benefit of it
    #not having 0.8 estimated df.  
    mod1 <- do.call("gls", args = list(form1,
                                       weights=varIdent(form=~1|Fgset2,
                                                        fixed=c(
                                                          FL=1
                                                        )),
                                       
                                       data=tmp))
    # mod1 <- lm(total~Fgset2,
    #            data=mods)
    summary(mod1)
    anova(mod1)
    plot(mod1)
    tmp$resid <- residuals(mod1,type="normalized")
    # residplot(mod1)
    ggplot(tmp,aes(x=Fgset2,y=resid))+
      geom_boxplot()
    if(x=="Carbon"&type1=="Root.0.30cm"){
      em <- emmeans(mod1,specs=~Fgset2,
                    options = list(tran = "log"), type = "response",
                    sigmaAdjust = TRUE)
    }else{
      em <- emmeans(mod1,specs=~Fgset2,
                    sigmaAdjust = TRUE)
    }
    pairs(em)
    cld <- CLD(em,Letters = LETTERS,reversed=TRUE)
    cld
    cld <- cld %>% select(Fgset2,.group)
    #
    out <- as.data.frame(em)
    #
    out1 <- left_join(out,cld)
    #
    out1$.group <- str_trim(out1$.group)
    colnames(out1)[2] <- "emmean"
    return(out1)
  }
  nuts <- unique(dat$nut)
  nuts <- set_names(nuts)
  out_table <- map_dfr(.x = nuts,
                       .f = tabs,.id = "Nutrient",type1=type1)
  return(out_table)
}
####
types <- unique(dat$pool)
types <- set_names(types)
tabs <- map_dfr(.x = types,
                .f = allstuff,.id = "pool")
#######
#test
testdat <- dat %>% 
  filter(nut=="Calcium") %>% 
  filter(pool=="soil")
mod1 <- gls(val~Fgset2,weights=varIdent(form=~1|Fgset2,
                                   fixed=c(FL=1)),
            data=testdat)
summary(mod1)
plot(mod1)
em_test <- emmeans(mod1,specs=~Fgset2,
              sigmaAdjust = TRUE)
out <- as.data.frame(em_test)
tabs_calc <- tabs %>% filter(pool=="Soil Difference") %>% 
  filter(Nutrient=="Calcium")
all(round(out$emmean,3)==round(tabs_calc$emmean,3))
#######
mono <- tabs %>% 
  filter(Fgset2 %in% c("G","F","L","GFL"))
#
####
library(grid)
library(gridExtra)
library(egg)
####
tabs$pool <- ifelse(tabs$pool=="soil","Soil Difference",tabs$pool)
####
fgfun <- function(x,pool1,title1){
  # x="Calcium"
  tmp <- tabs %>% 
    filter(Nutrient==x) %>% 
    filter(pool==pool1)
  #
  tmp$.group <- as.character(tmp$.group)
  tmp <- tmp %>% 
    group_by(pool) %>% 
    mutate(max = max(emmean+SE)*1.05) %>% 
    ungroup()
  #
  p1 <- ggplot(tmp,aes(x=Fgset2,y=emmean,fill=Fgset2))+
    geom_hline(yintercept = 0)+
    geom_bar(col='Black',stat="identity")+
    scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", 
                                 "#009E73", "#0072B2",
                                 "#D55E00", "#CC79A7"))+
    geom_errorbar(aes(ymax=emmean+SE,ymin=emmean-SE),col="Black",
                  width=0.7)+
    theme_classic(base_size = 10)+
    # geom_text(aes(y=min,label=adjusted_star),size=7,col="Black")+
    # ylab(expression(x,"(g"%.%"m)"^-2))+
    ggtitle(title1)+
    ylab(paste(x,pool1))+
    xlab("Functional Group Set")+
    guides(fill=FALSE)+
    theme(strip.text = element_blank(),
          plot.title = element_text(face="bold",family="Helvetica",size=8,
                                          margin=margin(0,0,0,0)),
                plot.margin = margin(t = 1,
                                     r = 1,
                                     l = 1,
                                     b = 1),
          strip.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y =element_text(size=8))+
    geom_text(aes(y=max,label=.group),size=2.5)
  # p1
  # p1 <- tag_facet(p1,
  #                 open = "", close = "",
  #                 vjust = 0.4,
  #                 tag_pool = letters)
  return(p1)
} 
fgfun(x = "Potassium",pool1 = "Abv",title1 = "a")
library(grid)
library(gridExtra)
nuts <- rep(unique(tabs$Nutrient)[c(5,3,4,1,2)],each=3)
nuts
pools0 <- rep(unique(tabs$pool)[c(2,3,1)],times=5)
pools0
title0 <- letters[1:15]
plots <- pmap(list(nuts,pools0,title0),
                   .f = fgfun)
#
list(nuts,pools0,title0)
#
plots[[1]] <- plots[[1]]+ylab(expression(paste("Aboveground Biomass ","(","g"%.%"m"^-2,")"))) 
plots[[1]]
plots[[2]] <- plots[[2]]+ylab(expression(paste("Belowground Biomass ","(","g"%.%"m"^-2,")")))
plots[[3]] <- plots[[3]]+ylab(expression(paste(Delta,"Soil Carbon "["(2017-1994)"]," (","g"%.%"m"^-2,")")))
plots[[3]]
#
plots[[4]] <- plots[[4]]+ylab(expression(paste("Shoot Nitrogen ","(","g"%.%"m"^-2,")")))
plots[[5]] <- plots[[5]]+ylab(expression(paste("Root Nitrogen ","(","g"%.%"m"^-2,")")))
plots[[6]] <- plots[[6]]+ylab(expression(paste(Delta,"Soil Nitrogen "["(2017-1994)"]," (","g"%.%"m"^-2,")")))
#
plots[[7]] <- plots[[7]]+ylab(expression(paste("Shoot Potassium ","(","g"%.%"m"^-2,")")))
plots[[8]] <- plots[[8]]+ylab(expression(paste("Root Potassium ","(","g"%.%"m"^-2,")")))
plots[[9]] <- plots[[9]]+ylab(expression(paste(Delta,"Soil Potassium "["(2017-1994)"]," (","g"%.%"m"^-2,")")))
#
plots[[10]] <- plots[[10]]+ylab(expression(paste("Shoot Calcium ","(","g"%.%"m"^-2,")")))
plots[[11]] <- plots[[11]]+ylab(expression(paste("Root Calcium ","(","g"%.%"m"^-2,")")))
plots[[12]] <- plots[[12]]+ylab(expression(paste(Delta,"Soil Calcium "["(2017-1994)"]," (","g"%.%"m"^-2,")")))
#
plots[[13]] <- plots[[13]]+ylab(expression(paste("Shoot Magnesium ","(","g"%.%"m"^-2,")")))
plots[[14]] <- plots[[14]]+ylab(expression(paste("Root Magnesium ","(","g"%.%"m"^-2,")")))
plots[[15]] <- plots[[15]]+ylab(expression(paste(Delta,"Soil Magnesium "["(2017-1994)"]," (","g"%.%"m"^-2,")")))
#
ft2=12
o1 <- exec(.fn = grid.arrange,grobs=plots,ncol=3,
           bottom=textGrob("Functional Group Set",
                           gp=gpar(cex=1,fontsize=ft2, 
                                   family="Helvitica")))
o1
ggsave(filename = "Figures_PRINT/FigureS8.pdf",
       dpi=600,
       plot = o1,
       height=240,width=183,unit="mm")
####
tabs2 <- tabs %>% filter(Fgset2 %in% c("G","L","F"))
tabs3 <- tabs2 %>% 
  filter(str_detect(.group, 'A'))
NROW(tabs3)
(NROW(tabs2)-NROW(tabs3))/NROW(tabs2)
####