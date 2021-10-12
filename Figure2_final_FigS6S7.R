##
library(tidyverse)
#custom functions
source("e120soils_functions.R")
#
e120 <- read.csv("SubmissionData/Dataset_S1.csv")
#
colnames(e120)
#get data in long format
plant <- e120 %>% 
  select(Plot,NumSp,Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm) %>% 
  pivot_longer(c(Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm),
               names_to = "nut",values_to = "val.g") %>% 
  separate(nut,into=c("nut","unit","type"),sep="_") %>% 
  select(-unit) %>%
  filter(nut!="Phosphorus")
unique(plant$type)
#
#colnames(e120)
#get biomass in long format
tmp <- e120 %>% select(Plot,NumSp,
                       AbvBioAnnProd.2015.2017,
                       Root030cm.2015.2017) %>% 
  pivot_longer(c(AbvBioAnnProd.2015.2017,Root030cm.2015.2017),
               names_to = "nut",values_to = "val.g") %>% 
  filter(NumSp!=0)
tmp$type <- ifelse(tmp$nut=="AbvBioAnnProd.2015.2017","Abv",
                   "Root.0.30cm")
tmp$nut <- "Biomass"
#
colnames(plant)
colnames(tmp)
#get new data with biomass as a nutrient 
plant <- rbind(plant,tmp)
#nutrient list 
nuts <- c("Calcium","Magnesium","Potassium","Nitrogen","Biomass")
#calculate means 
plant.m <- plant %>% 
  group_by(NumSp,nut,type) %>% 
  summarise(se=my.stand(val.g),
            val.g=mean(val.g)) %>% 
  filter(nut %in% nuts) %>% 
  arrange(type,val.g)
#get monoculture means 
mono_mean <- plant.m %>% 
  select(NumSp,type,nut,val.g) %>% 
  filter(NumSp==1) %>% 
  pivot_wider(names_from="NumSp",values_from="val.g")
colnames(mono_mean)[3] <- "mono"
#join together mono means and total plant data
plant <- left_join(plant,mono_mean)
#get percent increase
plant$diff <- pchange(plant$val.g,plant$mono)
########
write.csv(x = plant,file = "SubmissionData/Dataset_S6.csv",row.names = FALSE)
AppendixS3col <- colnames(plant)
write.csv(x = AppendixS3col,file = "Data/S3_DatasetS6.columns.csv",row.names = FALSE)

########
#get means 
plant.m1 <- plant %>% 
  group_by(NumSp,nut,type) %>% 
  summarise(se=my.stand(diff),
            diff=mean(diff)) %>% 
  filter(nut %in% nuts) %>% 
  arrange(type,diff)
#
plant$nut <- as.factor(plant$nut)
plant.m1$nut <- as.factor(plant.m1$nut)
#########
library(broom)
#get new data frame
newdat <- data.frame(NumSp=seq(from=1,to=16,by=0.2))
library(nlme)
#run models 
mods <- plant %>% 
  group_by(nut,type) %>% 
  nest() %>% 
  mutate(mod=map(data,~lm(diff~log(NumSp),
                          data=.)),
         tid=map(mod,glance),
         resid = map(mod,augment),
         fits = map(mod,augment,
                    newdata=newdat,
                    se_fit=TRUE))
#extract values 
mods1 <- mods %>%   
  unnest(tid) %>% 
  mutate(p.value.adj=p.adjust(p.value,"fdr")) %>% 
  # filter(p.value.adj<0.01) %>% 
  select(nut,type,r.squared:df.residual,p.value.adj)
#get fits 
fit <- mods %>% 
  unnest(fits) %>% 
  select(-c(data,mod,tid,resid))
colnames(fit)
colnames(fit)[4] <- "diff"
colnames(fit)[5] <- "se"
######
#get resids
resid <- mods %>% 
  unnest(resid) %>% 
  select(-c(data,mod,tid,fits))
#
residplots <- function(x){
  tmpdat <- resid %>% 
    filter(nut==x)
  p1 <- ggplot(tmpdat,aes(x=.fitted,y=.std.resid))+
    geom_point(shape=21)+
    facet_wrap(~type,scales="free")+
    ggtitle(x)+
    theme_bw()
  return(p1)
}
nuts <- unique(fit$nut)
residout <- map(nuts,residplots)
# residout
#
#some variables could be transformed, but the effects are
#extremely robust to a variety of transformations
#future analysis could test a variety of functional forms 
#to find the 'best' fit 
#function to test alternative models
altmods <- function(nut1,
                    x,y,
                    j,wt,form){
  # nut1="Potassium"
  # x="log(NumSp)"
  # y="diff"
  # wt=varExp(form=~fitted(.))
  # wt=varIdent(form=~1|as.factor(NumSp))
  # j="Abv"
  form1 <- as.formula(paste(y,"~",x,sep=""))
  
  tmpplant <- plant %>% 
    filter(nut==nut1) %>% 
    filter(type==j)
  mod1 <- lm(form1,data=tmpplant)
  mod2 <- do.call("gls", args = list(form1,
                                     weights=wt,
                                     data=tmpplant))
  print("lm")
  print(summary(mod1)$coef)
  print("gls")
  print(summary(mod2)$tTable)
  a <- summary(mod1)$coef[2,1]
  b <- summary(mod2)$tTable[2,1]
  out <- ((a-b)/((a+b)/2))*100
  print(paste("% Difference in effect size",round(out,3)))
  tmpplant$fit_lm <- fitted(mod1)
  tmpplant$fit_gls <- fitted(mod2)
  tmpplant$resid_lm <- residuals(mod1)
  tmpplant$resid_gls <- residuals(mod2,type="normalized")
  p1 <- ggplot(tmpplant,aes(x=fit_lm,y=resid_lm))+
    geom_point(shape=21)
  p2 <- ggplot(tmpplant,aes(x=fit_gls,y=resid_gls))+
    geom_point(shape=21)
  gridExtra::grid.arrange(p1,p2,ncol=2)
}
unique(plant$type)
#effects are robust 
# altmods(nut1 = "Potassium",x="log(NumSp)",y="diff",j = "Abv",wt =varExp(form=~fitted(.)))
# altmods(nut1 = "Potassium",x="log(NumSp)",y="diff",j = "Abv",wt=varIdent(form=~1|as.factor(NumSp)))
# altmods(nut1 = "Potassium",x="log(NumSp)",y="diff",j = "Root.0.30cm",wt =varExp(form=~fitted(.)))
# altmods(nut1 = "Potassium",x="log(NumSp)",y="diff",j = "Root.0.30cm",wt=varIdent(form=~1|as.factor(NumSp)))
# ######
levels(plant$nut)
plant$nut <- factor(plant$nut,levels(plant$nut)[c(5,1,3,2,4)])
plant.m1$nut <- factor(plant.m1$nut,levels(plant.m1$nut)[c(5,1,3,2,4)])
######
#get max and min
plant.m1 <- plant.m1 %>% 
  group_by(nut) %>% 
  mutate(max=max(diff)*1.05)
#tmp data for max and min
tmp2 <- plant.m1 %>% ungroup() %>% 
  select(nut,max) %>% distinct()
#
pfun <- function(x){
  # x="Potassium"
  #filter data sets 
  tmpdat <- plant.m1 %>% 
    filter(nut == x)
  fitdat <- fit %>% 
    filter(nut ==x )
  tmp2dat <- tmp2 %>% 
    filter(nut==x)
  #
  unique(tmpdat$type)
  #create a data set to jitter the monos
  jittermonmo.abv <- tmpdat %>% 
    filter(NumSp==1) %>% 
    filter(type=="Abv")
  jittermonmo.root <- tmpdat %>% 
    filter(NumSp==1) %>% 
    filter(type=="Root.0.30cm")
  tmpdat <- tmpdat %>% 
    filter(NumSp!=1)
  #generate plot 
  b1 <- ggplot(tmpdat,aes(x=log(NumSp),
                          y=diff,
                          fill=type,
                          col=type,
                          shape=type))+
    geom_errorbar(aes(ymax=diff+se,ymin=diff-se),
                  width=0.05,col="Black",
                  position = position_nudge(x = -0.1),
                  data=jittermonmo.abv)+
    geom_errorbar(aes(ymax=diff+se,ymin=diff-se),
                  width=0.05,col="Black",
                  position = position_nudge(x = 0.1),
                  data=jittermonmo.root)+
    geom_line(aes(y=diff),
              data=fitdat)+
    geom_line(aes(y=diff+se),
              data=fitdat,alpha=0.3)+
    geom_line(aes(y=diff-se),
              data=fitdat,alpha=0.3)+
    geom_errorbar(aes(ymax=diff+se,ymin=diff-se),
                  width=0.05,col="Black")+
    #####
  geom_point(size=3,col="Black")+
    geom_point(size=3,
               col="Black",
               aes(x=log(NumSp),
                   y=diff,
                   fill=type,
                   shape=type),
               position = position_nudge(x = -0.1),
               data=jittermonmo.abv)+
    geom_point(size=3,
               col="Black",
               aes(x=log(NumSp),
                   y=diff,
                   fill=type,
                   shape=type),
               position = position_nudge(x = 0.1),
               data=jittermonmo.root)+
    #####
  scale_shape_manual(values=c(21,23),name="",
                     labels=c("Shoots","Roots"))+
    scale_fill_brewer(name="",palette ="Set1",direction=-1,
                      labels=c("Shoots","Roots"))+
    scale_color_brewer(palette = "Set1",direction=-1)+
    guides(color="none")+
    theme_classic(base_size = 12)+
    scale_x_continuous(breaks=c(log(1),log(2),log(4),log(8),
                                log(16)),
                       labels=c("1","2",
                                "4","8",
                                "16"))+
    ylab("% Change Relative to Monoculture Mean")+
    theme(legend.position = "bottom",
          legend.key.height =  unit(0.01, "cm"),
          legend.key.width = unit(0.4,"cm"),
          legend.text = element_text(size=12),
          legend.box.spacing = unit(0,"cm"),
          legend.box.margin = margin(0,0,0,0,"cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          plot.margin = unit(c(0.2,0.14,0.1,0.1),"cm"))+
    xlab(expression(paste(
      "log"[e],"(Number of Plant Species)")))
  b1
  return(b1)
}
library(gridExtra)
library(grid)
#make each plot
p1 <- pfun("Nitrogen")+
  ggtitle("a")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(1,0,0,0)),
        plot.margin = margin(0,1.1,0,0))
p1
#possible source and credit 
#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
legendfunc<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend.1<-legendfunc(p1)
ft=12
leftjust <- 0.3
p1 <- p1+
  annotation_custom(textGrob(
    x = unit(leftjust, "npc"), 
    y = unit(0.9, "npc"),
    label="Nitrogen",
    gp=gpar(fontsize=ft,face="bold")))
p1
p1 <- p1+
  scale_y_continuous(breaks=seq(from=-20,160,by=20),
                     limits=c(-20,160),
                     labels=every_nth(seq(from=-20,160,by=20),2,inverse = FALSE))
  

#
p2 <- pfun("Potassium")+
  ggtitle("b")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(1,0,0,0)),
        plot.margin = margin(0,1.1,0,0))+
  annotation_custom(textGrob(
    x = unit(leftjust+0.1, "npc"), 
    y = unit(0.9, "npc"),
    label="Potassium",
    gp=gpar(fontsize=ft,face="bold")))
p2 <- p2+
  scale_y_continuous(breaks=seq(from=-50,400,by=50),
                     limits=c(-40,410),
                     labels=every_nth(seq(from=-50,400,by=50),2,inverse = FALSE))

p3 <- pfun("Calcium")+
  ggtitle("c")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(1,0,0,0)),
        plot.margin = margin(0,1.1,0,0))+
  annotation_custom(textGrob(
    x = unit(leftjust, "npc"), 
    y = unit(0.9, "npc"),
    label="Calcium",
    gp=gpar(fontsize=ft,face="bold")))
p3
p3 <- p3+
  scale_y_continuous(breaks=seq(from=-20,200,by=20),
                     limits=c(-30,220),
                     labels=every_nth(seq(from=-20,200,by=20),2,inverse = FALSE))

p4 <- pfun("Magnesium")+
  ggtitle("d")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(1,0,0,0)),
        plot.margin = margin(0,1.1,0,0))+
  annotation_custom(textGrob(
    x = unit(leftjust+0.12, "npc"), 
    y = unit(0.9, "npc"),
    label="Magnesium",
    gp=gpar(fontsize=ft,face="bold")))
p4 <- p4+
  scale_y_continuous(breaks=seq(from=-20,160,by=20),
                     limits=c(-20,170),
                     labels=every_nth(seq(from=-20,160,by=20),2,inverse = FALSE))

p5 <- pfun("Biomass")
p5 <- p5+
  theme_classic(base_size = 12)+
  ggtitle("e")+
  theme(legend.position = "bottom",
        legend.key.height =  unit(0.01, "cm"),
        legend.key.width = unit(0.4,"cm"),
        legend.text = element_text(size=12),
        legend.box.spacing = unit(0,"cm"),
        legend.box.margin = margin(0,0,0,0,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(1,0,0,0)),
        plot.margin = margin(0,2,0,0))+
  annotation_custom(textGrob(
    x = unit(0.35, "npc"), 
    y = unit(0.9, "npc"),
    label="Total\nBiomass",
    gp=gpar(fontsize=ft,face="bold")))

p5
p1 <- p1+guides(fill="none",shape="none")
p2 <- p2+guides(fill="none",shape="none")
p3 <- p3+guides(fill="none",shape="none")
p4 <- p4+guides(fill="none",shape="none")
p5 <- p5+guides(fill="none",shape="none")
#
ggsave(filename = "Figures_PRINT/Figure2.pdf",
       dpi=600,
       # device = cairo_ps,
       plot = grid.arrange(
         mylegend.1,arrangeGrob(p1,p2,p3,p4,ncol=4),
         heights=list(unit(4.5,"mm"),unit(75.5,"mm")),
         padding=unit(2,"mm"),
         left=textGrob("Relative Nutrient Content (%)",gp=gpar(cex=1,fontsize=12),rot=90),
         bottom=textGrob("Number of Plant Species (Log Scale)",
                         gp=gpar(cex=1,fontsize=12, family="Helvitica")))
       ,
       height = 80,width=178,unit="mm")
check <- grid.arrange(
  mylegend.1,arrangeGrob(p1,p2,p3,p4,ncol=4),
  heights=list(unit(4.5,"mm"),unit(75.5,"mm")),
  padding=unit(2,"mm"),
  left=textGrob("Relative Nutrient Content (%)",gp=gpar(cex=1,fontsize=12),rot=90),
  bottom=textGrob("Number of Plant Species (Log Scale)",
                  gp=gpar(cex=1,fontsize=12, family="Helvitica")))
saveRDS(check,"Figures_PRINT/Figure2.rds")
######
#get effect sizes for results 
tmp <- plant.m1 %>% 
  filter(type=="Abv") %>% 
  filter(nut %in% c("Calcium","Magnesium","Potassium")) %>% 
  filter(NumSp==16)
round(range(tmp$diff),-1)
tmp <- plant.m1 %>% 
  filter(type=="Root.0.30cm") %>% 
  filter(nut %in% c("Calcium","Magnesium","Potassium")) %>% 
  filter(NumSp==16)
round(range(tmp$diff),-1)
######
#make tables
library(flextable)
library(officer)
library(tidyverse)
colnames(mods1)
table <- mods1 %>% select(nut,type,statistic,p.value.adj,r.squared) %>% 
  arrange(type,nut,r.squared)
table$r.squared <- round(table$r.squared,2)
table$statistic <- round(table$statistic,3)
table$p.value <- round(table$p.value.adj,5)
table$p.value.adj <- NULL
#
colnames(table)
colnames(table) <- c("Nutrient","Shoots|Roots","F-Statistic",
                     "R^2","P-value")
#
colnames(table)
table$`Shoots|Roots` <- ifelse(table$`Shoots|Roots`=="Abv",
                               "Shoots","Roots")
table$`P-value` <- ifelse(table$`P-value`<0.0001,"<0.0001")

#
tmp <- flextable(table) %>% fontsize(size=12)
set_table_properties(tmp, width = 1, layout = "autofit")
tmp <- font(tmp,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 4", style = "Normal") %>%
  body_add_flextable(value = tmp)
print(doc, target = "Tables/SupplementalTable4_raw.docx")
#########
#these plots show the raw data 
#there is a nested function to plot each year
bothplot <- function(x){
  type1=x
  pfun_raw <- function(x,type1){
    # plant$type
    # x="Calcium";type1="Abv"
    tmp <- plant %>% 
      filter(nut==x) %>% 
      filter(type==type1)
    if(type1=="Abv"){
      col1="#377EB8";lab="Shoots";shp=21
    }else{col1="#E41A1C";lab="Roots";shp=23}
    b1 <- ggplot(tmp,aes(x=log(NumSp),
                         y=val.g,
                         fill=type,
                         col=type,
                         shape=type))+
      geom_point(size=3,col="Black",shape=shp)+
      scale_x_continuous(breaks=c(log(1),log(2),log(4),log(8),
                                  log(16)),
                         labels=c("1","2",
                                  "4","8",
                                  "16"))+
      scale_fill_manual(values=col1,name="",
                        labels=lab)+
      guides(color="none")+
      facet_wrap(~nut)+
      theme_classic(base_size = 12)+
      geom_smooth(method="lm",formula = "y~x",col="Black",fill="Black")+
      ylab("Nutrient Pool (gm2)")+
      theme(legend.position = "bottom",
            legend.key.height =  unit(0.01, "cm"),
            legend.key.width = unit(0.4,"cm"),
            legend.text = element_text(size=12),
            legend.box.spacing = unit(0,"cm"),
            legend.box.margin = margin(0,0,0,0,"cm"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            # strip.background = element_blank(),
            # strip.text.x = element_blank(),
            plot.margin = unit(c(0.2,0.14,0.1,0.1),"cm"))
    b1
    return(b1)
  }
  library(gridExtra)
  library(grid)
  library(grid)
  b1 <- pfun_raw("Nitrogen",type1)+
    ggtitle("a")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(1,0,0,0)),
          plot.margin = margin(0,1.1,0,0))
  b1
  legendfunc<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  mylegend.2<-legendfunc(b1)
  ft=12
  leftjust <- 0.3
  # b1 <- b1+
  #   annotation_custom(textGrob(
  #     x = unit(leftjust, "npc"), 
  #     y = unit(0.9, "npc"),
  #     label="Nitrogen",
  #     gp=gpar(fontsize=ft,face="bold")))
  # b1
  #
  b2 <- pfun_raw("Potassium",type1)+
    ggtitle("b")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(1,0,0,0)),
          plot.margin = margin(0,1.1,0,0))
  # annotation_custom(textGrob(
  #   x = unit(leftjust+0.1, "npc"), 
  #   y = unit(0.9, "npc"),
  #   label="Potassium",
  #   gp=gpar(fontsize=ft,face="bold")))
  b2
  b3 <- pfun_raw("Calcium",type1)+
    ggtitle("c")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(1,0,0,0)),
          plot.margin = margin(0,1.1,0,0))
  # annotation_custom(textGrob(
  #   x = unit(leftjust, "npc"), 
  #   y = unit(0.9, "npc"),
  #   label="Calcium",
  #   gp=gpar(fontsize=ft,face="bold")))
  b3
  b4 <- pfun_raw("Magnesium",type1)+
    ggtitle("d")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(1,0,0,0)),
          plot.margin = margin(0,1.1,0,0))
  # annotation_custom(textGrob(
  #   x = unit(leftjust+0.12, "npc"), 
  #   y = unit(0.9, "npc"),
  #   label="Magnesium",
  #   gp=gpar(fontsize=ft,face="bold")))
  b5 <- pfun_raw("Biomass",type1)
  b5 <- b5+
    theme_classic(base_size = 12)+
    ggtitle("e")+
    theme(legend.position = "bottom",
          legend.key.height =  unit(0.01, "cm"),
          legend.key.width = unit(0.4,"cm"),
          legend.text = element_text(size=12),
          legend.box.spacing = unit(0,"cm"),
          legend.box.margin = margin(0,0,0,0,"cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          # strip.background = element_blank(),
          # strip.text.x = element_blank(),
          plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(1,0,0,0)),
          plot.margin = margin(0,2,0,0))
  # annotation_custom(textGrob(
  #   x = unit(0.35, "npc"), 
  #   y = unit(0.8, "npc"),
  #   label="Total\nBiomass",
  #   gp=gpar(fontsize=ft,face="bold")))
  b5
  b1 <- b1+guides(fill="none",shape="none")
  b2 <- b2+guides(fill="none",shape="none")
  b3 <- b3+guides(fill="none",shape="none")
  b4 <- b4+guides(fill="none",shape="none")
  b5 <- b5+guides(fill="none",shape="none")
  #
  if(type1=="Abv"){
    fig=6
  }else{fig=7}
  name <- paste("Figures_PRINT/FigureS",fig,".pdf",sep="")
  
  ggsave(filename = name,
         dpi=600,
         plot = grid.arrange(
           mylegend.2,arrangeGrob(b1,b2,b3,b4,b5,ncol=2),
           heights=list(unit(4.5,"mm"),unit(75.5,"mm")),
           padding=unit(2,"mm"),
           left=textGrob(expression(paste("Biomass Pool ","(","g"%.%"m"^-2,")")),gp=gpar(cex=1,fontsize=12),rot=90),
           bottom=textGrob("Number of Plant Species (Log Scale)",
                           gp=gpar(cex=1,fontsize=12, family="Helvitica")))
         ,
         height = 160,width=120,unit="mm")
}
unique(plant$type)
bothplot("Abv")
bothplot("Root.0.30cm")

