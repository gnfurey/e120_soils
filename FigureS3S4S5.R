#
source("e120soils_functions.R")
library(tidyverse)
library(broom)
#get mgkg soil data
dat <- read.csv("SubmissionData/Dataset_S3.csv")
AppendixS3col <- colnames(dat)
write.csv(x = AppendixS3col,file = "Data/S3_DatasetS3.columns.csv",row.names = FALSE)

#get only surface depths 
dat <- dat %>% 
  filter(Depth =="0.20cm") %>% 
  filter(nut !="Organic.Matter")
str(dat)
dat$Year <- as.factor(dat$Year)
table(dat$nut)
colnames(dat)
#it's more useful to present these data as % as is more common in the literature 
dat$val <- ifelse(dat$nut %in% c(
                      "Carbon",
                      "Nitrogen"),
       dat$val/10000,dat$val)
#
#run a series of regressions 
unique(dat$nut)
mods <- dat %>%
  group_by(nut,Year) %>%
  filter(nut %in% c("Calcium",
                    "Magnesium",
                    "Potassium",
                    "Nitrogen",
                    "CEC",
                    "pH",
             "Carbon",
             "Phosphorus")) %>% 
  nest() %>%
  mutate(mod=map(data,~lm(val~log(NumSp),data=.)),
         mods.t = map(mod,glance),
         fits = map(mod,augment))
#extract data
mods.t <- mods %>% 
  unnest(mods.t) %>%
  arrange(desc(r.squared)) %>%
  mutate(p.value.adj=p.adjust(p.value,"fdr")) %>% 
  select(-c(data,mod))
#get fitted values
fits <- mods %>% 
  unnest(fits) %>% 
  select(-c(data,mod,mods.t))
#generate residual plots
residplots <- function(x){
  # x="Potassium"
  tmpdat <- fits %>% 
    filter(nut==x)
  p1 <- ggplot(tmpdat,aes(x=.fitted,y=.std.resid))+
    geom_point(shape=21)+
    facet_wrap(~Year,scales="free")+
    ggtitle(x)+
    theme_bw()
  p1
  return(p1)
}
#residual by fitted plots 
nuts <- unique(fits$nut)
residout <- map(nuts,residplots)
# residout#residuals show no major signs of concern
#check that stats match up
tmp1 <- mods.t %>% filter(nut=="Potassium") 
mod <- lm(val~log(NumSp),data=dat[
  dat$nut=="Potassium"&dat$Year==2017,
])
summary(mod)$r.squared==tmp1$r.squared[1]#must be TRUE 
###########
#some variables could be transformed, but the effects are
#extremely robust to a variety of transformations
#future analysis could test a variety of functional forms 
#to find the 'best' fit 
#function to test alternative models
##################
altmods <- function(nut1,
                    x,y,
                    j,wt,form,plot1){
  library(nlme)
  # nut1="Potassium"
  # x="log(NumSp)"
  # y="val"
  # # y="log(val)"
  # wt=varPower(form=~fitted(.))
  # wt=varIdent(form=~1|as.factor(NumSp))
  # j="2017"
  form1 <- as.formula(paste(y,"~",x,sep=""))
  form1
  tmpdat <- dat %>%
    filter(nut==nut1) %>%
    filter(Year==j)
  mod1 <- lm(form1,data=tmpdat)
  #boxcox
  # box1 <- MASS::boxcox(mod1)
  #
  mod2 <- do.call("gls", args = list(form1,
                                     weights=wt,
                                     data=tmpdat))
  print(paste("Nutrient = ",nut1))
  print("lm")
  print(summary(mod1)$coef)
  print("gls")
  print(summary(mod2)$tTable)
  a <- summary(mod1)$coef[2,1]
  b <- summary(mod2)$tTable[2,1]
  out <- ((a-b)/((a+b)/2))*100
  print(paste("% Difference in effect size",round(out,3)))
  tmpdat$fit_lm <- fitted(mod1)
  tmpdat$fit_gls <- fitted(mod2)
  tmpdat$resid_lm <- residuals(mod1)
  tmpdat$resid_gls <- residuals(mod2,type="normalized")
  p1 <- ggplot(tmpdat,aes(x=fit_lm,y=resid_lm))+
    geom_point(shape=21)
  p2 <- ggplot(tmpdat,aes(x=fit_gls,y=resid_gls))+
    geom_point(shape=21)
  p1
  p2
  if(plot1==TRUE){
    gridExtra::grid.arrange(p1,p2,ncol=2)
  }
}
#the p-value of the effect of log(NumSp) is extremely robust
#across various y transformations and error weights
#we chose to present the untransformed y-variable
# altmods(nut1 = "Potassium",plot1=TRUE,
#         x="log(NumSp)",y="val",j = "2017",wt =varExp(form=~fitted(.)))
# altmods(nut1 = "Potassium",plot1=TRUE,
#         x="log(NumSp)",y="val",j = "2017",wt=varIdent(form=~1|as.factor(NumSp)))
# altmods(nut1 = "Potassium",plot1=TRUE,
#         x="log(NumSp)",y="log(val)",j = "2017",wt =varExp(form=~fitted(.)))
# altmods(nut1 = "Potassium",plot1=TRUE,
#         x="log(NumSp)",y="log(val)",j = "2017",wt=varIdent(form=~1|as.factor(NumSp)))
#normal y different variance structure
# map(.x = nuts,x="log(NumSp)",.f = altmods,
#     plot1=FALSE,
#     y="val",j = "2017",wt=varIdent(form=~1|as.factor(NumSp)))
# #log transformed Y
# map(.x = nuts,x="log(NumSp)",
#     plot1=FALSE,
#     .f = altmods,y="log(val)",j = "2017",wt =varExp(form=~fitted(.)))
###########
#calculate mean and se 
dat.m <- dat %>% group_by(Year,NumSp,nut) %>%
  summarise(se=my.stand(val,na.rm = TRUE),
            val=mean(val,na.rm = TRUE))
#merge data together 
dat.m <- left_join(dat.m,mods.t)
############
#calculate percent increase 
pdat <- dat.m %>% 
  filter(NumSp %in% c(1,16)) %>% 
  select(nut,Year,NumSp,val) %>% 
  spread(NumSp,val)
#use function to get percent difference 
pdat$diff <- pchange(pdat$`16`,pdat$`1`)
pdat_2017 <- pdat %>% filter(Year==2017)
pdat_2017#
###############
#a function to generate plots for each nutrient 
plot_point <- function(x){
  # x="C"
  #new data frame for fits 
  newdat <- data.frame(NumSp=seq(from=1,to=16,by=0.2))
  #a new data set for each nutrient
  dat1 <- dat %>% filter(nut==x)
  #get fits 
  fit <- dat1 %>% 
    group_by(nut,Year) %>%
    nest() %>%
    mutate(mod=map(data,~lm(val~log(NumSp),data=.)),
           mods.t = map(mod,augment,
                        newdata=newdat,
                        se_fit=TRUE)) %>% 
    unnest(mods.t) %>% 
    select(-c(data,mod))
  #rename fits and se 
  colnames(fit)[4] <- "val"
  colnames(fit)[5] <- "se"
  #get max and min 
  tmp <- dat.m %>% filter(nut==x) 
  minval <- min(tmp$val)-max(tmp$se)*1.5
  maxval <- max(tmp$val)+max(tmp$se)*0.15
  #generate plot 
  p1 <- ggplot(tmp,aes(x=log(NumSp),
                       y=val,
                       fill=Year,
                       col=Year,
                       shape=Year))+
    geom_errorbar(aes(ymax=val+se,ymin=val-se),
                  width=0.1)+
    theme_bw(base_size = 12,base_family = "Helvetica")+
    xlab("Number of Species")+
    ylab("")+
    # scale_y_continuous(expand = expansion(0.04,0.04))+
    scale_x_continuous(expand = expansion(mult=c(0.05,0.05)),
                       breaks=c(log(1),log(2),log(4),log(8),log(16)),
                       labels = c("1","2","4","8","16"))+
    scale_fill_brewer(name="Depth_Year",palette = "Dark2")+
    scale_color_brewer(name="Depth_Year",palette = "Dark2")+
    scale_shape_manual(values=c(23,21),name="Depth_Year",
                       labels=c("1994","2017"))+
    guides(color="none")+
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.key.height =  unit(0.02, "cm"),
          legend.key.width = unit(0.4,"cm"),
          legend.text = element_text(size=15),
          legend.box.spacing = unit(0,"cm"),
          legend.box.margin = margin(0,0,0,0,"cm"),
          axis.title.x = element_blank(),
          # axis.text.x = element_text(),
          panel.grid = element_blank(),
          # title = element_blank(),
          plot.margin=margin(t = 0,r = 1.4,l = 0.8,b = 0))+
    geom_point(size=3,col="Black")
  p1
  #P did not depend on biodiversity so we do not show those lines 
  if(x!="Phosphorus"){
    p1 <- p1+
      geom_line(aes(x=log(NumSp),y=val),data=fit)+
      geom_line(aes(x=log(NumSp),y=val+se),alpha=0.3,data=fit)+
      geom_line(aes(x=log(NumSp),y=val-se),alpha=0.3,data=fit)
    p1}
  return(p1)
}
#
p1 <- plot_point("Carbon")
p1 <- p1+scale_y_continuous(breaks=seq(from=0.4,0.8,by=0.025),
                            limits=c(0.42,0.78),
                            labels=every_nth(seq(from=0.4,0.8,by=0.025),2,inverse = TRUE))
p1
# p1
#possible source and credit 
#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
legendfunc<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
#
mylegend.1<-legendfunc(p1)

# p2
p2 <- plot_point("Nitrogen")
p2 <- p2+
  scale_y_continuous(breaks=seq(from=0.04,0.065,by=0.0025),
                     limits=c(0.04,0.065),
                     labels=every_nth(seq(from=0.04,0.065,by=0.0025),2,inverse = TRUE))
  
p2
p3 <- plot_point("Potassium")
p3 <- p3+  scale_y_continuous(breaks=seq(from=25,55,by=5),
                        limits=c(24,55),
                        labels=every_nth(seq(from=25,55,by=5),2,inverse = TRUE))

p3
p4 <- plot_point("Calcium")
p4 <- p4+
  scale_y_continuous(breaks=seq(from=350,600,by=25),
                        limits=c(340,600),
                        labels=every_nth(seq(from=350,600,by=25),2,inverse = TRUE))

# p5
p5 <- plot_point("Magnesium")
p5 <- p5+
  scale_y_continuous(breaks=seq(from=50,90,by=5),
                     limits=c(50,94),
                     labels=every_nth(seq(from=50,90,by=5),2,inverse = TRUE))

p5
p6 <- plot_point("CEC")+scale_y_continuous(breaks=seq(from=2.2,3.9,by=0.1),
                                            limits=c(2.2,3.9),
                                            labels=every_nth(seq(from=2.2,3.9,by=0.1),4,inverse = TRUE))

p6
p7 <- plot_point("pH")+scale_y_continuous(breaks=seq(from=5.3,6.3,by=0.1),
                                           limits=c(5.3,6.3),
                                           labels=every_nth(seq(from=5.3,6.3,by=0.1),2,inverse = TRUE))

p7
# p8
# p8 <- plot_point("Phosphorus")
# p8
#
#a separate plot is needed for P 
phos <- plot_point("Phosphorus")+
  ggtitle("a")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        
        plot.margin = margin(2,2,2,2))+
  ylab(expression(paste("Soil Bray Phosphorus ","(mg kg"^-1,")")))
phos <- phos+
  scale_y_continuous(breaks=seq(from=30,50,by=2.5),
                      limits=c(29,53),
                      labels=every_nth(seq(from=30,50,by=2.5),2,inverse = TRUE))
# ggsave(plot=phos,
#        filename = "Figures/Figure1_P.pdf",
#        height=100,width=89,unit="mm")
#convert to Rdata file for export in another plotting script 
# saveRDS(phos, file = "P_numsp.rds")
###
ft <- 12
leftjust <- 0.2
labelfun <- function(plot0,title1,top,right,left,bottom,ylab1){
  plot0 <- plot0+theme(legend.position = "none")+
    ggtitle(title1)+
    #the code below can add text into the plot if needed
    # annotation_custom(textGrob(
    #   x = unit(leftjust, "npc"), 
    #   y = unit(0.9, "npc"),
    #   label="C",
    #   gp=gpar(fontsize=12,face="bold")))+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(t = top,
                               r = right,
                               l = left,
                               b = bottom))+
    ylab(ylab1)
  return(plot0)
}
p1 <- labelfun(plot0 = p1,
               title1 = "a",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("Carbon ","(%)")))
# 
p2 <- labelfun(plot0 = p2,
               title1 = "b",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("Nitrogen ","(%)")))
# p3
p3 <- labelfun(plot0 = p3,
               title1 = "c",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("Potassium ","(mg kg"^-1,")")))
p3
p4 <- labelfun(plot0 = p4,
               title1 = "d",
               top = 0,
               right = 1.2,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("Calcium ","(mg kg"^-1,")")))
p4
p5 <- labelfun(plot0 = p5,
               title1 = "e",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("Magnesium ","(mg kg"^-1,")")))
p6 <- labelfun(plot0 = p6,
               title1 = "f",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("CEC ","(","meq"%.%"100g"^-1,")")))
p6
p7 <- labelfun(plot0 = p7,
               title1 = "g",
               top = 0,
               right = 0.8,
               left = 0.8,
               bottom = 0,
               ylab1 = expression(paste("soil pH ","(1:1w)")))
p7
#
phos1 <- labelfun(plot0 = phos,
                  title1 = "h",
                  top = 0,
                  right = 1.2,
                  left = 0.8,
                  bottom = 0,
                  ylab1 = expression(paste("Phosphorus ","(mg kg"^-1,")")))
###### 
ft2 <- 12
#generate final plot
#this may break in the future as better methods develop 
#make each item a grob and then arrange those in a grid 
library(gridExtra)
library(grid)
ggsave(filename = "Figures_PRINT/FigureS3_ppm.pdf",
       dpi=600,
       plot = grid.arrange(
         mylegend.1,
         arrangeGrob(p1,p2,p3,p4,p5,p6,p7,phos1,nrow=2),
         heights=list(unit(7,"mm"),unit(93,"mm")),
         padding=unit(3,"mm"),
         # left=textGrob("Soil Chemistry (0-20 cm)",rot=90,
         #               gp=gpar(cex=1,fontsize=ft2, family="Helvetica")),
         bottom=textGrob("Number of Plant Species (Log Scale)",
                         gp=gpar(cex=1,fontsize=ft2, family="Helvetica"))),
       height=130,width=183,unit="mm")
###########
library(flextable)
library(officer)
library(tidyverse)
colnames(mods.t)
mods1 <- mods.t %>% select(nut,Year,statistic,p.value.adj,r.squared) %>% 
  arrange(Year,nut)
table <- mods1
table$r.squared <- round(table$r.squared,2)
table$statistic <- round(table$statistic,3)
table$p.value <- round(table$p.value.adj,5)
table$p.value.adj <- NULL
#
table <- table %>% filter(nut %ni% c("CEC","pH"))
colnames(table) <- c("Soil Variable","Year","F-Statistic",
                     "R^2","P-value")

#
tmp <- flextable(table) %>% fontsize(size=12)
set_table_properties(tmp, width = 1, layout = "autofit")
tmp <- font(tmp,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 3", style = "Normal") %>%
  body_add_flextable(value = tmp)
print(doc, target = "Tables/SupplementalTable3_ppm_raw.docx")
#########
#plot function to show all the data 
yearfun <- function(Year2){
  Year1 <- Year2
  plot_all <- function(x,Year1){
    # x="Nitrogen";Year1=2017
    newdat <- data.frame(NumSp=seq(from=1,to=16,by=0.2))
    dat1 <- dat %>% filter(nut==x) %>% 
      filter(NumSp!=0) %>% 
      filter(Year==Year1)
    fit <- dat1 %>% 
      group_by(nut,Year) %>%
      nest() %>%
      mutate(mod=map(data,~lm(val~log(NumSp),data=.)),
             mods.t = map(mod,augment,
                          newdata=newdat,
                          se_fit=TRUE)) %>% 
      unnest(mods.t) %>% 
      select(-c(data,mod))
    colnames(fit)[4] <- "val"
    colnames(fit)[5] <- "se"
    #
    if(Year1==2017){
      shp=21;col1="#d95f02"
    }else{shp=23;col1="#1b9e77"}
    p1 <- ggplot(dat1,aes(x=log(NumSp),y=val,fill=Year,col=Year,shape=Year))+
      theme_bw(base_size = 11.5,base_family = "Helvetica")+
      xlab("Number of Species")+
      ylab("")+
      scale_y_continuous(
                         n.breaks=5,
                         breaks = waiver())+
      scale_x_continuous(breaks=c(log(1),log(2),log(4),log(8),
                                  log(16)),
                         labels=c("1","2",
                                  "4","8",
                                  "16"))+
      scale_fill_manual(values=col1,name=Year1)+
      scale_color_manual(values=col1)+
      guides(color="none")+
      theme(legend.title = element_blank(),
            legend.position = "top",
            legend.key.height =  unit(0.02, "cm"),
            legend.key.width = unit(0.4,"cm"),
            legend.text = element_text(size=15),
            legend.box.spacing = unit(0,"cm"),
            legend.box.margin = margin(0,0,0,0,"cm"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=8),
            # axis.text.x = element_text(),
            panel.grid = element_blank(),
            # title = element_blank(),
            plot.margin=margin(t = 0,r = 1.4,l = 0.8,b = 0))+
      geom_point(size=2,shape=shp,col="Black")+
      geom_line(aes(y=val),data=fit,col="Black")+
      geom_line(aes(y=val+se),alpha=0.3,data=fit,col="Black")+
      geom_line(aes(y=val-se),alpha=0.3,data=fit,col="Black")
    p1
    return(p1)
  }
  b1 <- plot_all("Carbon",Year1)
  b1
  mylegend.2<-legendfunc(b1)
  ft <- 8
  leftjust <- 0.2
  b1 <- b1+theme(legend.position = "none")+
    ggtitle("a")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Carbon ","(mg kg"^-1,")")))
  # 
  b2 <- plot_all("Nitrogen",Year1)
  #
  b2 <- b2+theme(legend.position = "none")+
    ggtitle("b")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Nitrogen ","(mg kg"^-1,")")))
  
  #
  b3 <- plot_all("Potassium",Year1)
  b3 <- b3+theme(legend.position = "none")+
    ggtitle("c")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Potassium ","(mg kg"^-1,")")))
  #
  b4 <- plot_all("Calcium",Year1)
  b4 <- b4+theme(legend.position = "none")+
    ggtitle("d")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Calcium ","(mg kg"^-1,")")))
  #
  b5 <- plot_all("Magnesium",Year1)
  #
  b5 <- b5+theme(legend.position = "none")+
    ggtitle("e")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Magnesium ","(mg kg"^-1,")")))
  #
  b6 <- plot_all("CEC",Year1)
  #
  b6 <- b6+theme(legend.position = "none")+
    ggtitle("f")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("CEC ","(","meq"%.%"100g"^-1,")")))
  #
  b7 <- plot_all("pH",Year1)
  b7 <- b7+theme(legend.position = "none")+
    ggtitle("g")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("soil pH ","(1:1w)")))
  #
  b8 <- plot_all("Phosphorus",Year1)
  #
  b8 <- b8+theme(legend.position = "none")+
    ggtitle("h")+
    theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                    margin=margin(0,0,0,0)),
          
          plot.margin = margin(0,0.8,0,0))+
    ylab(expression(paste("Phosphorus ","(mg kg"^-1,")")))
  #
  ft2 <- 8
  #
  if(Year1==2017){
    fig=4
  }else{fig=5}
  
  name <- paste("Figures_PRINT/FigureS",fig,".pdf",sep="")
  ggsave(filename = name,
         plot = grid.arrange(
           mylegend.2,
           arrangeGrob(b1,b2,b3,b4,b5,b6,b7,b8,ncol=2,nrow=4),
           heights=list(unit(7,"mm"),unit(93,"mm")),
           padding=unit(3,"mm"),
           # left=textGrob("Soil Chemistry (0-20 cm)",rot=90,
           #               gp=gpar(cex=1,fontsize=ft2, family="Helvetica")),
           bottom=textGrob("Number of Plant Species (Log Scale)",
                           gp=gpar(cex=1,fontsize=ft2, family="Helvetica"))),
         height=130,width=120,unit="mm")
}
yearfun("1994")
yearfun("2017")

