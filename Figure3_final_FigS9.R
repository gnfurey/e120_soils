library(tidyverse)
# require(devtools)
#we used compact letter display before this update
# install_version("emmeans", version = "1.5.2", repos = "http://cran.us.r-project.org")
library(emmeans)
#custom functions
source("e120soils_functions.R")
#read in data
e120 <- read.csv("SubmissionData/Dataset_S1.csv")
#
table(e120$Fgset2)
#long format
colnames(e120)
dat2 <- e120 %>% select(Plot:Fgset2,Calcium_total.soil.plant.g.m2:
                          Potassium_total.soil.plant.g.m2) %>% 
  pivot_longer(Calcium_total.soil.plant.g.m2:Potassium_total.soil.plant.g.m2,names_to = "nut",
               values_to="total") %>% 
  separate(nut,into=c("nut","drop"),sep="_") %>% 
  select(-drop)
dat2$Fgset2 <- as.factor(dat2$Fgset2)
dat2$Fgset2 = factor(dat2$Fgset2,
                     levels(dat2$Fgset2)[c(3,1,7,4,2,6,5)])

#create a new dataset in case of errors 
fset <- dat2
fset$Fgset2 <- as.factor(fset$Fgset2)
fset$Fgset2 = factor(fset$Fgset2,
                     levels(fset$Fgset2)[c(3,1,7,4,2,6,5)])
#function to run statistics 
tabs <- function(x,y){
  # x="Phosphorus"
  mods <- fset %>% 
    filter(nut==x)
  library(nlme)
  library(emmeans)
  #the FL combination only has 5 points and it fails 
  #when using a separate weight. I gave it a fixed weight at the cost
  #of increasing its degrees of freedom with the benefit of it
  #not having 0.8 estimated df.  
  #source for using do.call because of some weird formula parsing
  #https://stackoverflow.com/questions/7666807/anova-test-fails-on-lme-fits-created-with-pasted-formula
  mod1 <- do.call("gls", args = list(total~Fgset2,
                                     weights=varIdent(form=~1|Fgset2,
                                                      fixed=c(
                                                        FL=1
                                                      )),
                                     
                                     data=mods))
  summary(mod1)
  anova(mod1)
  #get resid
  mods$resid <- residuals(mod1,type="normalized")
  # residplot(mod1)
  #plot resid 
  ggplot(mods,aes(x=Fgset2,y=resid))+
    geom_boxplot()
  #get least square means 
  em <- emmeans(mod1,specs=~Fgset2,
                # mode="df.error",
                sigmaAdjust = TRUE)
  em
  # pwpm(em)
  pairs(em)
  ##
  #run t-test different from zero
  test <- test(em, null = 0, side = "both")
  #blatantly ignore the suggestion not to use CLD display
  cld <- CLD(em,Letters = LETTERS,reversed=TRUE)
  #get letter groups 
  cld <- cld %>% select(Fgset2,.group)
  #
  out <- as.data.frame(em)
  test <- as.data.frame(test)
  #
  out1 <- left_join(out,test)
  out1 <- left_join(out1,cld)
  out1$star <- ifelse(out1$p.value<0.05,"*","")
  #remove whitespace 
  out1$.group <- str_trim(out1$.group)
  return(out1)
}
#get args 
nuts <- unique(fset$nut)
nuts
#setting names so that we get nice rownames
nuts <- set_names(nuts)
#run function using purrr
out_table <- map_dfr(.x = nuts,
                     .f = tabs,.id = "Nutrient")
#adjustt-test pvals
out_table$p.value_adjusted <- p.adjust(out_table$p.value,"fdr")
out_table$adjusted_star <- ifelse(out_table$p.value_adjusted<0.05,"*","")
##############
#change names to short nutrients 
# out_table$Nutrient <- get_shortnut(out_table$Nutrient)
#function to plot means
fgfun <- function(x){
  # x <- "Calcium"
  #get raw data 
  mods <- fset %>% 
    filter(nut==x)
  #get lsmeans 
  out1 <- out_table %>% 
    filter(Nutrient==x)
  # return(out1)
  # min <- min(mods$total)*1.4
  #
  if(x!="Phosphorus"){
  min <- 0
  # max <- max(mods$total)*1.05
  max <- max(out1$emmean+out1$SE)*1.05
  print(max)
  }else{
  min <- -4  
  }
  #generate plot 
  e120 <- read.csv("SubmissionData/Dataset_S1.csv")
  
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "soil_diffgm2_0.20cm_",replacement ="" )
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_ea",replacement ="" )
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_aa",replacement ="" )
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_LOI",replacement ="" )
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_1.1w",replacement ="" )
  colnames(e120) <- str_replace_all(string = colnames(e120),pattern = "_bray",replacement ="" )
  
  colnames(e120)[6:14] <- paste(colnames(e120)[6:14],"soil_soil",sep="_")
  colnames(e120)
  #
  datm <- e120 %>% select(Plot,NumSp,
                         Fgset,Fgset2,
                         Calcium_soil_soil:
                           Potassium_soil_soil,
                         Nitrogen_gm2_Abv:Calcium_gm2_Root.0.30cm) %>% 
    pivot_longer(Calcium_soil_soil:Calcium_gm2_Root.0.30cm,
                 names_to="nut",values_to="val") %>% 
    separate(nut,into = c("nut","drop","pool"),sep = "_") %>% 
    mutate(drop=NULL) %>% 
    filter(nut %in% c("Nitrogen",
                      "Phosphorus",
                      "Potassium","Calcium","Magnesium")) %>% 
    group_by(Fgset2,nut,pool) %>%
    summarise(val=mean(val)) %>% 
    filter(nut==x)
  #
  datm$Fgset2 <- as.factor(datm$Fgset2)
  datm$Fgset2 = factor(datm$Fgset2,
                       levels(datm$Fgset2)[c(3,1,7,4,2,6,5)])
  
  #
  datm <- left_join(datm,out1)
  #
  datm <- as.data.frame(datm)
  #
  p1 <- ggplot(datm,aes(x=Fgset2,
                  y=val,
                  fill=pool))+
    geom_hline(yintercept = 0)+
    geom_bar(col='Black',stat="identity")+
    scale_fill_manual(name="",
                      values = c("#999999", "#E69F00", "#56B4E9"),
                      labels = c("Shoots","Roots","Soil"))+
    geom_errorbar(inherit.aes = FALSE,
                  aes(x=Fgset2,ymax=emmean+SE,ymin=emmean-SE),
                  col="Black",width=0.3,data=out1)+
    theme_classic(base_size = 7)+
    # geom_text(aes(y=min,label=adjusted_star),size=7,col="Black")+
    ylab(x)+
    xlab("Functional Group Set")+
    ggtitle(unique(out1$title))+
    # guides(fill=FALSE)+
    theme(
      legend.position = "bottom",
      legend.key.height =  unit(0.02, "cm"),
      legend.key.width = unit(0.4,"cm"),
      legend.text = element_text(size=10),
      legend.box.spacing = unit(0,"cm"),
      legend.box.margin = margin(0,0,0,0,"cm"),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(face="bold",family="Helvetica",size=8,
                                margin=margin(1,0,0,0)),
      plot.margin = margin(0,2,0,2))+
    geom_point(aes(x=Fgset2,
                   y=emmean),
               data=out1,
               size=2,
               fill="Black",
               col="Black")
  # data=out1)
  #
  if(x!="Phosphorus"){
    datm <- datm %>%
      ungroup()
    datm <- as.data.frame(datm)
    out1 <- out1 %>% 
      ungroup()
    tmp00 <- out1 %>% 
      select(Fgset2,.group) %>% 
      distinct()
    p1 <- p1+
      geom_text(inherit.aes = FALSE,
      aes(x=Fgset2,y=max,label=.group,
          group=NULL),
                       size=2.4,
      col="Black",
      data=tmp00)
    p1
    # p1
    }
  p1
  return(p1)
} 
library(grid)
nuts <- c("Nitrogen",
          "Potassium",
          "Calcium",
          "Magnesium",
          "Phosphorus")
plots <- map(.x = nuts,.f = fgfun)
plots[[1]]
#
ht <- 0.8
leftjust <- 0.2
ft <- 8
#here we extract each plot as a list 
#then we annotate each plot with a label
plots[[1]] <- plots[[1]]+
  annotation_custom(textGrob(
    x = unit(leftjust, "npc"),
    y = unit(ht, "npc"),
    label="Nitrogen",
    gp=gpar(fontsize=ft,face="bold")))+
  ggtitle("a")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))
plots[[1]]
plots[[2]] <- plots[[2]]+ggtitle("b")+
  annotation_custom(textGrob(
    x = unit(leftjust, "npc"),
    y = unit(ht, "npc"),
    label="Potassium",
    gp=gpar(fontsize=ft,face="bold")))+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))
plots[[3]] <- plots[[3]]+ggtitle("c")+
  annotation_custom(textGrob(
    x = unit(leftjust, "npc"),
    y = unit(ht, "npc"),
    label="Calcium",
    gp=gpar(fontsize=ft,face="bold")))+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))
plots[[4]] <- plots[[4]]+ggtitle("d")+
  annotation_custom(textGrob(
    x = unit(0.25, "npc"),
    y = unit(ht, "npc"),
    label="Magnesium",
    gp=gpar(fontsize=ft,face="bold")))+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))

library(gridExtra)
library(grid)
ft2 <- 10
legendfunc<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
mylegend.1<-legendfunc(plots[[1]])
p1 <- plots[[1]]+theme(legend.position = "none")
p1 <- p1+
  scale_y_continuous(breaks=seq(from=-20,60,by=10),
                     limits=c(-11,63),
                     labels=every_nth(seq(from=-20,60,by=10),2,inverse = TRUE))


p2 <- plots[[2]]+theme(legend.position = "none")
p2
p2 <- p2+
  scale_y_continuous(breaks=seq(from=-2,12,by=2),
                     limits=c(-3.5,13),
                     labels=every_nth(seq(from=-2,12,by=2),2,inverse = FALSE))

p3 <- plots[[3]]+theme(legend.position = "none")
p3
p3 <- p3+ scale_y_continuous(breaks=seq(from=0,60,by=10),
                        limits=c(0,64),
                        labels=every_nth(seq(from=0,60,by=10),2,inverse = TRUE))

p4 <- plots[[4]]+theme(legend.position = "none")
p4
p4 <- p4+
  scale_y_continuous(breaks=seq(from=0,11.5,by=1.5),
                     limits=c(0,11.5),
                     labels=every_nth(seq(from=0,11.5,by=1.5),2,inverse = TRUE))

#run grid.arrange function across the list using extra args 
#
o1 <- grid.arrange(
  mylegend.1,
  arrangeGrob(p1,p2,p3,p4,ncol=2),
  heights=list(unit(7,"mm"),unit(93,"mm")),
  padding=unit(3,"mm"),
  bottom=textGrob("Functional Group Set",
                  gp=gpar(cex=1,fontsize=ft2, 
                          family="Helvitica")),
  left=
    textGrob(expression(paste(
      "Change in Nutrient Pool (Biomass + Soils) ","(g"%.%"m"^-2,")")),rot=90,
      gp=gpar(cex=1,fontsize=ft2, 
                       family="Helvitica")))
o1
ggsave(plot=o1,
       dpi=600,
       filename = "Figures_PRINT/Figure3.pdf",
       height=11,width=8.7,unit="cm")
####
saveRDS(o1,"Figures_PRINT/Figure3.rds")
####
#generate a separate plot for P 
phos <- plots[[5]]+
  # annotation_custom(textGrob(
  #   x = unit(.8, "npc"),
  #   y = unit(.1, "npc"),
  #   label="Phosphorus",
  #   gp=gpar(fontsize=15,face="bold")))+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0),
        axis.text = element_text(size=10),
        axis.title.y = element_text(size=10,angle = 90),
        axis.title.x = element_text(size=10),
        legend.position = "top")+
  ylab(expression(paste("Change in Phosphorus Pool (Biomass + Soils) ","(g"%.%"m"^-2,")")))
phos
# saveRDS(phos, file = "P_func.rds")

ggsave(plot=phos,
       dpi = 600,
       filename = "Figures_PRINT/FigureS9.pdf",
       height=100,width=89,unit="mm")
#
##############
library(flextable)
library(officer)
library(tidyverse)
#rename data set
table <- out_table
table$Fgset2 = factor(table$Fgset2,
                     levels(table$Fgset2)[c(2,5,1,4,7,6,3)])
table$Fgset2
colnames(table)
#round 
table$Fgroup
table <- table %>% 
  mutate(across(c(emmean:t.ratio),~round(.x,2)))
table$df <- round(table$df,0) 
table$p.value <- round(table$p.value,5)
table$p.value <- ifelse(table$p.value<0.0001,"<0.0001",table$p.value)
table$p.value_adjusted <- round(table$p.value_adjusted,5)
table$p.value_adjusted <- ifelse(table$p.value_adjusted<0.0001,"<0.0001",table$p.value_adjusted)
#
colnames(table)
table <- table %>% select(-c(star,p.value,p.value_adjusted,adjusted_star))
#
colnames(table)
#
colnames(table)[2] <- "Fgroup"
#
table$Fgroup
colnames(table)[3] <- "Mean"
colnames(table)[9] <- "Group"
tmp <- flextable(table) %>% fontsize(size=12)
#
table <- table %>% arrange(Nutrient,Fgroup)
table
table$t.ratio <- NULL
tmp <- flextable(table)
tmp <- font(tmp,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 5", style = "Normal") %>%
  body_add_flextable(value = tmp)
print(doc, target = "Tables/SupplementalTable5_raw.docx")

