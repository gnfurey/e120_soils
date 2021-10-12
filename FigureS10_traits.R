####
library(tidyverse)
dat1 <- read.csv("SubmissionData/Dataset_S7.csv")
source("e120soils_functions.R")
#
AppendixS3col <- colnames(dat1)
write.csv(x = AppendixS3col,file = "Data/S3_DatasetS7.columns.csv",row.names = FALSE)

####
#get functional group
dat1$func <- get_funcgroup(get_specidfrom5(dat1$Species))
####
#get means 
dat1m <- dat1 %>% 
  group_by(func,nut) %>% 
  summarise(se=my.stand(val),
            val=mean(val))
#
table(dat1$nut,dat1$func)
###########
#generate mean ± SE for all traits 
plotfun <- function(x){
  # x="Potassium"
  tmp <- dat1m %>% 
    filter(nut==x)
  tmp1 <- dat1 %>% 
    filter(nut==x)
  p1 <- ggplot(tmp,aes(x=func,y=val,fill=func))+
    # geom_boxplot()+ 
    scale_fill_manual(
      name="Functional Group",
      values=
        c("#7fc97f","#1b9e77","#7570b3","#d95f02"))+
    # geom_point(shape=21,data=dat1)+
    theme_classic(base_size = 12)+
    ylab(x)+
    xlab("Functional Group")+
    theme(legend.position = "top",
          axis.title.x = element_blank())+
    geom_bar(stat="identity",col="Black")+
    # geom_jitter(data=tmp1,shape=21,size=3,height = 0,
    #             width=0.2)
    geom_errorbar(aes(ymax=val+se,ymin=val-se),width=0.5)
  p1
  return(p1)
}
plots <- map(.x = unique(dat1$nut),.f = plotfun)
p1 <- plots[[1]]
p1
p2 <- plots[[2]]
p3 <- plots[[3]]
p4 <- plots[[4]]
p5 <- plots[[5]]
p6 <- plots[[6]]
labelfun <- function(plot0,title1,top,right,left,bottom,ylab1){
  plot0 <- plot0+theme(legend.position = "none")+
    ggtitle(title1)+
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
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = "Shoot Tissue Calcium (%)")
p1
p2 <- labelfun(plot0 = p2,
               title1 = "b",
               top = 0,
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = "Shoot Tissue Magnesium (%)")
p3 <- labelfun(plot0 = p3,
               title1 = "c",
               top = 0,
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = "Shoot Tissue Nitrogen (%)")
p4 <- labelfun(plot0 = p4,
               title1 = "d",
               top = 0,
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = "Shoot Tissue Phosphorus (%)")
p5 <- labelfun(plot0 = p5,
               title1 = "e",
               top = 0,
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = "Shoot Tissue Potassium (%)")
p6 <- labelfun(plot0 = p6,
               title1 = "f",
               top = 0,
               right = 1.5,
               left = 1.2,
               bottom = 0,
               ylab1 = expression(paste("Root Biomass ","(g"%.%"m"^-2,")")))
p6
#
library(gridExtra)
library(grid)
ft2=12
ggsave(filename = "Figures_PRINT/FigureS10.pdf",
       dpi=600,
       plot = grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2,
                           padding=unit(5,"mm"),
                           bottom=textGrob("Functional Group",
                                           gp=gpar(cex=1,fontsize=ft2, family="Helvetica"))),
       height=130,width=183,unit="mm")
####