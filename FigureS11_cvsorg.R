#
library(tidyverse)
#read in data
e120 <- read.csv("SubmissionData/Dataset_S1.csv")
#long format 
dat <- e120 %>% select(Plot,NumSp,
                       gm2_2017_0.20cm_Calcium_aa:gm2_2017_0.20cm_Potassium_aa,
                       gm2_1994_0.20cm_Calcium_aa:gm2_1994_0.20cm_Potassium_aa) %>% 
  pivot_longer(gm2_2017_0.20cm_Calcium_aa:gm2_1994_0.20cm_Potassium_aa,
               names_to="nut",values_to="val") %>% 
  separate(nut,into = c("type","Year","Depth","nut","extract"),sep = "_")
#subset variables
dat1 <- dat %>% 
  filter(nut %in% c("CEC","Carbon")) %>% 
  select(-extract) %>% 
  spread(nut,val)
#note for future reference
#there are many variables that may be correlated with CEC
#There is considerable theoretical and empirical evidence that soil organic carbon
#should increase CEC in a highly sandy soil
#plotting function  
yearfun <- function(x){
  # x=2017
  library(grid)
  #subset
  tmp <- dat1 %>% 
    filter(Year==x)
  #points
  if(x==2017){
    shp=21;col1="#d95f02"
  }else{shp=23;col1="#1b9e77"}
  #run major axis regression 
  sma1 <- smatr::sma(CEC~
                       Carbon,
                     data=tmp)
  print(summary(sma1))
  #
  r2 <- round(sma1$r2[[1]],2) 
  r2
  coef <- sma1$coef[[1]]
  a <- coef[1,1]
  b <- coef[2,1]
  #prove that I did indeed pass grade 7 math 
  tmp$fit <- tmp$Carbon*b+a
  #plot 
  print(min(tmp$Carbon))
  print(min(tmp$CEC))
  print(max(tmp$Carbon))
  print(max(tmp$CEC))
  p1 <- ggplot(tmp,aes(x=Carbon,
                      y=CEC,
                      fill=Year,
                      shape=Year))+
  # facet_wrap(~Year)+
    geom_point(size=4,shape=shp,col="Black")+
    scale_fill_manual(values=col1,name=x)+
    scale_color_manual(values=col1)+
    guides(shape="none")+
    scale_x_continuous(limits=c(580,5000))+
    scale_y_continuous(limits=c(1,7.5))+
  theme_bw(base_size = 12)+
    geom_line(aes(y=fit))+
  # geom_abline(slope=b,intercept = a,data=tmp)+
  xlab(expression(paste("Soil Organic Carbon ","(","g"%.%"m"^-2,")")))+
    ylab(expression(paste("Cation Exchange Capacity ","(","meq"%.%"100g"^-1,")")))+
    guides(fill="none")
  p1
  return(p1)
}
#run for each year
yearfun("1994")
yearfun("2017")
#add titles etc .
out1 <- yearfun("1994")+
  ggtitle("a")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))
out1
out2 <- yearfun("2017")+
  ggtitle("b")+
  theme(plot.title = element_text(face="bold",family="Helvetica",size=8,
                                  margin=margin(0,0,0,0)),
        plot.margin = margin(0,0.8,0,0))
out2
library(gridExtra)
ggsave(plot = grid.arrange(out1,out2,ncol=2),
  filename = "Figures_PRINT/FigureS11.pdf",
       height=130,width=183,unit="mm")
#
#over all points
# ggplot(dat1,aes(x=Carbon,y=CEC))+
#   geom_point()
#