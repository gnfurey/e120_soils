library(colorBlindness)
library(ggplot2)
p1 <- readRDS("Figures_PRINT/Figure1.rds")
cvdPlot(p1)
p2 <- readRDS("Figures_PRINT/Figure2.rds")
cvdPlot(p2)
p3 <- readRDS("Figures_PRINT/Figure3.rds")
cvdPlot(p3)
#test plot 4 not using ggplot
colvec <- c("#1b9e77","#7fc97f",
  "#d95f02","#7570b3"
)
tmp <- data.frame(x=c(1:4),y=c(4),col=colvec,treat=c("C4","C3","Legume","Forb"))
p4 <- ggplot(tmp,aes(x=x,y=y,fill=treat))+
  geom_bar(stat="Identity",col="Black")+
  scale_fill_manual(values = colvec)
p4
cvdPlot(p4)
