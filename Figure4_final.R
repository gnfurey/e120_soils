#create a new data set 
dat.w <- read.csv("SubmissionData/Dataset_S2.csv")
AppendixS3col <- colnames(dat.w)
write.csv(x = AppendixS3col,file = "Data/S3_DatasetS2.columns.csv",row.names = FALSE)

newdat <- dat.w
#resid
residplots <- function(mod){
  foo <- data.frame(fitted=fitted(mod),resid=residuals(mod))
  p1 <- ggplot(foo,aes(x=fitted(mod),y=residuals(mod)))+
    geom_point(shape=21)+
    theme_bw()
  return(p1)
}
#
#
#run regression
mymod0 <- lm(RootBiomass~1,data=newdat)
mymod1 <- lm(RootBiomass~Potassium,data=newdat)
mymod2 <- lm(RootBiomass~Nitrogen,data=newdat)
mymod3 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat)
summary(mymod3)#fig 4 caption
mymod3.noc3 <- lm(RootBiomass~Potassium+Nitrogen,data=newdat[newdat$func!="C3",])
summary(mymod3.noc3)#fig 4 caption
#
AIC(mymod0,mymod1,mymod2,mymod3)
#backwards selection
mymod.step <- step(mymod3)
summary(mymod.step)
#run some tests - results are robust 
car::Anova(mymod3,type="III")
shapiro.test(residuals(mymod3))
lmtest::bptest(mymod3)
car::ncvTest(mymod3)
#
#create vectors 
z=newdat$RootBiomass
y=newdat$Nitrogen
x=newdat$Potassium
#get regression data set
regvar <- data.frame(x=x,y=y,z=z)
#run
mymod <- lm(z~x+y,data=regvar)
summary(mymod)
#
# car::scatter3d(x=x,y=z,z=y)
# s3d <- scatterplot3d::scatterplot3d(x=x,y=y,z=z,
#                              zlim=c(min(z)*0.9,max(z)*1.1),
#                              highlight.3d = TRUE,
#                              scale.y=1, angle=-130)
#get fits
fit <- fitted(mymod)
#get resids
res <- residuals(mymod,type = "response")
#set gridlines
grid.lines = 20
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
#get matrix
xy <- expand.grid( x = x.pred, y = y.pred)
#get predictions
z.pred <- matrix(predict(mymod, newdata = xy),
                 nrow = grid.lines, ncol = grid.lines)
#get colors 
newdat$col <- ifelse(newdat$func=="C4","#1b9e77",
                     ifelse(newdat$func=="L","#d95f02",
                            ifelse(newdat$func=="F","#7570b3","#7fc97f")))
######
#see ?rgl
library(rgl)
clear3d()
persp3d(x = x.pred,
        y = y.pred,
        z = z.pred,
        axes=TRUE,
        box=TRUE,
        add=TRUE,
        ticks=TRUE,
        ylab="",
        zlab="",
        zlim=c(min(z)*0.9,max(z)*1.1),
        front="lines",back="lines",
        aspect=TRUE)
aspect3d(1,1,1)
persp3d(x = x.pred,
        add=TRUE,
        box=TRUE,
        y = y.pred,
        z = z.pred,
        xlab="",
        ylab="",
        zlab="",
        color = "grey",
        zlim=c(min(z)*0.9,max(z)*1.1),
        alpha=0.4)
aspect3d(1,1,1)
##We adjusted this point to avoid overplotting
#Andropogon and Sorghastrum have similar values 
# newdat[newdat$Species=="Andge","Nitrogen"] <- 0.77#was .709
# newdat[newdat$Species=="Andge","Potassium"] <- 0.91# was .956
x
y
x[3] <- 0.91
y[3] <- 0.77


plot3d(x=x, y=y, z=z,size=2,
       lit=TRUE,
       ambient="black",
       smooth=TRUE,
       shininess=100,
       type="s",col=newdat$col,
       xlab="",
       add=TRUE,
       ylab="",
       zlab="")
aspect3d(1,1,1)
#
zstart <- z
zfinish <- z-res
zboth <- c(rbind(zstart,zfinish))
segments3d(x=rep(x,2,each=2), y=rep(y,2,each=2),
           z=zboth, lwd=5,
           # col="red",
           xlab="",
           add=TRUE,
           ylab="",
           zlab="")
aspect3d(1,1,1)
bbox3d(color = c("black", "black"),
       # yat = c(0.5,1,1.5,2),
       xlen = 0, ylen = 0, zlen = 0,
       # specular = "grey",
       expand=1,
       front="lines",
       back="lines",
       # yunit = "pretty",
       # xunit=""
       # ylim=c(0,2),
       # xat=c(0,0.5,2),
       # yat=c(0,2),
       # zat=c(0,2),
       # zlim=c(0,2),
       edge="z-+",
       edge="y+-", alpha = 1)
aspect3d(1,1,1)
#
########
axes3d(edges = c("x+-"),cex=2,
       # at=c(0.7,1.4,2.1),
       labels=c(0.8,"",1.2,"",1.6,"",2),
       tick=TRUE,
       nticks=7)
axes3d(edges = c("y+-"),cex=2,
       # at=c(0.6,1.2,1.8)),
       labels=c(0.6,"",1,"",1.4,"",1.6,"",2),
       tick=TRUE,
       nticks=7)
axes3d(edges = c("z+-"),cex=2,
       # at=c(170,370,670)),
       labels=c(200,
                "","",400,
                "","",600,
                "","",800),
       tick=TRUE,
       nticks=10)
axes3d(edges = c("z-+"),cex=2,
       # at=c(170,370,670)),
       labels=c(200,
                "","",400,
                "","",600,
                "","",800),
       tick=TRUE,
       nticks=10)
library(webshot)
# rgl.snapshot("Figures_Print/Figure4_raw.png")
snapshot3d("Figures_Print/Figure4_raw.png",webshot = FALSE,
           height=2000,width=2000*1.61)
#
newdat$col <- ifelse(newdat$func=="C4","#1b9e77",
                     ifelse(newdat$func=="L","#d95f02",
                            ifelse(newdat$func=="F","#7570b3","#7fc97f")))
clear3d()
par3d(windowRect = c(100, 100, 612, 612))
legend3d("center",
         legend = c('C4 Grass',
                    "C3 Grass",
                    'Legume',
                    'Forb'),
         cex=5,
         pch = 16, col = c("#1b9e77","#7fc97f",
                           "#d95f02","#7570b3"
                           ))
snapshot3d("Figures_Print/Legend.png",
           height=2000,width=2000,
           webshot = FALSE)
#,
# mtext3d("Plant % K", edge="x+-", line = 4,cex=2)
# mtext3d("Root Biomass", edge="z+-", line = 3,cex=2)
# mtext3d("Plant % N", edge="y+-", line = 4,cex=2)
#
##########
library(flextable)
library(officer)
library(tidyverse)
source("e120soils_functions.R")
supptab <- dat.w %>% select(Species,func) %>% distinct()
supptab$Family <- get_family(get_specidfrom5(supptab$Species))
supptab$Species <- get_species(get_specidfrom5(supptab$Species))
#
tmp <- flextable(supptab) %>% fontsize(size=12)
set_table_properties(tmp, width = 1, layout = "autofit")
tmp <- font(tmp,fontname = "Times")
doc <- read_docx() %>% 
  body_add_par(value = "Supplemental Table 7", style = "Normal") %>%
  body_add_flextable(value = tmp)
print(doc, target = "Tables/SupplementalTable7_raw.docx")
