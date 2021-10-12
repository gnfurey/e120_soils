#functions

#function to get shortnames 
get_shortnut <- function(x){
  # x="Calcium"
  dat <- read.csv("Data/nutmatch.csv",stringsAsFactors = FALSE)
  out <- dat[match(x,dat$Nutrient),"Short"]
  out
  return(out)
}
#general metadata for ccesr (user beware!)
path <- "Data/ccesrPlantSpeciesData.csv"
#get family
get_family<-function(Specid) {
  taxdat <-read.csv(path)
  taxdat[match(Specid, taxdat$Specid),"Family"]
}
#get functional group
get_funcgroup_noc3 <- function(Specid){
  # Specid="Amoca"
  taxdat <-read.csv(path)
  fg <- taxdat[match(Specid, taxdat$Specid),"FunctionalGroup"]
  fg <- as.character(fg)
  fg <- ifelse(fg=="C3","Grass",
               ifelse(fg=="C4","Grass",
                      ifelse(fg=="F","Forb",
                             ifelse(fg=="L","Legume",
                                    fg))))
  return(fg)
}
#
#fiveIDs
get_specidfrom5<-function(idnum) {
  taxdat <- read.csv(path)
  taxdat[match(idnum, taxdat$X5Lspecid),"Specid"]
}
#grassforb 
get_grassforb<-function(Specid) {
  taxdat <-read.csv(path)
  fg <- taxdat[match(Specid, taxdat$Specid),"FunctionalGroup"]
  ifelse(fg=="C3","Mono",
         ifelse(fg=="C4","Mono",
                ifelse(fg=="S","Mono",
                       ifelse(fg=="F","Dicot",
                              ifelse(fg=="L","Dicot",
                                     ifelse(fg=="W","Woody","mistake"))))))
}
#specid
get_specid<-function(Species) {
  taxdat <-read.csv(path)
  taxdat[match(Species, taxdat$Species),"Specid"]
}
#functiongroup
get_funcgroup<-function(Specid) {
  taxdat <-read.csv(path)
  fg <- taxdat[match(Specid, taxdat$Specid),"FunctionalGroup"]
  fg
}
#https://stackoverflow.com/questions/2676554/in-r-how-to-find-the-standard-error-of-the-mean
my.stand <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
#species 
get_species<-function(Specid) {
  taxdat <-read.csv(path)
  taxdat[match(Specid, taxdat$Specid),"Species"]
} 
#invert %in% 
"%ni%" <- Negate("%in%")
#
#https://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-rhttps://stackoverflow.com/questions/34533472/insert-blanks-into-a-vector-for-e-g-minor-tick-labels-in-r
every_nth <- function(x, nth, empty = TRUE, inverse = FALSE){
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}
#percentchange 
pchange <- function(new,old){
  x <- ((new-old)/abs(old))*100
  return(x)
}
