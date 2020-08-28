# 
## Programa MONITORA-ICMBio
## Analise de dados
## Elildo Carvalho Jr @ ICMBio/CENAP
## Partes do script escritas por Jorge Ahumada @ Conservation International
#

## ----Load libraries------------------------
library(TeachingDemos)
library(lubridate)
library(unmarked)
#library(ggplot2)
library(dplyr)
library(chron)
library(vegan)
#library(activity)
#library(ggmap)
library(here)


## ----Source this file--------
source(here("bin", "ahumada_codes.R"))

## ----Load data and do some fixes-------
#dataRBG <- f.readin.fix.data("Wild_ID_RBG_2016.csv")
dataRBG <- f.readin.fix.data(here("data", "Wild_ID_RBG_2017.csv"))

levels(dataRBG$bin)[levels(dataRBG$bin)=="Dasypus sp"] <- "Dasypus unknown" # rename by name
levels(dataRBG$bin)[levels(dataRBG$bin)=="Psophia viridis"] <- "Psophia unkown" # rename by name

dataRBG$Photo.time <- as.POSIXct(dataRBG$Photo.time, origin = dataRBG$Photo.Date, tz="GMT") # fixing format

## ----Filter independent records-------

dataRBG <- f.separate.events(dataRBG, 60) # Group by events that are 60 minutes apart
dataRBG<- distinct(dataRBG, bin, grp, .keep_all = TRUE) # filter retaining only independent records (distinctive combinations of bin and grp)

## ----Keep only species of interest -------

prepare.icmbio <- function(x) {
  # Select species subset for analysis
  species.list <- c("Cuniculus paca", "Dasyprocta prymnolopha", "Dasypus unknown","Leopardus pardalis", "Leopardus wiedii",
                    "Mazama americana", "Mazama nemorivaga", "Mitu tuberosum", "Myrmecophaga tridactyla", "Panthera onca",
                    "Pecari tajacu", "Penelope superciliaris", "Priodontes maximus", "Psophia unknown", 
                    "Puma concolor", "Tamandua tetradactyla", "Tapirus terrestris", "Tayassu pecari", 
                    "Tinamus tao") # create list of species to be included on analysis
  df <- subset(x, bin %in% species.list) # keeping only species in species.list
  df$bin <- factor(df$bin) # to remove excluded species from factor levels otherwise they will end up as zeros in the paMatrix
  df$Camera.Trap.Name <- factor(df$Camera.Trap.Name)
  assign("df", df, envir = globalenv())
}

prepare.icmbio(dataRBG)
dataRBG <- df


## ----Split sampling into two halves

start <- min(dataRBG$Camera.Start.Date)
end <- max(dataRBG$Camera.End.Date)
end-start

# check effort per site
df <- distinct(dataRBG, Camera.Trap.Name, Camera.Start.Date, Camera.End.Date, .keep_all=F)
df$effort <- as.numeric(with(df, Camera.End.Date - Camera.Start.Date))
hist(df$effort)

# split dataset into 1st and 2nd halves
first.half <- subset(dataRBG, Photo.Date < start+((end-start)/2))
second.half <- subset(dataRBG, Photo.Date > start+((end-start)/2))


# create dataframe with camera.trap.names and lat-long
df1 <- distinct(dataRBG, Camera.Trap.Name, Latitude, Longitude)

# fill df1 for effort
df1$eff.1 <- NA
df1$eff.2 <- NA

for(i in 1:nrow(df1)) {
  start <- min(dataRBG$Camera.Start.Date)
  end <- max(dataRBG$Camera.End.Date)
  mid <- start+((end-start)/2)
  cam <- subset(dataRBG, Camera.Trap.Name == df1[i,1])
  eff.1 <- as.numeric(mid-min(cam$Camera.Start.Date))
  eff.2 <- as.numeric(max(cam$Camera.End.Date)-mid)
  df1[i,4] <- eff.1
  df1[i,5] <- eff.1
}

# fill df1 for species trapping rates
species <- sort(unique(dataRBG$bin))

for (j in 1:length(species)) {
  for (i in 1:nrow(df1)) {
  cam <- subset(dataRBG, Camera.Trap.Name == df1[i,1])
  cam.1st <- subset(cam, Photo.Date < mid)
  cam.2nd <- subset(cam, Photo.Date > mid)
  cam.1st.sp.i <- subset(cam.1st, bin == species[j])
  cam.2nd.sp.i <- subset(cam.2nd, bin == species[j])
  df1[i,5+j] <- (nrow(cam.2nd.sp.i)/df1[i,5]) - (nrow(cam.1st.sp.i)/df1[i,4]) # delta TR
}}
colnames(df1) <- c(names(df1[1:5]), as.vector(species))
colnames(df1) <- make.names(colnames(df1), unique=TRUE) # insert dots in colnames
#View(df1)

hist(df1$`Cuniculus paca`)
hist(df1$`Dasyprocta prymnolopha`)


covars <- read.csv("/home/elildojr/Documents/r/Analysis code/team_gurupi_covars.csv")
with(covars, plot(altitude.gps, altitude.srtm))

# Join covars into df1
df2 <- left_join(df1, covars, by="Camera.Trap.Name")

barplot(df2$Dasyprocta.prymnolopha)
with(subset(df2, df2$Dasyprocta.prymnolopha != 0), plot(df2$altitude.gps, df2$Dasyprocta.prymnolopha, xlab="Elevation", ylab="deltaTR", pch=16))

barplot(df2$Tapirus.terrestris)
with(subset(df2, df2$Tapirus.terrestris != 0), plot(altitude.gps, Tapirus.terrestris, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Mazama.americana != 0), plot(altitude.gps, Mazama.americana, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Mazama.nemorivaga != 0), plot(altitude.gps, Mazama.nemorivaga, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Cuniculus.paca != 0), plot(altitude.gps, Cuniculus.paca, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Myrmecophaga.tridactyla != 0), plot(altitude.gps, Myrmecophaga.tridactyla, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Psophia.unknown != 0), plot(altitude.gps, Psophia.unknown, xlab="Elevation", ylab="deltaTR", pch=16))

with(subset(df2, df2$Mitu.tuberosum != 0), plot(altitude.gps, Mitu.tuberosum, xlab="Elevation", ylab="deltaTR", pch=16))



##--------NEXT STEPS------------

# Get distance to water or HAND covars
# Try for other sites and years: this was already at the end of the dry season: ranges already contracted by then?
# To work for other sites the covars file must of course have the corresponding Camera.Trap.Names...





