#
## Script to estimate animal densities withouth individual recognition
## written by Elildo Carvalho Jr @ ICMBio/CENAP
## data reading and fixing script written by Jorge Ahumada @ Conservation International
#

## ----Load libraries------------------------
#library(TeachingDemos)
#library(lubridate)
#library(unmarked)
library(ggplot2)
library(dplyr)
library(chron)
#library(vegan)
library(activity)
#library(ggmap)
library(here)
library(openair)

source(here("bin", "ahumada_codes.R"))
source(here("bin", "lifepulse.R"))

x <- f.readin.fix.data(here("data", "Wild_ID_RBG_2016.csv"))

# Run function lifepulse

lifepulse(x, "Dasyprocta prymnolopha", "min")
dasyprocta.min <- min60[1441:2880,]
with(dasyprocta.min, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Number of Dasyprocta records per minute in a 24-h window", ylab="Mean number of individuals / min"), )

fit0 <- lm(dasyprocta.hour$Number.of.Animals ~ dasyprocta.hour$date)
plot.ts(dasyprocta.hour$Number.of.Animals)
lines(dasyprocta.hour$date, predict(fit0), col = "red")

# Plotting fitted trend
lifepulse(x, "Dasyprocta prymnolopha", "day")
dasyprocta.day <- min60
with(dasyprocta.day, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Number of Dasyprocta records per hour in a 7-day window", ylab="Mean number of individuals / hour"), )
fit0 <- lm(dasyprocta.day$Number.of.Animals ~ dasyprocta.day$date) # fit can be used to extract slopes and use in the regression
lines(dasyprocta.day$date, predict(fit0), col = "red")

with(dasyprocta.day, scatter.smooth(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Number of Dasyprocta records per hour in a 7-day window", ylab="Mean number of individuals / hour"))
fit0 <- lm(dasyprocta.day$Number.of.Animals ~ dasyprocta.day$date) # fit can be used to extract slopes and use in the regression
lines(dasyprocta.day$date, predict(fit0), col = "red")


lifepulse(gurupi2017, "Dasyprocta prymnolopha", "day")
dasyprocta.day <- min60[1:60,]
with(dasyprocta.day, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Number of Dasyprocta records per day in a 60-day window", ylab="Mean number of individuals / day"), )

# Test with ttr package for time series

library("TTR")
dasyprocta.hourSMA200 <- SMA(dasyprocta.hour$Number.of.Animals,n=200)
plot.ts(dasyprocta.hourSMA200)

dasyproctaforecasts <- HoltWinters(dasyprocta.hour, beta=FALSE, gamma=FALSE)
plot(dasyproctaforecasts)

lifepulse(gurupi2017, "Tapirus terrestris", "day")
tapirus.day <- min60[1:60,]
with(tapirus.day, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Number of Tapirus records per day in a 60-day window", ylab="Mean number of individuals / day"), )

# in separate cameras
camera.names <- sort(unique(gurupi2017$Camera.Trap.Name))
lifepulse(subset(gurupi2017, Camera.Trap.Name == camera.names[1]), "Dasyprocta prymnolopha", "2 day")
dasyprocta.day.cam02 <- min60
lifepulse(subset(gurupi2017, Camera.Trap.Name == camera.names[2]), "Dasyprocta prymnolopha", "2 day")
dasyprocta.day.cam03 <- min60
lifepulse(subset(gurupi2017, Camera.Trap.Name == camera.names[3]), "Dasyprocta prymnolopha", "2 day")
dasyprocta.day.cam04 <- min60

par(mfrow=c(3,1))
with(dasyprocta.day.cam02[1:31,], plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Dasyprocta, CT-RBG-1-02", ylab="Mean number of individuals / day"), )
  fit1 <- lm(dasyprocta.day.cam02[1:31,]$Number.of.Animals ~ dasyprocta.day.cam02[1:31,]$date) # fit can be used to extract slopes and use in the regression
  lines(dasyprocta.day.cam02[1:31,]$date, predict(fit1), col = "red")
with(dasyprocta.day.cam03, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Dasyprocta, CT-RBG-1-03", ylab="Mean number of individuals / day"), )
  fit2 <- lm(dasyprocta.day.cam03$Number.of.Animals ~ dasyprocta.day.cam03$date) # fit can be used to extract slopes and use in the regression
  lines(dasyprocta.day.cam03$date, predict(fit2), col = "red")
with(dasyprocta.day.cam04, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, main="Dasyprocta, CT-RBG-1-04", ylab="Mean number of individuals / day"), )
  fit3 <- lm(dasyprocta.day.cam04$Number.of.Animals ~ dasyprocta.day.cam04$date) # fit can be used to extract slopes and use in the regression
  lines(dasyprocta.day.cam04$date, predict(fit3), col = "red")

# extract linear trend for each camera:
camera.names <- sort(unique(gurupi2017$Camera.Trap.Name))
trend <- rep(NA, length(camera.names))    
data.trends <- data.frame(camera.names,trend)
  # fill object with density estimates
  for(i in 1:nrow(data.trends))    # criando contador
  {
    lifepulse(subset(gurupi2017, Camera.Trap.Name == data.trends[i,"camera.names"]), "Dasyprocta prymnolopha", "4 day")
    data.temp <- min60[complete.cases(min60), ]
    fit <- lm(data.temp$Number.of.Animals ~ data.temp$date) # fit can be used to extract slopes and use in the regression
    data.trends[i,"trend"] <- fit[["coefficients"]][["data.temp$date"]]
  }

setwd("/home/elildojr/Documents/r/timedens")
write.csv(data.trends, "data.trends.csv")
test <- read.csv("data_trends.csv", header=T)

hist(test$trend)
with(test, plot(elevation, trend))
with(test, plot(distance.to.river, log(trend))
modtest <- with(test, glm(trend~distance.to.river+elevation))
summary(modtest)

lifepulse(subset(x, Camera.Trap.Name == "CT-TDM-1-08"), "Dasyprocta unknown", "hour")

lifepulse(subset(x, Camera.Trap.Name == "CT-TDM-1-08"), "Tapirus terrestris", "hour")

lifepulse(subset(x, Camera.Trap.Name == "CT-TDM-1-08"), "Crypturellus unknown", "hour")
pauxi <<- min60
with(pauxi, plot(date, Number.of.Animals, type="l", cex=.5, ylab="Mean number of individuals / hour"))

with(t.pecari, plot(date, Number.of.Animals, type="l", cex=.5, ylab="Mean number of individuals / hour"))


lifepulse(x, "Dasyprocta unknown", "day")
lifepulse(x, "Dasyprocta unknown", "12 hour")

lifepulse(x, "Tapirus terrestris", "day")
lifepulse(x, "Tapirus terrestris", "12 hour")

lifepulse(x, "Peccaries", "day")
lifepulse(x, "Peccaries", "12 hour")

lifepulse(x, "Crypturellus unknown", "day")
lifepulse(x, "Crypturellus unknown", "12 hour")

lifepulse(x, "Cuniculus paca", "day")

lifepulse(x, "Panthera onca", "day")


lifepulse(x, "Crax fasciolata", "week")
lifepulse(x, "Crax fasciolata", "day")
lifepulse(x, "Psophia viridis")

