
# Function lifepulse
# by Elildo Carvalho Jr @ ICMBio/CENAP

lifepulse <- function(x,y,z) # x = data, y = species, z = avg.time
{
  # subset species and time window (only times with all cameras operating simultaneously)
  subset.species <- subset(x, bin == y)
  #subset.species <- subset(subset.species, Photo.Date >= max(subset.species$Camera.Start.Date) & Photo.Date <= min(subset.species$Camera.End.Date))
  #discrete.times <- (as.numeric(min(subset.species$Camera.End.Date)-max(subset.species$Camera.Start.Date)))*24*60*60
  #*length(unique(subset.time.window$Camera.Trap.Name))
  
  subset.species$date <- subset.species$td.photo # create "date" column because openair demands it
  subset.species <- subset.species[,c("date","Number.of.Animals")] # just to clean dataframe
  df1 <- data.frame(table(subset.species$date, subset.species$Number.of.Animals))
  df1 <- df1[,1:2]
  names(df1) <- c("date", "Number.of.Animals")
  df1$date <- as.POSIXct(df1$date)
  df1$Number.of.Animals <- as.numeric(df1$Number.of.Animals)
  
  # Improving: filling empty cases in time series
  #join this with another data frame that has all the time intervals you are expecting.
  
  #First, construct data frame using a difference of time intervals (1 sec in this case):
  #dif <- difftime(min(x$Camera.End.Date), max(x$Camera.Start.Date), units="secs")
  #timespan <- seq(1:dif)
  dif <- difftime(min(x$td.photo), max(x$td.photo), units="secs")
  timespan <- seq(1:dif)
  df <- as.POSIXlt(timespan, origin = max(x$Camera.Start.Date))
  df <- data.frame(df,NA)
  names(df) <- c("date", "Number.of.Animals")
  df$Number.of.Animals <- as.numeric(df$Number.of.Animals)
  head(df)
  
  #Join the two:
  result <- dplyr::full_join(df, df1)
  result[is.na(result)] <- 0
  result <<- result
  
  library(openair)
  min60 <-  timeAverage(result, avg.time = z, statistic = "mean")
  min60[is.na(min60)] <- 0
  min60$Number.of.Animals[min60$Number.of.Animals > 1] <- NA # removing crazy high values, CHECK THEM LATER
  with(min60, plot(date, Number.of.Animals, type="l", pch = 16, col="black", cex=.5, ylab="Mean number of individuals / hour"))
  #timePlot(result, pollutant = "Number.of.Animals")
  #min60[is.na(min60)] <- 0 # for filling the data.trends object, perhaps later change for complete.cases in fit
  min60 <<- min60
}  
