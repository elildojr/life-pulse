# timedens
# Function to estimate animal densities withouth individual recognition
# by Elildo Carvalho Jr @ ICMBio/CENAP

timedens <- function(x,y,z,w) # x = data, y = species, effective detection distance, angle
{
  # create empty object with 1000 rows:
  dens <- rep(NA, 1000)
  
  subset.species <- subset(x, bin == y) # subset species
  time.window <- as.numeric((max(subset.species$td.photo)-min(subset.species$td.photo))*24*60*60) # total length of survey in seconds
  area.sampled <- 60*(pi*(z^2)*w/360)*1e-06 # area sampled by 60 cameras in km2. For individual cameras it is given by the circular sector with radius z and angle w
  animal.records <-as.numeric(subset.species$Number.of.Animals) # for groups the number in EACH picture
  empty.records <- rep(0,time.window-length(animal.records))
  state.space <- c(animal.records, empty.records) # the collection of all observed densities
  # I could have estimated densities for every second of survey in sequence, but it was easier to do as in the three previous lines
  
  # THIS IS THE FORMULA FOR DENSITY AS A TIME-AVERAGE:
  # density <- sum(state.space)/(area.sampled*time.window)
  # Instead of using the formula above I sample the state.space to to obtain confidence intervals. Perhaps there is an easier way, this is a first draft:
  
  # fill object dens with density estimates
  for(i in 1:length(dens))    # criando contador
  {
    sample.state.space <- sample(state.space, 3600) # samples of 3600 seconds (equivalent to 1-hour)
    density.i <- sum(sample.state.space)/(area.sampled*3600) # observed density in sample
    dens[i] <- density.i
  }
  
  # Graphical outputs
  
  par(mfrow=c(1,2))
  
  # Plot the collection of observed density estimates
  xrange <- 1:length(dens)
  plot(xrange, dens, type="n", ylim=c(0,(max(dens)*1.2)), xlab="State-space sample", ylab="Density estimate" , main = "Observed density estimates", cex=0.5)
  lines(xrange, dens, col = "steelblue", lwd=0.5)
  
  # density plot
  dens.dens <- density(dens)
  plot(dens.dens, col = "steelblue", main = "PDF population density estimates", cex=0.5)
  
  # Return mean and confidence intervals
  meandens <- mean(dens)
  ci.dens <- round(quantile(dens, probs=c(0.025, 0.975), na.rm=T), 4) 
  ret <- list("Time-averaged density (ind/km2)" = meandens, "95% CI" = ci.dens)
  return(ret)
}

