options(repr.plot.width=4, repr.plot.height=3)

# Function to plot maps
mapPlot <- function(x, loc, title = "", lims = NULL) {
  if (is.null(lims)) {
    lims <- c(min(x), max(x))
  }
  
  dataToPlot = data.frame(X = loc[, 1], Y = loc[, 2], Value = x)
  mapPoints <- ggmap(map) +
    geom_point(
      data = dataToPlot,
      aes(x = X, y = Y, color = Value),
      alpha = 0.9,
      # size = 0.05
      size = 1.5,
      shape = 15
    ) +
    scale_color_gradientn(colors = rev(rainbow(7)), limits = lims) +
    labs(title = title, fill = "") +
    guides(color = guide_colorbar(barwidth = 1)) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
  print(mapPoints)
  
}

# Function to estimate the empirical extremogram between two locations
computeExtremalCoeff <- function(pos1, pos2, database, local_quantiles) {
  theta <- length(which((database[, pos1] > local_quantiles[pos1] & !is.na(database[,pos1])) &
                          (database[, pos2] > local_quantiles[pos2] & !is.na(database[,pos2]))
  )) /
    length(which((database[, pos1] > local_quantiles[pos1] & !is.na(database[, pos1]))))
  return(theta)
}

computeTemporalExtremo <- function(timeLag, dataBaseRescaled, timesDists, timePairs) {
  theta <- length(which((dataBaseRescaled[timePairs[timesDists == timeLag,1],] > 0) & (dataBaseRescaled[timePairs[timesDists == timeLag,2],] > 0))) /
    length(which((dataBaseRescaled[timePairs[timesDists == timeLag,1],] > 0)))
  return(theta)
}

fit_spatial_vario <- function(theta, estimates){
  # print(theta)
  lambda = (theta[3])
  beta = (theta[2])
  alpha = theta[1]
  # nugget = theta[6]
  nugget = 0
  
  angle <- theta[4]
  scaling <- theta[5]
  
  if(alpha < 0 | alpha >= 2 | beta >= 2 | lambda < 0 | nugget < 0){return(1e50)}
  if(angle < -pi/2 | angle >= pi/2 | scaling <= 0){return(1e50)}
  H <- matrix(c(cos(angle), -sin(angle), scaling * sin(angle), scaling * cos(angle)), 2, 2)  
  
  windShiftedLocations <- matrix(0, length(estimates$dist) ,2)
  windShiftedLocations[,1] <- estimates$dist * cos(estimates$angle * pi / 180)
  windShiftedLocations[,2] <- estimates$dist * sin(estimates$angle * pi / 180) 
  
  distXY <- sqrt((windShiftedLocations[,1])^2 + (windShiftedLocations[,2])^2)
  
  windShiftedLocationsW <- windShiftedLocations
  windShiftedLocationsW[,1] <- (H[1,1] * windShiftedLocations[,1] + H[1,2] * windShiftedLocations[,2])
  windShiftedLocationsW[,2] <- (H[2,1] * windShiftedLocations[,1] + H[2,2] * windShiftedLocations[,2])
  
  distXY_anisotropic <- sqrt((windShiftedLocationsW[,1])^2 + (windShiftedLocationsW[,2])^2)
  gamma <- rep(0, length(distXY_anisotropic))
  gamma[distXY_anisotropic > 0] <- nugget + ((1 + (sqrt(distXY_anisotropic[distXY_anisotropic > 0]^2 / lambda))^alpha )^(beta / alpha) - 1) / (2^(beta / alpha) - 1)
  # gamma[distXY_anisotropic > 0] <- nugget + beta * (1 - exp((-sqrt(distXY_anisotropic[distXY_anisotropic > 0]^2 / lambda))^alpha) )
  
  extremogram <- 2 * (1 - pnorm(sqrt(gamma / 2)))
  
  mls <- sum(((extremogram - estimates$extremogram)^2)[distXY < 1000])
  # mls <- sum((gamma[distXY < 1000] -  estimates$vario[distXY < 1000])^2 )
  
  return(mls)
}


objective_function_temporal <- function(theta, mean_extremo){
  lambda = (theta[2])
  beta = (theta[1])
  alpha = theta[3]
  
  distXTemp <- (0:7)
  
  gamma <- rep(0, length(distXTemp))
  gamma <- ((1 + (sqrt(distXTemp^2 / lambda))^alpha )^(beta / alpha) - 1) / (2^(beta / alpha) - 1)
  
  extremogram <- 2 * (1 - pnorm(sqrt(gamma / 2)))
  
  mls <- sum((extremogram -  mean_extremo)^2)
  
  return(mls)
}

