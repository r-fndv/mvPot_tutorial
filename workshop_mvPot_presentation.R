# Load dataset
load("DATA_TRAINING_SUBSET.RData")

# Load set of utility functions (e.g. plotting maps)
source("utils.R")

# Number of core tu use for computation (not necessary to have more than one to reproduce example here)
nCores <- 1

# Laod package for parallelisation and create cluster instance
if(nCores > 1){
  library(parallel)
  cl <- makeCluster(nCores,type="FORK")
  clusterSetRNGStream(cl)
} else {
  cl <- NULL
}

# Set to true if you want to run heavy computations (tractable on laptop but too long for the workshop)
run_heavy_computations <- FALSE

# Load package for manipulation of coordinates
library(sp)

# Create coordinates object
coords <- as.data.frame(loc.sub)
names(coords) <-  c("X", "Y")
coordinates(coords) <- c("X", "Y")

# Define system of coordinates
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Project lon/lat coordinates to UTM (m) system of coordinates
projected_coordinates <-  spTransform(coords, CRS("+proj=utm +zone=37 ellps=WGS84"))

# Rescale to have new coordinates in km
projected_coordinates@coords <-  projected_coordinates@coords / 1000

#Load ggplot2 for plot
library(ggplot2)

# Set to true if you want to produce maps
map_plotting <- FALSE

if(map_plotting == TRUE){
  # Load package for map plotting
  library(ggmap)
  
  # Specify here your onwn google_api_key
  api <- readLines("google_api_key")
  register_google(key = api)
  
  # Downlaod map
  map <- get_map(location = c(mean(loc.sub[,1]),mean(loc.sub[,2])), maptype = "toner-2010", zoom = 7)
  
  #Plot one sample from the dataset
  mapPlot(anom.training.sub[1,], loc.sub, title = "Temperature Anomaly")
}

####################################
######## DATA EXPLORATION ##########
####################################

# Data set summary
print("Data set summary:")
summary(as.vector(anom.training.sub))

# Compute local mean and 0.95quantile
quantile_anomaly <- rep(NA, length(anom.training.sub[1,]))
mean_anomaly <- rep(NA, length(anom.training.sub[1,]))

# Set quantile level
quantileLevel <- 0.95
for(i in 1:length(quantile_anomaly)){
  quantile_anomaly[i] <- quantile(anom.training.sub[,i], quantileLevel, na.rm = TRUE)
  mean_anomaly[i] <- mean(anom.training.sub[,i], na.rm = TRUE)
}

print("Local mean:")
summary(mean_anomaly)
print("Local empirical 0.95 quantile:")
summary(quantile_anomaly)

#Plot results
if(map_plotting == TRUE){
  mapPlot(mean_anomaly, loc.sub, title = "Local Mean")
}

if(map_plotting == TRUE){
  mapPlot(quantile_anomaly, loc.sub, title = "Local 0.95 quantile")
}

####################################
######## TAIL MARGINAL MODEL #######
####################################

# Parameters for generalized Pareto distributions
scales <- rep(NA, length(anom.training.sub[1,]))
shapes <- rep(NA, length(anom.training.sub[1,]))

#List of fitted objects to check results
list_fit <- list()

# Local GP fit above the 0.95 quantile
for(i in 1:length(scales)){
  print(paste(i, "/", length(scales)))
  
  # Likelihood based GPD fit
  fitGP <- evd::fpot(
    anom.training.sub[, i],
    threshold = quantile_anomaly[i],
    std.err = FALSE,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  
  # Store estimation results
  scales[i] <- fitGP$estimate[1]
  shapes[i] <- fitGP$estimate[2]
  list_fit[[i]] <- fitGP
}

#Check results
summary(shapes)
summary(scales)

# Location to plot 
posistion_to_plot <- 500

#Plot fit quality
plot(list_fit[[posistion_to_plot]], which = 2)

# Plot maximum likelhood tail parameter estimates
if(map_plotting == TRUE){
  mapPlot(shapes, loc.sub, title = "Estimated tail indexes (xi)")
}

# Plot maximum likelhood scale estimates
if(map_plotting == TRUE){
  mapPlot(scales, loc.sub, title = "Estimated scale parameters (sigma)")
}

####################################
###### MARGINAL NORMALIZATION  #####
####################################

# Create object for normalized data
normalized_database <- matrix(NA, nrow = nrow(anom.training.sub), ncol = ncol(anom.training.sub))

# Standarize to unit Frechet or generalized Pareto with unit tail and scale and zero location
frechet = FALSE

for(i in 1:length(scales)){
  # Compute local empirical CDF
  empiricalCdf <- ecdf(anom.training.sub[,i])
  
  if(frechet == TRUE){
    # Use empirical cdf below the threshold
    normalized_database[(anom.training.sub[,i] <= quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] <- -1 / log(empiricalCdf(anom.training.sub[(anom.training.sub[,i] <= quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i]))
    # Use estimated GP distribution above the threshold
    normalized_database[(anom.training.sub[,i] > quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] <-
      -1 / log(1 - (1 - quantileLevel) * (1 + shapes[i]*(anom.training.sub[(anom.training.sub[,i] > quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] - quantile_anomaly[i]) / scales[i])^(-1 / shapes[i]))
  } else {
    # Use empirical cdf below the threshold
    normalized_database[(anom.training.sub[,i] <= quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] <- 1 / (1 - empiricalCdf(anom.training.sub[(anom.training.sub[,i] <= quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i])) - 1
    # Use estimated GP distribution above the threshold
    normalized_database[(anom.training.sub[,i] > quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] <- 
      1 /  ((1 - quantileLevel) * (1 + shapes[i]*(anom.training.sub[(anom.training.sub[,i] > quantile_anomaly[i] & !is.na(anom.training.sub[,i])),i] - quantile_anomaly[i]) / scales[i])^(-1 / shapes[i])) - 1
  }
}

#Plot one sample
if(map_plotting == TRUE){
  sample_position <- 1
  mapPlot(normalized_database[sample_position,], loc.sub)
}

#####################################################
###### ESTIMATE DEPENDENCE CENSORED LIKELIHOOD  #####
#####################################################

# Define threshold for censoring on rescaled data
if(frechet == TRUE){
  rescaled_threshold <- evd::qgev(0.95,1,1,1)
} else {
  rescaled_threshold <- evd::qgpd(0.95,0,1,1)
}

#Set seed for subset selection; full censored likelihood is really computationnally heavy and
# a vector of size 1254 is not tractable on a laptop
set.seed("294")

# Define subset size
subset_location_size <- 25
# Random subset selection
index_for_estimation <- sample(size = subset_location_size, x = 1:length(quantile_anomaly), replace = FALSE)

# Subset of locations
locations_for_estimation <- projected_coordinates@coords[index_for_estimation,]

# Subset of data
database_for_estimation <- normalized_database[,index_for_estimation]

# Plot one sample of the database subset
if(map_plotting == TRUE){
  sample_position <- 500
  mapPlot(database_for_estimation[sample_position,],loc.sub[index_for_estimation,],
          lims = c(min(normalized_database[sample_position,], na.rm = TRUE), max(normalized_database[sample_position,], na.rm = TRUE)),
          title = "Location subset")
}

# Find samples without missing data
observations_without_na <- apply(database_for_estimation,1, function(x){any(is.na(x))})

# Subset database to keep observations without missing data
database_for_estimation <- database_for_estimation[!observations_without_na,]

# Find observations for wich at least one component is exceeding a threshold
database_excesses <- database_for_estimation[apply(database_for_estimation, 1, function(x){any(x > rescaled_threshold)}),]

# Transform matrix of observation to list
list_database_excesses <- lapply(1:length(database_excesses[,1]), function(i){database_excesses[i,]})

# Define variogram model - Here the model by Schlather and Moreva
vario <- function(h, alpha = 1, beta = 1, lambda = 1, A = matrix(c(1,0,0,1),2,2)){
  if(norm(h,type = "2") > 0){
    ((1 + (sqrt(norm(A %*%h,type = "2")^2 / lambda))^alpha )^(beta / alpha) - 1) / (2^(beta / alpha) - 1)
  } else {0}
}

# Create objective function to be minimized
objectiveFunctionCensored = function(parameter, excesses, loc, vario, rescaled_threshold, nCores, cl){
  # Print parameters to follow the evolution of the optimition
  print(parameter)
  
  # Enforce boundaries on the variogram parameters
  if(parameter[1] < 0 | parameter[1] > 2 | parameter[2] > 2 | parameter[3] < 0){return(1e50)}
  if(parameter[4] < -pi/2 | parameter[4] > pi/2 | parameter[5] < 0){return(1e50)}
  
  # Set the variogram parameters
  A <- matrix(c(cos(parameter[4]), -sin(parameter[4]), parameter[5] * sin(parameter[4]),parameter[5] * cos(parameter[4])), 2, 2)
  variogram_model <- function(h){
    vario(h, alpha = parameter[1], beta = parameter[2], lambda = parameter[3], A = A)
  }
  
  # Compute consored likelihood from mvPot package
  mvPot::censoredLikelihoodBR(excesses, loc, variogram_model, rescaled_threshold, p = 4999, vec = mvPot::genVecQMC(4999, length(loc[,1]))$genVec, nCores = nCores, cl = cl)
}

if(run_heavy_computations == TRUE){
  # Set seed for quasi Monte Carlo estimate and `ensure' computation time
  set.seed(56892734)
  
  # Optimize the objective function to estimate variogram parameters using censored likelihood takes about 3 hours on one core 
  ref_time <- proc.time()
  estimates_censored_likelihood <- optim(par = c(1.7726726, -4.7176393, 9325.3082619, 1.2818892, 0.7136902),
                                         fn = objectiveFunctionCensored,
                                         excesses = list_database_excesses,
                                         loc = as.data.frame(locations_for_estimation),
                                         vario = vario,
                                         nCores = nCores,
                                         cl = cl,
                                         rescaled_threshold = rep(rescaled_threshold, subset_location_size),
                                         control = list(maxit = 100000, trace = 3),
                                         method = "Nelder-Mead")
  final_time <- proc.time() - ref_time
  
  # Save results
  save(estimates_censored_likelihood, file = "estimated_dependence.RData")
} else {
  load("estimated_dependence.RData")
}

########################################
############ MODEL CHECKING  ###########
########################################

# Compute all possible pairwise combinations
pairs_combination_indexes <- expand.grid(1:subset_location_size,1:subset_location_size)

# Create vector of empirical extremogram
empirical_extremogram <- rep(NA, length(pairs_combination_indexes[,1]))

# Compute empirical extremogram for each pair of locations
for(k in 1:length(empirical_extremogram)){
  empirical_extremogram[k] <- computeExtremalCoeff(pairs_combination_indexes[k,1], pairs_combination_indexes[k,2], anom.training.sub[,index_for_estimation], quantile_anomaly[index_for_estimation])
}

# Compute vector of distance between each pairs
distances_between_pairs <- sqrt((projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,1]],1] - projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,2]],1])^2 +
                                  (projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,1]],2] - projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,2]],2])^2)

# Compute orientation of each pairs
matrixOfAngles <- atan((projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,1]],2] - projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,2]],2]) /
                         (projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,1]],1] - projected_coordinates@coords[index_for_estimation[pairs_combination_indexes[,2]],1])) / pi * 180
matrixOfAngles[is.nan(matrixOfAngles)] <- 0

# Store results in data frame for easier plotting
estimates <- data.frame(extremogram = empirical_extremogram[!is.nan(empirical_extremogram)], dist = distances_between_pairs[!is.nan(empirical_extremogram)], angle = matrixOfAngles[!is.nan(empirical_extremogram)],
                        x = distances_between_pairs[!is.nan(empirical_extremogram)] * cos(matrixOfAngles[!is.nan(empirical_extremogram)]), y = distances_between_pairs[!is.nan(empirical_extremogram)] * sin(matrixOfAngles[!is.nan(empirical_extremogram)]))

# Plot empirical extremogram as function of distance
p_Extremogram_Dist <- ggplot(data = estimates, aes(x = dist, y = extremogram, color = angle)) + geom_point(alpha = 0.5) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7))) + ylim(0,1)
print(p_Extremogram_Dist)

# Plot empirical extremogram as a map
p_Extremogram_Map <- ggplot(data = estimates, aes(x = x, y = y, color = extremogram)) + geom_point(alpha = 0.5) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7)), limits = c(0,1))
print(p_Extremogram_Map)

# Define variogram with estimated parameters
variogram_model_censored <- function(h){
  # Anisotrpy matrix
  A <- matrix(c(cos(estimates_censored_likelihood$par[4]), -sin(estimates_censored_likelihood$par[4]), estimates_censored_likelihood$par[5] * sin(estimates_censored_likelihood$par[4]),
                estimates_censored_likelihood$par[5] * cos(estimates_censored_likelihood$par[4])), 2, 2)
  #Variogram function
  vario(h, alpha = estimates_censored_likelihood$par[1], beta = estimates_censored_likelihood$par[2], lambda = estimates_censored_likelihood$par[3], A = A)
}

# Compute theoretical extremogram from variogram
extremogram_model_censored <- rep(0,length(distances_between_pairs))
for(i in 1:length(distances_between_pairs)){
  # Vector of coordinates difference between pairs
  coordinates_difference <- c(distances_between_pairs[i] * cos(matrixOfAngles[i]), distances_between_pairs[i] * sin(matrixOfAngles[i]))
  # Compute extremogram between pairs
  extremogram_model_censored[i] <- 2 * (1 - pnorm(sqrt(variogram_model_censored(coordinates_difference) / 2)))
}

# Plot empirical variogram cloud along with estimated model
pExtremogram <- ggplot(data = estimates, aes(x = dist, y = extremogram, color = angle)) + geom_point(alpha = 0.2) + geom_smooth(se = FALSE) + 
  geom_line(data = data.frame(extremogram = extremogram_model_censored, dist = distances_between_pairs, angle = matrixOfAngles), mapping =  aes(x = dist, y = extremogram, color = angle)) +
  guides(color = guide_colorbar(barheight = 30, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7))) + ylim(0,1)
print(pExtremogram)

################################################
###### ESTIMATE DEPENDENCE GRADIENT SCORE  #####
################################################

# Defined weighting function
weighting_function <- function(x, u){
  x * (1 - exp(-(sum(x) / u - 1)))
}

# Partial derivative of weighting function
weighting_function_derivative <- function(x, u){
  (1 - exp(-(sum(x) / u - 1))) + (x / (u)) * exp( - (sum(x) / u - 1))
}

# Define objective function for gradient score
objective_function_gradient_score = function(parameter, excesses, locations_coordinates, vario, weighting_function, weighting_function_derivative, threshold, nCores, cl = NULL){
  
  # Print parameter to follow optimization evolution
  print(parameter)
  
  #Define the variogram corresponding to input parameters
  if(parameter[1] < 0 | parameter[1] > 2 | parameter[2] > 2 | parameter[3] < 0){return(1e50)}
  if(parameter[4] < -pi/2 | parameter[4] > pi/2 | parameter[5] < 0){return(1e50)}
  A <- matrix(c(cos(parameter[4]), -sin(parameter[4]), parameter[5] * sin(parameter[4]),parameter[5] * cos(parameter[4])), 2, 2)
  variogram_model <- function(h){
    vario(h, alpha = parameter[1], beta = parameter[2], lambda = parameter[3], A = A)
  }
  
  #Compute score
  mvPot::scoreEstimation(excesses, locations_coordinates, variogram_model, weighting_function, weighting_function_derivative, u = threshold, nCores = nCores, cl = cl)
}

# Compute spatial accumulation for each rescaled observations
sums <- apply(database_for_estimation,1, sum)

#Defined threshold for exceedance - We chose 0.9 quantile to obtain about the same number of observations as for censored likelihood
threshold <- quantile(sums, 0.9)

# Create matrix of exceedances for spatial accumulation
spatial_accumulation_excesses <- database_for_estimation[sums > threshold,]
list_spatial_accumulation_excesses <- lapply(1:length(spatial_accumulation_excesses[,1]), function(i){spatial_accumulation_excesses[i,]})

#Estimate the parameter - Quick to run on single core
estimates_gradient_score <- optim(par = c(1.9686318,   -2.5746466, 3917.9238753,    0.5421206,    0.5714497),
                                  fn = objective_function_gradient_score,
                                  excesses = list_spatial_accumulation_excesses,
                                  locations_coordinates = as.data.frame(locations_for_estimation),
                                  vario = vario,
                                  weighting_function = weighting_function,
                                  weighting_function_derivative = weighting_function_derivative,
                                  threshold = threshold,
                                  nCores = nCores,
                                  cl = cl,
                                  control = list(maxit = 10000, trace = 1),
                                  method = "Nelder-Mead")


save(estimates_gradient_score, file = "gradient_score_estimate.RData")

# Define variogram with estimated parameters
variogram_model_score <- function(h){
  # Anisotrpy matrix
  A <- matrix(c(cos(estimates_gradient_score$par[4]), -sin(estimates_gradient_score$par[4]), estimates_gradient_score$par[5] * sin(estimates_gradient_score$par[4]),
                estimates_gradient_score$par[5] * cos(estimates_gradient_score$par[4])), 2, 2)
  # Variogram function
  vario(h, alpha = estimates_gradient_score$par[1], beta = estimates_gradient_score$par[2], lambda = estimates_gradient_score$par[3], A = A) 
}

# Compute theoretical extremogram from variogram
extremogram_model_score <- rep(0,length(distances_between_pairs))
for(i in 1:length(distances_between_pairs)){
  # Vector of coordinates difference between pairs
  coordinates_difference <- c(distances_between_pairs[i] * cos(matrixOfAngles[i]), distances_between_pairs[i] * sin(matrixOfAngles[i]))
  # Compute extremogram between pairs
  extremogram_model_score[i] <- 2 * (1 - pnorm(sqrt(variogram_model_score(coordinates_difference) / 2)))
}

# Plot empirical variogram cloud along with estimated model
pExtremogram <- ggplot(data = estimates, aes(x = dist, y = extremogram, color = angle)) + geom_point(alpha = 0.2) + geom_smooth(se = FALSE) + 
  geom_line(data = data.frame(extremogram = extremogram_model_score, dist = distances_between_pairs, angle = matrixOfAngles), mapping =  aes(x = dist, y = extremogram, color = angle)) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7))) + ylim(0,1)
print(pExtremogram)

############################################
###### SIMULATIONS FROM FITTED MODELS  #####
############################################

# Produce simulation with same seed for comparison
if(run_heavy_computations == TRUE){
  set.seed(1093)
  simulation_censored <- mvPot::simulPareto(n = 1, loc = as.data.frame(projected_coordinates@coords), vario = variogram_model_censored, nCores = 1)
  set.seed(1093)
  simulation_score <- mvPot::simulPareto(n = 1, loc = as.data.frame(projected_coordinates@coords), vario = variogram_model_score, nCores = 1)
  save(simulation_censored, simulation_score, file = "simulations.RData")
} else {
  load("simulations.RData")
}

# Marginal back-transform
if(frechet == TRUE){
  simulation_censored_uniform <- evd::pgev(simulation_censored[[1]], 1,1,1)
  simulation_score_uniform <- evd::pgev(simulation_score[[1]], 1,1,1)
} else{
  simulation_censored_uniform <- evd::pgpd(simulation_censored[[1]], 0,1,1)
  simulation_score_uniform <- evd::pgpd(simulation_score[[1]], 0,1,1)
}

simulation_censored_rescaled <- rep(NA, length(simulation_censored[[1]]))
simulation_score_rescaled <- rep(NA, length(simulation_score[[1]]))

for(i in 1:length(simulation_score_rescaled)){
  if(simulation_censored_uniform[i] < 0.95){simulation_censored_rescaled[i] <- quantile(anom.training.sub[,i], probs =  simulation_censored_uniform[i], na.rm = TRUE)} else {
    simulation_censored_rescaled[i] <- evd::qgpd((simulation_censored_uniform[i] - 0.95) / 0.05, loc = quantile_anomaly[i], scale = scales[i], shape = shapes[i])
  }
  if(simulation_score_uniform[i] < 0.95){simulation_score_rescaled[i] <- quantile(anom.training.sub[,i], probs =  simulation_score_uniform[i], na.rm = TRUE)} else {
    simulation_score_rescaled[i] <- evd::qgpd((simulation_score_uniform[i] - 0.95) / 0.05, loc = quantile_anomaly[i], scale = scales[i], shape = shapes[i])
  }
}

print("Summary censored likelihood")
summary(simulation_censored_rescaled)
print("Summary gradient scoring rule")
summary(simulation_score_rescaled)

# Plot simulations
plot_limits <- c(min(c(simulation_censored_rescaled, simulation_score_rescaled)), max(c(simulation_censored_rescaled, simulation_score_rescaled)))
if(map_plotting == TRUE){
  mapPlot(simulation_censored_rescaled, loc.sub, title = "Censored likelihood", lims = plot_limits)
}

if(map_plotting == TRUE){
  mapPlot(simulation_score_rescaled, loc.sub, title = "Gradient Scoring Rule", lims = plot_limits)
}
