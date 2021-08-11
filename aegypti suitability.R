library('seegSDM')
library('ncdf4')

set.seed(207)

maxTSI <- 1269295.1
fmaxTSI <- 1276854.5

city_names <- c(
  'ahmadabad',
  'bangalore',
  'bangkok',
  'beijing',
  'bogota',
  'buenos_aires',
  'cairo',
  'chengdu',
  'chennai',
  'chongqing',
  'dar_es_salaam',
  'delhi',
  'dhaka',
  'guangzhou',
  'ho_chi_minh',
  'hyderabad',
  'istanbul',
  'jakarta',
  'karachi',
  'kinshasa',
  'kolkata',
  'lagos',
  'lahore',
  'lima',
  'london',
  'los_angeles',
  'luanda',
  'manila',
  'mexico_city',
  'moscow',
  'mumbai',
  'nanjing',
  'new_york',
  'osaka',
  'paris',
  'rio_de_janeiro',
  'seoul',
  'shanghai',
  'tehran',
  'tianjin',
  'tokyo'
)

aegypti <- read.csv("./aegypti.csv")

cleanData <- function(city,
                      period = 'Present',
                      na = 1, np = 0.1, mu = 100000) {
  ncFile <- paste('./covariates_cities/covariates_', city, '.nc', sep='')
  if (period == 'Future') {
    ncFile <- paste('./fcovariates_cities/fcovariates_', city, '.nc', sep='')
  }

  df <- nc_open(ncFile)
  
  Tmean <- ncvar_get(df, 'Tmean')
  TSI <- ncvar_get(df, 'TSI')
  TSI <- pmax(TSI, 0) / max(maxTSI, fmaxTSI)
  rainmax <- ncvar_get(df, 'RAINMAX')
  rainmin <- ncvar_get(df, 'RAINMIN')
  rh <- ncvar_get(df, 'RH_MIN')
  urb <- ncvar_get(df, 'MOD_FRC_URB2D')
  long <- ncvar_get(df, 'XLONG')
  lat <- ncvar_get(df, 'XLAT')
  lu <- ncvar_get(df, 'LU_INDEX')

  template <- createTemplate(as.vector((lu)), dim(TSI)[2], dim(TSI)[1], xmn = min(long), xmx = max(long), ymn = min(lat), ymx = max(lat))
  
  tsi <- createCov(as.vector(TSI), template)
  Tmean <- createCov(as.vector(Tmean), template)
  rainmax <- createCov(as.vector(rainmax), template)
  rainmin <- createCov(as.vector(rainmin), template)
  rh <- createCov(as.vector(rh), template)
  urb <- createCov(as.vector(urb), template)

  tsi <- flip(tsi, direction = 'y')
  Tmean <- flip(Tmean, direction = 'y')
  rainmax <- flip(rainmax, direction = 'y')
  rainmin <- flip(rainmin, direction = 'y')
  rh <- flip(rh, direction = 'y')
  urb <- flip(urb, direction = 'y')

  covs <- createCovs(Tmean, tsi, rainmax, rainmin, rh, urb)
  occ <- dplyr::filter(aegypti, Longitude < max(long),
                Longitude >= min(long),
                Latitude <= max(lat),
                Latitude >= min(lat))

  threshold <- -25
  consensus <- tsi
  consensus <- consensus * 200 - 100
  abs_consensus <- (1 - (consensus + 100) / 200)
  values(abs_consensus) <- pmax(values(abs_consensus), 0)
  pres_consensus <- calc(consensus,
                         fun = function(x) {
                           ifelse(x <= threshold,
                                  0,
                                  (x + 100) / 200)
                         })

  if (nrow(occ) == 0) {
    p_pres <- matrix(ncol = 2, nrow = 0)
    p_abs <- p_pres
    colnames(p_abs) <- c('Longitude', 'Latitude')
    colnames(p_pres) <- c('Longitude', 'Latitude')
    p_abs <- data.frame(p_abs)
    p_pres <- data.frame(p_pres)
  }
  else {
    p_abs <- createBackground(consensus = abs_consensus, occ = occ, n = ceiling(nrow(occ)*na), mu = mu)
    p_pres <- tryCatch(createBackground(consensus = pres_consensus, occ = occ, n = ceiling(nrow(occ)*np), mu = mu),
             error = function(cond) {
               p_pres <- matrix(ncol = 2, nrow = 0)
               colnames(p_pres) <- c('Longitude', 'Latitude')
               p_pres <- data.frame(p_pres)
               return(p_pres)
             }
             )
  }

  data <- createData(occ, p_abs, p_pres, covs)

  nc_close(df)
  return(c(covs, data))
}

createTemplate <- function(lu, nrow, ncol, xmn, xmx, ymn, ymx) {
  r <- raster(nrow = nrow, ncol = ncol, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  values(r) <- lu
  r[r[] == 17] <- NA
  r <- r * 0
  return(r)
}

createBackground <- function(consensus, occ, n, mu) {
  background <- bgDistance(n = n,
                           distance = mu,
                           points = SpatialPoints(occ[, c(3, 4)], consensus@crs),
                           raster = consensus,
                           prob = TRUE,
                           replace = TRUE,
                           spatial = FALSE)
  
  colnames(background) <- c('Longitude', 'Latitude')
  background <- data.frame(background)
  return(background)
}

createCov <- function(cov, template) {
  r <- template
  values(r) <- cov
  r <- mask(r, template)
  return(r)
}

createCovs <- function(Tmean, tsi, max_rain, min_rain, rh, urb) {
  b <- brick(Tmean, tsi, max_rain, min_rain, rh, urb)
  names(b) <- c('Tmean', 'tsi', 'max_rain', 'min_rain', 'rh', 'urb')
  return(b)
}

createData <- function(occurrence, p_abs, p_pres, covariates) {
  # combine the occurrence and background records
  dat <- rbind(cbind(PA = rep(1, nrow(occurrence)),
                     type = rep('occurrence', nrow(occurrence)),
                     occurrence[, c('Longitude', 'Latitude')]),
               cbind(PA = rep(1, nrow(p_pres)),
                     type = rep('presence', nrow(p_pres)),
                     p_pres[, c('Longitude', 'Latitude')]),
               cbind(PA = rep(0, nrow(p_abs)),
                     type = rep('absence', nrow(p_abs)),
                     p_abs[, c('Longitude', 'Latitude')])
               )
  
  # extract covariate values for each data point
  dat_covs <- extract(covariates, dat[, c('Longitude', 'Latitude')])
  
  # combine covariates with the other info
  dat_all <- cbind(dat, dat_covs)
  
  return(dat_all)
}

plot_model <- function(city, case) {
  df_present <- nc_open(paste('./covariates_cities/covariates_', city, '.nc', sep=''))
  df_future <- nc_open(paste('./fcovariates_cities/fcovariates_', city, '.nc', sep=''))
  
  long <- ncvar_get(df_present, 'XLONG')
  lat <- ncvar_get(df_present, 'XLAT')
  plu <- ncvar_get(df_present, 'LU_INDEX')
  flu <- ncvar_get(df_future, 'LU_INDEX')
  
  ptemplate <- createTemplate(as.vector(plu), dim(long)[2], dim(long)[1], xmn = min(long), xmx = max(long), ymn = min(lat), ymx = max(lat))
  ftemplate <- createTemplate(as.vector(flu), dim(long)[2], dim(long)[1], xmn = min(long), xmx = max(long), ymn = min(lat), ymx = max(lat))

  ptemplate <- flip(ptemplate)
  ftemplate <- flip(ftemplate)
  
  pmean <- nc_open(paste('./', case, '/present/present_mean_', city, '.nc', sep = ''))
  fmean <- nc_open(paste('./', case, '/future/future_mean_', city, '.nc', sep = ''))

  present <- createCov(t(ncvar_get(pmean, 'mean')), ptemplate)
  future <- createCov(t(ncvar_get(fmean, 'mean')), ftemplate)
  
  color1 <- colorRampPalette(c("blue", "yellow", "red"))(255)
  color2 <- colorRampPalette(c("blue", "white", "red"))(255)
  
  png(paste('./', case, '/plots/', city, '.png', sep = ''), height = 640, width = 960)
  
  par(mfrow=c(2, 2), mar=c(3, 3, 3, 3))
  plot(present, zlim = c(0, 1), main='Present')
  
  plot(future, zlim = c(0, 1), main='Future')
  plot(future - present, zlim = c(-1, 1), main='Difference', col=color2)
  dev.off()
  
  nc_close(pmean)
  nc_close(fmean)
  nc_close(df_present)
  nc_close(df_future)
}

ensembleData <- function(par) {
  na <- par$na
  np <- par$np
  mu <- par$mu
  data <- data.frame()
  for (city in city_names) {
    city_data <- cleanData(city, 'Present', na, np, ndef, mu)
    
    city_data_frame <- data.frame(PA = city_data$PA,
                                  type = city_data$type,
                                  Longitude = city_data$Longitude,
                                  Latitude = city_data$Latitude,
                                  tsi = city_data$tsi,
                                  max_rain = city_data$max_rain,
                                  min_rain = city_data$min_rain,
                                  rh = city_data$rh)
    data <- rbind(data, city_data_frame)
  }
  return(data)
}

trainEnsembleModel <- function(parallel) {
  na <- c(1, 2, 4, 6)
  np <- c(0, 0.01, 0.05, 0.075, 0.1)
  mu <- c(150000, 170000, 200000)
  pars <- expand.grid(na = na, np = np, mu = mu)
  sub <- function(i, pars) {
    pars[i, ]
  }
  par_list <- lapply(1:nrow(pars), sub, pars)

  sfInit(cpus = 5, parallel = TRUE)
  sfLibrary(seegSDM)
  sfLibrary(ncdf4)
  sfLibrary(raster)

  sfExportAll()
  
  data_list <- sfLapply(par_list, ensembleData)
  model_list <- sfLapply(data_list, runBRT, gbm.x = 5:8, gbm.y = 1, gbm.coords = 2:3)
  stat_list <- sfLapply(model_list, getStats, pwd = FALSE)
  
  sfStop()
  
  return(c(data_list, model_list, stat_list))
}

covsChange <- function(city) {
  p <- cleanData(city, 'Present')
  f <- cleanData(city, 'Future')
  if (is.null(p) || is.null(f)) return(NULL)
  
  color <- colorRampPalette(c("blue", "white", "red"))(255)
  
  changeTmean <- f[[1]]$Tmean - p[[1]]$Tmean
  changetsi <- f[[1]]$tsi - p[[1]]$tsi
  changeMaxRain <- f[[1]]$max_rain - p[[1]]$max_rain
  changeMinRain <- f[[1]]$min_rain - p[[1]]$min_rain
  changeRH <- f[[1]]$rh - p[[1]]$rh
  changeUrb <- f[[1]]$urb - p[[1]]$urb
  
  getzlim <- function(r) {
    end <- max(abs(r@data@min), abs(r@data@max))
    return(c(-abs(end), abs(end)))
  }
  
  png(paste('./covsChange/covsChange_', city, '.png', sep = ''), height = 960, width = 1280)
  
  par(mfrow = c(3, 2))
  plot(changeTmean, zlim = getzlim(changeTmean), main = 'Change in Annual Mean Temperature (°C)', col = color)
  plot(changetsi, zlim = getzlim(changetsi), main = 'Change in Temperature Suitability Index', col = color)
  plot(changeMaxRain, zlim = getzlim(changeMaxRain), main = 'Change in Annual Monthly Maximum Total Rainfall (mm)', col = color)
  plot(changeMinRain, zlim = getzlim(changeMinRain), main = 'Change in Annual Monthly Minimum Total Rainfall (mm)', col = color)
  plot(changeRH, zlim = getzlim(changeRH), main = 'Change in Relative Humidity (%)', col = color)
  plot(changeUrb, zlim = getzlim(changeUrb), main = 'Change in Urbanization', col = color)
  
  dev.off()
}

genEnsemblePreds <- function(city) {
  print(city)
  present <- cleanData(city, 'Present')
  future <- cleanData(city, 'Future')
  
  p <- lapply(model_list, makePreds, present)
  f <- lapply(model_list, makePreds, future)
  
  p <- combinePreds(brick(p), parallel = TRUE, ncore = 5)
  f <- combinePreds(brick(f), parallel = TRUE, ncore = 5)
  
  p$uncertainty <- p[[4]]-p[[3]]
  f$uncertainty <- f[[4]]-f[[3]]
  
  writeRaster(p$mean, paste('./present_mean_', city, '.nc', sep = ''), format = 'CDF', overwrite = TRUE)
  writeRaster(p$uncertainty, paste('./present__uncertainty_', city, '.nc', sep = ''), format = 'CDF', overwrite = TRUE)
  writeRaster(f$mean, paste('./future_mean_', city, '.nc', sep = ''), format = 'CDF', overwrite = TRUE)
  writeRaster(f$uncertainty, paste('./future_uncertainty', city, '.nc', sep = ''), format = 'CDF', overwrite = TRUE)
}