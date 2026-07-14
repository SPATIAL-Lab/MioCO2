parseFranks = function(d, ages, condense = TRUE, fixA = FALSE){
  # Parse values from sheet, obs, parameters, constants

  # Drop aggregate estimate rows based on missing stomatal density
  d = d[!is.na(d$Dab), ]
  
  data = list()
  
  # Data obs 
  data.names = c("d13Cp", "Dab", "GCLab", "GCWab", "Dad", "GCLad", "GCWad")
  data.sd = c(0.2, 1.5e5, 3e-7, 1.5e-7, 1.5e5, 3e-7, 1.5e-7)
  
  for(i in seq_along(data.names)){
    ci = match(data.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = data.sd[i]
    d.sub[d.sub[, 2] == 0, 2] = data.sd[i]
    data[[i]] = d.sub
  }
  
  # Free parameters
  mp.names = c("A0", "CiCa0", "gb", "s1", "s2", "s3",
               "s4", "s5")
  mp.sd = c(0.25, 0.05, 0.05, 0.001, 0.05, 0.01, 0.01, 0.001)
  
  for(i in seq_along(mp.names)){
    ci = match(mp.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = mp.sd[i]
    if(mp.names[i] == "s1"){
      # Special case if Pl is measured directly 
      d.sub[d.sub[, 2] == 0, 2] = 0.001
    } else{
      d.sub[d.sub[, 2] == 0, 2] = mp.sd[i]
    }
    data[[i + 7]] = d.sub
  }
  
  # Fixed parameters
  c.names = c("CO2_0", "PL_GCL")
  for(i in seq_along(c.names)){
    ci = match(c.names[i], names(d))
    if(is.na(ci) & c.names[i] == "PL_GCL"){
      d.sub = rep(0.5, nrow(d))
    } else{
      d.sub = d[, ci]
    }
    data[[i + 15]] = d.sub
  }
  
  names(data) = c(data.names, mp.names, c.names)
  
  ## s3 values of 1 are common but impossible, substitute
  data$s3[data$s3[, 1] == 1, 1] = 0.97
  
  # Abaxial adaxial indicies
  data$ind.ab = which(data$Dad$Dad == 0)
  data$ind.ad = which(data$Dad$Dad != 0)
  
  # Transform CiCa0
  ## Temporary insertion
  data$CiCa0$CiCa0[is.na(data$CiCa0$CiCa0)] = 0.65
  data$CO2_0[is.na(data$CO2_0)] = 380
  
  data$Ci0 = data$CiCa0 * data$CO2_0
  names(data$Ci0) = c("Ci0", "eCi0")
  data = data[!(names(data) %in% c("CiCa0", "CO2_0"))]
  
  # Age index
  data$ai.plant = rep(0)
  for(i in seq_along(d$age_mean)){
    data$ai.plant[i] = which.min(abs(d$age_mean[i] - ages))
  }
  
  # Condense sites and taxa
  if(condense){
    ## Taxa index
    ### Collect unknowns and assign dummy names to them
    d$genus[d$genus == "unknown"] = NA
    d$genus[is.na(d$genus)] = seq(sum(is.na(d$genus)))
    
    genera = unique(d$genus)
    data$gi = match(d$genus, genera)
    
    ## Condense species level parameters
    ### First occurrence of each species
    fo = match(genera, d$genus)
    data$gb = data$gb[fo, ]
    if(fixA == FALSE){
      data$A0 = data$A0[fo, ]
    }
    data$Ci0 = data$Ci0[fo, ]
    data$s5 = data$s5[fo, ]
    data$s4 = data$s4[fo, ]
    data$s3 = data$s3[fo, ]
    data$s2 = data$s2[fo, ]
  } else{
    data$gi = data$ai.plant = seq_len(nrow(d))
  }
  
  data$n.gen = length(unique(data$gi))
  
  return(data)
}

inits = function() {
  list("ca.s" = runif(length(d13Ca[, 1]), 0.5, 2))  
}  

prepMod = function(data){
  # Read base model
  basemod = readLines("code/models/forwardFranksMultiAbAd.R")
  
  # Find lines
  ad.fl = grep("# Adaxial species", basemod)
  ab.fl = grep("# Abaxial species", basemod)
  ad.ll = ab.fl - 1
  ab.ll = grep("# Taxon priors", basemod)
  
  # Remove unneeded code
  if(length(data$ind.ab) == 0){
    mod = basemod[-c(ab.fl:ab.ll)]
  } else if(length(data$ind.ad) == 0){
    mod = basemod[-c(ad.fl:ad.ll)]
  } else{
    mod = basemod
  }
  
  # Write
  writeLines(mod, file.path(tempdir(), "model.txt"))
}

prepMod.fixA = function(data){
  # Read base model
  basemod = readLines("code/models/forwardFranksMultiAbAd_fixA.R")
  
  # Find lines
  ad.fl = grep("# Adaxial species", basemod)
  ab.fl = grep("# Abaxial species", basemod)
  ad.ll = ab.fl - 1
  ab.ll = grep("# Taxon priors", basemod)
  
  # Remove unneeded code
  if(length(data$ind.ab) == 0){
    mod = basemod[-c(ab.fl:ab.ll)]
  } else if(length(data$ind.ad) == 0){
    mod = basemod[-c(ad.fl:ad.ll)]
  } else{
    mod = basemod
  }
  
  # Write
  writeLines(mod, file.path(tempdir(), "model.txt"))
}


parseStan = function(d, condense = TRUE, fixA = FALSE){
  # Parse values from sheet, obs, parameters, constants
  
  # Drop aggregate estimate rows based on missing stomatal density
  d = d[!is.na(d$Dab), ]
  
  data = list()
  
  # Data obs 
  data.names = c("d13Cp", "Dab", "GCLab", "GCWab", "Dad", "GCLad", "GCWad")
  data.sd = c(0.2, 1.5e5, 3e-7, 1.5e-7, 1.5e5, 3e-7, 1.5e-7)
  
  for(i in seq_along(data.names)){
    ci = match(data.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = data.sd[i]
    d.sub[d.sub[, 2] == 0, 2] = data.sd[i]
    data[[i]] = as.numeric(d.sub[1, ])
  }
  
  # Free parameters
  mp.names = c("d13Ca", "A0", "CiCa0", "gb", "s1", "s2", "s3",
               "s4", "s5")
  mp.sd = c(0.5, 0.25, 0.05, 0.05, 0.001, 0.05, 0.01, 0.01, 0.001)
  
  for(i in seq_along(mp.names)){
    ci = match(mp.names[i], names(d))
    d.sub = d[, ci:(ci + 1)]
    d.sub[is.na(d.sub[, 2]), 2] = mp.sd[i]
    if(mp.names[i] == "s1"){
      # Special case if Pl is measured directly 
      d.sub[d.sub[, 2] == 0, 2] = 0.001
    } else{
      d.sub[d.sub[, 2] == 0, 2] = mp.sd[i]
    }
    data[[i + 7]] = as.numeric(d.sub[1, ])
  }
  
  # Fixed parameters
  c.names = c("CO2_0", "b", "gamma", "PL_GCL")
  for(i in seq_along(c.names)){
    ci = match(c.names[i], names(d))
    if(is.na(ci) & c.names[i] == "PL_GCL"){
      d.sum = rep(0.5, nrow(d))
    } else{
      d.sub = d[, ci]
    }
    data[[i + 16]] = d.sub
  }
  
  c.names[3] = "gam"
  names(data) = c(data.names, mp.names, c.names)
  
  # Transform CiCa0
  data$Ci0 = data$CiCa0 * data$CO2_0
  names(data$Ci0) = c("Ci0", "eCi0")
  data = data[!(names(data) %in% c("CiCa0", "CO2_0"))]
  
  return(data)
}

# Add time series to plot w 2 prob density envelopes
tsdens = function(d, col = "black"){
  # Check dimensions of d
  if(ncol(d) != 6){stop("d cols should be should be time, 5%, 25%, 50%, 75%, 95% CI")}
  
  base.rgb = col2rgb(col)
  cols = c(rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25), 
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25),
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 1))
  
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 2], rev(d[, 6])), 
          col = cols[1], border = NA)
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 3], rev(d[, 5])), 
          col = cols[2], border = NA)
  lines(d[, 1], d[, 4], col = cols[3], lwd = 2)
}

# Make probability envelope plot
tsplot = function(a, d, v, col = "black", ylim = NA, add = FALSE){
  if(!inherits(d, "rjags")){
    stop("d must be rjags object")
  }
  d = d$BUGSoutput$sims.list
  
  if(!(v %in% names(d))){
    stop("v not a valid parameter")
  }
  d = d[[v]]
  
  if(length(a) != ncol(d)){
    stop("dimension mismatch")
  }
  d = apply(d, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  d = cbind(a, t(d))
  
  if(!add){
    if(is.na(sum(ylim))){
      ylim = range(d[, 2:6])
    }
    
    xlim = range(d[, 1])
    plot(0, 0, type = "n", xlim = rev(xlim), ylim = ylim, xlab = "Age (Mya)",
         ylab = v)
  }
  
  tsdens(d, col)
}

# Add point estimates to plot
pointplot = function(a, d, v, col = "black"){
  if(!inherits(d, "rjags")){
    stop("d must be rjags object")
  }
  d = d$BUGSoutput$sims.list
  
  if(!(v %in% names(d))){
    stop("v not a valid parameter")
  }
  d = d[[v]]
  
  d = t(apply(d, 2, quantile, probs = c(0.025, 0.5, 0.975)))
  
  if(length(a) != nrow(d)){
    stop("dimension mismatch")
  }
  
  arrows(a, d[, 1], a, d[, 3], 0.05, 90, 3, col)
  points(a, d[, 2], pch = 21, col = col, bg = "white")
}

parseSoil = function(d, ages, condense = TRUE){
  # Parse values from sheet, obs, parameters, constants
  
  # Current data structure
  ## datum
  ## soil
  
  # Long-term solution
  ## datum
  ## soil
  ## age
  ## site
  
  data.pass = list()
  
  # d13Cc
  ci = match("d13C_cc", names(d))
  d.sub = d[, ci:(ci + 1)]
  d.sub[, 2] = d.sub[, 2] / 2
  data.pass$d13Cc.obs = d.sub

  # d13Co
  ## Need to accommodate other OM types
  ci = match("d13Com_occluded", names(d))
  d.sub = d[, ci:(ci + 1)]
  d.sub[, 2] = d.sub[, 2] / 2
  data.pass$d13Co.obs = d.sub
  
  # MAT
  ci = match("tempC", names(d))
  d.sub = d[, ci:(ci + 1)]
  d.sub[, 2] = d.sub[, 2] / 2
  data.pass$temp.obs.soil = d.sub
  
  # MAP
#  ci = match("d13C_cc", names(d))
#  d.sub = d[, ci:(ci + 1)]
#  d.sub[, 2] = d.sub[, 2] / 2
#  data$MAP.obs = d.sub

  # lat
  data.pass$lat.soil = d$lat
  
  # number of data
  data.pass$n.obs.soil = nrow(d)
  
  # Age index
  data.pass$ai.soil = rep(0)
  for(i in seq_along(d$age_mean)){
    data.pass$ai.soil[i] = which.min(abs(d$age_mean[i] - ages))
  }
  
  ## Sites index
  sites = unique(d$lat)
  data.pass$si.soil = rep(0)
  for(i in seq_along(sites)){
    data.pass$si.soil[d$lat == sites[i]] = i
  }
  data.pass$n.sites.soil = length(sites)
  
  # Condense sites, add condense ages
  if(condense){
    ## Condense site level parameters
    ### First occurrence of each strat level
    sites = unique(data.pass$si.soil)
    fo = match(sites, data.pass$si.soil)
    data.pass$lat.soil = d$lat[fo]
  }
  
  return(data.pass)
}

parsePhyto = function(d, ages){

  tempC.lb <- 0
  tempC.ub <- 65
  sal.lb <- 15
  sal.ub <- 60
  
  # Step increments for sal (ppt) temp (degrees C) and press (bar)
  t.inc <- 0.25
  s.inc <- 0.25
  
  # Ranges of variables over which to evaluate
  tempC.vr <- seq(tempC.lb, tempC.ub, by=t.inc)
  sal.vr <- seq(sal.lb, sal.ub, by=s.inc)
  
  # Initiate arrays
  temp.vr <- c(1:length(tempC.vr))
  base2Darray <- c(1:(length(tempC.vr)*length(sal.vr)))
  dim(base2Darray) <- c(length(tempC.vr), length(sal.vr))
  Ksw_sta <- base2Darray
  K0a <- base2Darray
  
  # Constant (cm^3 bar mol^-1 K^-1)
  R <- 83.131
  
  # Calculate 2D array for K0 and Ksw (temp and sal dependent)
  for (i in 1:length(tempC.vr)){
    for (j in 1:length(sal.vr)){
      temp.vr[i] <- tempC.vr[i]+273.15
      Ksw_sta[i,j] <- exp(148.96502-13847.26/temp.vr[i]-23.6521*(log(temp.vr[i]))+(118.67/temp.vr[i]-5.977+1.0495*(log(temp.vr[i])))*(sal.vr[j]^0.5)-0.01615*sal.vr[j])
      K0a[i,j] <- exp(9345.17/temp.vr[i]-60.2409+23.3585*(log(temp.vr[i]/100))+sal.vr[j]*(0.023517-0.00023656*temp.vr[i]+0.0047036*((temp.vr[i]/100)^2)))
    }
  }
  
  n.temp <- nrow(K0a)
  n.sal <- ncol(K0a)

  # Size-based transfer function - Henderiks and Pagani 07 data
  lith.m = 1.859568
  lith.b = 0.3544168

  num <- function(x) suppressWarnings(as.numeric(x))
  
  mean_or_range <- function(x, xmin=NULL, xmax=NULL, default=NA){
    x <- num(x)
    
    if (!is.null(xmin) & !is.null(xmax)){
      x.range <- (num(xmin)+num(xmax))/2
      x <- ifelse(!is.na(x), x, x.range)
    }
    
    ifelse(!is.na(x), x, default)
  }
  
  sd_or_range <- function(x2s, xmin=NULL, xmax=NULL, default){
    x <- num(x2s)/2
    
    if (!is.null(xmin) & !is.null(xmax)){
      x.range <- (num(xmax)-num(xmin))/4
      x <- ifelse(!is.na(x), x, x.range)
    }
    
    ifelse(!is.na(x) & x > 0, x, default)
  }

  # Pull out the intermediate-sheet variables needed by the discrete PSM
  sample.id <- as.character(d$sample)
  sample.id <- ifelse(is.na(sample.id) | sample.id == "", paste0("sample_", seq_along(sample.id)), sample.id)
  lat <- num(d$lat)
  lon <- num(d$lon)
  lon <- ((lon + 180) %% 360) - 180
  age.mean <- mean_or_range(d$age_mean, d$age_min, d$age_max)
  
  d13Cmarker.mean <- mean_or_range(d$d13Corg_mean, d$d13Corg_min, d$d13Corg_max)
  d13Cmarker.sd <- sd_or_range(d$d13Corg_2s, d$d13Corg_min, d$d13Corg_max, 1)
  
  tempC.mean <- mean_or_range(d$tempC, d$tempC_min, d$tempC_max, 25)
  tempC.sd <- sd_or_range(d$tempC_2s, d$tempC_min, d$tempC_max, 5)
  
  po4.mean <- mean_or_range(d$po4_mean, d$po4_min, d$po4_max, 1.5)
  po4.sd <- sd_or_range(d$po4_2s, d$po4_min, d$po4_max, 0.25)
  
  lith.mean <- mean_or_range(d$lith_mean, d$lith_min, d$lith_max)
  lith.sd <- sd_or_range(d$lith_2s, d$lith_min, d$lith_max, 0.5)
  
  prox.in <- data.frame(sample = sample.id,
                        lat = lat,
                        lon = lon,
                        age = age.mean,
                        po4.prior = po4.mean,
                        po4.prior.sd = po4.sd,
                        tempC.prior = tempC.mean,
                        tempC.prior.sd = tempC.sd,
                        d13Cmarker.data = d13Cmarker.mean,
                        d13Cmarker.data.sd = d13Cmarker.sd,
                        len.lith.data = lith.mean,
                        len.lith.data.sd = lith.sd,
                        stringsAsFactors = FALSE)
  
  prox.in <- prox.in[complete.cases(prox.in[,c("sample", "lat", "lon", "age", 
                                               "po4.prior", "po4.prior.sd",
                                               "tempC.prior", "tempC.prior.sd", 
                                               "d13Cmarker.data",
                                               "d13Cmarker.data.sd")]),]
  prox.in = prox.in[prox.in$age < 25 & prox.in$age > 5, ]
  n.obs <- nrow(prox.in)
  
  ## Proxy data input
  temp.obs = data.frame("temp.m" = prox.in$tempC.prior,
                        "temp.sd" = prox.in$tempC.prior.sd)
  
  po4.obs = data.frame("po4.m" = prox.in$po4.prior,
                       "po4.sd" = prox.in$po4.prior.sd)
  
  d13Cmarker.obs = data.frame("d13Cmarker.m" = prox.in$d13Cmarker.data, 
                              "d13Cmarker.sd" = prox.in$d13Cmarker.data.sd)
  
  # Build compact observation vectors. Marker observations are required; lith observations are used only where present.
  lith.keep <- which(!is.na(prox.in$len.lith.data) & !is.na(prox.in$len.lith.data.sd) & 
                       prox.in$len.lith.data.sd > 0)
  
  n.lith <- length(lith.keep)
  
  # If there are no lith observations, pass one effectively uninformative row so JAGS can compile.
  if (n.lith < 1){
    n.lith <- 1
    lith.obs = data.frame("lith.obs" = lith.m * 1.5 + lith.b,
                          "lith.obs.sd" = 1e6)
    li = 1
  } else{
    lith.obs = data.frame("lith.obs" = prox.in$len.lith.data[lith.keep],
                          "lith.obs.sd" = prox.in$len.lith.data.sd[lith.keep])
    li = lith.keep
  }

  # Parameters to save as output
  parms <- c("tempC", "sal", "pCO2", "d13C.co2", "d13Ca", "po4", "rm", "b",
             "eps.p", "d13Cmarker", "len.lith", "coeff.po4", "coeff.rm")
  
  ## Age index
  ai = numeric(nrow(prox.in))
  for(i in seq_along(prox.in$age)){
    ai[i] = which.min(abs(prox.in$age[i] - ages))
  }
  
  ## Site index
  sites = unique(prox.in$lat + prox.in$lon)
  si = match(prox.in$lat + prox.in$lon, sites)
  
  # Select data to pass to jags
  data.pass <- list("n.obs.alk" = n.obs,
                    "n.lith" = n.lith,
                    "n.sites.alk" = length(sites),
                    "n.temp" = n.temp,
                    "n.sal" = n.sal,
                    "ai.alk" = ai,
                    "si.alk" = si,
                    "lith.m" = lith.m,
                    "lith.b" = lith.b,
                    "K0a" = K0a,
                    "Ksw_sta" = Ksw_sta,
                    "sal.lb" = sal.lb,
                    "tempC.lb" = tempC.lb,
                    "t.inc" = t.inc,
                    "s.inc" = s.inc,
                    "d13C.obs.alk" = d13Cmarker.obs,
                    "temp.obs.alk" = temp.obs,
                    "po4.obs" = po4.obs,
                    "lith.obs" = lith.obs,
                    "li" = li)
  
  return(data.pass)
}

parseAtmos = function(d, ages){
  
  d13Ca.obs = data.frame("d13Ca" = numeric(length(ages)),
                         "d13Ca.sd" = numeric(length(ages)))
  for(i in seq_along(ages)){
    d13Ca.obs[i, 1] = d$d13Ca_50[which.min(abs(ages[i] - d$age))]
    d13Ca.obs[i, 2] = d$sd[which.min(abs(ages[i] - d$age))]
  }
  
  return(list("d13Ca.obs" = d13Ca.obs))
}