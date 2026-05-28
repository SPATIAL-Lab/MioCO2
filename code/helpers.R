parseFranks = function(d, condense = TRUE, fixA = FALSE){
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
    data[[i + 16]] = d.sub
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
  
  # Condense sites and taxa
  if(condense){
    ## Sites index
    ci = match("ai", names(d))
    strats = unique(d[, ci])
    data$level = match(d[, ci], strats)
    
    ## Condense site level parameters
    ### First occurrence of each strat level
    fo = match(strats, d[, ci])
    data$d13Ca = data$d13Ca[fo, ]
    
    ## Taxa index
    ### Collect unknowns and assign dummy names to them
    d$genus[d$genus == "unknown"] = NA
    d$genus[is.na(d$genus)] = seq(sum(is.na(d$genus)))
    
    genera = unique(d$genus)
    data$genus = match(d$genus, genera)
    
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
    data$genus = data$level = seq_len(nrow(d))
  }
  
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
