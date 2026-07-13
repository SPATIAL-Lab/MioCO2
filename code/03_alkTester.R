library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Proxy file ----
prox.file <- "data/proxyData/phyto_Intermediate_combined.xlsx"
prox.sheet <- "data4PSM"
size.file <- "data/HP07_sizedata.csv"

# Generate look up tables for equilibrium constants
############################################################################################
# Set upper and lower STP bounds for equil constant array
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
############################################################################################

# Determine size-based transfer functions using linear regression and Henderiks and Pagani 07 data
############################################################################################
hp07 <- read.csv(size.file)
lm.lith <- lm(lith.size ~ cell.r, data=hp07)
lith.sum <- summary(lm.lith)

lith.m <- lith.sum$coefficients[2,1]
lith.b <- lith.sum$coefficients[1,1]
############################################################################################

# Load discrete proxy data for reconstruction from intermediate workbook
############################################################################################
prox.raw <- read.xlsx(prox.file, sheet = prox.sheet, startRow = 4, 
                      na.strings = c("", "NA", "#N/A"))
names(prox.raw) <- trimws(names(prox.raw))

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

state_prior <- function(x, x.sd, state.index, n.state, min.sd=NULL){
  x.m <- rep(NA, n.state)
  x.s <- rep(NA, n.state)
  
  for (j in 1:n.state){
    ii <- which(state.index == j & !is.na(x) & !is.na(x.sd) & x.sd > 0)
    w <- 1/x.sd[ii]^2
    x.m[j] <- sum(x[ii]*w)/sum(w)
    x.s[j] <- sqrt(1/sum(w))
  }
  
  if (!is.null(min.sd)) x.s <- pmax(x.s, min.sd)
  list(mean=x.m, sd=x.s)
}

# Pull out the intermediate-sheet variables needed by the discrete PSM
sample.id <- as.character(prox.raw$sample)
sample.id <- ifelse(is.na(sample.id) | sample.id == "", paste0("sample_", seq_along(sample.id)), sample.id)
doi <- as.character(prox.raw$doi)
doi <- ifelse(is.na(doi), "", doi)
lat <- num(prox.raw$lat)
lon <- num(prox.raw$lon)
lon <- ((lon + 180) %% 360) - 180
age.mean <- mean_or_range(prox.raw$age_mean, prox.raw$age_min, prox.raw$age_max)

d13Cmarker.mean <- mean_or_range(prox.raw$d13Corg_mean, prox.raw$d13Corg_min, prox.raw$d13Corg_max)
d13Cmarker.sd <- sd_or_range(prox.raw$d13Corg_2s, prox.raw$d13Corg_min, prox.raw$d13Corg_max, 1)

tempC.mean <- mean_or_range(prox.raw$tempC, prox.raw$tempC_min, prox.raw$tempC_max, 25)
tempC.sd <- sd_or_range(prox.raw$tempC_2s, prox.raw$tempC_min, prox.raw$tempC_max, 5)

po4.mean <- mean_or_range(prox.raw$po4_mean, prox.raw$po4_min, prox.raw$po4_max, 1.5)
po4.sd <- sd_or_range(prox.raw$po4_2s, prox.raw$po4_min, prox.raw$po4_max, 0.25)

lith.mean <- mean_or_range(prox.raw$lith_mean, prox.raw$lith_min, prox.raw$lith_max)
lith.sd <- sd_or_range(prox.raw$lith_2s, prox.raw$lith_min, prox.raw$lith_max, 0.5)

prox.in <- data.frame(sample=sample.id,
                      doi=doi,
                      lat=lat,
                      lon=lon,
                      age=age.mean,
                      po4.prior=po4.mean,
                      po4.prior.sd=po4.sd,
                      tempC.prior=tempC.mean,
                      tempC.prior.sd=tempC.sd,
                      d13Cmarker.data=d13Cmarker.mean,
                      d13Cmarker.data.sd=d13Cmarker.sd,
                      len.lith.data=lith.mean,
                      len.lith.data.sd=lith.sd,
                      stringsAsFactors=FALSE)

prox.in <- prox.in[complete.cases(prox.in[,c("sample", "doi", "lat", "lon", "age", "po4.prior", "po4.prior.sd",
                                             "tempC.prior", "tempC.prior.sd", "d13Cmarker.data",
                                             "d13Cmarker.data.sd")]),]
prox.in = prox.in[prox.in$age < 25 & prox.in$age > 5, ]
n.obs <- nrow(prox.in)

## Add d13Ca
d13Ca = read.csv("data/d13Ca_Cenozoic.csv")
d13Ca$sd = (d13Ca$d13Ca_97p5 - d13Ca$d13Ca_2p5) / 2

## Proxy data input
prox.in$d13Ca.obs = data.frame("d13Ca.m" = d13Ca$d13Ca_50[round(prox.in$age * 10, 0)],
                              "d13Ca.sd" = d13Ca$sd[round(prox.in$age * 10, 0)])

prox.in$temp.obs = data.frame("temp.m" = prox.in$tempC.prior,
                          "temp.sd" = prox.in$tempC.prior.sd)

prox.in$po4.obs = data.frame("po4.m" = prox.in$po4.prior,
                             "po4.sd" = prox.in$po4.prior.sd)

prox.in$d13Cmarker.obs = data.frame("d13Cmarker.m" = prox.in$d13Cmarker.data, 
                                    "d13Cmarker.sd" = prox.in$d13Cmarker.data.sd)

# Build compact observation vectors. Marker observations are required; lith observations are used only where present.
lith.keep <- which(!is.na(prox.in$len.lith.data) & !is.na(prox.in$len.lith.data.sd) & prox.in$len.lith.data.sd > 0)

n.lith <- length(lith.keep)

# If there are no lith observations, pass one effectively uninformative row so JAGS can compile.
if (n.lith < 1){
  n.lith <- 1
  lith.data <- lith.m*1.5 + lith.b
  lith.data.sd <- 1e6
  lith.state.index <- 1
} else {
  lith.data <- prox.in$len.lith.data[lith.keep]
  lith.data.sd <- prox.in$len.lith.data.sd[lith.keep]
  lith.state.index <- lith.keep
}

# Collapse environmental priors to state level, while retaining every observation in the likelihood.
#temp.state <- state_prior(prox.in$tempC.prior, prox.in$tempC.prior.sd, prox.in$state.index, n.state, min.sd=1)
#po4.state <- state_prior(prox.in$po4.prior, prox.in$po4.prior.sd, prox.in$state.index, n.state, min.sd=0.05)

#state.df <- aggregate(cbind(age, lat, lon) ~ state.index + state.id + age.bin + lat.bin + lon.bin,
#                      data=prox.in, FUN=mean)
#state.df <- state.df[order(state.df$state.index),]
#state.df$n.obs <- as.numeric(table(prox.in$state.index)[as.character(state.df$state.index)])

############################################################################################


# Select data to pass to jags
############################################################################################
data.pass <- list("n.obs" = n.obs,
                  "n.lith" = n.lith,
                  "n.temp" = n.temp,
                  "n.sal" = n.sal,
                  "lith.m" = lith.m,
                  "lith.b" = lith.b,
                  "K0a" = K0a,
                  "Ksw_sta" = Ksw_sta,
                  "sal.lb" = sal.lb,
                  "tempC.lb" = tempC.lb,
                  "t.inc" = t.inc,
                  "s.inc" = s.inc,
                  "d13Cmarker.obs" = prox.in$d13Cmarker.obs,
                  "d13Ca.obs" = prox.in$d13Ca.obs,
                  "temp.obs" = prox.in$temp.obs,
                  "po4.obs" = prox.in$po4.obs,
                  "len.lith.data" = lith.data,
                  "len.lith.data.sd" = lith.data.sd,
                  "lith.state.index" = lith.state.index)
############################################################################################


# Parameters to save as output
############################################################################################
parms <- c("tempC", "sal", "pCO2", "d13C.co2", "d13Ca", "po4", "rm", "b",
           "eps.p", "d13Cmarker", "len.lith", "coeff.po4", "coeff.rm")
############################################################################################


# Run the inversion using jags
############################################################################################
post.phyto <- jags.parallel(data = data.pass, model.file = "code/models/phyto.R", 
                         parameters.to.save = parms, inits = NULL, n.chains = 3, 
                         n.iter = 2e5, n.burnin = 1e5, n.thin = 100)
############################################################################################

View(post.phyto$BUGSoutput$summary)

plot(prox.in$age, post.phyto$BUGSoutput$median$pCO2)
plot(data.pass$d13Ca.obs[, 1], post.phyto$BUGSoutput$median$d13C.co2)
plot(data.pass$d13Cmarker.obs[, 1] - data.pass$d13Ca.obs[ai, 1], post.phyto$BUGSoutput$median$pCO2)
plot(post.phyto$BUGSoutput$median$rm, post.phyto$BUGSoutput$median$pCO2)

# Timeseries inversion

## Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

## d13Ca
d13Ca.obs = data.frame("d13Ca" = numeric(length(ages)),
                       "d13Ca.sd" = numeric(length(ages)))
for(i in seq_along(ages)){
  d13Ca.obs[i, 1] = d13Ca$d13Ca_50[which.min(abs(ages[i] - d13Ca$age))]
  d13Ca.obs[i, 2] = d13Ca$sd[which.min(abs(ages[i] - d13Ca$age))]
}

## Age index
ai = numeric(nrow(prox.in))
for(i in seq_along(prox.in$age)){
  ai[i] = which.min(abs(prox.in$age[i] - ages))
}

## Site index
sites = unique(prox.in$lat + prox.in$lon)
si = match(prox.in$lat + prox.in$lon, sites)

# Select data to pass to jags
############################################################################################
data.pass <- list("n.obs" = n.obs,
                  "n.lith" = n.lith,
                  "n.steps" = length(ages),
                  "n.sites" = length(sites),
                  "n.temp" = n.temp,
                  "n.sal" = n.sal,
                  "ai" = ai,
                  "dt" = stepsize,
                  "si" = si,
                  "lith.m" = lith.m,
                  "lith.b" = lith.b,
                  "K0a" = K0a,
                  "Ksw_sta" = Ksw_sta,
                  "sal.lb" = sal.lb,
                  "tempC.lb" = tempC.lb,
                  "t.inc" = t.inc,
                  "s.inc" = s.inc,
                  "d13Cmarker.obs" = data.frame(prox.in$d13Cmarker.data, prox.in$d13Cmarker.data.sd),
                  "d13Ca.obs" = d13Ca.obs,
                  "temp.obs" = prox.in$temp.obs,
                  "po4.obs" = prox.in$po4.obs,
                  "len.lith.data" = lith.data,
                  "len.lith.data.sd" = lith.data.sd,
                  "lith.state.index" = lith.state.index)
############################################################################################

# Run the inversion using jags
############################################################################################
post.phyto.ts <- jags.parallel(data = data.pass, model.file = "code/models/phyto_ts.R", 
                         parameters.to.save = parms, inits = NULL, n.chains = 3, 
                         n.iter = 1e4, n.burnin = 5e3, n.thin = 1)
############################################################################################

View(post.phyto.ts$BUGSoutput$summary)

tsplot(ages, post.phyto.ts, "pCO2")
pointplot(prox.in$age, post.phyto, "pCO2")

plot(ages, post.phyto.ts$BUGSoutput$median$rm[1,] , ylim = c(1, 5), type = "l")
for(i in 2:n.sites){
  lines(ages, post.phyto.ts$BUGSoutput$median$rm[i,])
}

plot(ages, post.phyto.ts$BUGSoutput$median$po4[1,] , ylim = c(0.1, 2), type = "l")
for(i in 2:n.sites){
  lines(ages, post.phyto.ts$BUGSoutput$median$po4[i,])
}

plot(data.pass$d13Cmarker.obs[, 1] - data.pass$d13Ca.obs[ai, 1], post.phyto.ts$BUGSoutput$median$pCO2[ai])
plot(data.pass$d13Cmarker.data - data.pass$d13Ca.obs[ai, 1], post.phyto$BUGSoutput$median$pCO2)

save(post.phyto, post.phyto.ts, file = "bigout/phyto.rda")
