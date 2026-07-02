# Phytoplankton forward PSM driver for discrete paleo reconstruction, mui calibration integrated in jags model
# Incorporates both modern observational data and GIG data to calibrate mui = f(radius, po4)
# Dustin T. Harper
############################################################################################

# Load libraries
############################################################################################
library(rjags)
library(R2jags)
library(readxl)

############################################################################################


# File paths
############################################################################################
cal.file <- "caldata_culture.csv"
size.file <- "HP07_sizedata.csv"
gig.file <- "timeseriesGIG.csv"
prox.file <- "phyto_Intermediate_combined.xlsx"
prox.sheet <- "data4PSM"

if (!file.exists(prox.file)) stop("Proxy intermediate file not found: ", prox.file)
############################################################################################


# Multi linear regression model for prior slopes and intercepts to calculate mu(i)
############################################################################################
# Read in instantaneous growth rate (mu,i) culture calibration data from Aloisi et al. (2015)
# with calculated mean radius from measured coccosphere using Henderiks and Pagani (2007) transfer functions
cal.df <- read.csv(cal.file)
cal.df <- subset(cal.df, !is.na(radius) & !is.na(po4) & !is.na(mui) & radius <= 10 & po4 <= 2)
cal.df$radius <- cal.df$radius*1e-6

# Generate multiple linear regression model mu,i (as a function of [PO4] and mean radius)
igr.model <- lm(mui ~ po4 + radius, data=cal.df)
igr.model.sum <- summary(igr.model)

# Load prior coefficients and SEs for mui = f(po4, rm) from multi linear regression
po4.co.lr <- igr.model.sum$coefficients[2,1]
po4.se.lr <- igr.model.sum$coefficients[2,2]
r.co.lr <- igr.model.sum$coefficients[3,1]
r.se.lr <- igr.model.sum$coefficients[3,2]
y.int.lr <- igr.model.sum$coefficients[1,1]
y.int.se.lr <- igr.model.sum$coefficients[1,2]
############################################################################################


# Determine size-based transfer functions using linear regression and Henderiks and Pagani 07 data
############################################################################################
hp07 <- read.csv(size.file)
lm.lith <- lm(lith.size ~ cell.r, data=hp07)
lith.sum <- summary(lm.lith)

lith.m <- lith.sum$coefficients[2,1]
lith.b <- lith.sum$coefficients[1,1]
############################################################################################


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


# Load GIG calibration data to evaluate against
############################################################################################
prox.in.gig <- read.csv(gig.file)
prox.in.gig <- prox.in.gig[,1:12]
names(prox.in.gig) <- c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data",
                        "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")
prox.in.gig <- prox.in.gig[complete.cases(prox.in.gig[,c("site", "age", "po4.prior", "d13Cmarker.data", "d13Cmarker.data.sd", "d13Cpf.data",
                                                         "d13Cpf.data.sd", "len.lith.data", "len.lith.data.sd", "Uk.data", "Uk.data.sd", "iceco2.data")]),]

# Site index GIG calibration data
prox.in.gig <- prox.in.gig[order(prox.in.gig$site),]
prox.in.gig <- transform(prox.in.gig, site.index=as.numeric(factor(site)))
site.index.gig <- prox.in.gig$site.index

# Set prior distributions for GIG calibration data

# Temperature (degrees C)
tempC.m.gig <- 22
tempC.p.gig <- 1/5^2

# Salinity (ppt)
sal.m.gig <- 35
sal.p.gig <- 1/0.5^2

# pCO2 (uatm)
pco2.m.gig <- 250
pco2.p.gig <- 1/40^2

# d13C of aqueous CO2 (per mille)
d13C.co2.m.gig <- -8
d13C.co2.p.gig <- 1/1^2

# Concentration of phosphate (PO4; umol/kg)
po4.m.gig <- unique(prox.in.gig$po4.prior)
po4.p.gig <- 1/0.25^2

# Mean cell radius (m)
rm.m.gig <- 1.5e-6
rm.p.gig <- 1/(0.5e-6)^2
############################################################################################


# Load discrete proxy data for reconstruction from intermediate workbook
############################################################################################
# State bins define the spatial/temporal resolution of each discrete inversion state
state.age.width <- 5
state.lat.width <- 20
state.lon.width <- 20
state.age.origin <- 0
state.lat.origin <- -90
state.lon.origin <- -180

prox.raw <- readxl::read_excel(prox.file, sheet=prox.sheet, skip=3, na=c("", "NA", "#N/A"))
prox.raw <- as.data.frame(prox.raw)
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

# Define discrete state bins using user-selected temporal and spatial integration windows
prox.in$age.bin <- floor((prox.in$age-state.age.origin)/state.age.width)
prox.in$lat.bin <- floor((prox.in$lat-state.lat.origin)/state.lat.width)
prox.in$lon.bin <- floor((prox.in$lon-state.lon.origin)/state.lon.width)
prox.in$state.id <- paste(prox.in$age.bin, prox.in$lat.bin, prox.in$lon.bin, sep="_")
prox.in$state.index <- as.numeric(factor(prox.in$state.id))
prox.in <- prox.in[order(prox.in$state.index, prox.in$age, prox.in$lat, prox.in$lon),]

state.levels <- levels(factor(prox.in$state.id))
n.state <- length(state.levels)
n.obs <- nrow(prox.in)

# Build compact observation vectors. Marker observations are required; lith observations are used only where present.
marker.keep <- which(!is.na(prox.in$d13Cmarker.data) & !is.na(prox.in$d13Cmarker.data.sd) & prox.in$d13Cmarker.data.sd > 0)
lith.keep <- which(!is.na(prox.in$len.lith.data) & !is.na(prox.in$len.lith.data.sd) & prox.in$len.lith.data.sd > 0)

n.marker <- length(marker.keep)
n.lith <- length(lith.keep)

if (n.marker < 1) stop("No d13C marker observations available after filtering")

# If there are no lith observations, pass one effectively uninformative row so JAGS can compile.
if (n.lith < 1){
  n.lith <- 1
  lith.data <- lith.m*1.5 + lith.b
  lith.data.sd <- 1e6
  lith.state.index <- 1
} else {
  lith.data <- prox.in$len.lith.data[lith.keep]
  lith.data.sd <- prox.in$len.lith.data.sd[lith.keep]
  lith.state.index <- prox.in$state.index[lith.keep]
}

# Collapse environmental priors to state level, while retaining every observation in the likelihood.
temp.state <- state_prior(prox.in$tempC.prior, prox.in$tempC.prior.sd, prox.in$state.index, n.state, min.sd=1)
po4.state <- state_prior(prox.in$po4.prior, prox.in$po4.prior.sd, prox.in$state.index, n.state, min.sd=0.05)

state.df <- aggregate(cbind(age, lat, lon) ~ state.index + state.id + age.bin + lat.bin + lon.bin,
                      data=prox.in, FUN=mean)
state.df <- state.df[order(state.df$state.index),]
state.df$n.obs <- as.numeric(table(prox.in$state.index)[as.character(state.df$state.index)])

# Dimensions
n.cd <- nrow(cal.df)
n.gig <- nrow(prox.in.gig)
n.prox <- n.state

# Set prior distributions for discrete proxy states
# Temperature (degrees C)
tempC.m <- temp.state$mean
tempC.p <- 1/temp.state$sd^2

# Salinity (ppt)
sal.m <- rep(35, n.state)
sal.p <- rep(1/2^2, n.state)

# pCO2 (uatm)
pco2.m <- rep(450, n.state)
pco2.p <- rep(1/250^2, n.state)
pco2.lb <- rep(50, n.state)
pco2.ub <- rep(1500, n.state)

# d13C of aqueous CO2 (per mille)
d13C.co2.m <- rep(-8, n.state)
d13C.co2.p <- rep(1/1^2, n.state)

# Concentration of phosphate (PO4; umol/kg)
po4.m.cd <- 1.5
po4.m <- po4.state$mean
po4.p <- 1/po4.state$sd^2

# Mean cell radius (m)
rm.m <- rep(1.5e-6, n.state)
rm.p <- rep(1/(0.5e-6)^2, n.state)


############################################################################################


# Select data to pass to jags
############################################################################################
data.pass <- list("n.cd" = n.cd,
                  "n.gig" = n.gig,
                  "n.prox" = n.prox,
                  "n.state" = n.state,
                  "n.marker" = n.marker,
                  "n.lith" = n.lith,
                  "n.temp" = n.temp,
                  "n.sal" = n.sal,
                  "po4.co.lr" = po4.co.lr,
                  "po4.se.lr" = po4.se.lr,
                  "r.co.lr" = r.co.lr,
                  "r.se.lr" = r.se.lr,
                  "y.int.lr" = y.int.lr,
                  "y.int.se.lr" = y.int.se.lr,
                  "lith.m" = lith.m,
                  "lith.b" = lith.b,
                  "radius.cd" = cal.df$radius,
                  "po4.cd.data" = cal.df$po4,
                  "mui.cd.data" = cal.df$mui,
                  "K0a" = K0a,
                  "Ksw_sta" = Ksw_sta,
                  "sal.lb" = sal.lb,
                  "tempC.lb" = tempC.lb,
                  "t.inc" = t.inc,
                  "s.inc" = s.inc,
                  "d13Cmarker.data.gig" = prox.in.gig$d13Cmarker.data,
                  "d13Cmarker.data.sd.gig" = prox.in.gig$d13Cmarker.data.sd,
                  "d13Cpf.data.gig" = prox.in.gig$d13Cpf.data,
                  "d13Cpf.data.sd.gig" = prox.in.gig$d13Cpf.data.sd,
                  "len.lith.data.gig" = prox.in.gig$len.lith.data,
                  "len.lith.data.sd.gig" = prox.in.gig$len.lith.data.sd,
                  "Uk.data.gig" = prox.in.gig$Uk.data,
                  "Uk.data.sd.gig" = prox.in.gig$Uk.data.sd,
                  "iceco2.data.gig" = prox.in.gig$iceco2.data,
                  "site.index.gig" = site.index.gig,
                  "tempC.m.gig" = tempC.m.gig,
                  "tempC.p.gig" = tempC.p.gig,
                  "sal.m.gig" = sal.m.gig,
                  "sal.p.gig" = sal.p.gig,
                  "pco2.m.gig" = pco2.m.gig,
                  "pco2.p.gig" = pco2.p.gig,
                  "d13C.co2.m.gig" = d13C.co2.m.gig,
                  "d13C.co2.p.gig" = d13C.co2.p.gig,
                  "po4.m.gig" = po4.m.gig,
                  "po4.p.gig" = po4.p.gig,
                  "rm.m.gig" = rm.m.gig,
                  "rm.p.gig" = rm.p.gig,
                  "d13Cmarker.data" = prox.in$d13Cmarker.data[marker.keep],
                  "d13Cmarker.data.sd" = prox.in$d13Cmarker.data.sd[marker.keep],
                  "marker.state.index" = prox.in$state.index[marker.keep],
                  "len.lith.data" = lith.data,
                  "len.lith.data.sd" = lith.data.sd,
                  "lith.state.index" = lith.state.index,
                  "tempC.m" = tempC.m,
                  "tempC.p" = tempC.p,
                  "sal.m" = sal.m,
                  "sal.p" = sal.p,
                  "pco2.m" = pco2.m,
                  "pco2.p" = pco2.p,
                  "pco2.lb" = pco2.lb,
                  "pco2.ub" = pco2.ub,
                  "d13C.co2.m" = d13C.co2.m,
                  "d13C.co2.p" = d13C.co2.p,
                  "po4.m" = po4.m,
                  "po4.m.cd" = po4.m.cd,
                  "po4.p" = po4.p,
                  "rm.m" = rm.m,
                  "rm.p" = rm.p)
############################################################################################


# Parameters to save as output
############################################################################################
parms <- c("tempC", "sal", "pco2", "pco2.gig", "d13C.co2", "po4", "rm", "b", "b.gig",
           "eps.p", "eps.p.gig", "d13Cmarker", "len.lith", "coeff.po4", "coeff.rm",
           "mui.y.int", "sigma.mui.cd", "eps.f", "eps.d", "eps.bob", "P.c")
############################################################################################


# Run the inversion using jags
############################################################################################
inv.out <- jags.parallel(data=data.pass, model.file="discrete_phyto_psm.R", parameters.to.save=parms,
                         inits=NULL, n.chains=3, n.iter=1e3,
                         n.burnin=5e2, n.thin=1)
############################################################################################


View(inv.out$BUGSoutput$summary)

