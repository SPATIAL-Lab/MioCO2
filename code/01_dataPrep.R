library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Get and parse proxy file names ----
## All files
pf = list.files("data/proxyData/", full.names = TRUE)

## Plants
pf.p = pf[grep("stomata", pf)]

## Soils
pf.s = pf[grep("paleosol", pf)]

## Foraminifera
pf.f = pf[grep("boron", pf)]

## Alkenones
pf.a = pf[grep("phyto", pf)]

# Compile plant data ----
## First file
pd.p = read.xlsx(pf.p[1], 1, startRow = 4)

## Others
if(length(pf.p) > 1){
  for(i in pf.p[2:length(pf.p)]){
    pd.p.a = read.xlsx(i, 1, startRow = 4)
    pd.p = rbind(pd.p, pd.p.a)
  }
}

## Get d13Ca data
d13Ca = read.csv("data/d13Ca_Cenozoic.csv")
d13Ca$sd = (d13Ca$d13Ca_97p5 - d13Ca$d13Ca_2p5) / 2

## Age index in 100kyr bins
pd.p$ai = round(pd.p$age_mean * 10, 0)

## Add d13Ca to dataset
pd.p$d13Ca = d13Ca$d13Ca_50[pd.p$ai]
pd.p$ed13Ca = d13Ca$sd[pd.p$ai]

## Parse
pd.p.l = parseFranks(pd.p)

# Run inversion - works
## Parameters to save
parms = c("Pl", "l", "amax.scale", "D", "gc.scale", "ca", "meso.scale",
          "Ci0_m", "A0_m", "d13Ca_m", "A", "D13C", "gcop")

post = jags.parallel(pd.p.l, NULL, parms, "code/models/forwardFranksMultiAbAd.R", 
                     n.chains = 3, n.iter = 1e5, n.burnin = 1e4, n.thin = 1e2)

plot(d13Ca$age[unique(pd.p$ai)], post$BUGSoutput$median$ca)


## Looking at some of the other fields, comments to be addressed:
## - What to do with missing uncertainty values where measured quantities are
## reported, e.g., eDab?
## - Why are adaxial values only reported for some samples of a taxon and what
## to do about it?
## - Does "measured (n=6)" mean six obs on a single leaf, or six leaves?