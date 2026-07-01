library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Get and parse proxy files ----
## All files
pf = list.files("data/proxyData/", full.names = TRUE)

## Soils
pf.s = pf[grep("paleosol", pf)]

## Read file
pd.s = read.xlsx(pf.s[1], 1, startRow = 4)

# Per sample inversion
## Age index
ages = unique(pd.s$age_mean)
pd.s$ai = rep(0)
for(i in seq_along(ages)){
  pd.s$ai[pd.s$age_mean == ages[i]] = i
}

## Sites index
sites = unique(pd.s$lat)
pd.s$si = rep(0)
for(i in seq_along(sites)){
  pd.s$si[pd.s$lat == sites[i]] = i
}

## Parse
pd.s.l = parseSoil(pd.s)
pd.s.l$ages = ages
pd.s.l$sites = sites

## Add d13Ca
d13Ca = read.csv("data/d13Ca_Cenozoic.csv")
d13Ca$sd = (d13Ca$d13Ca_97p5 - d13Ca$d13Ca_2p5) / 2

## Age index in 100kyr bins
pd.s.l$d13Ca.obs = data.frame("d13Ca.m" = d13Ca$d13Ca_50[round(pd.s.l$ages * 10, 0)],
                            "d13Ca.sd" = d13Ca$sd[round(pd.s.l$ages * 10, 0)])

## Run it
## Parameters to save
parms = c("pCO2", "d13Ca", "MAT", "MAP")

post = jags.parallel(pd.s.l, NULL, parms, "code/models/soil.R", 
                     n.chains = 3, n.iter = 1e5, n.burnin = 1e4, n.thin = 10)

View(post$BUGSoutput$summary)
plot(ages, post$BUGSoutput$median$pCO2)
i = sample(seq_along(ages), 1)
plot(post$BUGSoutput$sims.list$MAP[, i], post$BUGSoutput$sims.list$pCO2[, i])
plot(post$BUGSoutput$sims.list$MAT[, i], post$BUGSoutput$sims.list$pCO2[, i])
plot(post$BUGSoutput$sims.list$d13Ca[, i], post$BUGSoutput$sims.list$pCO2[, i])
plot(pd.s.l$d13Cc.obs$d13C_cc - pd.s.l$d13Co.obs$d13Com_occluded, 
     post$BUGSoutput$median$pCO2)

# Run timeseries inversion ----
## Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

## Re-populate data
## Parse
pd.s.l = parseSoil(pd.s)
pd.s.l$sites = sites

pd.s.l$d13Ca.obs = data.frame("d13Ca" = numeric(length(ages)), 
                          "ed13Ca" = numeric(length(ages)))
for(i in seq_along(ages)){
  pd.s.l$d13Ca.obs[i, 1] = d13Ca$d13Ca_50[which.min(abs(ages[i] - d13Ca$age))]
  pd.s.l$d13Ca.obs[i, 2] = d13Ca$sd[which.min(abs(ages[i] - d13Ca$age))]
}

## Re-populate age index
for(i in seq_along(pd.s$age_mean)){
  pd.s.l$ai[i] = which.min(abs(pd.s$age_mean[i] - ages))
}

## Add timeseries info
pd.s.l$dt = stepsize
pd.s.l$nstep = length(ages)

parms = c("pCO2", "d13Ca", "MAT", "MAP", "PCQ_to", "PCQ_pf", "tsc",
          "f_R", "pore", "D13Cr")

## Run it
post.ts = jags.parallel(pd.s.l, NULL, parms, "code/models/soil_ts.R", 
                        n.chains = 3, n.iter = 1e3, n.burnin = 1e2)

View(post.ts$BUGSoutput$summary)

tsplot(ages, post.ts, "pCO2")
pointplot(pd.s$age_mean, post, "pCO2")
tsplot(ages, post.ts, "MAT")
points(pd.s$age_mean, pd.s$tempC)
tsplot(ages, post.ts, "d13Ca")
points(d13Ca$age, d13Ca$d13Ca_50)
tsplot(ages, post.ts, "MAP")
tsplot(ages, post.ts, "PCQ_to")
tsplot(ages, post.ts, "PCQ_pf")
tsplot(ages, post.ts, "pore")

# Current data structure
## datum
## soil

# Long-term solution
## datum
## soil
## age
## site

