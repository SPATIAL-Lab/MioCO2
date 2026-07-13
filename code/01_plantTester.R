library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Get and parse proxy files ----
## All files
pf = list.files("data/proxyData/", full.names = TRUE)

## Plants
pf.p = pf[grep("stomata", pf)]

## First file
pd.p = read.xlsx(pf.p[1], 1, startRow = 4)

## Others
if(length(pf.p) > 1){
  for(i in pf.p[2:length(pf.p)]){
    pd.p.a = read.xlsx(i, 1, startRow = 4)
    pd.p = rbind(pd.p, pd.p.a)
  }
}

# Data cleaning and additions
## Get d13Ca data
d13Ca = read.csv("data/d13Ca_Cenozoic.csv")
d13Ca$sd = (d13Ca$d13Ca_97p5 - d13Ca$d13Ca_2p5) / 2
d13Ca$age = as.integer(d13Ca$a * 100, 0) / 100

## Age index in 100kyr bins
pd.p$ai = round(pd.p$age_mean * 10, 0)

## Add d13Ca to dataset
pd.p$d13Ca.obs = d13Ca$d13Ca_50[pd.p$ai]
pd.p$ed13Ca = d13Ca$sd[pd.p$ai]

## Parse
pd.p.l = parseFranks(pd.p)

## Site ages for plotting
siteAges = d13Ca$age[unique(pd.p$ai)]

# Run independent sample inversion ----
## Parameters to save
parms = c("pCO2", "d13Ca")

post.plants = jags.parallel(pd.p.l, NULL, parms, "code/models/plant.R", 
                     n.chains = 3, n.iter = 3e5, n.burnin = 1e5, n.thin = 1e2)

View(post.plants$BUGSoutput$summary)
plot(d13Ca$age[unique(pd.p$ai)], post.plants$BUGSoutput$median$pCO2)

# Run timeseries inversion ----
## Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

## Parse
pd.p.l = parseFranks(pd.p)

## Re-populate age index
for(i in seq_along(pd.p$age_mean)){
  pd.p$ai[i] = which.min(abs(pd.p$age_mean[i] - ages))
}

pd.p.l$level = pd.p$ai

## Re-populate d13Ca
pd.p.l$d13Ca.obs = data.frame("d13Ca" = numeric(length(ages)), 
                          "ed13Ca" = numeric(length(ages)))
for(i in seq_along(ages)){
  pd.p.l$d13Ca.obs[i, 1] = d13Ca$d13Ca_50[which.min(abs(ages[i] - d13Ca$age))]
  pd.p.l$d13Ca.obs[i, 2] = d13Ca$sd[which.min(abs(ages[i] - d13Ca$age))]
}

## Add timeseries info
pd.p.l$dt = stepsize
pd.p.l$nstep = length(ages)

## run it
parms = c("pCO2", "d13Ca")
post.plants.ts = jags.parallel(pd.p.l, NULL, parms, "code/models/plant_ts.R", 
                     n.chains = 3, n.iter = 1e5, n.burnin = 5e4, n.thin = 50)

View(post.plants.ts$BUGSoutput$summary)
tsplot(ages, post.plants.ts, "pCO2")
pointplot(siteAges, post.plants, "pCO2")

save(post.plants, post.plants.ts, file = "bigout/plants.rda")