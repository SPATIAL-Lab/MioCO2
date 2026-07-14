library(openxlsx)
source("code/helpers.R")
library(R2jags)

load("bigout/phyto.rda")
load("bigout/plants.rda")
load("bigout/soils.rda")

# Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

# Initialize data.pass
data.pass = list("dt" = stepsize)
data.pass$n.steps = length(ages)

# Add atmospheric data
d13Ca = read.csv("data/d13Ca_Cenozoic.csv")
d13Ca$sd = (d13Ca$d13Ca_97p5 - d13Ca$d13Ca_2p5) / 2

data.pass = append(data.pass, parseAtmos(d13Ca, ages))

# Alkenones
# Load proxy data for reconstruction from intermediate workbook
prox.raw <- read.xlsx("data/proxyData/phyto_Intermediate_combined.xlsx", 
                      sheet = "data4PSM", startRow = 4, 
                      na.strings = c("", "NA", "#N/A"))
names(prox.raw) <- trimws(names(prox.raw))

# Testing...
prox.raw = prox.raw[prox.raw$sample != "phytoplankton_Bolton_2016_13", ]

data.pass = append(data.pass, parsePhyto(prox.raw, ages))

# Soils
pf = list.files("data/proxyData/", full.names = TRUE)
pf.s = pf[grep("paleosol", pf)]
pd.s = read.xlsx(pf.s[1], 1, startRow = 4)

data.pass = append(data.pass, parseSoil(pd.s, ages))

# Plants
pf = list.files("data/proxyData/", full.names = TRUE)
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


data.pass = append(data.pass, parseFranks(pd.p, ages))

parms = c("pCO2", "d13Ca", "MAT", "MAP", "S_z", "b", "pCO2.eps.ac")

post.ts = jags.parallel(data.pass, NULL, parms, "code/models/joint_ts.R", 
                        n.chains = 3, n.iter = 6e3, n.burnin = 3e3)

data.pass.natmos = data.pass
data.pass.natmos$d13Ca.obs[, 2] = rep(1e3)
post.ts.natmos = jags.parallel(data.pass.natmos, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 3, n.iter = 6e3, n.burnin = 3e3)

data.pass.nalk = data.pass
data.pass.nalk$d13C.obs.alk[, 2] = rep(1e3)
post.ts.nalk = jags.parallel(data.pass.nalk, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 3, n.iter = 6e3, n.burnin = 3e3)

data.pass.nsoil = data.pass
data.pass.nsoil$d13Cc.obs[, 2] = rep(1e3)
post.ts.nsoil = jags.parallel(data.pass.nsoil, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 3, n.iter = 6e3, n.burnin = 3e3)

data.pass.nplant = data.pass
data.pass.nplant$d13Cp[, 2] = rep(1e3)
post.ts.nplant = jags.parallel(data.pass.nplant, NULL, parms, "code/models/joint_ts.R", 
                              n.chains = 3, n.iter = 6e3, n.burnin = 3e3)

View(post.ts$BUGSoutput$summary)

tsplot(ages, post.ts, "pCO2", ylim = c(150, 1500))

tsplot(ages, post.phyto.ts, "pCO2", col = "lightblue", add = TRUE)
tsplot(ages, post.plants.ts, "pCO2", col = "green3", add = TRUE)
tsplot(ages, post.soils.ts, "pCO2", col = "brown", add = TRUE)

pointplot(ages[unique(data.pass$ai.plant)], post.plants, "pCO2", col = "green3")
pointplot(ages[data.pass$ai.alk], post.phyto, "pCO2", col = "lightblue")

points(ages[data.pass$ai.alk], post.phyto$BUGSoutput$median$pCO2, pch = 21,
       bg = data.pass$si.alk)



points(ages[data.pass$ai.alk], post.phyto$BUGSoutput$median$pCO2, pch = 21, bg = "lightblue")
points(ages[unique(data.pass$ai.plant)], post.plants$BUGSoutput$median$pCO2, 
       pch = 21, bg = "green3")
points(ages[unique(data.pass$ai.soil)], post.soils$BUGSoutput$median$pCO2, 
       pch = 21, bg = "brown")


tsplot(ages, post.ts, "MAT", ylim = c(8, 35))
points(prox.in$age, prox.in$tempC.prior, col = "lightblue")
points(pd.s$age_mean, pd.s$tempC, col = "brown")
