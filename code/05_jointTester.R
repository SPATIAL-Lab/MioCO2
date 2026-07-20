library(openxlsx)
source("code/helpers.R")
library(R2jags)

load("bigout/phyto.rda")
load("bigout/plants.rda")
load("bigout/soils.rda")

# Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

# Load data
data.pass = loadTS(ages)

parms = c("pCO2", "d13Ca", "MAT", "MAP", "S_z", "b", "pCO2.v")

post.ts = jags.parallel(data.pass, NULL, parms, "code/models/joint_ts.R", 
                        n.chains = 5, n.iter = 1e5, n.burnin = 5e4, 
                        n.thin = 50)

data.pass.natmos = data.pass
data.pass.natmos$d13Ca.obs[, 2] = rep(1e3)
post.ts.natmos = jags.parallel(data.pass.natmos, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 5, n.iter = 2e4, n.burnin = 1e4, 
                           n.thin = 10)

data.pass.nalk = data.pass
data.pass.nalk$d13C.obs.alk[, 2] = rep(1e3)
post.ts.nalk = jags.parallel(data.pass.nalk, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 5, n.iter = 2e4, n.burnin = 1e4, 
                           n.thin = 10)

data.pass.nsoil = data.pass
data.pass.nsoil$d13Cc.obs[, 2] = rep(1e3)
post.ts.nsoil = jags.parallel(data.pass.nsoil, NULL, parms, "code/models/joint_ts.R", 
                           n.chains = 5, n.iter = 2e4, n.burnin = 1e4, 
                           n.thin = 10)

data.pass.nplant = data.pass
data.pass.nplant$d13Cp[, 2] = rep(1e3)
post.ts.nplant = jags.parallel(data.pass.nplant, NULL, parms, "code/models/joint_ts.R", 
                              n.chains = 5, n.iter = 2e4, n.burnin = 1e4, 
                              n.thin = 10)

View(post.ts$BUGSoutput$summary)

save(post.ts, file = "bigout/post.rda")
save(post.ts.natmos, file = "bigout/post_natmos.rda")
save(post.ts.nalk, file = "bigout/post_nalk.rda")
save(post.ts.nsoil, file = "bigout/post_nsoil.rda")
save(post.ts.nplant, file = "bigout/post_nplant.rda")


tsplot(ages, post.ts, "pCO2")
tsplot(ages, post.ts.natmos, "pCO2")
tsplot(ages, post.ts.nalk, "pCO2")
tsplot(ages, post.ts.nsoil, "pCO2")
tsplot(ages, post.ts.nplant, "pCO2")

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
