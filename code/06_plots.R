library(terra)
library(R2jags)
library(MASS)

# MCMC output
lf = list.files("bigout/", full.names = TRUE)
for(i in lf){
  load(i)
}

# Age vector
stepsize = 0.1
ages = seq(25 - stepsize / 2, 5 - stepsize / 2, by = -stepsize)

# Load data
data.pass = loadTS(ages)

## Need to run code in loadTS to get lats/longs
prox.raw = prox.raw[prox.raw$age_mean > 5, ]
prox.raw = prox.raw[prox.raw$age_mean < 25, ]

# Geo coords
ll.alk = prox.raw[, c("lon", "lat")]
ll.alk = unique(ll.alk)
alk.shp = vect(ll.alk, geom = c("lon", "lat"), crs = "WGS84")

ll.plant = pd.p[, c("lon", "lat")]
ll.plant = unique(ll.plant)
plant.shp = vect(ll.plant, geom = c("lon", "lat"), crs = "WGS84")

ll.soil = pd.s[, c("lon", "lat")]
ll.soil = unique(ll.soil)
soil.shp = vect(ll.soil, geom = c("lon", "lat"), crs = "WGS84")

# Continents
cont = vect("data/continents")

# Map 
png("out/map.png", width = 7, height = 4, units = "in", res = 600)
par(mar = rep(1, 4))
plot(cont, col = "grey80", axes = FALSE)
points(ll.alk, pch = 21, bg = "lightblue", cex = 1.5)
points(ll.plant, pch = 21, bg = "green3", cex = 1.5)
points(ll.soil, pch = 21, bg = "brown", cex = 1.5)
legend(-180, 0, c("Alkenones", "Plants", "Soils"), bg = "white", pch = 21,
       pt.bg = c("lightblue", "green3", "brown"), box.lty = 0, cex = 0.75)
dev.off()

# Individual proxies
png("out/proxies.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", col = "white", ylim = c(100, 1500),
       ylab = expression("[CO"[2]*"] (ppm)"))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey90")
pointplot(ages[data.pass$ai.alk], post.phyto, "pCO2", col = "lightblue")
pointplot(ages[data.pass$ai.soil], post.soils, "pCO2", col = "brown")
pointplot(ages[unique(data.pass$ai.plant)], post.plants, "pCO2", col = "green3")
dev.off()

# Reconstruction
ceno = read.csv("../CenoCO2/out/100kyrCO2.csv")
ceno[, 2:6] = exp(ceno[, 2:6])

png("out/curves.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", col = "white", ylab = expression("[CO"[2]*"] (ppm)"))
tsdens(ceno, col = "gold")
tsplot(ages, post.ts, "pCO2", add = TRUE)
legend("topright", legend = c("CenoCO2PIP", "This work"), lwd = 2, 
       col = c("gold", "black"), bty = "n")
dev.off()

# Sensitivity test - d13Ca
png("out/noAtm.png", width = 6, height = 4, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", ylim = c(250, 600), ylab = expression("[CO"[2]*"] (ppm)"))
tsplot(ages, post.ts.natmos, "pCO2", add = TRUE, "blue")
legend("topright", legend = c(expression("No "*delta^{13}*"C"[atm]), "Control"), lwd = 2, 
       col = c("blue", "black"), bty = "n")
dev.off()

# Sensitivity test - proxies
png("out/sens.png", width = 12, height = 4, units = "in", res = 600)
layout(matrix(c(1, 2, 3), nrow = 1))

par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", ylim = c(250, 500), ylab = expression("[CO"[2]*"] (ppm)"))
tsplot(ages, post.ts.nsoil, "pCO2", add = TRUE, "brown")
legend("topright", legend = c("No soils", "Control"), lwd = 2, 
       col = c("brown", "black"), bty = "n")

par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", ylim = c(200, 500), ylab = expression("[CO"[2]*"] (ppm)"))
tsplot(ages, post.ts.nplant, "pCO2", add = TRUE, "green3")
legend("topright", legend = c("No plants", "Control"), lwd = 2, 
       col = c("green3", "black"), bty = "n")

par(mar = c(5, 5, 1, 1))
tsplot(ages, post.ts, "pCO2", ylim = c(250, 800), ylab = expression("[CO"[2]*"] (ppm)"))
tsplot(ages, post.ts.nalk, "pCO2", add = TRUE, "lightblue")
legend("topright", legend = c("No alkenones", "Control"), lwd = 2, 
       col = c("lightblue", "black"), bty = "n")
dev.off()

# Mechanisms
png("out/mech.png", width = 8, height = 4, units = "in", res = 600)
layout(matrix(c(1, 2), nrow = 1))

par(mar = c(5, 5, 1, 1))
plot(post.soils.ts$BUGSoutput$median$S_z, post.ts$BUGSoutput$median$S_z,
     xlab = "S(z), soil inversion (ppm)", ylab = "S(z), joint inversion (ppm)",
     pch = 21, bg = "brown", cex = 1.5)
abline(0, 1)

par(mar = c(5, 5, 1, 1))
plot(post.phyto.ts$BUGSoutput$median$b[-13], post.ts$BUGSoutput$median$b,
     xlab = "b, alkenone inversion (uM)", ylab = "b, joint inversion (uM)",
     pch = 21, bg = "lightblue", cex = 1.5)
abline(0, 1)
dev.off()

# Climate sensitivity
dbls = as.vector(post.ts$BUGSoutput$sims.list$pCO2)
dbls = log(dbls / 200) / log(2)

MAT = as.vector(post.ts$BUGSoutput$sims.list$MAT)

kd = kde2d(dbls, MAT, n = 50)

dbls.med = as.vector(log(post.ts$BUGSoutput$median$pCO2 / 200) / log(2))
mod = lm(post.ts$BUGSoutput$median$MAT ~ dbls.med)

png("out/ClimSens.png", width = 5, height = 5, units = "in", res = 600)
par(mar = c(5, 5, 1, 1))

image(kd, xlab = expression("CO"[2]*" doublings"), ylab = "MAT",
      xlim = c(0.4, 1.2), ylim = c(9, 25))
abline(mod, lw = 2)
points(dbls.med, post.ts$BUGSoutput$median$MAT, pch = 21, bg = "white")

box()
dev.off()
