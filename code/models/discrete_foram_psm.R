model{
  
  ################################################################################
  # likelihood function
  ################################################################################
  
  for (i in 1:n_d11B){
    d11B.data[i] ~ dnorm(d11Bcarb[i], d11B.p[i])
    d11B.p[i] <- 1 / d11B.data.sigma[i]^2
  }
  
  for (i in 1:n_mgca){
    mgca.data[i] ~ dnorm(mgcacarb, mgca.p[i])
    mgca.p[i] <- 1 / mgca.data.sigma[i]^2
  }
  
  for (i in 1:n_d18O){
    d18O.data[i] ~ dnorm(d18Ocarb, d18O.p[i])
    d18O.p[i] <- 1 / d18O.data.sigma[i]^2
  }
  
  ################################################################################
  # environmental priors
  ################################################################################
  
  tempC ~ dnorm(tempC.m, 1 / tempC.sigma^2)T(tempC.min, tempC.max)
  pH ~ dnorm(pH.m, 1 / pH.sigma^2)T(pH.min, pH.max)
  sal ~ dnorm(sal.m, 1 / sal.sigma^2)T(sal.min, sal.max)
  dic ~ dnorm(dic.m, 1 / dic.sigma^2)T(dic.min, dic.max)
  press ~ dnorm(press.m, 1 / press.sigma^2)T(press.min, press.max)
  
  xca ~ dnorm(xca.m, 1 / xca.sigma^2)T(xca.min, xca.max)
  xmg ~ dnorm(xmg.m, 1 / xmg.sigma^2)T(xmg.min, xmg.max)
  xso4 ~ dnorm(xso4.m, 1 / xso4.sigma^2)T(xso4.min, xso4.max)
  d11Bsw ~ dnorm(d11Bsw.m, 1 / d11Bsw.sigma^2)T(d11Bsw.min, d11Bsw.max)
  d18Osw ~ dnorm(d18Osw.m, 1 / d18Osw.sigma^2)T(d18Osw.min, d18Osw.max)
  
  ################################################################################
  # proxy system priors
  ################################################################################
  
  alpha ~ dnorm(alpha.m, 1 / alpha.sigma^2)T(alpha.min, alpha.max)
  epsilon <- (alpha - 1) * 1000
  
  sw.sens ~ dnorm(sw.sens.m, 1 / sw.sens.sigma^2)T(sw.sens.min, sw.sens.max)
  
  for (g in 1:n_cal){
    m.vital[g] ~ dnorm(m.vital.m[g], 1 / m.vital.sigma[g]^2)
    c.vital[g] ~ dnorm(c.vital.m[g], 1 / c.vital.sigma[g]^2)
    c.final[g] <- c.vital[g] + c.correction[g]
  }
  
  Hp ~ dnorm(Hp.m, 1 / Hp.sigma^2)T(Hp.min, Hp.max)
  A ~ dnorm(A.m, 1 / A.sigma^2)T(A.min, A.max)
  
  d18O_pHcorr ~ dnorm(d18O_pHcorr.m, 1 / d18O_pHcorr.sigma^2)T(d18O_pHcorr.min, d18O_pHcorr.max)
  salcorrco ~ dnorm(salcorrco.m, 1 / salcorrco.sigma^2)T(salcorrco.min, salcorrco.max)
  pHcorrco ~ dnorm(pHcorrco.m, 1 / pHcorrco.sigma^2)T(pHcorrco.min, pHcorrco.max)
  
  indexop ~ dnorm(indexop.m, 1 / indexop.sigma^2)T(indexop.min, indexop.max)
  seccal.d18O ~ dnorm(seccal.d18O.m, 1 / seccal.d18O.sigma^2)
  
  A_Daeron ~ dnorm(17.57, 1 / 0.43^2)T(16.71, 18.43)
  B_Daeron ~ dnorm(29.89, 1 / 0.06^2)T(29.77, 30.01)
  
  ################################################################################
  # constants
  ################################################################################
  
  R <- 83.131
  
  # modern
  xcam <- 10.2821
  xmgm <- 52.8171
  xso4m <- 28.24
  mgcaswm <- xmgm / xcam
  
  s1_ca <- 5 / 1000
  s1_mg <- 17 / 1000
  s1_so4 <- 208 / 1000
  
  s2_ca <- 157 / 1000
  s2_mg <- 420 / 1000
  s2_so4 <- 176 / 1000
  
  sspc_ca <- 185 / 1000
  sspc_mg <- 518 / 1000
  sspc_so4 <- 106 / 1000
  
  ################################################################################
  # carbonate chemistry
  ################################################################################
  
  temp <- tempC + 273.15
  mgcasw <- xmg / xca
  
  K0 <- exp(9345.17 / temp - 60.2409 + 23.3585 * log(temp / 100) + sal * 
    (0.023517 - 0.00023656 * temp + 0.0047036 * ((temp / 100)^2)))
  
  Ks1m_st <- exp(2.83655 - 2307.1266 / temp - 1.5529413 * log(temp) -
    ((0.20760841 + 4.0484 / temp) * sqrt(sal)) + 0.0846834 * sal - 0.00654208 * sal^1.5 +
    log(1 - (0.001005 * sal)))
  
  Ks2m_st <- exp(-9.226508 - 3351.6106 / temp - 0.2005743 * log(temp) -
    ((0.106901773 + 23.9722 / temp) * sqrt(sal)) + 0.1130822 * sal - 0.00846934 * sal^1.5 +
    log(1 - (0.001005 * sal)))
  
  logKsspcm_st <- -171.9065 - 0.077993 * temp + 2839.319 / temp + 71.595 * (log(temp) / log(10)) +
    (-0.77712 + 0.0028426 * temp + 178.34 / temp) * sal^0.5 - 0.07711 * sal + 0.0041249 * sal^1.5
  
  Ksspcm_st <- 10^logKsspcm_st
  
  lnKsB_st <- ((-8966.9 - 2890.53 * sal^0.5 - 77.942 * sal + 1.728 * sal^1.5 - 0.0996 * sal^2) / temp) +
    148.0248 + 137.1942 * sal^0.5 + 1.62142 * sal - (24.4344 + 25.085 * sal^0.5 + 0.2474 * sal) * log(temp) +
    0.053105 * sal^0.5 * temp
  
  KsB_st <- exp(lnKsB_st)
  
  Ksw_st <- exp(148.96502 - 13847.26 / temp - 23.6521 * log(temp) +
    (118.67 / temp - 5.977 + 1.0495 * log(temp)) * sal^0.5 - 0.01615 * sal)
  
  delV1 <- -25.50 + 0.1271 * tempC
  delV2 <- -15.82 + -0.0219 * tempC
  delVspc <- -48.76 + 0.5304 * tempC
  delVB <- -29.48 + 0.1622 * tempC + (2.608 / 1000) * tempC^2
  delVw <- -25.60 + 0.2324 * tempC + (-3.6246 / 1000) * tempC^2
  
  delk1 <- (-3.08 / 1000) + (0.0877 / 1000) * tempC
  delk2 <- (1.13 / 1000) + (-0.1475 / 1000) * tempC
  delkspc <- (-11.76 / 1000) + (0.3692 / 1000) * tempC
  delkB <- -2.84 / 1000
  delkw <- (-5.13 / 1000) + (0.0794 / 1000) * tempC
  
  Ks1m <- exp(-((delV1 / (R * temp)) * press) + ((0.5 * delk1) / (R * temp)) * press^2) * Ks1m_st
  
  Ks2m <- exp(-((delV2 / (R * temp)) * press) + ((0.5 * delk2) / (R * temp)) * press^2) * Ks2m_st
  
  Ksspcm <- exp(-((delVspc / (R * temp)) * press) + ((0.5 * delkspc) / (R * temp)) * press^2) * Ksspcm_st
  
  KsB <- exp(-((delVB / (R * temp)) * press) + ((0.5 * delkB) / (R * temp)) * press^2) * KsB_st
  
  Ksw <- exp(-((delVw / (R * temp)) * press) + ((0.5 * delkw) / (R * temp)) * press^2) * Ksw_st
  
  Ks1 <- Ks1m * (1 + (s1_ca * (xca / xcam - 1) + s1_mg * (xmg / xmgm - 1) + s1_so4 * (xso4 / xso4m - 1)))
  
  Ks2 <- Ks2m * (1 + (s2_ca * (xca / xcam - 1) + s2_mg * (xmg / xmgm - 1) + s2_so4 * (xso4 / xso4m - 1)))
  
  Ksspc <- Ksspcm * (1 + (sspc_ca * (xca / xcam - 1) + sspc_mg * (xmg / xmgm - 1) + sspc_so4 * (xso4 / xso4m - 1)))
  
  ################################################################################
  # carbonate system
  ################################################################################
  
  hyd <- 10^(-pH)
  co2 <- (dic * 1e-3) / (1 + (Ks1 / hyd) + ((Ks1 * Ks2) / hyd^2))
  fco2 <- co2 / K0
  pco2_mol <- fco2 / 0.9968
  pco2 <- pco2_mol * 1e6
  
  ################################################################################
  # proxy forward models
  ################################################################################
  
  pKsB <- -(log(KsB) / log(10))
  t1 <- 10^(pKsB - pH)
  
  d11Bb <- ((t1 * epsilon) - (t1 * d11Bsw) - d11Bsw) /
    (-((t1 * alpha) + 1))
  
  for (i in 1:n_d11B){
    d11Bcarb[i] <- m.vital[cal_i[i]] * d11Bb + c.final[cal_i[i]]
  }
  
  d18Osw.sc <- d18Osw + sw.sens * (sal - 35)
  
  alpha.ccw <- exp((A_Daeron * 10^3 * (1 / temp - 1 / 297.7) + B_Daeron) / 10^3)
  epsilon.ccw <- (alpha.ccw - 1) * 10^3
  
  d18Osmow <- d18Osw.sc + epsilon.ccw
  d18Opdb <- 0.97001 * d18Osmow - 29.99
  
  d18Odiag <- d18Opdb - (pH - 8) * (d18O_pHcorr * 10)
  d18Ocarb <- d18Odiag * (1 - (indexop / 100)) + seccal.d18O * (indexop / 100)
  
  mgcacarb <- ((mgcasw / mgcaswm)^Hp) *
    exp((salcorrco / 100) * (sal - 35) + A * tempC - pHcorrco * (pH - 8.05))
}
