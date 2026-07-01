model{

  # Data model ----
  for(i in 1:nd){
    d13Cc.obs[i, 1] ~ dnorm(d13Cc[i], d13Cc.pre[i])
    d13Cc.pre[i] = 1 / d13Cc.obs[i, 2] ^ 2
  }
  
#  for(i in 1:length(d18Oc.ai)){
#    d18Oc.obs[i, 1] ~ dnorm(d18Oc[d18Oc.ai[i]], d18Oc.pre[i])
#    d18Oc.pre[i] = 1 / d18Oc.obs[i, 2] ^ 2
#  }
  
  for(i in 1:nd){
    d13Co.obs[i, 1] ~ dnorm(d13Co[i], d13Co.pre[i])
    d13Co.pre[i] = 1 / d13Co.obs[i, 2] ^ 2
  }
  
  for(i in 1:length(nd)){
    MAT.obs[i, 1] ~ dnorm(MAT[ai[i]], MAT.pre[i])
    MAT.pre[i] = 1 / MAT.obs[i, 2] ^ 2
  }
  
  for(i in 1:length(ages)){
    d13Ca.obs[i, 1] ~ dnorm(d13Ca[i], d13Ca.pre[i])
    d13Ca.pre[i] = 1 / d13Ca.obs[i, 2] ^ 2
  }
  
#  for(i in 1:length(MAP.ai)){
#    MAP.obs[i, 1] ~ dnorm(MAP[MAP.ai[i]], MAP.pre[i])
#    MAP.pre[i] = 1 / MAP.obs[i, 2] ^ 2
#  }
  
#  for(i in 1:length(MS.ai)){
#    MS.obs[i, 1] ~ dnorm(MS[MS.ai[i]], MS.pre[i])
#    MS.pre[i] = 1 / MS.obs[i, 2] ^ 2
#  }

#  for(i in 1:length(D47c.ai)){
#    D47c.obs[i, 1] ~ dnorm(D47c[D47c.ai[i]], D47c.pre[i])
#    D47c.pre[i] = 1 / D47c.obs[i, 2] ^ 2
#  }
  
  # Process model ----
  for(i in 1:nd){
    ## Atmospheric humidity
    PPCQ[i] = MAP[ai[i]] * PCQ_pf[ai[i]] 
    h_m[i] = min(0.95, 0.25 + 0.7 * (PPCQ[i] / 900))
    ha[i] ~ dbeta(h_m[i] * 100 / (1 - h_m[i]), 100) # PCQ atmospheric humidity
    
    ## Depth to carbonate formation based on Retallack (2005) data, meters
    z.min[i] = MAP[ai[i]] * 0.0925 + 13.4
    z.thick[i] = abs(PPCQ[i] - MAP[ai[i]] / 4) * 0.74 + 17.3
    z.mean[i] = z.min[i] + z.thick[i] / 2
    z.beta[i] = z.mean[i] / (22 ^ 2)
    z.alpha[i] = z.mean[i] * z.beta[i]
    z[i] ~ dgamma(z.alpha[i], z.beta[i])
    # z[i] = (0.093 * MAP[i] + 13.12)
    z_m[i] = z[i] / 100

    ## Soil temperatures at depth z
    Tsoil[i] = MAT[ai[i]] + (PCQ_to[ai[i]] * sin(2 * 3.141593 * tsc[ai[i]] - 
                                                   z[i] / d)) / exp(z[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    Tair_PCQ[i] = MAT[ai[i]] + PCQ_to[ai[i]] * sin(2 * 3.141593 * tsc[ai[i]])
    
    ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
    PET_PCQ_D.1[i] = ifelse(ha[i] < 0.5, 
                       0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * 
                         (23.885 * Rs[si[i]] + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                       0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * 
                         (23.885 * Rs[si[i]] + 50))
    PET_PCQ_D[i] = max(PET_PCQ_D.1[i], 0.01)
    PET_PCQ[i] = PET_PCQ_D[i] * 90
    PET_D.1[i] = ifelse(ha[i] < 0.5, 
                            0.013 * (MAT[ai[i]] / (MAT[ai[i]] + 15)) * 
                          (23.885 * Rs[si[i]] + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                            0.013 * (MAT[ai[i]] / (MAT[ai[i]] + 15)) * 
                          (23.885 * Rs[si[i]] + 50))
    PET_D[i] = max(PET_D.1[i], 0.01)
    PET[i] = PET_D[i] * 365
    
    ## AET in mm/quarter from Budyko curve - Pike (1964)
    AET_var[i] ~ dgamma(1 / 0.2 ^ 2, 1 / 0.2 ^ 2) # noise parameter - Gentine (2012)
    AET_PCQ[i] = PPCQ[i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[i])) * 
                                                  AET_var[i])) ^ 2)))

    ## Average rooting depth
    AI[i] = PET[i] / MAP[ai[i]]
    L[i] = ifelse(AI[i] < 1.4, (-2 * AI[i]^2 + 2.5 * AI[i] + 1) * 100, 60)
    
    ## Carbon isotopes ----
    ### Free air porosity
    FAP.1[i] = min((pore[ai[i]] - ((PPCQ[i] - AET_PCQ[i]) / 
                                     (L[i] * 10 * pore[ai[i]]))), 
                   pore[ai[i]] - 0.05)
    FAP[i] = max(FAP.1[i], 0.01)
    
    ### Soil respiration rate 
    R_PCQ_D_m1[i] = 1.25 * exp(0.05452 * Tair_PCQ[i]) * PPCQ[i] / (127.77 + PPCQ[i])
    R_PCQ_D_m[i] = R_PCQ_D_m1[i] * f_R[ai[i]] # (gC/m2/d)
    R_beta[i] = R_PCQ_D_m[i] / (R_PCQ_D_m[i] * 0.15) ^ 2
    R_alpha[i] = R_PCQ_D_m[i] * R_beta[i]
    R_PCQ_D[i] ~ dgamma(R_alpha[i], R_beta[i])
    
    ### Convert to molC/cm3/s
    R_PCQ_D.1[i] = R_PCQ_D[i] / (12.01 * 100 ^ 2)  # from gC/m2/d to molC/cm2/d
    R_PCQ_S[i] = R_PCQ_D.1[i] / (24 * 3600)  # molC/ cm2 / s
    R_PCQ_S_0[i]= R_PCQ_S[i] / (L[i] * pore[ai[i]]) # Quade et al. (2007)
    
    ### CO2 diffusion
    Dair[i] = 0.1369 * (Tsoil.K[i] / 273.15) ^ 1.958
    DIFC[i] = FAP[i] * tort[si[i]] * Dair[i]
    
    ### S(z)
    k[i] = L[i] / (2 * log(2)) # Respiration characteristic production depth (cm) - Quade (2007)
    S_z_mol[i] = k[i] ^ 2 * R_PCQ_S_0[i] / DIFC[i] * (1 - exp(-z[i] / k[i])) # (mol/cm3)
    S_z[i] = S_z_mol[i] * (0.08206 * Tsoil.K[i] * 10^9) # ppmv

    ### d13C of soil-respired CO2
    DD13_water[i] = 25.09 - 1.2 * (MAP[ai[i]] + 975) / 
      (27.2 + 0.04 * (MAP[ai[i]] + 975))
    D13C_plant[i] = (28.26 * 0.22 * (pCO2[ai[i]] + 23.9)) / 
      (28.26 + 0.22 * (pCO2[ai[i]] + 23.9)) - DD13_water[i] # schubert & Jahren (2015)
    D13C_off[i] ~ dnorm(0, 1 / 2 ^ 2) # Noise term
    d13Cr[i] = d13Ca[ai[i]] - D13C_plant[i] + D13C_off[i]
#    d13Cr[i] = d13Ca[ai[i]] + D13Cr[ai[i]]
    d13Co[i] = d13Cr[i] + SOM.frac[si[i]]
    
    ### d13C of pedogenic carbonate
    d13Cs[i] = (pCO2[ai[i]] * d13Ca[ai[i]] + S_z[i] * 
                  (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[ai[i]])
    d13Cc[i] = ((1 + (11.98 - 0.12 * Tsoil[i]) / 1e3) * (d13Cs[i] + 1e3)) - 1e3
    
    ## Oxygen isotopes ----
    ### Rainfall isotopes
#    Tair_OOS[i] = (4 * MAT[i] - Tair_PCQ[i]) / 3
#    d18p[i] ~ dnorm(-15 + 0.58 * (Tair_OOS[i] * (1 - PCQ_pf[i]) * (1 - spre[i]) + Tair_PCQ[i] * PCQ_pf[i]), 1 / 1 ^ 2) # Precipitation d18O, ppt
#    R18p[i] = (d18p[i] / 1e3 + 1) * R18.VSMOW
    
    ### Equilibrium fractionation (Horita and Wesolowski 1994)
#    alpha18.eq[i] = 1 / exp(((1.137e6 / (Tsoil.K[i] ^ 2) - 0.4156e3/Tsoil.K[i] - 2.0667) /1e3))
    
    ### Atmospheric water vapor isotopes
#    R18a[i] = R18p[i] * alpha18.eq[i]
    
    ### Soil evaporation from AET
#    E1[i] = ETR[i] * AET_PCQ[i]
#    E[i] = max(E1[i], 1) # minimum of 1 mm
#    E_s[i] = E[i] / (1e3 * 90 * 24 * 3600) # soil evaporation rate in m/sec
    
    ### Water vapor diffusivity
#    es[i] = (0.611 * exp(17.502 * Tsoil[i] / (Tsoil[i] + 240.97))) * 1e3 # saturated water vapor pressure from Tetens formula
#    N.sat[i] = 0.01802 * es[i] / (Rgas * Tsoil.K[i]) # saturated water vapor concentration at a given temperature
#    Dv.soil[i] = Dv.air * tort * (pore[i] - 0.05) # effective diffusivity of water vapor in soil (m2/s)
#    z.bar[i] = N.sat[i] * Dv.soil[i] / (E_s[i] * rho) # penetration depth (m)
#    z.ef1[i] = (1 - ha[i]) * z.bar[i] # the thickness of the water vapor phase region (m)
#    z.ef[i] = max(z.ef1[i], 1e-10)
    
    ### Liquid water diffusivity (m2/s) (Easteal 1984)
#    Dl[i] = exp(1.6766 + 1.6817 * (1e3 / Tsoil.K[i]) - 0.5773 * (1e3 / Tsoil.K[i]) ^ 2) * 1e-9 
#    Dl.soil[i] = Dl[i] * pore[i] * tort # effective diffusivity of liquid water (m2/s)
#    z.hat[i] = Dl.soil[i] / E_s[i] # the decay length (mean penetration depth)
    
    ### The evaporation front
#    h.ef[i] = ha[i] + z.ef[i] / z.bar[i] # humidity at the evaporation front
#    R18.ef[i] = (alpha18.diff * R18p[i] * (z.ef[i] / z.bar[i]) + 
#                   ha[i] * R18a[i]) / (h.ef[i] * alpha18.eq[i]) # isotopic composition at the evaporation front
    
    ### Isotope composition of soil water at depth z
#    hs[i] = min(ha[i] + z_m[i] / z.bar[i], 1)
#    z.f[i] = (pore[i] / a.theta) * log(z_m[i] / z.ef[i]) # the modified depth function
#    R18s[i] = ifelse(z_m[i] <= z.ef[i], 
#                      (alpha18.diff * R18p[i] * z_m[i] / z.bar[i] + ha[i] * R18a[i]) / 
#                        (hs[i] * alpha18.eq[i]),
#                      (R18.ef[i] - R18p[i]) * exp(-z.f[i] / z.hat[i]) + R18p[i])
#    d18s[i] = ((R18s[i] / R18.VSMOW) - 1) * 1e3
    
    ### Isotope composition of soil carbonate
#    alpha18_c_w_eq[i] = exp((1.61e4 / Tsoil.K[i] - 24.6) / 1e3) # Wostbrock (2020)
#    R18c[i] = R18s[i] * alpha18_c_w_eq[i]
#    d18Oc[i] = (R18c[i] / R18.VPDB - 1) * 1e3
#    D47c[i] = 0.0391e6 / Tsoil.K[i] ^ 2 + 0.154 # Andersen (2021)
  }
  
  # Site level time-dependent ----
  # here cheating since there's only one record 
  for(i in 1:length(ages)){
    ## Primary ----
    MAT[i] ~ dunif(4, 25) # mean annual temperature
    PCQ_to[i] ~ dunif(5, 20)
    MAP[i] ~ dunif(1e2, 1e3) # mean annual precipitation, mm
    PCQ_pf[i] ~ dunif(0.02, 0.25) # PCQ precipitation fraction
    
    ## Secondary soil ----
    tsc[i] ~ dunif(0, 0.5)
    f_R[i] ~ dbeta(2, 16) # ratio of PCQ to mean annual respiration rate
#    ETR[i] ~ dbeta(0.06 * 1e3 / 0.94, 1e3) # Soil evaporation / AET
#    spre[i] ~ dbeta(27, 22)
    pore[i] ~ dunif(0.45, 0.54) # soil porosity
#    D13Cr[i] ~ dunif(-22, -10) # photosynthetic discrimination
  }

  # Global time-dependent ----  
  for(i in 1:length(ages)){
    ## Primary environmental ----
    pCO2[i] = ca.s[i] * 1e3
    ca.s[i] ~ dunif(0.1, 2)
    d13Ca[i] ~ dunif(-8, -3)
  }
  
  # Site dependent ----
  for(i in 1:length(sites)){
    Ra[i] = 42.608 - 0.3538 * abs(lat[i]) # total radiation at the top of the atmosphere
    Rs[i] = Ra[i] * 0.16 * sqrt(12) # daily temperature range assumed to be 12
    tort[i] ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
    SOM.frac[i] ~ dunif(0.5, 1.5)
  }

  # Constants ----
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  ### Isotope ratio constants
  R13.VPDB = 0.011237
#  R18.VSMOW = 0.0020052
#  R18.VPDB = 0.0020672
#  alpha18.diff = 1.028489
#  a.theta = 0.05 # rate of increase of water content with depth (m-1) (Barnes and Allison, 1983)
#  Rgas = 8.314462 # gas constant
#  rho = 1e3 # liquid water density (kg/m3)
#  Dv.air = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)
}
  