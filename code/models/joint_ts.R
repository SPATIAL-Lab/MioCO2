model {
  
  # Likelihoods ----
  ## Alkanones
  for(i in 1:n.obs.alk){
    d13C.obs.alk[i, 1] ~ dnorm(d13C.alk[i], 1 / d13C.obs.alk[i, 2] ^ 2)
    temp.obs.alk[i, 1] ~ dnorm(tempC.alk[si.alk[i], ai.alk[i]], 
                               1 / temp.obs.alk[i, 2] ^ 2)
    po4.obs[i, 1] ~ dnorm(po4[si.alk[i], ai.alk[i]], 1 / po4.obs[i, 2] ^ 2)
  }
  
  ## Coco liths
  for(i in 1:n.lith){
    lith.obs[i, 1] ~ dnorm(len.lith[li[i]], 1 / lith.obs[i, 2] ^ 2)
  }
  
  ## Soil carbonate
  for(i in 1:n.obs.soil){
    d13Cc.obs[i, 1] ~ dnorm(d13Cc[i], 1 / d13Cc.obs[i, 2] ^ 2)
    d13Co.obs[i, 1] ~ dnorm(d13Co[i], 1 / d13Co.obs[i, 2] ^ 2)
    temp.obs.soil[i, 1] ~ dnorm(temp.soil[ai.soil[i]], 
                                1 / temp.obs.soil[i, 2] ^ 2)
  }
  
  ## Adaxial plants
  for(i in ind.ad){
    d13Cp[i, 1] ~ dnorm(d13C_m[i], 1 / d13Cp[i, 2] ^ 2)
    Dab[i, 1] ~ dnorm(D.ab[i], 1 / Dab[i, 2] ^ 2)
    GCLab[i, 1] ~ dnorm(Pl.ab[i] / s1_m[i], 1 / GCLab[i, 2] ^ 2)
    GCWab[i, 1] ~ dnorm(l.ab[i] / s2_m[gi[i]], 1 / GCWab[i, 2] ^ 2)
    Dad[i, 1] ~ dnorm(D.ad[i], 1 / Dad[i, 2] ^ 2)
    GCLad[i, 1] ~ dnorm(Pl.ad[i] / s1_m[i], 1 / GCLad[i, 2] ^ 2)
    GCWad[i, 1] ~ dnorm(l.ad[i] / s2_m[gi[i]], 1 / GCWad[i, 2] ^ 2)
  }
  
  ## Abaxial plants
  for(i in ind.ab){
    d13Cp[i, 1] ~ dnorm(d13C_m[i], 1 / d13Cp[i, 2] ^ 2)
    Dab[i, 1] ~ dnorm(D[i], 1 / Dab[i, 2] ^ 2)
    GCLab[i, 1] ~ dnorm(Pl[i] / s1_m[i], 1 / GCLab[i, 2] ^ 2)
    GCWab[i, 1] ~ dnorm(l[i] / s2_m[gi[i]], 1 / GCWab[i, 2] ^ 2)
  }

  ## Global atmospheric d13C
  for(i in 1:n.steps){
    d13Ca.obs[i, 1] ~ dnorm(d13Ca[i], 1 / d13Ca.obs[i, 2] ^ 2)
  }
  
  # PSMs ----
  ## Alkanones
  for(i in 1:n.obs.alk){
    temp.alk[i] = tempC.alk[si.alk[i], ai.alk[i]] + 273.15
    
    ### Pull K0 and Ksw from driver look up tables
    t.index[i] = min(max(round((tempC.alk[si.alk[i], ai.alk[i]] - tempC.lb) / t.inc) + 1, 1), 
                     n.temp)
    s.index[i] = min(max(round((sal[si.alk[i], ai.alk[i]] - sal.lb) / s.inc) + 1, 1), 
                     n.sal)
    K0[i] = K0a[t.index[i], s.index[i]]
    Ksw.noP[i] = Ksw_sta[t.index[i], s.index[i]]
    
    ### Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem)
    fco2[i] = pCO2[ai.alk[i]] * 0.9968
    co2[i] = fco2[i] * K0[i] * 1e-3 # mol/m^3 (uM)
    
    ### Calculate length of the coccolith from mean radius - add uncertaintiy here?
    len.lith[i] = lith.m * rm[si.alk[i], ai.alk[i]] + lith.b
    
    ### Calculate cell carbon content
    cell.vol.um[i] = 4 / 3 * pi * rm[si.alk[i], ai.alk[i]] ^ 3 # in um^3
    gam.c[i] = 1.46e-14 * cell.vol.um[i]
    
    ### Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) 
    ### from SST using eqn. 8 of Rau et al. (1996)
    DTi[i] = 5.019e-6 * exp(-(19510 / (R.gc * temp.alk[i])))
    DT[i] = DTi[i] * (0.9508 - 7.389e-4 * tempC.alk[si.alk[i], ai.alk[i]])
    
    ### Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 
    ### and 7 of Rau et al. (1996); rate coeffs. follow Zhang et al. (2020)
    k.one[i] = 8718 * (exp(-(62800 / (R.gc * temp.alk[i]))))
    k.two[i] = (680.5 - 4.72 * sal[si.alk[i], ai.alk[i]]) * 1e8 * 
      (exp(-(69400 / (R.gc * temp.alk[i]))))
    hydrox[i] = Ksw.noP[i] / hyd.const
    k.pr[i] = k.one[i] * hydrox[i] + k.two[i]
    rk[i] = sqrt(DT[i] / k.pr[i])
    
    ### Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.err[i] ~ dnorm(0, 1 / mui.y.var)
    mui[i] = coeff.po4 * po4[si.alk[i], ai.alk[i]] + coeff.rm * rm[si.alk[i], ai.alk[i]] + 
      mui.y.int + mui.err[i]

    ### Calculate Qs (the co2 flux into the cell per unit surface area of the cell 
    ### membrane)
    Qr[i] = mui[i] / log(2) * gam.c[i]
    Qs[i] = Qr[i] / (4 * pi * ((rm[si.alk[i], ai.alk[i]] * 1e-6) ^ 2))
    
    ### Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 
    ### 15 of Rau et al (1996)
    b[i] = ((eps.f - eps.d) * Qs[i] * 
              ((rm[si.alk[i], ai.alk[i]] * 1e-6 / 
                  ((1 + rm[si.alk[i], ai.alk[i]] * 1e-6 / rk[i]) * DT[i])) 
               + 1 / P.c)) * 1e3
    
    ### Calculate eps.p from [CO2](aq), eps.f and b
    eps.p[i] = eps.f - (b[i] / (co2[i] * 1e3))
    
    ### Calculate d13C.co2 from d13Ca
    eps.co2.aq.g[i] = -373 / temp.alk[i] + 0.19
    alpha.co2.aq.g[i] = 1 + eps.co2.aq.g[i] / 1e3
    d13C.co2[i] = d13Ca[ai.alk[i]] * alpha.co2.aq.g[i] + eps.co2.aq.g[i]
    
    ### Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
    d13C.biomass[i] = d13C.co2[ai.alk[i]] / (1 + eps.p[i] / 1e3) - eps.p[i]
    
    ### Calculate d13C.marker from d13C.biomass (eqn 22)
    d13C.alk[i] = d13C.biomass[i] / (1 + eps.bob / 1e3) - eps.bob
  }

  ## Soils
  for(i in 1:n.obs.soil){
    ## Atmospheric humidity
    PPCQ[i] = MAP[ai.soil[i]] * PCQ_pf[ai.soil[i]] 
    h_m[i] = min(0.95, 0.25 + 0.7 * (PPCQ[i] / 900))
    ha[i] ~ dbeta(h_m[i] * 100 / (1 - h_m[i]), 100) # PCQ atmospheric humidity
    
    ## Depth to carbonate formation based on Retallack (2005) data, meters
    z.min[i] = MAP[ai.soil[i]] * 0.0925 + 13.4
    z.thick[i] = abs(PPCQ[i] - MAP[ai.soil[i]] / 4) * 0.74 + 17.3
    z.mean[i] = z.min[i] + z.thick[i] / 2
    z.beta[i] = z.mean[i] / (22 ^ 2)
    z.alpha[i] = z.mean[i] * z.beta[i]
    z[i] ~ dgamma(z.alpha[i], z.beta[i])
    z_m[i] = z[i] / 100
    
    ## Soil temperatures at depth z
    Tsoil[i] = temp.soil[ai.soil[i]] + 
      (PCQ_to[ai.soil[i]] * sin(2 * pi * tsc[ai.soil[i]] - z[i] / d)) / 
      exp(z[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    Tair_PCQ[i] = temp.soil[ai.soil[i]] + PCQ_to[ai.soil[i]] * 
      sin(2 * pi * tsc[ai.soil[i]])
    
    ## Potential Evapotranspiration - Hargreaves and Samani (1982) and Turc (1961)
    PET_PCQ_D.1[i] = ifelse(ha[i] < 0.5, 
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * 
                              (23.885 * Rs[si.soil[i]] + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                            0.013 * (Tair_PCQ[i] / (Tair_PCQ[i] + 15)) * 
                              (23.885 * Rs[si.soil[i]] + 50))
    PET_PCQ_D[i] = max(PET_PCQ_D.1[i], 0.01)
    PET_PCQ[i] = PET_PCQ_D[i] * 90
    PET_D.1[i] = ifelse(ha[i] < 0.5, 
                        0.013 * (temp.soil[ai.soil[i]] / (temp.soil[ai.soil[i]] + 15)) * 
                          (23.885 * Rs[si.soil[i]] + 50) * (1 + ((0.5 - ha[i]) / 0.7)),
                        0.013 * (temp.soil[ai.soil[i]] / (temp.soil[ai.soil[i]] + 15)) * 
                          (23.885 * Rs[si.soil[i]] + 50))
    PET_D[i] = max(PET_D.1[i], 0.01)
    PET[i] = PET_D[i] * 365
    
    ## AET in mm/quarter from Budyko curve - Pike (1964)
    AET_var[i] ~ dgamma(1 / 0.2 ^ 2, 1 / 0.2 ^ 2) # noise parameter - Gentine (2012)
    AET_PCQ[i] = PPCQ[i] * (1 / (sqrt(1 + (1 / ((PET_PCQ[i] / (PPCQ[i])) * 
                                                  AET_var[i])) ^ 2)))
    
    ## Average rooting depth
    AI[i] = PET[i] / MAP[ai.soil[i]]
    L[i] = ifelse(AI[i] < 1.4, (-2 * AI[i] ^ 2 + 2.5 * AI[i] + 1) * 100, 60)
    
    ## Carbon isotopes ----
    ### Free air porosity
    FAP.1[i] = min((pore[ai.soil[i]] - ((PPCQ[i] - AET_PCQ[i]) / 
                                     (L[i] * 10 * pore[ai.soil[i]]))), 
                   pore[ai.soil[i]] - 0.05)
    FAP[i] = max(FAP.1[i], 0.01)
    
    ### Soil respiration rate 
    R_PCQ_D_m1[i] = 1.25 * exp(0.05452 * Tair_PCQ[i]) * PPCQ[i] / (127.77 + PPCQ[i])
    R_PCQ_D_m[i] = R_PCQ_D_m1[i] * f_R[ai.soil[i]] # (gC/m2/d)
    R_beta[i] = R_PCQ_D_m[i] / (R_PCQ_D_m[i] * 0.15) ^ 2
    R_alpha[i] = R_PCQ_D_m[i] * R_beta[i]
    R_PCQ_D[i] ~ dgamma(R_alpha[i], R_beta[i])
    
    ### Convert to molC/cm3/s
    R_PCQ_D.1[i] = R_PCQ_D[i] / (12.01 * 100 ^ 2)  # from gC/m2/d to molC/cm2/d
    R_PCQ_S[i] = R_PCQ_D.1[i] / (24 * 3600)  # molC/ cm2 / s
    R_PCQ_S_0[i]= R_PCQ_S[i] / (L[i] * pore[ai.soil[i]]) # Quade et al. (2007)
    
    ### CO2 diffusion
    Dair[i] = 0.1369 * (Tsoil.K[i] / 273.15) ^ 1.958
    DIFC[i] = FAP[i] * tort[si.soil[i]] * Dair[i]
    
    ### S(z)
    k[i] = L[i] / (2 * log(2)) # Respiration characteristic production depth (cm) - Quade (2007)
    S_z_mol[i] = k[i] ^ 2 * R_PCQ_S_0[i] / DIFC[i] * (1 - exp(-z[i] / k[i])) # (mol/cm3)
    S_z[i] = S_z_mol[i] * (0.08206 * Tsoil.K[i] * 1e9) # ppmv
    
    ### d13C of soil-respired CO2
    DD13_water[i] = 25.09 - 1.2 * (MAP[ai.soil[i]] + 975) / 
      (27.2 + 0.04 * (MAP[ai.soil[i]] + 975))
    D13C_plant[i] = (28.26 * 0.22 * (pCO2[ai.soil[i]] + 23.9)) / 
      (28.26 + 0.22 * (pCO2[ai.soil[i]] + 23.9)) - DD13_water[i] # schubert & Jahren (2015)
    D13C_off[i] ~ dnorm(0, 1 / 2 ^ 2) # Noise term
    d13Cr[i] = d13Ca[ai.soil[i]] - D13C_plant[i] + D13C_off[i]
    d13Co[i] = d13Cr[i] + SOM.frac[si.soil[i]]
    
    ### d13C of pedogenic carbonate
    d13Cs[i] = (pCO2[ai.soil[i]] * d13Ca[ai.soil[i]] + S_z[i] * 
                  (1.0044 * d13Cr[i] + 4.4))/(S_z[i] + pCO2[ai.soil[i]])
    d13Cc[i] = ((1 + (11.98 - 0.12 * Tsoil[i]) / 1e3) * (d13Cs[i] + 1e3)) - 1e3
  }
  
  ## Adaxial plants
  for(i in ind.ad){
    ## Pl to obs scaling
    s1_m[i] ~ dgamma(s1[i, 1] * s1.beta[i], s1.beta[i])
    s1.beta[i] = s1[i, 1] / s1[i, 2] ^ 2
    
    # Franks model
    d13C_m[i] = d13Ca[ai.plant[i]] - D13C[i]
    D13C[i] = eps.a + (eps.b - eps.a) * ci[i] / pCO2[ai.plant[i]]
    ci[i] = pCO2[ai.plant[i]] - A[i] / gcop[i] - 1 / meso.scale[gi[i]]
    
    ## Based on data I've seen A should have noise added
    A[i] = (-q.b[i] - sqrt(q.b[i] ^ 2 - 4 * q.a[i] * q.c[i])) / (2 * q.a[i])
    
    q.a[i] = 1 / gcop[i] * (Ci0_m[gi[i]] - gamma)
    q.b[i] = gamma * (-2 * A0_m[gi[i]] / gcop[i] + pCO2[ai.plant[i]] - 
                        1 / meso.scale[gi[i]] - 
                        2 * Ci0_m[gi[i]] + 2 * gamma) + 
      Ci0_m[gi[i]] * (-A0_m[gi[i]] / gcop[i] - pCO2[ai.plant[i]] + 
                           1 / meso.scale[gi[i]])
    q.c[i] = A0_m[gi[i]] * (gamma * (2 * pCO2[ai.plant[i]] - 
                                          2 / meso.scale[gi[i]] - 
                                          Ci0_m[gi[i]] - 2 * gamma) +
                                 Ci0_m[gi[i]] * 
                                 (pCO2[ai.plant[i]] - 1 / meso.scale[gi[i]]))
    
    # Stomatal conductance ----
    # Individual level
    gcop[i] = (1 / gcop.g.ab[i] + 1 / gb_m[gi[i]]) ^ -1 + 
      (1 / gcop.g.ad[i] + 1 / gb_m[gi[i]]) ^ -1
    
    gcop.g.ab[i] = gcmax.ab[i] * gc.scale[gi[i]]
    gcop.g.ad[i] = gcmax.ad[i] * gc.scale[gi[i]]
    
    gcmax.ab[i] = (d.v * D.ab[i] * amax.ab[i]) / 
      (l.ab[i] + ((pi / 2) * sqrt(amax.ab[i] / pi))) / 1.6
    gcmax.ad[i] = (d.v * D.ad[i] * amax.ad[i]) / 
      (l.ad[i] + ((pi / 2) * sqrt(amax.ad[i] / pi))) / 1.6
    
    amax.ab[i] = SA.ab[i] * amax.scale[gi[i]] / D.ab[i] 
    amax.ad[i] = SA.ad[i] * amax.scale[gi[i]] / D.ad[i] 
    
    # Stomatal geometry calculations ----
    # Individual level
    ## Stomatal density
    D.ab[i] = SA.ab[i] / (pi * (Pl.ab[i] / 2) ^ 2)
    D.ad[i] = SA.ad[i] / (pi * (Pl.ad[i] / 2) ^ 2)
    
    ## Stomatal pore area
    SA.ab[i] = SA_gc.ab[i] * PL_GCL[i] ^ 2
    SA.ad[i] = SA_gc.ad[i] * PL_GCL[i] ^ 2
    
    SA_gc.ab[i] ~ dbeta(1.5, 20)
    SA_gc.ad[i] ~ dbeta(1.5, 20)
    
    ## Pore length
    Pl.ab[i] = gcl_m.ab[i] * PL_GCL[i] * 1e-6
    Pl.ad[i] = gcl_m.ad[i] * PL_GCL[i] * 1e-6
    
    gcl_m.ab[i] ~ dgamma(2, 0.08)
    gcl_m.ad[i] ~ dgamma(2, 0.08)
    
    ## Pore depth   
    l.ab[i] ~ dunif(1e-6, 1e-4)
    l.ad[i] ~ dunif(1e-6, 1e-4)
  }
  
  ## Abaxial plants
  for(i in ind.ab){
    ## Pl to obs scaling
    s1_m[i] ~ dgamma(s1[i, 1] * s1.beta[i], s1.beta[i])
    s1.beta[i] = s1[i, 1] / s1[i, 2] ^ 2
    
    # Franks model
    d13C_m[i] = d13Ca[ai.plant[i]] - D13C[i]
    D13C[i] = eps.a + (eps.b - eps.a) * ci[i] / pCO2[ai.plant[i]]
    ci[i] = pCO2[ai.plant[i]] - A[i] / gcop[i] - 1 / meso.scale[gi[i]]
    
    ## Based on data I've seen A should have noise added
    A[i] = (-q.b[i] - sqrt(q.b[i] ^ 2 - 4 * q.a[i] * q.c[i])) / (2 * q.a[i])
    
    q.a[i] = 1 / gcop[i] * (Ci0_m[gi[i]] - gamma)
    q.b[i] = gamma * (-2 * A0_m[gi[i]] / gcop[i] + pCO2[ai.plant[i]] - 
                        1 / meso.scale[gi[i]] - 
                        2 * Ci0_m[gi[i]] + 2 * gamma) + 
      Ci0_m[gi[i]] * (-A0_m[gi[i]] / gcop[i] - pCO2[ai.plant[i]] + 
                           1 / meso.scale[gi[i]])
    q.c[i] = A0_m[gi[i]] * (gamma * (2 * pCO2[ai.plant[i]] - 
                                          2 / meso.scale[gi[i]] - 
                                          Ci0_m[gi[i]] - 2 * gamma) +
                                 Ci0_m[gi[i]] * 
                                 (pCO2[ai.plant[i]] - 1 / meso.scale[gi[i]]))
    
    # Stomatal conductance ----
    # Individual level
    gcop[i] = (1 / gcop.g[i] + 1 / gb_m[gi[i]]) ^ -1
    gcop.g[i] = gcmax[i] * gc.scale[gi[i]]
    gcmax[i] = (d.v * D[i] * amax[i]) / (l[i] + ((pi / 2) * sqrt(amax[i] / pi))) / 1.6
    amax[i] = SA[i] * amax.scale[gi[i]] / D[i] 
    
    # Stomatal geometry calculations ----
    # Individual level
    ## Stomatal density
    D[i] = SA[i] / (pi * (Pl[i] / 2) ^ 2)
    
    ## Stomatal pore area
    SA[i] = SA_gc[i] * PL_GCL[i] ^ 2
    SA_gc[i] ~ dbeta(1.5, 20)
    
    ## Pore length
    Pl[i] = gcl_m[i] * PL_GCL[i] * 1e-6
    gcl_m[i] ~ dgamma(2, 0.08)
    
    ## Pore depth   
    l[i] ~ dunif(1e-6, 1e-4)
  }
  
  # Environmental model ----
  ## Site-dependent alkenone
  for(i in 1:n.sites.alk){
    for(j in 2:n.steps){
      ### Temperature (degrees C)
      tempC.alk[i, j] = MAT[j] + toff.alk[i, j]
      toff.alk[i, j] = toff.alk[i, j - 1] + toff.eps.alk[i, j]
      toff.eps.alk[i, j] ~ dnorm(toff.eps.alk[i, j - 1] * (toff.phi.alk ^ dt), 
                                 toff.pc.alk)
      
      ### Salinity (ppt)
      sal[i, j] ~ dgamma(sal[i, j - 1] ^ 2 / sal.v, sal[i, j - 1] / sal.v)
      
      ### Concentration of phosphate (PO4; umol/kg)
      po4[i, j] ~ dgamma(po4[i, j - 1] ^ 2 / po4.v, po4[i, j - 1] / po4.v)T(0.1, 2)
      
      ### Mean cell radius (m)
      rm[i, j] ~ dgamma(rm[i, j - 1] ^ 2 / rm.v, rm[i, j - 1] / rm.v)T(1, 5)
    }
    
    ### Initial conditions
    toff.alk[i, 1] ~ dnorm(10, 1 / 4 ^ 2)
    toff.eps.alk[i, 1] = 0
    
    sal[i, 1] ~ dgamma(35 ^ 2 / 4, 35 / 4)
    
    po4[i, 1] ~ dunif(0.1, 1.5)
    
    rm[i, 1] ~ dgamma(1.5 ^ 2 / 0.25, 1.5 / 0.25)
  }
  
  ### Priors
  toff.pc.alk = toff.tau.alk * ((1 - toff.phi.alk ^ 2) / (1 - toff.phi.alk ^ (2 * dt)))
  toff.tau.alk ~ dgamma(10, 0.1)
  toff.phi.alk ~ dbeta(2, 5)
  
  sal.v ~ dgamma(10, 100)
  
  po4.v ~ dgamma(10, 500)
  
  rm.v ~ dgamma(10, 1000)
  
  ## Time-dependent soils
  for(i in 2:n.steps){
    temp.soil[i] = MAT[i] + toff.soil[i]
    toff.soil[i] = toff.soil[i - 1] + toff.eps.soil[i]
    toff.eps.soil[i] ~ dnorm(toff.eps.soil[i - 1] * (toff.phi.soil ^ dt), 
                               toff.pc.soil)
      
    PCQ_to[i] = PCQ_to[i - 1] + PCQ_to.eps[i]
    PCQ_to.eps[i] ~ dnorm(PCQ_to.eps[i - 1] * (PCQ_to.phi ^ dt), PCQ_to.pc)
    
    MAP[i] = MAP[i - 1] * (1 + MAP.eps[i])
    MAP.eps[i] ~ dnorm(MAP.eps[i - 1] * (MAP.phi ^ dt), MAP.pc)T(-0.99,)
    
    PCQ_pf[i] ~ dbeta(PCQ_pf[i - 1] * PCQ_pf.v, (1 - PCQ_pf[i - 1]) * PCQ_pf.v)

    tsc[i] = tsc[i - 1] + tsc.eps[i]
    tsc.eps[i] ~ dnorm(tsc.eps[i - 1] * (tsc.phi ^ dt), tsc.pc)
    
    f_R[i] ~ dbeta(f_R[i - 1] * f_R.v, (1 - f_R[i - 1]) * f_R.v)

    pore[i] ~ dbeta(pore[i - 1] * pore.v, (1 - pore[i - 1]) * pore.v)    
  }
  
  ### Initial conditions
  toff.soil[1] ~ dnorm(10, 1 / 4 ^ 2)
  toff.eps.soil[1] = 0
  
  PCQ_to[1] ~ dunif(7, 15) # PCQ temperature offset, C
  PCQ_to.eps[1] = 0
  
  MAP[1] ~ dunif(1e2, 1e3) # mean annual precipitation, mm
  MAP.eps[1] = 0
  
  PCQ_pf[1] ~ dunif(0.02, 0.25) # PCQ precipitation fraction

  tsc[1] ~ dunif(0, 0.5) # seasonal offset of PCQ for thermal diffusion
  tsc.eps[1] = 0
  
  f_R[1] ~ dbeta(2, 16) # ratio of PCQ to mean annual respiration rate

  pore[1] ~ dunif(0.45, 0.54) # soil porosity

  ### Priors
  toff.pc.soil = toff.tau.soil * ((1 - toff.phi.soil ^ 2) / (1 - toff.phi.soil ^ (2 * dt)))
  toff.tau.soil ~ dgamma(10, 0.1)
  toff.phi.soil ~ dbeta(2, 5)
  
  PCQ_to.pc = PCQ_to.tau * ((1 - PCQ_to.phi ^ 2) / (1 - PCQ_to.phi ^ (2 * dt)))
  PCQ_to.tau ~ dgamma(10, 1e-1)
  PCQ_to.phi ~ dbeta(2, 5)
  
  MAP.pc = MAP.tau * ((1 - MAP.phi ^ 2) / (1 - MAP.phi ^ (2 * dt)))
  MAP.tau ~ dgamma(10, 1e-2) # percentage
  MAP.phi ~ dbeta(2, 5)
  
  PCQ_pf.v ~ dgamma(10, 1e-3) # percentage

  tsc.pc = tsc.tau * ((1 - tsc.phi ^ 2) / (1 - tsc.phi ^ (2 * dt)))
  tsc.tau ~ dgamma(10, 1e-6)
  tsc.phi ~ dbeta(2, 5)
  
  f_R.v ~ dgamma(10, 1e-3) # percentage

  pore.v ~ dgamma(1, 1e-4) # was 1e-2

  ## Site-dependent soils
  for(i in 1:n.sites.soil){
    Ra[i] = 42.608 - 0.3538 * abs(lat.soil[i]) # total radiation at the top of the atmosphere
    Rs[i] = Ra[i] * 0.16 * sqrt(12) # daily temperature range assumed to be 12
    tort[i] ~ dbeta(0.7 * 100 / 0.3, 100) # soil tortuosity
    SOM.frac[i] ~ dunif(0.5, 1.5)
  }
  
  ## Taxon-dependent plant
  # Taxon priors ----
  for(i in 1:n.gen){
    s2_m[i] ~ dgamma(s2[i, 1] * s2.beta[i], s2.beta[i])
    s2.beta[i] = s2[i, 1] / s2[i, 2] ^ 2
    
    amax.scale[i] ~ dbeta(s3[i, 1] * amax.v[i], (1 - s3[i, 1]) * amax.v[i]) I (0.0001, 0.9999) # aka s3
    amax.v[i] = (s3[i, 1] * (1 - s3[i, 1])) / s3[i, 2] ^ 2 - 1
    
    gc.scale[i] ~ dbeta(s4[i, 1] * gc.v[i], (1 - s4[i, 1]) * gc.v[i]) # aka s4
    gc.v[i] = (s4[i, 1] * (1 - s4[i, 1])) / s4[i, 2] ^ 2 - 1
    
    meso.scale[i] ~ dbeta(s5[i, 1] * meso.v[i], (1 - s5[i, 1]) * meso.v[i]) # aka s5
    meso.v[i] = (s5[i, 1] * (1 - s5[i, 1])) / s5[i, 2] ^ 2 - 1
    
    Ci0_m[i] ~ dgamma(Ci0[i, 1] * Ci0.beta[i], Ci0.beta[i])
    Ci0.beta[i] = Ci0[i, 1] / Ci0[i, 2] ^ 2
    
    A0_m[i] ~ dgamma(A0[i, 1] * A0.beta[i], A0.beta[i])
    A0.beta[i] = A0[i, 1] / A0[i, 2] ^ 2
    
    gb_m[i] ~ dgamma(gb[i, 1] * gb.beta[i], gb.beta[i])
    gb.beta[i] = gb[i, 1] / gb[i, 2] ^ 2
  }
  
  ## Global
  for(i in 2:n.steps){
    ### pCO2 (uatm)
#    pCO2[i] ~ dgamma(pCO2[i - 1] ^ 2 / pCO2.v, pCO2[i - 1] / pCO2.v)
    pCO2[i] = pCO2[i-1] + pCO2.eps[i]
    pCO2.eps[i] ~ dnorm(pCO2.eps[i-1] * pCO2.eps.ac, pCO2.pre) 
    
    ### d13Ca
    d13Ca[i] = d13Ca[i - 1] + d13Ca.eps[i]
    d13Ca.eps[i] ~ dnorm(d13Ca.eps[i - 1] * (d13Ca.phi ^ dt), d13Ca.pc)
    
    ### Global MAT in degrees C
    MAT[i] = MAT[i - 1] + MAT.eps[i]
    MAT.eps[i] ~ dnorm(MAT.eps[i - 1] * MAT.phi ^ dt, MAT.pc)
  }
  
  ### Initial conditions
  pCO2[1] = pCO2.s * 1e3
  pCO2.s ~ dunif(0.25, 1)
  pCO2.eps[1] = 0
  
  d13Ca.eps[1] = 0
  d13Ca[1] ~ dunif(-8, -3)
  
  MAT[1] ~ dunif(10, 20)
  MAT.eps[1] = 0
  
  ### Priors
#  pCO2.v ~ dgamma(10, 1e-2)
  pCO2.eps.ac ~ dunif(0.5, 0.99)
  pCO2.pre ~ dgamma(1, 0.1)
  
  d13Ca.pc = d13Ca.tau * ((1 - d13Ca.phi ^ 2) / (1 - d13Ca.phi ^ (2 * dt)))
  d13Ca.phi ~ dbeta(5, 2)
  d13Ca.tau ~ dgamma(5, 1e-2)
  
  MAT.pc = MAT.tau * ((1 - MAT.phi ^ 2) / (1 - MAT.phi ^ (2 * dt)))
  MAT.tau ~ dgamma(10, 1)
  MAT.phi ~ dbeta(2, 5)
  
  ############################################################################################
  # Shared constants and calibration parameters
  ############################################################################################
  
  # eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰), 
  # post-calibration
  eps.f ~ dnorm(26.6, 1 / 0.3 ^ 2)
  
  # epsilon of diffusive transport of CO2(aq) in water, post-calibration
  eps.d ~ dnorm(0.656, 1 / 0.1 ^ 2)
  
  # epsilon of biomass/biomarker (C isotope fractionation b/w algal biomass and 
  # biomarkers, or bob [b over b]), post-calibration
  eps.bob ~ dnorm(4.4, 1 / 0.2 ^ 2)
  
  # cell wall permeability to CO2(aq) in m/s, post-calibration
  P.c ~ dgamma(4.539e-5 ^ 2 / 1.4e-6 ^ 2, 4.539e-5 / 1.4e-6 ^ 2)
  
  # gas constant in J / K*mol
  R.gc = 8.3143
  
  # pH value for calculating rk - held constant here; varying this has almost 
  # zero effect on the model
  hyd.const = 10 ^ (-8)
  
  # Coefficient for mu(i) - PO4 multi linear regression, post-calibration
  coeff.po4 = 1.998e-6
  
  # Coefficient for radius - PO4 multi linear regression, post-calibration
  coeff.rm = -20e-6
  
  # Coefficient for y intercept for mu(i) multi linear regression
  mui.y.int = 4.011e-5
  
  # Approximate error variance for mu(i) calibration
  mui.y.var = 2.5e-12 
  
  # Diffusivity of CO2 in air
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3))
  
  # Isotope ratio of VPDB
  R13.VPDB = 0.011237
  
  pi = 3.14159265
  
  eps.a = 4.4
  
  d.v = 0.000940096
  
  eps.b = 30
  
  gamma = 40
  
}
