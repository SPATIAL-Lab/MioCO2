model {
  
  ############################################################################################
  # Discrete proxy data likelihoods
  ############################################################################################
  # The observation data are still evaluated row by row, but each observation points to a
  # spatial/temporal state. The expensive PSM is evaluated once per state below.
  for(i in 1:n.obs){
    d13Cmarker.data[i] ~ dnorm(d13Cmarker[i], 1 / d13Cmarker.data.sd[i] ^ 2)

    d13Ca.obs[i, 1] ~ dnorm(d13Ca[i], 1 / d13Ca.obs[i, 2] ^ 2)
  }
  
  for(i in 1:n.lith){
    len.lith.data[i] ~ dnorm(len.lith[lith.state.index[i]], 
                             1 / len.lith.data.sd[i] ^ 2)
  }
  ############################################################################################
  
  
  ############################################################################################
  # Proxy System Model - discrete proxy states
  ############################################################################################
  for (i in 1:n.obs){
    # Use temperature directly from the intermediate-file state prior
    temp[i] = tempC[i] + 273.15
    
    # Pull K0 and Ksw from driver look up tables
    t.index[i] = min(max(round((tempC[i]-tempC.lb)/t.inc)+1, 1), n.temp)
    s.index[i] = min(max(round((sal[i]-sal.lb)/s.inc)+1, 1), n.sal)
    K0[i] = K0a[t.index[i], s.index[i]]
    Ksw.noP[i] = Ksw_sta[t.index[i], s.index[i]]
    
    # Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem)
    fco2[i] = pCO2[i] * 0.9968
    co2[i] = fco2[i] * K0[i] * 1e-3 # mol/m^3 (uM)
    
    # Calculate length of the coccolith from mean radius (rm)
    len.lith[i] = lith.m * (rm[i] * 1e6) + lith.b # in um
    
    # Calculate cell carbon content
    cell.vol.um[i] = 4 / 3 * 3.141593 * ((rm[i] * 1e6) ^ 3) # in um^3
    gam.c[i] = 1.46e-14 * cell.vol.um[i]
    
    # Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
    DTi[i] = 5.019e-6 * exp(-(19510 / (R.gc * temp[i])))
    DT[i] = DTi[i] * (0.9508 - 7.389e-4 * tempC[i])
    
    # Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996); rate coeffs. follow Zhang et al. (2020)
    k.one[i] = 8718 * (exp(-(62800 / (R.gc * temp[i]))))
    k.two[i] = (680.5 - 4.72 * sal[i]) * 1e8 * (exp(-(69400 / (R.gc * temp[i]))))
    hydrox[i] = Ksw.noP[i] / hyd.const
    k.pr[i] = k.one[i] * hydrox[i] + k.two[i]
    rk[i] = sqrt(DT[i] / k.pr[i])
    
    # Calculate foram d13C from env parm d13C.co2 (retained as a derived value, not a discrete-data likelihood)
    #eps.aog[i] = (-373 / temp[i]) + 0.19 # fractionation b/w CO2 (aq) and CO2 (g)
    #d13C.co2g[i] = ((d13C.co2[i] + 1000) / (eps.aog[i] / 1000 + 1)) - 1000
    #d13Cpf[i] = ((11.98 - 0.12 * tempC[i]) / 1000 + 1)*(d13C.co2g[i] + 1000) - 1000
    
    # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.err[i] ~ dnorm(0, 1 / mui.y.var)
    mui[i] = coeff.po4 * po4[i] + coeff.rm * rm[i] + mui.y.int + mui.err[i]
#    mui[i] = max(mui.raw[i], 1*10^-6)
    
    # Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
    Qr[i] = mui[i] / log(2) * gam.c[i]
    Qs[i] = Qr[i] / (4 * 3.141593 * (rm[i] ^ 2))
    
    # Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996)
    b[i] = ((eps.f - eps.d) * Qs[i] * ((rm[i] / ((1 + rm[i] / rk[i]) * DT[i])) 
                                       + 1 / P.c)) * 1e3
    
    # Calculate eps.p from [CO2](aq), eps.f and b
    eps.p[i] = eps.f - (b[i] / (co2[i] * 1e3))
    
    # Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
    d13C.biomass[i] = ((d13C.co2[i] + 1e3) / ((eps.p[i] / 1e3) + 1)) - 1e3
    
    # Calculate d13C.marker from d13C.biomass (eqn 22)
    d13Cmarker[i] = ((d13C.biomass[i] + 1e3) / ((eps.bob/1e3) + 1)) - 1e3
    
    # Calculate d13Ca from d13C.co2
    eps.co2.aq.g[i] = -373 / temp[i] + 0.19
    alpha.co2.aq.g[i] = 1 + eps.co2.aq.g[i] / 1e3
    d13Ca[i] = (d13C.co2[i] + 1e3) / alpha.co2.aq.g[i] - 1e3
    
  }
  ############################################################################################
  
  
  ############################################################################################
  # Environmental prior distributions - discrete proxy states
  ############################################################################################
  for (i in 1:n.obs){
    # Temperature (degrees C)
    tempC[i] ~ dnorm(tempC.m[i], tempC.p[i])
    
    # Salinity (ppt)
    sal[i] ~ dgamma(sal.m[i] ^ 2 / sal.v[i], sal.m[i] / sal.v[i])
    
    # pCO2 (uatm)
    pCO2[i] = ca.s[i] * 1e3
    ca.s[i] ~ dunif(0.1, 2)  
    
    # d13C of aqueous CO2 (per mille)
    d13C.co2[i] ~ dnorm(d13C.co2.m[i], d13C.co2.p[i])
    
    # Concentration of phosphate (PO4; umol/kg)
    po4[i] ~ dgamma(po4.m[i] ^ 2 / po4.v[i], po4.m[i] / po4.v[i])
    
    # Mean cell radius (m)
    rm[i] ~ dgamma(rm.m[i] ^ 2 / rm.v[i], rm.m[i] / rm.v[i])
  }
  ############################################################################################
 
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
  
  # Uk'37 temperature calibration (Conte et al., 2006; sediment - AnnO linear model) 
  # with 1.1C = reported se of estimation
  Uk.sl = 29.876
  Uk.int = 1.334
  Uk.cal.se = 1.1
  
  # gas constant in J / K*mol
  R.gc = 8.3143
  
  # pH value for calculating rk - held constant here; varying this has almost 
  # zero effect on the model
  pH.const = 8
  hyd.const = 10 ^ (-pH.const)
  
  # Coefficient for mu(i) - PO4 multi linear regression, post-calibration
  coeff.po4 = 1.998e-6
  
  # Coefficient for radius - PO4 multi linear regression, post-calibration
  coeff.rm = -20
  
  # Coefficient for y intercept for mu(i) multi linear regression
  mui.y.int = 4.011e-5
  
  # Approximate error variance for mu(i) calibration
  mui.y.var = 2.5e-12 
  
  
}
