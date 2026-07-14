model {
  
  ############################################################################################
  # Shared constants and calibration parameters
  ############################################################################################
  
  # eps.f = max fractionation by RuBisCO at infinite CO2 (generally 25 to 28‰)
  eps.f ~ dnorm(27, 1/0.3^2)
  
  # epsilon of diffusive transport of CO2(aq) in water
  eps.d ~ dnorm(0.7, 1 / 0.1 ^ 2)
  
  # epsilon of biomass/biomarker (C isotope fractionation b/w algal biomass and biomarkers, or bob [b over b])
  eps.bob ~ dnorm(4.80, 1 / 0.2 ^ 2)
  
  # cell wall permeability to CO2(aq) in m/s taken from Zhang et al. (2020)
  P.c ~ dgamma(5.09e-5 ^ 2 / 0.16e-5 ^ 2, 5.09e-5 / 0.16e-5 ^ 2)
  
  # Uk'37 temperature calibration (Conte et al., 2006; sediment - AnnO linear model) with 1.1C = reported se of estimation
  Uk.sl = 29.876
  Uk.int = 1.334
  Uk.cal.se = 1.1
  
  # gas constant in J / K*mol
  R.gc = 8.3143
  
  # pH value for calculating rk - held constant here; varying this has almost zero effect on the model
  pH.const = 8
  hyd.const = 10 ^ (-pH.const)
  
  # Coefficient for mu(i) - PO4 multi linear regression
  coeff.po4 ~ dnorm(po4.co.lr, 1/(po4.se.lr)^2)
  
  # Coefficient for radius - PO4 multi linear regression
  coeff.rm ~ dnorm(r.co.lr, 1/(r.se.lr)^2)
  
  # Coefficient for y intercept for mu(i) multi linear regression
  mui.y.int ~ dnorm(y.int.lr, 1/(y.int.se.lr)^2)
  ############################################################################################
  
  
  ############################################################################################
  # Culture calibration data
  ############################################################################################
  sigma.mui.cd ~ dunif(0, 1)
  tau.mui.cd = 1/sigma.mui.cd^2
  
  for (i in 1:n.cd){
    # Likelihood function
    mui.cd.data[i] ~ dnorm(mui.cd[i], tau.mui.cd)
    radius.cd[i] ~ dnorm(rm.cd[i], 1/(0.25*10^-6)^2)
    po4.cd.data[i] ~ dnorm(po4.cd[i], 1/0.05^2)
    
    # Prior distributions
    rm.cd[i] ~ dnorm(1.5*10^-6, 1/(0.5*10^-6)^2)T(0,5*10^-6)
    po4.cd[i] ~ dnorm(po4.m.cd, 1/0.25^2)T(0,2)
    
    # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.raw.cd[i] = coeff.po4*po4.cd[i] + coeff.rm*rm.cd[i] + mui.y.int
    mui.cd[i] = max(mui.raw.cd[i], 1*10^-6)
  }
  ############################################################################################
  
  
  ############################################################################################
  # GIG calibration data likelihoods
  ############################################################################################
  for (i in 1:n.gig){
    d13Cmarker.data.gig[i] ~ dnorm(d13Cmarker.gig[i], d13Cmarker.p.gig[i])
    d13Cmarker.p.gig[i] = 1/d13Cmarker.data.sd.gig[i]^2
    
    d13Cpf.data.gig[i] ~ dnorm(d13Cpf.gig[i], d13Cpf.p.gig[i])
    d13Cpf.p.gig[i] = 1/d13Cpf.data.sd.gig[i]^2
    
    len.lith.data.gig[i] ~ dnorm(len.lith.gig[i], len.lith.p.gig[i])
    len.lith.p.gig[i] = 1/len.lith.data.sd.gig[i]^2
    
    Uk.data.gig[i] ~ dnorm(Uk.gig[i], Uk.p.gig[i])
    Uk.p.gig[i] = 1/Uk.data.sd.gig[i]^2
    
    iceco2.data.gig[i] ~ dnorm(pco2.gig[i], 1/6^2)
  }
  ############################################################################################
  
  
  ############################################################################################
  # Discrete proxy data likelihoods
  ############################################################################################
  # The observation data are still evaluated row by row, but each observation points to a
  # spatial/temporal state. The expensive PSM is evaluated once per state below.
  for (i in 1:n.marker){
    d13Cmarker.data[i] ~ dnorm(d13Cmarker[marker.state.index[i]], d13Cmarker.p[i])
    d13Cmarker.p[i] = 1/d13Cmarker.data.sd[i]^2
  }
  
  for (i in 1:n.lith){
    len.lith.data[i] ~ dnorm(len.lith[lith.state.index[i]], len.lith.p[i])
    len.lith.p[i] = 1/len.lith.data.sd[i]^2
  }
  ############################################################################################
  
  
  ############################################################################################
  # Proxy System Model - GIG calibration data
  ############################################################################################
  for (i in 1:n.gig){
    # Calculate Uk'37 from temperature
    Uk.gig[i] = (tempC.gig[i] + Uk.int)/Uk.sl
    temp.gig[i] = tempC.gig[i] + 273.15
    
    # Pull K0 and Ksw from driver look up tables
    t.index.gig[i] = min(max(round((tempC.gig[i]-tempC.lb)/t.inc)+1, 1), n.temp)
    s.index.gig[i] = min(max(round((sal.gig[i]-sal.lb)/s.inc)+1, 1), n.sal)
    K0.gig[i] = K0a[t.index.gig[i], s.index.gig[i]]
    Ksw.noP.gig[i] = Ksw_sta[t.index.gig[i], s.index.gig[i]]
    
    # Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem)
    fco2.gig[i] = pco2.gig[i]*0.9968
    co2.gig[i] = fco2.gig[i]*K0.gig[i]*10^-3 # mol/m^3 (uM)
    
    # Calculate length of the coccolith from mean radius (rm)
    len.lith.gig[i] = lith.m*(rm.gig[i]*10^6) + lith.b # in um
    
    # Calculate cell carbon content
    cell.vol.gig[i] = 4/3*3.141593*((rm.gig[i])^3) # in m^3
    cell.vol.um.gig[i] = 4/3*3.141593*((rm.gig[i]*10^6)^3) # in um^3
    gam.c.gig[i] = 1.46*10^-14*cell.vol.um.gig[i]
    
    # Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
    DTi.gig[i] = 5.019*(10^-6)*exp(-(19510/(R.gc*temp.gig[i])))
    DT.gig[i] = DTi.gig[i]*(0.9508 - 7.389*(10^-4)*tempC.gig[i])
    
    # Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996); rate coeffs. follow Zhang et al. (2020)
    k.one.gig[i] = 8718*(exp(-(62800/(R.gc*temp.gig[i]))))
    k.two.gig[i] = (680.5-4.72*sal.gig[i])*10^8*(exp(-(69400/(R.gc*temp.gig[i]))))
    hydrox.gig[i] = Ksw.noP.gig[i] / hyd.const
    k.pr.gig[i] = k.one.gig[i]*hydrox.gig[i] + k.two.gig[i]
    rk.gig[i] = sqrt(DT.gig[i]/k.pr.gig[i])
    
    # Calculate foram d13C from env parm d13C.co2 (rearranged eqns 19 and 20)
    eps.aog.gig[i] = (-373/temp.gig[i]) + 0.19 # fractionation b/w CO2 (aq) and CO2 (g)
    d13C.co2g.gig[i] = ((d13C.co2.gig[i] + 1000) / (eps.aog.gig[i]/1000 + 1)) - 1000
    d13Cpf.gig[i] = ((11.98 - 0.12*tempC.gig[i])/1000 + 1)*(d13C.co2g.gig[i] + 1000) - 1000
    
    # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.raw.gig[i] = coeff.po4*po4.gig[i] + coeff.rm*rm.gig[i] + mui.y.int
    mui.gig[i] = max(mui.raw.gig[i], 1*10^-6)
    
    # Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
    Qr.gig[i] = mui.gig[i]/log(2) * gam.c.gig[i]
    Qs.gig[i] = Qr.gig[i] / (4*3.141593*(rm.gig[i]^2))
    
    # Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996)
    b.gig[i] = ((eps.f - eps.d) * Qs.gig[i] * ((rm.gig[i]/((1 + rm.gig[i]/rk.gig[i])*DT.gig[i])) + 1/P.c)) * 10^3
    
    # Calculate eps.p from [CO2](aq), eps.f and b
    eps.p.gig[i] = eps.f - (b.gig[i]/(co2.gig[i]*10^3))
    
    # Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
    d13C.biomass.gig[i] = ((d13C.co2.gig[i] + 1000) / ((eps.p.gig[i] / 1000) + 1)) - 1000
    
    # Calculate d13C.marker from d13C.biomass (eqn 22)
    d13Cmarker.gig[i] = ((d13C.biomass.gig[i] + 1000) / ((eps.bob/1000) + 1)) - 1000
  }
  ############################################################################################
  
  
  ############################################################################################
  # Environmental prior distributions - GIG calibration data
  ############################################################################################
  for (i in 1:n.gig){
    # Temperature (degrees C)
    tempC.gig[i] ~ dnorm(tempC.m.gig, tempC.p.gig)
    
    # Salinity (ppt)
    sal.gig[i] ~ dnorm(sal.m.gig, sal.p.gig)T(0,)
    
    # pCO2 (uatm)
    pco2.gig[i] ~ dnorm(pco2.m.gig, pco2.p.gig)T(125,375)
    
    # d13C of aqueous CO2 (per mille)
    d13C.co2.gig[i] ~ dnorm(d13C.co2.m.gig, d13C.co2.p.gig)
    
    # Concentration of phosphate (PO4; umol/kg)
    po4.gig[i] ~ dnorm(po4.m.gig[site.index.gig[i]], po4.p.gig)T(0,2)
    
    # Mean cell radius (m)
    rm.gig[i] ~ dnorm(rm.m.gig, rm.p.gig)T(0,5*10^-6)
  }
  ############################################################################################
  
  
  ############################################################################################
  # Proxy System Model - discrete proxy states
  ############################################################################################
  for (i in 1:n.state){
    # Use temperature directly from the intermediate-file state prior
    temp[i] = tempC[i] + 273.15
    
    # Pull K0 and Ksw from driver look up tables
    t.index[i] = min(max(round((tempC[i]-tempC.lb)/t.inc)+1, 1), n.temp)
    s.index[i] = min(max(round((sal[i]-sal.lb)/s.inc)+1, 1), n.sal)
    K0[i] = K0a[t.index[i], s.index[i]]
    Ksw.noP[i] = Ksw_sta[t.index[i], s.index[i]]
    
    # Calculate aqueous CO2 from atmospheric CO2 with Henry's law (carb chem)
    fco2[i] = pco2[i]*0.9968
    co2[i] = fco2[i]*K0[i]*10^-3 # mol/m^3 (uM)
    
    # Calculate length of the coccolith from mean radius (rm)
    len.lith[i] = lith.m*(rm[i]*10^6) + lith.b # in um
    
    # Calculate cell carbon content
    cell.vol[i] = 4/3*3.141593*((rm[i])^3) # in m^3
    cell.vol.um[i] = 4/3*3.141593*((rm[i]*10^6)^3) # in um^3
    gam.c[i] = 1.46*10^-14*cell.vol.um[i]
    
    # Calculate DT (temperature-sensitive diffusivity of C02(aq)in seawater) from SST using eqn. 8 of Rau et al. (1996)
    DTi[i] = 5.019*(10^-6)*exp(-(19510/(R.gc*temp[i])))
    DT[i] = DTi[i]*(0.9508 - 7.389*(10^-4)*tempC[i])
    
    # Calculate rK (reacto-diffusive length) from SSS and SST using eqns. 6 and 7 of Rau et al. (1996); rate coeffs. follow Zhang et al. (2020)
    k.one[i] = 8718*(exp(-(62800/(R.gc*temp[i]))))
    k.two[i] = (680.5-4.72*sal[i])*10^8*(exp(-(69400/(R.gc*temp[i]))))
    hydrox[i] = Ksw.noP[i] / hyd.const
    k.pr[i] = k.one[i]*hydrox[i] + k.two[i]
    rk[i] = sqrt(DT[i]/k.pr[i])
    
    # Calculate foram d13C from env parm d13C.co2 (retained as a derived value, not a discrete-data likelihood)
    eps.aog[i] = (-373/temp[i]) + 0.19 # fractionation b/w CO2 (aq) and CO2 (g)
    d13C.co2g[i] = ((d13C.co2[i] + 1000) / (eps.aog[i]/1000 + 1)) - 1000
    d13Cpf[i] = ((11.98 - 0.12*tempC[i])/1000 + 1)*(d13C.co2g[i] + 1000) - 1000
    
    # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.raw[i] = coeff.po4*po4[i] + coeff.rm*rm[i] + mui.y.int
    mui[i] = max(mui.raw[i], 1*10^-6)
    
    # Calculate Qs (the co2 flux into the cell per unit surface area of the cell membrane)
    Qr[i] = mui[i]/log(2) * gam.c[i]
    Qs[i] = Qr[i] / (4*3.141593*(rm[i]^2))
    
    # Calculate b in uM from Qs, rmean, rK, DT, Pc, eps.f and eps.d using eqn. 15 of Rau et al (1996)
    b[i] = ((eps.f - eps.d) * Qs[i] * ((rm[i]/((1 + rm[i]/rk[i])*DT[i])) + 1/P.c)) * 10^3
    
    # Calculate eps.p from [CO2](aq), eps.f and b
    eps.p[i] = eps.f - (b[i]/(co2[i]*10^3))
    
    # Calculate d13C.biomass from d13C.co2 and epsilon.p (eqn 2)
    d13C.biomass[i] = ((d13C.co2[i] + 1000) / ((eps.p[i] / 1000) + 1)) - 1000
    
    # Calculate d13C.marker from d13C.biomass (eqn 22)
    d13Cmarker[i] = ((d13C.biomass[i] + 1000) / ((eps.bob/1000) + 1)) - 1000
  }
  ############################################################################################
  
  
  ############################################################################################
  # Environmental prior distributions - discrete proxy states
  ############################################################################################
  for (i in 1:n.state){
    # Temperature (degrees C)
    tempC[i] ~ dnorm(tempC.m[i], tempC.p[i])
    
    # Salinity (ppt)
    sal[i] ~ dnorm(sal.m[i], sal.p[i])T(0,)
    
    # pCO2 (uatm)
    pco2[i] ~ dnorm(pco2.m[i], pco2.p[i])T(pco2.lb[i], pco2.ub[i])
    
    # d13C of aqueous CO2 (per mille)
    d13C.co2[i] ~ dnorm(d13C.co2.m[i], d13C.co2.p[i])
    
    # Concentration of phosphate (PO4; umol/kg)
    po4[i] ~ dnorm(po4.m[i], po4.p[i])T(0,2)
    
    # Mean cell radius (m)
    rm[i] ~ dnorm(rm.m[i], rm.p[i])T(0,5*10^-6)
  }
  ############################################################################################
  
}
