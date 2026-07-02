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
  for (i in 1:n.cd){
    # Likelihood function
    mui.cd.data[i] ~ dnorm(mui.cd[i], 1 / 2.6e-12) # variance of lm residuals

    # Calculate instantaneous growth rate (mu,i) from [PO4] and rmean
    mui.cd[i] = coeff.po4*po4.cd[i] + coeff.rm*rm.cd[i] + mui.y.int
#    mui.cd[i] = max(mui.raw.cd[i], 1e-6)

    # Data w/ error
    rm.cd[i] ~ dnorm(radius.cd[i], 1 / 0.25e-6 ^ 2)
    po4.cd[i] ~ dnorm(po4.cd.data[i], 1 / 0.05 ^ 2)
    
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
    mui.gig[i] = coeff.po4*po4.gig[i] + coeff.rm*rm.gig[i] + mui.y.int
    #mui.gig[i] = max(mui.raw.gig[i], 1*10^-6)
    
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
  
}
