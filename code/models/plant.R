model {
 
  # Atmosphere
  for(i in 1:length(d13Ca.obs[, 1])){
    d13Ca.obs[i, 1] ~ dnorm(d13Ca[i], 1 / d13Ca.obs[i, 2] ^ 2)
  }
  
  # Adaxial species --- how to deal with zero length???
  for(i in ind.ad){
    # Likelihood ----
    d13Cp[i, 1] ~ dnorm(d13C_m[i], 1 / d13Cp[i, 2] ^ 2)
    Dab[i, 1] ~ dnorm(D.ab[i], 1 / Dab[i, 2] ^ 2)
    GCLab[i, 1] ~ dnorm(Pl.ab[i] / s1_m[i], 1 / GCLab[i, 2] ^ 2)
    GCWab[i, 1] ~ dnorm(l.ab[i] / s2_m[gi[i]], 1 / GCWab[i, 2] ^ 2)
    Dad[i, 1] ~ dnorm(D.ad[i], 1 / Dad[i, 2] ^ 2)
    GCLad[i, 1] ~ dnorm(Pl.ad[i] / s1_m[i], 1 / GCLad[i, 2] ^ 2)
    GCWad[i, 1] ~ dnorm(l.ad[i] / s2_m[gi[i]], 1 / GCWad[i, 2] ^ 2)
    
    ## Pl to obs scaling
    s1_m[i] ~ dgamma(s1[i, 1] * s1.beta[i], s1.beta[i])
    s1.beta[i] = s1[i, 1] / s1[i, 2] ^ 2
    
    # Franks model
    d13C_m[i] = d13Ca[ai.plant[i]] - D13C[i]
    D13C[i] = a + (b - a) * ci[i] / pCO2[ai.plant[i]]
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
    # Individual ai.plant
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
    # Individual ai.plant
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
  
  # Abaxial species --- how to deal with zero length???
  for(i in ind.ab){
    # Likelihood ----
    d13Cp[i, 1] ~ dnorm(d13C_m[i], 1 / d13Cp[i, 2] ^ 2)
    Dab[i, 1] ~ dnorm(D[i], 1 / Dab[i, 2] ^ 2)
    GCLab[i, 1] ~ dnorm(Pl[i] / s1_m[i], 1 / GCLab[i, 2] ^ 2)
    GCWab[i, 1] ~ dnorm(l[i] / s2_m[gi[i]], 1 / GCWab[i, 2] ^ 2)

    ## Pl to obs scaling
    s1_m[i] ~ dgamma(s1[i, 1] * s1.beta[i], s1.beta[i])
    s1.beta[i] = s1[i, 1] / s1[i, 2] ^ 2
    
    # Franks model
    d13C_m[i] = d13Ca[ai.plant[i]] - D13C[i]
    D13C[i] = a + (b - a) * ci[i] / pCO2[ai.plant[i]]
    ci[i] = pCO2[ai.plant[i]] - A[i] / gcop[i] - 1 / meso.scale[gi[i]]
    
    ## Based on data I've seen A should have noise added
    A[i] ~ dgamma(A.1[i] ^ 2 / 0.25, A.1[i] / 0.25)
    A.1[i] = (-q.b[i] - sqrt(q.b[i] ^ 2 - 4 * q.a[i] * q.c[i])) / (2 * q.a[i])
    
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
    # Individual ai.plant
    gcop[i] = (1 / gcop.g[i] + 1 / gb_m[gi[i]]) ^ -1
    gcop.g[i] = gcmax[i] * gc.scale[gi[i]]
    gcmax[i] = (d.v * D[i] * amax[i]) / (l[i] + ((pi / 2) * sqrt(amax[i] / pi))) / 1.6
    amax[i] = SA[i] * amax.scale[gi[i]] / D[i] 
    
    # Stomatal geometry calculations ----
    # Individual ai.plant
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

  # Locality priors ----
  for(i in 1:length(d13Ca.obs[, 1])){
    pCO2[i] = pCO2.s[i] * 1e3
    pCO2.s[i] ~ dunif(0.1, 2)
    d13Ca[i] ~ dunif(-8, -3)
  }
  
  # Constants ----
  pi = 3.14159265
  a = 4.4
  d.v = 0.000940096
  b = 30
  gamma = 40
}
