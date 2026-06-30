library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Get and parse proxy files ----
## All files
pf = list.files("data/proxyData/", full.names = TRUE)

## Alkenones
pf.a = pf[grep("phyto", pf)]

