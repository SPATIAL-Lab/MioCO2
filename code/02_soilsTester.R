library(openxlsx)
source("code/helpers.R")
library(R2jags)

# Get and parse proxy files ----
## All files
pf = list.files("data/proxyData/", full.names = TRUE)

## Soils
pf.s = pf[grep("paleosol", pf)]

## Read file
pd.s = read.xlsx(pf.s[1], 1, startRow = 4)

# Current data structure
## datum
## soil

# Long-term solution
## datum
## soil
## age
## site

