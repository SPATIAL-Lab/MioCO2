

################################################################################
# Discrete foram PSM inversion driver
################################################################################

library(rjags)
library(R2jags)
library(openxlsx)

################################################################################
# paths
################################################################################

model_file <- "discrete_foram_psm.R"

proxy_file <- file.path("boron_Intermediate_combined.xlsx")
template_sheet <- "data4PSM"

################################################################################
# proxy data
################################################################################

prox_raw <- read.xlsx(proxy_file, sheet = template_sheet, startRow = 5, colNames = FALSE)

proxy_df <- data.frame(
  doi = as.character(prox_raw[[2]]), # B
  lat = as.numeric(prox_raw[[6]]), # F
  lon = as.numeric(prox_raw[[7]]), # G
  age = as.numeric(prox_raw[[8]]), # H
  species = trimws(as.character(prox_raw[[13]])), # M
  method = toupper(trimws(as.character(prox_raw[[14]]))), # P
  d11B = as.numeric(prox_raw[[15]]), # Q
  d11B.2sd = as.numeric(prox_raw[[16]]), # R
  mgca = as.numeric(prox_raw[[17]]), # S
  mgca.2sd = as.numeric(prox_raw[[18]]), # T
  stringsAsFactors = FALSE
)

proxy_df$method <- toupper(trimws(proxy_df$method))
proxy_df$method <- gsub("[^A-Z]", "", proxy_df$method)
proxy_df$method[grepl("MCICPMS", proxy_df$method)] <- "MCICPMS"
proxy_df$method[grepl("TIMS", proxy_df$method)] <- "TIMS"

clean.d11B <- proxy_df[complete.cases(proxy_df[, c("d11B", "d11B.2sd", "species", "method")]), ]
clean.mgca <- proxy_df[complete.cases(proxy_df[, c("mgca", "mgca.2sd")]), ]

d11B.data <- clean.d11B$d11B
d11B.data.sigma <- clean.d11B$d11B.2sd / 2

mgca.data <- clean.mgca$mgca
mgca.data.sigma <- clean.mgca$mgca.2sd / 2

# The current combined boron intermediate sheet does not provide d18O
# This placeholder keeps the existing d18O likelihood from constraining the model
d18O.data <- c(NA_real_)
d18O.data.sigma <- c(1)

n_d11B <- length(d11B.data)
n_mgca <- length(mgca.data)
n_d18O <- length(d18O.data)

################################################################################
# boron vital-effect calibration by species and method
################################################################################

d11B.species <- clean.d11B$species
d11B.method <- clean.d11B$method

species_group <- rep(NA_character_, length(d11B.species))

species_group[d11B.species %in% c(
  "Trilobatus trilobus",
  "Globigerinoides trilobus",
  "Globigerinoides sacculifer",
  "T. trilobus"
)] <- "T_trilobus"

species_group[d11B.species %in% c(
  "Globoturborotalita pseudopraebulloides",
  "G. praebulloides"
)] <- "G_pseudopraebulloides"

species_group[d11B.species %in% c(
  "Ciperoella angulisuturalis"
)] <- "C_angulisuturalis"

species_group[d11B.species %in% c(
  "Globigerinoides ruber white (sensu stricto)",
  "Globigerinoides ruber white"
)] <- "G_ruber_white"

cal <- data.frame(
  species_group = c(
    "T_trilobus",
    "T_trilobus",
    "G_pseudopraebulloides",
    "G_pseudopraebulloides",
    "C_angulisuturalis",
    "C_angulisuturalis",
    "G_ruber_white",
    "G_ruber_white"
  ),
  method = c(
    "TIMS",
    "MCICPMS",
    "TIMS",
    "MCICPMS",
    "TIMS",
    "MCICPMS",
    "TIMS",
    "MCICPMS"
  ),
  m.vital.m = c(
    0.73,
    0.833,
    0.62,
    0.62,
    0.62,
    0.62,
    0.60,
    0.60
  ),
  m.vital.sigma = c(
    0.04,
    0.05,
    0.055,
    0.055,
    0.055,
    0.055,
    0.10,
    0.10
  ),
  c.vital.m = c(
    6.42,
    2.69,
    9.52,
    9.52,
    9.52,
    9.52,
    8.87,
    8.87
  ),
  c.vital.sigma = c(
    0.82,
    1.00,
    1.01,
    1.01,
    1.01,
    1.01,
    1.50,
    1.50
  ),
  c.correction = c(
    -2.8,
    0,
    -4.7,
    -4.7,
    -3.76,
    -3.76,
    0,
    0
  ),
  stringsAsFactors = FALSE
)

sample_cal_key <- paste(species_group, d11B.method, sep = "__")
cal$key <- paste(cal$species_group, cal$method, sep = "__")

keep_cal <- !is.na(species_group) & sample_cal_key %in% cal$key

d11B.data <- d11B.data[keep_cal]
d11B.data.sigma <- d11B.data.sigma[keep_cal]
d11B.species <- d11B.species[keep_cal]
d11B.method <- d11B.method[keep_cal]
species_group <- species_group[keep_cal]
sample_cal_key <- sample_cal_key[keep_cal]

cal_keys_used <- unique(sample_cal_key)
cal_used <- cal[match(cal_keys_used, cal$key), , drop = FALSE]

cal_i <- match(sample_cal_key, cal_used$key)
n_cal <- nrow(cal_used)
n_d11B <- length(d11B.data)

m.vital.m <- cal_used$m.vital.m
m.vital.sigma <- cal_used$m.vital.sigma
c.vital.m <- cal_used$c.vital.m
c.vital.sigma <- cal_used$c.vital.sigma
c.correction <- cal_used$c.correction

boron_cal_samples <- data.frame(
  species = d11B.species,
  method = d11B.method,
  species_group = species_group,
  cal_i = cal_i,
  stringsAsFactors = FALSE
)

print(boron_cal_samples)
print(cal_used)

################################################################################
# prior distributions
################################################################################

# Environmental state
tempC.m <- 25
tempC.sigma <- 3
tempC.min <- 0
tempC.max <- 45

pH.m <- 7.8
pH.sigma <- 0.05
pH.min <- 7.0
pH.max <- 8.5

sal.m <- 35
sal.sigma <- 1
sal.min <- 25
sal.max <- 45

dic.m <- 2.0
dic.sigma <- 0.3
dic.min <- 0.5
dic.max <- 5.0

press.m <- 5
press.sigma <- 1
press.min <- 0
press.max <- 20

# Seawater chemistry approx Miocene values from Zeebe & Tyrrell (2019)
xca.m <- 12
xca.sigma <- 0.5
xca.min <- 1
xca.max <- 40

xmg.m <- 45
xmg.sigma <- 0.5
xmg.min <- 1
xmg.max <- 100

xso4.m <- 25
xso4.sigma <- 0.5
xso4.min <- 1
xso4.max <- 40

d11Bsw.m <- 39.0 # approx Miocene value from Ring et al. (2025)
d11Bsw.sigma <- 0.5
d11Bsw.min <- d11Bsw.m - 2 * d11Bsw.sigma
d11Bsw.max <- d11Bsw.m + 2 * d11Bsw.sigma

d18Osw.m <- -0.5 # Miocene approx
d18Osw.sigma <- 0.2
d18Osw.min <- d18Osw.m - 2 * d18Osw.sigma
d18Osw.max <- d18Osw.m + 2 * d18Osw.sigma

# Scalar proxy-system parameters
alpha.m <- 1.0272
alpha.sigma <- 0.0003
alpha.min <- alpha.m - 2 * alpha.sigma
alpha.max <- alpha.m + 2 * alpha.sigma

sw.sens.m <- 0.558
sw.sens.sigma <- 0.03
sw.sens.min <- sw.sens.m - 2 * sw.sens.sigma
sw.sens.max <- sw.sens.m + 2 * sw.sens.sigma

# Mg/Ca calibration
A.m <- 0.061
A.sigma <- 0.005
A.min <- A.m - 2 * A.sigma
A.max <- A.m + 2 * A.sigma

Hp.m <- 0.74
Hp.sigma <- 0.05
Hp.min <- Hp.m - 2 * Hp.sigma
Hp.max <- Hp.m + 2 * Hp.sigma

salcorrco.m <- 4.2
salcorrco.sigma <- 0.4
salcorrco.min <- salcorrco.m - 2 * salcorrco.sigma
salcorrco.max <- salcorrco.m + 2 * salcorrco.sigma

# Mr/Ca pH correction (sensitvity per 1.0 pH unit)
pHcorrco.m <- 0.003
pHcorrco.sigma <- 0.001
pHcorrco.min <- 0
pHcorrco.max <- 0.02

# d18O pH correction
d18O_pHcorr.m <- -0.089
d18O_pHcorr.sigma <- 0.02
d18O_pHcorr.min <- d18O_pHcorr.m - 2 * d18O_pHcorr.sigma
d18O_pHcorr.max <- 0

# Diagenetic overprint
indexop.m <- 10
indexop.sigma <- 5
indexop.min <- 0
indexop.max <- 100

seccal.d18O.m <- 0.85
seccal.d18O.sigma <- 0.1

################################################################################
# data passed to JAGS
################################################################################

data2pass <- list(
  d11B.data = d11B.data,
  d11B.data.sigma = d11B.data.sigma,
  mgca.data = mgca.data,
  mgca.data.sigma = mgca.data.sigma,
  d18O.data = d18O.data,
  d18O.data.sigma = d18O.data.sigma,
  n_d11B = n_d11B,
  n_mgca = n_mgca,
  n_d18O = n_d18O,

  tempC.m = tempC.m,
  tempC.sigma = tempC.sigma,
  tempC.min = tempC.min,
  tempC.max = tempC.max,
  pH.m = pH.m,
  pH.sigma = pH.sigma,
  pH.min = pH.min,
  pH.max = pH.max,
  sal.m = sal.m,
  sal.sigma = sal.sigma,
  sal.min = sal.min,
  sal.max = sal.max,
  dic.m = dic.m,
  dic.sigma = dic.sigma,
  dic.min = dic.min,
  dic.max = dic.max,
  press.m = press.m,
  press.sigma = press.sigma,
  press.min = press.min,
  press.max = press.max,

  xca.m = xca.m,
  xca.sigma = xca.sigma,
  xca.min = xca.min,
  xca.max = xca.max,
  xmg.m = xmg.m,
  xmg.sigma = xmg.sigma,
  xmg.min = xmg.min,
  xmg.max = xmg.max,
  xso4.m = xso4.m,
  xso4.sigma = xso4.sigma,
  xso4.min = xso4.min,
  xso4.max = xso4.max,
  d11Bsw.m = d11Bsw.m,
  d11Bsw.sigma = d11Bsw.sigma,
  d11Bsw.min = d11Bsw.min,
  d11Bsw.max = d11Bsw.max,
  d18Osw.m = d18Osw.m,
  d18Osw.sigma = d18Osw.sigma,
  d18Osw.min = d18Osw.min,
  d18Osw.max = d18Osw.max,

  alpha.m = alpha.m,
  alpha.sigma = alpha.sigma,
  alpha.min = alpha.min,
  alpha.max = alpha.max,
  sw.sens.m = sw.sens.m,
  sw.sens.sigma = sw.sens.sigma,
  sw.sens.min = sw.sens.min,
  sw.sens.max = sw.sens.max,

  m.vital.m = m.vital.m,
  m.vital.sigma = m.vital.sigma,
  c.vital.m = c.vital.m,
  c.vital.sigma = c.vital.sigma,
  c.correction = c.correction,
  cal_i = cal_i,
  n_cal = n_cal,
  
  A.m = A.m,
  A.sigma = A.sigma,
  A.min = A.min,
  A.max = A.max,
  Hp.m = Hp.m,
  Hp.sigma = Hp.sigma,
  Hp.min = Hp.min,
  Hp.max = Hp.max,
  salcorrco.m = salcorrco.m,
  salcorrco.sigma = salcorrco.sigma,
  salcorrco.min = salcorrco.min,
  salcorrco.max = salcorrco.max,
  pHcorrco.m = pHcorrco.m,
  pHcorrco.sigma = pHcorrco.sigma,
  pHcorrco.min = pHcorrco.min,
  pHcorrco.max = pHcorrco.max,

  d18O_pHcorr.m = d18O_pHcorr.m,
  d18O_pHcorr.sigma = d18O_pHcorr.sigma,
  d18O_pHcorr.min = d18O_pHcorr.min,
  d18O_pHcorr.max = d18O_pHcorr.max,

  indexop.m = indexop.m,
  indexop.sigma = indexop.sigma,
  indexop.min = indexop.min,
  indexop.max = indexop.max,
  seccal.d18O.m = seccal.d18O.m,
  seccal.d18O.sigma = seccal.d18O.sigma
)

################################################################################
# parameters to monitor
################################################################################

parms2save <- c(
  "tempC",
  "pH",
  "sal",
  "dic",
  "press",
  "xca",
  "xmg",
  "xso4",
  "d11Bsw",
  "d18Osw",
  "d18Osw.sc",
  "pco2",
  "d11Bcarb",
  "d18Ocarb",
  "mgcacarb",
  "alpha",
  "sw.sens",
  "m.vital",
  "c.vital",
  "c.final",
  "A",
  "Hp",
  "salcorrco",
  "pHcorrco",
  "d18O_pHcorr",
  "indexop",
  "seccal.d18O"
)

################################################################################
# run model
################################################################################

set.seed(1)

jout <- jags.parallel(
  model.file = model_file,
  parameters.to.save = parms2save,
  data = data2pass,
  inits = NULL,
  n.chains = 3,
  n.iter = 5000,
  n.burnin = 2500,
  n.thin = 2
)

print(jout)

################################################################################
# posterior object
################################################################################

sims.list <- jout$BUGSoutput$sims.list



