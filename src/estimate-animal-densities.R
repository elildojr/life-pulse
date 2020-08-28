#
## Script to estimate animal densities based on time-averages withouth individual recognition
## written by Elildo Carvalho Jr @ ICMBio/CENAP
## data reading and fixing script written by Jorge Ahumada @ Conservation International
#

# ---- Load libraries ---------------
library(here)

# ----Source file with Ahumada codes (for reading data)
source(here("bin", "ahumada_codes.R"))
source(here("bin", "timedens.R"))

# ----Read data-----------------------------
x <- f.readin.fix.data(here("data", "Wild_ID_RBG_2016.csv"))


# ---- Use function timedens to estimate densities based on time-averages

# User must specify x, y, z, w:
# x = dataset
# y = species binomial name (column "bin" created by f.readin.fix.data)
# z = effective detection distance. See table 2 in Rowcliffe et al. 2011 doi: 10.1111/j.2041-210X.2011.00094.x for a reference
# w = camera detection angle (e.g., 52 degrees in Bushnell trophy cam), see also Rowcliffe et al. 2011
# Note: The function assumes 60 cameras operating simultaneously from the first to the last day of survey
# to correct for number of cameras change the formula for area.sampled in function


# Test
timedens(x, "Dasyprocta prymnolopha", 3, 50)
timedens(x, "Mitu tuberosum", 3, 50)
timedens(x, "Psophia unknown", 3, 50)
