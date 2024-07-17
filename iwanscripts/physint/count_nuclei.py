from rate_conversions import count_nuclei, hl_to_mass
import math

ELECTRON_MASS = 511e3 # electron volts
TE_130_MASS   = 127.60 * 1.67e-27 # ~130u
TE_130_ABUNDANCE  = 0.3408

SNOT_MASS = 782e3
SNOT_LOADING = 0.005

print count_nuclei(SNOT_MASS, SNOT_LOADING, TE_130_ABUNDANCE, TE_130_MASS)
