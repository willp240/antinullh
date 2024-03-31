import math

ELECTRON_MASS = 511e3 # electron volts
TE_130_MASS   = 127.60 * 1.67e-27 # ~130u
TE_130_ABUNDANCE  = 0.3408

SNOT_MASS = 782e3
SNOT_LOADING = 0.005

# bb with light neutrino exchange
def mass_to_hl(bb_mass, phase_space_fac, nuc_mat_el):
    return (ELECTRON_MASS/bb_mass) ** 2 /phase_space_fac/(nuc_mat_el ** 2)
    
def hl_to_mass(hl, phase_space_fac, nuc_mat_el):
    return 1000 * ELECTRON_MASS/math.sqrt(phase_space_fac * hl) / nuc_mat_el

def count_nuclei(detector_mass, loading, abundance, iso_mass):
    return detector_mass * loading * abundance/iso_mass

def counts_to_hl(counts, live_time, sig_eff,
                 detector_mass, loading, abundance, iso_mass):
    n_nuc = count_nuclei(detector_mass, loading, abundance, iso_mass)
    return math.log(2)/counts * n_nuc * sig_eff * live_time

def expected_bb_events(bb_mass, phase_space_fac, nuc_mat_el, detector_mass,
                       loading, abundance, iso_mass):
    return count_nuclei(detector_mass, loading, abundance, iso_mass) \
           * math.log(2)/mass_to_hl(bb_mass, phase_space_fac, nuc_mat_el)




    


    
