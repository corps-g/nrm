# This example shows k-inf as a function of burnup for several 
# cladding configurations.
import scipy as sp
import pickle

from atf_models import rho

def kinf(p, B, T_F, T_C, C_B) :
    return 1.0 / (1.0 - rho(p, B, T_F, T_C, C_B))

# burnup, GWd/MTU
B = sp.linspace(0, 50, 1000)
# fuel temperature, K
T_F = 900.0
# coolant temperature, K
T_C = 580.0
# boron concentration, ppm
C_B = 0.0

# Zr-4
k_inf_Zr4 = kinf({}, B, T_F, T_C, C_B)
# FeCrAl - 100 micron
k_inf_FeCrAl_100 = kinf({'t_fecral': 100.0}, B, T_F, T_C, C_B)
# FeCrAl - 300 micron
k_inf_FeCrAl_300 = kinf({'t_fecral': 300.0}, B, T_F, T_C, C_B)
# SiC - 100 micron
k_inf_SiC_100 = kinf({'t_sic': 100.0}, B, T_F, T_C, C_B)
# SiC - 300 micron
k_inf_SiC_300 = kinf({'t_sic': 300.0}, B, T_F, T_C, C_B)

# Reactivity defect
rho_FeCrAl_100 = (k_inf_FeCrAl_100-k_inf_Zr4) / k_inf_FeCrAl_100*1e5
rho_FeCrAl_300 = (k_inf_FeCrAl_300-k_inf_Zr4) / k_inf_FeCrAl_300*1e5
rho_SiC_100 = (k_inf_SiC_100-k_inf_Zr4) / k_inf_FeCrAl_100*1e5
rho_SiC_300 = (k_inf_SiC_300-k_inf_Zr4) / k_inf_FeCrAl_300*1e5

# save for plotting
pickle.dump({'B': B, 
             'k_inf_Zr4': k_inf_Zr4, 
             'k_inf_FeCrAl_100': k_inf_FeCrAl_100, 
             'k_inf_FeCrAl_300': k_inf_FeCrAl_300,
             'k_inf_SiC_100': k_inf_SiC_100, 
             'k_inf_SiC_300': k_inf_SiC_300, 
             'rho_FeCrAl_100': rho_FeCrAl_100, 
             'rho_FeCrAl_300': rho_FeCrAl_300, 
             'rho_SiC_100': rho_SiC_100,              
             'rho_SiC_300': rho_SiC_300}, open('example_1.p', 'wb'))
