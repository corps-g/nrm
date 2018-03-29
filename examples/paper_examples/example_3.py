# This example computes temperatures and power-peaking factors for 
# various thicknesses of FeCrAl and SiC.

import scipy as sp
import pickle 
import time

from nrm import NRM
from atf_models import rho, m2, k_cladding

def run() :

    # parameter dictionary
    p = {}                 
    p['number_batches'] = 3       
    p['leakage_penalty'] = 0.04
    p['assembly_width'] = 21.5036
    p['assembly_power'] = 3.4/193
    p['active_height'] = 366.0
    p['fuel_radius'] = 0.4096
    p['cladding_inner_radius'] = 0.4180
    p['cladding_outer_radius'] = 0.4750
    p['number_pins'] = 264
    p['power_share'] = 'reactivity'
        
    T_F = 900*sp.ones(p['number_batches'])    # batch fuel temperatures (K)
    T_C = 580*sp.ones(p['number_batches'])    # batch moderator temperatures (K)
    
    num_thick = 20
    #thick = sp.linspace(0.0, 500, num_thick)
    thick = sp.logspace(-1, sp.log10(5*10**2), num_thick)
    T_F_FeCrAl = sp.zeros((3, num_thick))
    T_C_FeCrAl = sp.zeros((3, num_thick))
    T_F_SiC = sp.zeros((3, num_thick))
    T_C_SiC = sp.zeros((3, num_thick))
    PPF_FeCrAl = sp.zeros((3, num_thick))
    PPF_SiC = sp.zeros((3, num_thick))
    
    solver = NRM(p, rho=rho, m2=m2, k_cladding=k_cladding)

    
    for i in range(num_thick) :    
        p['t_fecral'] = thick[i]
        p['t_sic'] = 0.0
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        T_F_FeCrAl[:, i] = T_F[:]
        T_C_FeCrAl[:, i] = T_C[:]
        PPF_FeCrAl[:, i] = ppf[:]
    
        p['t_fecral'] = 0.0
        p['t_sic'] = thick[i]
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        T_F_SiC[:, i] = T_F[:]
        T_C_SiC[:, i] = T_C[:]
        PPF_SiC[:, i] = ppf[:]
    
    pickle.dump({'thick': thick,
                 'T_F_FeCrAl': T_F_FeCrAl,
                 'T_C_FeCrAl': T_C_FeCrAl,
                 'PPF_FeCrAl': PPF_FeCrAl,
                 'T_F_SiC': T_F_SiC,
                 'T_C_SiC': T_C_SiC,
                 'PPF_SiC': PPF_SiC}, open('example_3.p', 'wb'))

run()
