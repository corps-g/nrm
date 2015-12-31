# This example computes cycle lengths for 
# various thicknesses of FeCrAl and SiC.

import numpy as np
import cPickle as pickle 
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
    
    T_F = 900*np.ones(p['number_batches'])    # batch fuel temperatures (K)
    T_C = 580*np.ones(p['number_batches'])    # batch moderator temperatures (K)
    
    num_thick = 20
    thick = np.linspace(0.0, 500, num_thick)
    B_c_FeCrAl = np.zeros((2, num_thick))
    B_c_SiC = np.zeros((2, num_thick))
    
    solver = NRM(p, rho=rho, m2=m2, k_cladding=k_cladding)
    
    te = time.time()
    for i in range(num_thick):    
    
        print 'thickness = %f micron' % thick[i]
    
        # FeCrAl
        p['power_share'] = 'equal'
        
        p['t_fecral'] = thick[i]
        p['t_sic'] = 0.0
        p['power_share'] = 'equal'
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        B_c_FeCrAl[0, i] = B[-1]/len(T_F)
        p['power_share'] = 'reactivity'
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        B_c_FeCrAl[1, i] = B[-1]/len(T_F)
    
        p['t_sic'] = thick[i]
        p['t_fecral'] = 0.0
        p['power_share'] = 'equal'
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        B_c_SiC[0, i] = B[-1]/len(T_F)
        p['power_share'] = 'reactivity'
        B, ppf, T_F, T_C = solver.solve(T_F, T_C)
        B_c_SiC[1, i] = B[-1]/len(T_F)
    te = time.time() - te
    
    print "elapsed time: %f seconds" % te
    
    pickle.dump({'thick':thick, 'B_c_FeCrAl': B_c_FeCrAl, 'B_c_SiC': B_c_SiC},\
                open('example_2.p', 'w'))

run()
