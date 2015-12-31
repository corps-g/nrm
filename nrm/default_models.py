# Defines the default data models for neutronics and thermal-hydraulics

import scipy as sp
atan2 = sp.arctan2 
sqrt = sp.sqrt 

def rho(p, B, T_F, T_C, C_B) :
    #return 0.2 - B*0.2/30.
    return 0.265701782943972 + \
           0.0173604082545277*atan2(0.0173604082545277, abs(B)**0.5) - \
           2.64479010825198e-5*T_F - \
           0.00319948157566795*B - 5.72926698137041e-6*B*T_C  
           
def m2(p, B, T_F, T_C, C_B) :
    #return 60.0 - B*.2
    return 1.41350127802683*T_C + 0.00112134680090497*B**2 + \
           (397008.794613071 - 0.696210782691167*C_B)/T_C - \
           1439.98856519259 - 0.000249491110525838*B*T_C
           
def k_fuel(p, B, T_K) : 
    """ Calculate fuel conductivity (W/m-K)
    
    From J.D. Hales et al. (2013) "Bison Theory" (NFIR model)
    """
    # kelvin to celsius
    T = T_K - 273.15;
    # thermal recovery function
    rf = 0.5*(1.0+sp.tanh((T-900.0)/150.0))
    # phonon contribution at start of thermal recovery [Hales eq. 8.14]
    kps = 1.0/(0.09592+0.00614*B-0.000014*B**2+(0.00025-0.00000181*B)*T)
    # phonon contribution at the end of thermal recovery [Hales eq. 8.15]
    kpend = 1.0/(0.09592+0.0026*B+(0.00025-0.00000027*B)*T)
    # unirradiated material at 95% th. density [Hales eq. 8.17]
    kel = 0.0132*sp.exp((0.00188)*T)
    k = (1.0-rf)*kps + rf*kpend + kel
    return k
    
def k_cladding(p, B, T) :
    """ Returns thermal conductivity (W/m) for Zr-4 at temperature T (K)
    
    Reference: INL/EXT-14-32591, Revision 1
    """
    return 7.511 + 2.088e-2*T - 1.45e-5*T**2 + 7.668e-9*T**3
    
def h_gap(p, B, T) :
    """ Returns the gap conductance (W/m^2-K) as a function of burnup.  

    This very simple model is based on Figure 8.24 of Todreas and Kazimi.  
    A more realistic model would require direct coupling to a fuel performance
    code or a correlation built from output of such a code. 
    """
    if B < 30.0 :
        h = 0.6 + 0.5*B/30.0
    else :
        h = 1.1
    return 10000.0*h