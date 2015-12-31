# Defines data models for neutronics and thermal-hydraulics analysis of 
# fuel with protective layers of FeCrAl or SiC on a Zircaloy substrate
from nrm.default_models import k_fuel, h_gap
import numpy as np
atan2 = np.arctan2

def rho(p, B, T_F, T_C, C_B) :
    t1 = p.get('t_fecral', 0.0)
    t2 = p.get('t_sic', 0.0)
    a = rho_n(t1, t2, B, T_F, T_C, C_B)
    b = drho_TF(t1, t2, B, T_F, T_C, C_B) * (T_F - 900)
    c = drho_TC(t1, t2, B, T_F, T_C, C_B) * (T_C - 580)
    d = drho_CB(t1, t2, B, T_F, T_C, C_B) * (C_B - 900)
    return a + b + c + d
    
def rho_n(t1, t2, B, T_F, T_C, C_B):
    return 0.189191355597358 + 1.03127898873563e-5*t2 + \
        0.0182701657558425*atan2(0.047863097316982, B) + \
        1.87120440460553e-8*t1**2 - 0.000101201929816802*t1 - \
        0.00692731147074213*B
        
def drho_TF (t1, t2, B, T_F, T_C, C_B):
    return 3.67779316136406e-10*t2 + 1.88453108088694e-11*t1*B + \
    2.99036198934836e-6*atan2(1.03432017102404e-9, 4.32782583441201e-12*B**2) \
    - 2.3479749123139e-5 - 1.55088064463373e-9*t1 - \
    2.00083642922382e-5*atan2(4.32782583441201e-12*B, 3.53879022031715e-10) - \
    4.57451854075048e-7*atan2(2.87730200198639e-11, 4.32782583441201e-12*B**2)
    
def drho_TC(t1, t2, B, T_F, T_C, C_B):
    return 5.17451736346039e-8*t2 + 2.58256031392406e-8*t1 + \
    1.73019057111078e-5*atan2(B, 0.0604702101669828) - 9.50595239549672e-5 - \
    0.000252658620269585*atan2(2.20870099307213 + B + 0.00225477340696242*t2, \
    25.1901924655218)
    
def drho_CB (t1, t2, B, T_F, T_C, C_B):
    return 4.49383643308143e-10*t2 + 2.16366537735464e-10*t1*B + \
    3.88963558342896e-6*atan2(3.93055887273536e-9*B, 1.56134738600496e-8) - \
    6.344933749376e-5 - 4.19488057273303e-7*B - 3.85361926498037e-9*B**2
    
def m2(p, B, T_F, T_C, C_B) :
    t1 = p.get('t_fecral', 0.0)
    t2 = p.get('t_sic', 0.0)
    a = m2_n(t1, t2, B, T_F, T_C, C_B)
    b = dm2_TF(t1, t2, B, T_F, T_C, C_B) * (T_F - 900)
    c = dm2_TC(t1, t2, B, T_F, T_C, C_B) * (T_C - 580)
    d = dm2_CB(t1, t2, B, T_F, T_C, C_B) * (C_B - 900)
    return a + b + c + d
    
def m2_n(t1, t2, B, T_F, T_C, C_B):
    return 19.0993846779234*atan2(t1, -3792.00816694545 - t1) - \
    0.000854117522725537*t2 - 0.0167254413395886*B*atan2(t1,\
    -1394.83489391658) - 1.70177760223332*atan2(B +\
    0.139743374206404/(B - 0.0167254413395886) + atan2(t2, \
    0.000854117522725537 - t1 - t2), 25.4410176338263)   
    
def dm2_TF(t1, t2, B, T_F, T_C, C_B):
    return 1.14195849992072e-7*t1 + 2.04954693314212e-8*t2 + \
    6.18196432362702e-8*B**2 + 5.37036840825969e-5*atan2(B, \
    atan2(6.18196432362702e-8*B**2, 2.04954693314212e-8)) - \
    0.000910107267275084 - 1.73400777561809e-5*atan2(0.104278409419016*t1, \
    -10.6773814299913 - B)

def dm2_TC(t1, t2, B, T_F, T_C, C_B):
    return 0.220865758013557 + 6.6568349721182e-8*t1*B + 0.00714439416183893*\
    atan2(16.5298645830687, B) + 0.00104797566471284*\
    atan2(7.25098475175339e-9, B) + 7.25098475175339e-9*t1**2 - \
    4.38990606667827e-6*t2 - 2.48220306084701e-5*t1 - 0.000279275167686063*B

def dm2_CB(t1, t2, B, T_F, T_C, C_B):
    return 5.02371172919458e-6*B + 2.36342474601454e-7*t1 + \
    2.74313384828021e-8*t2 - 0.00106746273951085 - 0.000169489148631781*\
    atan2(6.61292547375856, B + atan2(B - 0.158892639202784, \
    0.00106746273951085*t1 - 0.0711784200386005 - 9.22151673718882e-5*B*t1**2))
    
def k_cladding(p, B, T) :
    """ Returns the cladding thermal conductivity given the temperature T.

    Rather than explicitly model a multi-layer cladding, a simple 
    volume-averaging of the conductivities is used instead.
    """
    k_i = k_zr4(T)
    t_fecral = p.get('t_fecral', 0.0)
    t_sic = p.get('t_sic', 0.0)
    assert(t_fecral >= 0.0 and t_sic >= 0.0)
    if t_fecral > 0.0 :
        t, k_o = t_fecral, k_fecral(T)
    else :
        t, k_o = t_sic, k_sic(T)
    # zr-4 volume fraction
    r_co, r_ci = p['cladding_outer_radius'], p['cladding_inner_radius']
    v_i = ((r_co-t)**2-r_ci**2)/r_co**2;
    # coating volume fraction
    v_o = (r_co**2-(r_ci+t)**2)/r_co**2;
    # volume-weighted conductivity
    return v_i*k_i + v_o*k_o;

def k_zr4(T) :
    """ Returns thermal conductivity (W/m) for Zr-4 at temperature T (K)
    
    Reference: INL/EXT-14-32591, Revision 1
    """
    return 7.511 + 2.088e-2*T - 1.45e-5*T**2 + 7.668e-9*T**3
       
def k_fecral(T) :
    """ Returns thermal conductivity (W/m) for FeAlCr at temperature T (K).
    
    Note: Valid for 300 K <= T <= 1773 K.
 
    Reference: INL/EXT-14-32591, Revision 1
    """
    return 2.53282 + 3.2532e-2*T - 2.2e-5*T**2 + 8.5645e-9*T**3
    
def k_sic(T):
    """ Returns thermal conductivity (W/m) for SiC at temperature T (K).
    
    Note: Valid for 200 K <= T <= 1968 K.
    
    Reference: INL/EXT-14-32591, Revision 1
    """
    return 194.776655 - 0.36061185*T + 3.3084327e-4*T**2 \
          - 1.4600630e-7*T**3 + 2.4758791e-11 * T**4
          
          
