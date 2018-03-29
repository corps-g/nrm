"""
Defines data models for neutronics and thermal-hydraulics analysis of 
a standard WH-style, 17x17 assembly with no burnable absorber.  Note,
all dimensions come from 
   K. Smith, et al., “Benchmarks for Quantifying Fuel Reactivity Depletion 
   Uncertainty,” Electric Power Research Institute (EPRI), Palo Alto, 
   CA, Technical Report Number 1022909, (2011).
"""

from nrm.default_models import k_fuel, h_gap
import numpy as np
import os

inp_template = \
"""TTL * STD
PDE 104.5 'KWL'
PRE 155
TFU={0:.2f} TMO={1:.2f} BOR={2:.2f}
PWR 17 1.2598 21.5036
FUE 1,10.34/{3:.2f}
PIN 1 0.4096 0.4180 0.4750 /'1' 'AIR','CAN'
PIN 2 0.5610 0.6120 / 'MOD' 'BOX'
PIN 3 0.5610 0.6120 / 'MOD' 'BOX'
LPI
2
1 1
1 1 1
3 1 1 3
1 1 1 1 1
1 1 1 1 1 3
3 1 1 3 1 1 1
1 1 1 1 1 1 1 1
1 1 1 1 1 1 1 1 1
DEP -60
STA
COE ,, -60 
TFU {4:.2f}
TMO {5:.2f}
BOR {6:.2f}
STA
END
"""

class CASMO4:
    
    def __init__(self, T_F0, T_C0, C_B0, e = 4.0, degree=1, run=True):
        """ Initialize the CASMO4 model generator.
        
        Given historical fuel temperature, coolant temperature, and 
        boron concentration, and given the fuel enrichment, produces
        and runs CASMO4, reads the log file, and first-order Taylor
        expansion to evaluate rho and M2 as functions of T_F, T_C, and C_B.
        
        Inputs:
            T_F0 : float
                Historical fuel temperature
            T_C0 : float
                Historical coolant temperature   
            C_B0 : float
                Historical boron concentration
               e : float
                Fuel U-235 enrichment (mass %)
          degree : int
                Degree of polynomial used to represent B dependence of
                rho, M2, and their derivatives with respect to T_F, T_C, and
                B_C.
             run : bool
                If True, the CASMO4 calculation will be performed.  This is
                the default.  By setting it to false, the existing CASMO4
                output will be used if it exists.  
        """        
        
        assert 0 < e <= 5
        
        delta_T_F = 50.0
        delta_T_C = 20.0
        delta_C_B = 50.0
        
        T_F = T_F0 + delta_T_F
        T_C = T_C0 + delta_T_C
        C_B = C_B0 + delta_C_B
        
        # Make input
        with open('casmo.inp', 'w') as f:
            
            f.write(inp_template.format(T_F0, T_C0, C_B0, e, T_F, T_C, C_B))
        
        # Run CASMO4, assuming in your path, and overwrite existing files
        try:
            if run or not os.path.isfile('casmo.log'):
                os.system('casmo4 -k casmo.inp')
        except:
            print("running CASMO4 was unsuccessful!")
            
        # Read the output
        with open('casmo.log', 'r') as f:      
            lines = f.readlines()      
            starts = {'base': 0, 'T_F': 0, 'T_C': 0, 'C_B': 0}
            for i in range(len(lines)):
                if '** C A S M O - 4E SUMMARY **' in lines[i]:
                    starts['base'] = i + 6
                elif 'COE' in lines[i] and 'TFU' in lines[i]:
                    starts['T_F'] = i + 1
                elif 'COE' in lines[i] and 'TMO' in lines[i]:
                    starts['T_C'] = i + 1               
                elif 'COE' in lines[i] and 'BOR' in lines[i]:
                    starts['C_B'] = i + 1    
            data = {}
            for key in starts:
                data[key] = {}
                data[key]['B'] = []
                data[key]['rho'] = []
                data[key]['M2'] = []
                i = starts[key]
                while True:
                    words = lines[i][40:].split()
                    if len(words) != 8:
                        break
                    data[key]['B'].append(float(words[0]))
                    kinf = float(words[1])
                    data[key]['rho'].append((kinf-1)/kinf)
                    data[key]['M2'].append(float(words[3]))
                    i += 1
            self.data = data
            
        assert len(data['base']['B']) == len(data['T_F']['B'])
        assert len(data['base']['B']) == len(data['T_C']['B'])
        assert len(data['base']['B']) == len(data['C_B']['B'])

        # Process the data
        # assume rho(BU, T_F, T_C, C_B) = f(BU) + (T_F-T_F0)*drho/dT_F + ...
        
        B = np.array(data['base']['B'])
        rho = np.array(data['base']['rho'])
        M2 = np.array(data['base']['M2'])
        drho_dT_F = (np.array(np.array(data['T_F']['rho'])) - rho)/delta_T_F
        drho_dT_C = (np.array(np.array(data['T_C']['rho'])) - rho)/delta_T_C
        drho_dC_B = (np.array(np.array(data['C_B']['rho'])) - rho)/delta_C_B
        dM2_dT_F = (np.array(np.array(data['T_F']['M2'])) - M2)/delta_T_F
        dM2_dT_C = (np.array(np.array(data['T_C']['M2'])) - M2)/delta_T_C
        dM2_dC_B = (np.array(np.array(data['C_B']['M2'])) - M2)/delta_C_B
        
        self.coefs = {}
        
        # Fit rho and derivatives to polynomials in BU
        self.coefs['rho'] = np.polyfit(B, rho, degree)
        self.coefs['drho_dT_F'] = np.polyfit(B, drho_dT_F, degree)
        self.coefs['drho_dT_C'] = np.polyfit(B, drho_dT_C, degree)
        self.coefs['drho_dC_B'] = np.polyfit(B, drho_dC_B, degree)
  
        # Fit M2 and derivatives to polynomials in BU
        self.coefs['M2'] = np.polyfit(B, M2, degree)
        self.coefs['dM2_dT_F'] = np.polyfit(B, dM2_dT_F, degree)
        self.coefs['dM2_dT_C'] = np.polyfit(B, dM2_dT_C, degree)
        self.coefs['dM2_dC_B'] = np.polyfit(B, dM2_dC_B, degree)
        
        self.T_F0 = T_F0
        self.T_C0 = T_C0
        self.C_B0 = C_B0
        
    def rho(self, p, B, T_F, T_C, C_B) :   
        """ Return reactivity.
        
            Inputs:
                   p : dict
                    Parameters for NRM.  Not used here.
                   B : float or array-like
                    Burnup
                 T_F : float
                    Instantaneous fuel temperature   
                 T_C : float
                    Instantaneous coolant temperature   
                 C_B : float
                    Instantaneous boron concentration
        """
        dT_F = T_F - self.T_F0 
        dT_C = T_C - self.T_C0
        dC_B = C_B - self.C_B0
        
        rho_0 = np.polyval(self.coefs['rho'], B)
        drho_dT_F = np.polyval(self.coefs['drho_dT_F'], B)
        drho_dT_C = np.polyval(self.coefs['drho_dT_C'], B)
        drho_dC_B = np.polyval(self.coefs['drho_dC_B'], B)
        return rho_0 + dT_F*drho_dT_F + dT_C*drho_dT_C + dC_B*drho_dC_B

    def m2(self, p, B, T_F, T_C, C_B) :   
        """ Computes the migration area for given conditions.
        
            Inputs:
                   p : dict
                    Parameters for NRM.  Not used here.
                   B : float or array-like
                    Burnup
                 T_F : float
                    Instantaneous fuel temperature   
                 T_C : float
                    Instantaneous coolant temperature   
                 C_B : float
                    Instantaneous boron concentration
                    
            Returns:
                 migration area : float or array-like
        """
        dT_F = T_F - self.T_F0 
        dT_C = T_C - self.T_C0
        dC_B = C_B - self.C_B0
        
        M2_0 = np.polyval(self.coefs['M2'], B)
        dM2_dT_F = np.polyval(self.coefs['dM2_dT_F'], B)
        dM2_dT_C = np.polyval(self.coefs['dM2_dT_C'], B)
        dM2_dC_B = np.polyval(self.coefs['dM2_dC_B'], B)
        return M2_0 + dT_F*dM2_dT_F + dT_C*dM2_dT_C + dC_B*dM2_dC_B    
        
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    c = CASMO4(900.0, 580.0, 900.0, 4.0, degree=1, run=False)
 
    # fuel temperature coefficient
    B = np.linspace(0, 60)
    rho_0 = c.rho(None, B, 900, 580, 900)
    rho_1 = c.rho(None, B, 1000, 580, 900)
    plt.plot(B, (rho_1-rho_0)/100*1e5)
    plt.xlabel('burnup (MWd/kg)')
    plt.ylabel(r'$\alpha_F$ (pcm/K)')
    