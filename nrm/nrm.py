import scipy as sp
from scipy.optimize import root, fsolve
from . import default_models

class NRM(object) :
    """ 
    Implementation of a coupled-physics, nonlinear reactivity method.
    
    Attributes
    ----------
    p : dict
        User input parameters
    rho : function
        Reactivity model     
    m2 : function
        Migration area model
        
    """
    
    ###########################################################################
    # PUBLIC INTERFACE                                                        #
    ###########################################################################
    
    def __init__(self, p, **kwargs) :
        """ Initialize the NRM solver.
        
        FINISH ME.
        """
        self.p = p
        self.rho = kwargs.get('rho', default_models.rho)
        self.m2 = kwargs.get('m2', default_models.m2)
        self.k_fuel = kwargs.get('k_fuel', default_models.k_fuel) 
        self.k_cladding = kwargs.get('k_cladding', default_models.k_cladding) 
        self.h_gap = kwargs.get('h_gap', default_models.h_gap) 

        # other parameters
        self.mdot = kwargs.get('mdot', 88.0)
        self.cp = kwargs.get('cp', 5500.0)
        self.h = kwargs.get('h', 30000.0)
        self.T_inlet = kwargs.get('T_inlet', 560.0)
        
    def solve(self, T_F, T_C) :
        """ Solve the coupled neutronics and thermal-hydraulics equations.
        
        Parameters
        ----------
        T_F : array_like
            Initial guess for batch-wise fuel temperatures

        Returns
        -------
        T_F : array_like
            Converged, batch-wise fuel temperatures
     
        FINISH ME.
        """

        T_F_old = 1.0*T_F
        T_C_old = 1.0*T_C
        for iteration in range(100):
            B, ppf = self.compute_cycle(T_F_old, T_C_old)
            T_F, T_C = self.compute_temperatures(B, ppf, T_F_old, T_C_old, 0.0)
            e = max([max(T_F-T_F_old), max(T_C-T_C_old)])
            if e < 0.01:
                break            
            T_F_old = 1.0*T_F
            T_C_old = 1.0*T_C
        return B, ppf, T_F, T_C
    
    def compute_cycle(self, T_F, T_C) :
        """ Compute batch-wise EOC burnups and cycle-averaged PPF.
        """
        
        if self.p['power_share'] == 'equal' :
            B, ppf = self._compute_cycle_equal_weighted_power(T_F, T_C)
        else :
            B, ppf = self._compute_cycle_reactivity_weighted_power(T_F, T_C)
  
        return B, ppf
    
    def compute_temperatures(self, B, ppf, T_F_old, T_C_old, C_B) :
        """ Compute batch-wise, cycle-averaged temperatures.
        """
        # assembly powers (W)
        q = self.p['assembly_power'] * ppf * 1.0e9
        # fuel, cladding, and gap radii, m
        r_f  = self.p['fuel_radius']/100.0
        r_co = self.p['cladding_outer_radius']/100.0
        r_ci = self.p['cladding_inner_radius']/100.0
        r_g  = 0.5*(r_ci+r_f)
        # active fuel height, m
        H = self.p['active_height']/100.0
        # extrapolation length, m
        #delta_H = 0.0
 
        
        T_F = 1*T_F_old
        T_C = 1*T_C_old        
        
        for i in range(len(T_F)) :
            # average and max assembly linear heat generation rate, W/m
            qp_avg = q[i]/H          
            #qp_max = (sp.pi/2.0/(H+2.0*delta_H)) / \
            #         sp.sin(sp.pi/(2.0*H+4.0*delta_H))
            # axially-averaged moderator temperature
            T_C[i] = self.T_inlet + 0.5*q[i]/self.mdot/self.cp
            # cladding thermal conductivity, W/m.K, assuming temperature is 
            # equal to the 0.12*T_F + 0.88*T_C, which is consistent with the 
            # lattice physics modeling
            T_c = 0.88*T_C[i]+0.12*T_F[i]
            k_c = self.k_cladding(self.p, B[i], T_c)
            # average gap temperature, K, and conductance, W/m^2.K
            T_g = 0.5*(T_c+T_F[i])
            h_g = self.h_gap(self.p, B[i], T_g)
            # effective non-fuel thermal conductivity, W/m.K
            k_NF = 1.0 / ((0.5/sp.pi) * (1.0/h_g/r_g + \
                                         sp.log(r_co/r_ci)/k_c + 1.0/r_co/self.h))
            # fuel thermal conductivity, W/m.K
            k_F = self.k_fuel(self.p, B[i], T_F[i])
            # average pin linear heat generation rate, W/m
            qp_pin = qp_avg / self.p['number_pins']
            # axially- and radially-averaged fuel temperature
            T_F[i] = qp_pin*(0.125/sp.pi/k_F + 1.0/k_NF) + T_C[i]
        return T_F, T_C
    
    ###########################################################################
    # PRIVATE IMPLEMENTATION                                                  #
    ###########################################################################

    def _compute_cycle_equal_weighted_power(self, T_F, T_C) :
        """ Computes cycle burnups and peaking factors assuming 
            equal batch powers.
        """
        N = len(T_F)
        rho_L = self.p['leakage_penalty']       
        
        # shorten function call by eliminating p and boron dependence
        rho = lambda b, t_f, t_c: self.rho(self.p, b, t_f, t_c, 0.0)   
        
        # equal power sharing implies equal temperatures--using the average
        T_Fa, T_Ca = sp.mean(T_F), sp.mean(T_C)
        
        # linearize the reactivity, i.e., rho ~ rho_0 + AB.  (this may fail
        # if poison is still dominant at 10 GWd/MTU)
        B_a, B_b = 10.0, 20.0
        rho_a, rho_b = rho(B_a, T_Fa, T_Ca), rho(B_b, T_Fa, T_Ca)
        A = (rho_b-rho_a)/(B_b-B_a)
        rho_0 = rho_a - A*B_a
        
        # then B_s and B_c are *approximately*
        B_s = (rho_L - rho_0)/A
        B_c = 2.0*B_s/(len(T_F)+1)
        
        # solve f(B_c) = mean(rho)-rho_L = 0 via scipy's root finder
        f = lambda B : sp.mean(rho(B*sp.arange(1, N+1), T_Fa, T_Ca)) - rho_L
        B_c = root(f, B_c).x[0]

        # compute batch-wise, EOC burnups and associated peaking factors
        B = B_c * sp.arange(1, N+1)
        ppf = sp.ones(N)
        return B, ppf

    def _compute_cycle_reactivity_weighted_power(self, T_F, T_C) :
        """ Computes cycle burnups and peaking factors using 
            reactivity-weighted powers.
        """

        # get approximate burnups using equal-power approximation
        B, ppf = self._compute_cycle_equal_weighted_power(T_F, T_C)
        
        # solve f(B_c) = sum(w_i*rho_i) - rho_L = 0 via scipy's root finder
        f = lambda B : \
            self._reactivity_weighted_power_residual(B, T_F, T_C)            
        B = fsolve(f, B)

        # compute the batch power peaking factors
        f_bar = self._cycle_average_power_fraction(B, T_F, T_C, 0.0)
        ppf = len(T_F)*f_bar/sum(f_bar)
        return B, ppf
    
    def _power_fraction(self, B, T_F, T_C, C_B) :
        """ Computes the batch power fraction.
        """
        theta = 1.0 + (1./4.)*self.p['assembly_width']**2 / \
                self.m2(self.p, B, T_F, T_C, C_B);
        N = self.p['number_batches']
        f = (1.0/N) / (1.0-theta*self.rho(self.p, B, T_F, T_C, C_B))
        return f

    def _cycle_average_power_fraction(self, B, T_F, T_C, C_B) :
        """ Computes cycle-averaged power fractions.
        """
        f_bar = sp.zeros(self.p['number_batches'])
        BB = sp.array([0] + [i for i in B])
        for i in range(self.p['number_batches']) :
            #  evaluate f for fixed T's, etc. for any B
            fun = lambda x : self._power_fraction(x, T_F[i], T_C[i], C_B);
            # integrate f over the cycle burnup range and divide by range
            sol = sp.integrate.quad(fun, BB[i], BB[i+1], epsrel=1.0e-4)
            f_bar[i] = sol[0]
            f_bar[i] = f_bar[i] / (BB[i+1]-BB[i])
        return f_bar
    
    def _reactivity_weighted_power_residual(self, B, T_F, T_C) :
        """ Computes the nonlinear residual.
        
        The batch cycle burnups satisfy the nonlinear equation
        
                      B(1) - B(N)*f_bar(1)
                      B(2) - B(N)*f_bar(2)
              r(B) =     ...                     = 0
                      B(N-1) - B(N)*f_bar(N-1)
                      rho_bar - rho_L
        
        where r(B) is the residual.  It is assumed boron is at zero ppm
        for this calculation.
        """
        # both the eoc and cycle-averaged power fractions must be 
        # normalized.  
        f = self._power_fraction(B, T_F, T_C, 0.0)
        f = f / sum(f);         
        f_bar = self._cycle_average_power_fraction(B, T_F, T_C, 0.0)  
        f_bar = f_bar / sum(f_bar)
        # compute batch and core reactivities
        rho_batch = self.rho(self.p, B, T_F, T_C, 0.0)
        rho_bar = sum(f*rho_batch)
        # compute batch burnups
        B_computed = 0*B
        for i in range(len(T_F)) :
            B_computed[i] = B[-1] * sum(f_bar[0:i+1])
        # compute the residual 
        r = 0*B;
        r[:-1] = B[:-1] - B_computed[:-1]
        r[-1] = rho_bar - self.p['leakage_penalty']
        return r