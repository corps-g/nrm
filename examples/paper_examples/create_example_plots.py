# Although MATLAB and Octave can be used to make nice plots, it's hard to 
# do it for both at the same time due to some apparent API inconsistencies. 
# Rather, all the output for the examples was saved as MATLAB ".mat" files and
# are read here by SciPy for plotting with matplotlib.

import matplotlib.pyplot as plt
import matplotlib as mpl
# Many of these are set to reproduce MPL v1x-style plots that were in the
# published paper (for better or worse).
mpl.rcParams['figure.figsize'] = [8.0, 6.0]
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['legend.numpoints'] = 1
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['legend.edgecolor'] = 'inherit'
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['lines.dashed_pattern'] = [6, 6]
mpl.rcParams['lines.dashdot_pattern'] = [3, 5, 1, 5]
mpl.rcParams['lines.dotted_pattern'] = [1, 3]
mpl.rcParams['lines.scale_dashes'] = False

import pickle
import numpy as np

dash_dot_dot = [8, 4, 2, 4, 2, 4]
long_dash = [16, 4]

# example 1 -------------------------------------------------------------------#
plt.clf()
data = pickle.load(open('example_1.p','rb'))
B = data['B']
k_inf_Zr4 = data['k_inf_Zr4']
k_inf_FeCrAl_100 = data['k_inf_FeCrAl_100']
k_inf_FeCrAl_300 = data['k_inf_FeCrAl_300']
k_inf_SiC_100 = data['k_inf_SiC_100']
k_inf_SiC_300 = data['k_inf_SiC_300']
rho_FeCrAl_100 = data['rho_FeCrAl_100']
rho_SiC_100 = data['rho_SiC_100']
rho_FeCrAl_300 = data['rho_FeCrAl_300']
rho_SiC_300 = data['rho_SiC_300']

plt.figure(1) 
plt.plot(B, k_inf_Zr4,        'k-')
plt.plot(B, k_inf_FeCrAl_100, 'k--')
plt.plot(B, k_inf_FeCrAl_300, 'k-.')
plt.plot(B, k_inf_SiC_100,    'k:')
plt.plot(B, k_inf_SiC_300,    'k-', dashes=dash_dot_dot)
plt.xlabel('burnup (GWd/MTU)')
plt.ylabel('$k_{\infty}$')
plt.legend(['Zr-4', 
            'FeCrAl, 100 $\mu$m', 
            'FeCrAl, 300 $\mu$m',
            'SiC, 100 $\mu$m',
            'SiC, 300 $\mu$m'], loc=0)
plt.savefig('example_1a.pdf')

plt.figure(2) 
plt.plot(B, rho_FeCrAl_100, 'k-',  
         B, rho_FeCrAl_300, 'k--', 
         B, rho_SiC_100,    'k-.',
         B, rho_SiC_300,    'k:')         
plt.xlabel('burnup (GWd/MTU)')
plt.ylabel('reactivity defect (pcm)')
plt.legend(['FeCrAl, 100 $\mu$m', 
            'FeCrAl, 300 $\mu$m',
            'SiC, 100 $\mu$m',
            'SiC, 300 $\mu$m'],  loc=(0.02,0.2))
plt.savefig('example_1b.pdf')

# example 2 -------------------------------------------------------------------#
data = pickle.load(open('example_2.p', 'rb'))
t = data['thick']
B_c_FeCrAl = data['B_c_FeCrAl']
B_c_SiC = data['B_c_SiC']
plt.figure(3) 
plt.plot(t, B_c_FeCrAl[0], 'k-', 
         t, B_c_FeCrAl[1], 'k--',
         t, B_c_SiC[0], 'k-.',
         t, B_c_SiC[1], 'k:')
plt.xlabel('thickness ($\mu$m)')
plt.ylabel('cycle length (GWd/MTU)')
plt.legend(['FeCrAl, equal', 'FeCrAl, unequal', 
            'SiC, equal', 'SiC, unequal'], loc=(0.02,0.05))
plt.savefig('example_2.pdf')

# example 3 -------------------------------------------------------------------#
data = pickle.load(open('example_3.p', 'rb'))
t = data['thick']
T_F = np.zeros((2, 3, len(t)))
T_C = np.zeros((2, 3, len(t)))
PPF = np.zeros((2, 3, len(t)))
T_F[0] = data['T_F_FeCrAl']
T_F[1] = data['T_F_SiC']
T_C[0] = data['T_C_FeCrAl']
T_C[1] = data['T_C_SiC']
PPF[0] = data['PPF_FeCrAl']
PPF[1] = data['PPF_SiC']

plt.figure(4) 
plt.plot(t, T_F[0][0], 'k-', 
         t, T_F[0][1], 'k--', 
         t, T_F[0][2], 'k-.',
         t, T_F[1][0], 'k:')
plt.plot(t, T_F[1][1], 'k-', dashes=dash_dot_dot)
plt.plot(t, T_F[1][2], 'k-', dashes=long_dash)
plt.xlabel('thickness ($\mu$m)')
plt.ylabel('fuel temperature (K)')
plt.legend(['FeCrAl, batch 1', 'FeCrAl, batch 2', 'FeCrAl, batch 3',
            'SiC, batch 1', 'SiC, batch 2', 'SiC, batch 3'], loc=(0.02,0.4))
plt.savefig('example_3a.pdf')

plt.figure(5) 
plt.plot(t, T_C[0][0], 'k-', 
         t, T_C[0][1], 'k--', 
         t, T_C[0][2], 'k-.',
         t, T_C[1][0], 'k:')
plt.plot(t, T_C[1][1], 'k-', dashes=dash_dot_dot)
plt.plot(t, T_C[1][2], 'k-', dashes=long_dash)
plt.xlabel('thickness ($\mu$m)')
plt.ylabel('coolant temperature (K)')
plt.legend(['FeCrAl, batch 1', 'FeCrAl, batch 2', 'FeCrAl, batch 3',
            'SiC, batch 1', 'SiC, batch 2', 'SiC, batch 3'], loc=(0.02,0.4))
plt.savefig('example_3b.pdf')

plt.figure(6) 
plt.plot(t, PPF[0][0], 'k-', 
         t, PPF[0][1], 'k--', 
         t, PPF[0][2], 'k-.',
         t, PPF[1][0], 'k:')
plt.plot(t, PPF[1][1], 'k-', dashes=dash_dot_dot)
plt.plot(t, PPF[1][2], 'k-', dashes=long_dash)
plt.xlabel('thickness ($\mu$m)')
plt.ylabel('power-peaking factor')
plt.legend(['FeCrAl, batch 1', 'FeCrAl, batch 2', 'FeCrAl, batch 3',
            'SiC, batch 1', 'SiC, batch 2', 'SiC, batch 3'], loc=(0.02,0.4))
plt.savefig('example_3c.pdf')

