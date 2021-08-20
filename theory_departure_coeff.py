#!/usr/bin/python3
# Filename: bn.py
# Aim: to calculate the bn

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl

def main(args):

#==============================================================================
# Check functions from atom.py
# --------------------------------
#    n_u = args.N_u
#    n_l = args.N_l
#    f_ul = rrl.atom.oscillator_strength(n_u,n_l,method=args.method)
#    g_ul = rrl.atom.gaunt_factor(n_u,n_l,method=args.method)
#    A_ul0= rrl.atom.Einstein_Coefficient_A(n_u,n_l,method=args.method)
#    A_ul1= rrl.atom.Einstein_Coefficient_A(n_u,n_l)
#    lt = rrl.atom.energy_level_lifetime(n_u)
#    alpha = rrl.atom.recombination_rate(args.N,args.Te,method=args.method)

# -----------------------------------------------------------------------------
# Print out results
# -----------------------------------------------------------------------------
#    print('N_u: {:d};\nN_l: {:d};\nMethod: {:d}'.format(n_u,n_l, args.method))
#    print('Oscillator Strength: f = {:f}'.format(f_ul))
#    print('Gaunt Factor:        g = {:f}'.format(g_ul))
#    print('Einstein A:          A = {:f}'.format(A_ul0))
#    print('Einstein A default:  A = {:f}'.format(A_ul1))
#    print('EnergyLeve Lifetime:     {:f}'.format(lt))
#    print('Recombination rate: {:e}'.format(alpha)) 
#==============================================================================
    n_arr = np.arange(args.N_min,args.N_max+1)
    bn = rrl.bn.departure_coefficients(n_min=args.N_min,n_max=args.N_max,Te=args.Te,Ne=args.Ne)
    plt.plot(n_arr,bn)
    plt.show()
    return 0
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Te',   type=float, default=10000, \
            help='Electron Temperature in K, default 10000 K')
    parser.add_argument('--Ne',   type=float, default=10000, \
            help='Electron Density in cm-3, default 10000 cm-3')
    parser.add_argument('--N',    type=int,   default=40, \
            help='Primary quantum level, default 40')
    parser.add_argument('--N_min', type=int,   default=40, \
            help='lower limit of Primary quantum level, default 10')
    parser.add_argument('--N_max', type=int,   default=200, \
            help='upper limit Primary quantum level, default 507')
    parser.add_argument('--N_l',  type=int,   default=10, \
            help='lower level, default 10')
    parser.add_argument('--N_u',  type=int,   default=11, \
            help='upper level, default 11')
    parser.add_argument('--case', type=int,   default=2, \
            help='1 for Case A, 2 for Case B, default 2')
    #parser.add_argument('--method', type=int,   default=1,   help='method')
    parser.add_argument('--method', type=str,   default='S59', \
            help='method')

    args = parser.parse_args()

    start_time = time.time()
    main(args)
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
