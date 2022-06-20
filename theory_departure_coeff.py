#!/usr/bin/python3

# Filename: bn.py
# Aim: to calculate the bn

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl

def plot(n_arr, bn, beta_n, labels, format='png', issave=False, isshow=True):
#==============================================================================
    fig = plt.figure(figsize=(10,6))
#==============================================================================
    ax1 = plt.subplot(211)
    for i, lab in enumerate(labels):
        ax1.plot(n_arr,bn[i],label=lab)

    ax1.set_title('Departure Coefficient')
    ax1.set_xlabel('Principal energy level $n$')
    ax1.set_ylabel('$b_n$')
    ax1.legend()
#==============================================================================
    ax2 = plt.subplot(212)
    for i, lab in enumerate(labels):
        ax2.plot(n_arr,beta_n[i],label=lab)

    ax2.set_title('Amplification Factor')
    ax2.set_xlabel('Principal energy level $n$')
    ax2.set_ylabel(r'$\beta_n$')
    ax2.legend()
#==============================================================================
    plt.tight_layout()
    if issave:
        plt.savefig('departure_coeff.png',format=format,dpi=300,bbox_inches='tight')
    if isshow:
        plt.show()
    plt.close(fig)
    return 0
#==============================================================================

def main(args):

    n_arr = np.arange(args.N_min,args.N_max+1)
    Te = args.Te

    if args.Ne is None:
        Ne_list = [10, 1e2, 1e3, 1e4]
    else:
        Ne_list = [args.Ne]

    labels = []
    bn_list = []
    beta_n_list = []

    for Ne in Ne_list:
        labels.append('Te:{:.0f}K; Ne:{:.0E}'.format(Te,Ne)+'cm$^{-3}$;')
        bn = rrl.bn.departure_coefficients(n_min=args.N_min, n_max=args.N_max, \
                                           n_cut=args.N_cut, Te=Te, Ne=Ne)
        beta_n = rrl.bn.amplification_factor(bn,Te)
        bn_list.append(bn)
        beta_n_list.append(beta_n)
    plot(n_arr,bn_list,beta_n_list,labels,issave=args.save)
    return 0
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--Te',   type=float, default=1e4, \
            help='Electron Temperature in K')
    parser.add_argument('--Ne',   type=float, \
            help='Electron Density in cm-3')
    parser.add_argument('--N_cut', type=int,   default=800, \
            help='cut of Primary quantum level, default 800')
    parser.add_argument('--N_min', type=int,   default=20, \
            help='lower limit of Primary quantum level, default 30')
    parser.add_argument('--N_max', type=int,   default=300, \
            help='upper limit Primary quantum level, default 300')
    parser.add_argument('--case', type=int,   default=2, \
            help='1 for Case A, 2 for Case B, default 2')
    parser.add_argument('--method', type=str,   default='S59', \
            help='method')
    parser.add_argument('--save', action='store_true',help='set to save figure')

    args = parser.parse_args()

    start_time = time.time()
    main(args)
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
