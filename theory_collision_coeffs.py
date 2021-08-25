#!/usr/bin/python3
# Filename: bn.py
# Aim: to calculate the bn

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl

def plot(n_arr, C_lu, C_ul, C_ni, C_in, labels, format='png',issave=False,isshow=True):
# =======================================================
    fig = plt.figure(figsize=(10,6))
# =======================================================
    ax1 = plt.subplot(221)
# --------------------------------
    y1_a = np.empty(n_arr.shape)
    y1_b = np.empty(n_arr.shape)
    y1_c = np.empty(n_arr.shape)
    
    for l, lab in enumerate(labels):
        for i,n in enumerate(n_arr):
            try:
                y1_a[i] = C_lu[i,i+1,l]
            except IndexError:
                y1_a[i] = np.nan
            try:
                y1_b[i] = C_lu[i,i+2,l]
            except IndexError:
                y1_b[i] = np.nan
            try:
                y1_c[i] = C_lu[i,i+3,l]
            except IndexError:
                y1_c[i] = np.nan

# --------------------------------
        ax1.plot(n_arr,y1_a,'-', label='$\Delta = 1$, C$_{n,n+1}$, '+lab)
        ax1.plot(n_arr,y1_b,'-.',label='$\Delta = 2$, C$_{n,n+2}$, '+lab)
        ax1.plot(n_arr,y1_c,'--',label='$\Delta = 3$, C$_{n,n+3}$, '+lab)
    
    ax1.set_title('Rate Coefficient of Collision excitation')
    ax1.set_xlabel('Primary level $n$')
    ax1.set_ylabel('Rate ($cm^{3} s^{-1}$)')
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.legend()
 # =======================================================
    ax2 = plt.subplot(222)
# --------------------------------
    y2_a = np.empty(n_arr.shape)
    y2_b = np.empty(n_arr.shape)
    y2_c = np.empty(n_arr.shape)
    
    for l, lab in enumerate(labels):
        for i,n in enumerate(n_arr):
            try:
                y2_a[i] = C_ul[i,i-1,l]
            except IndexError:
                y2_a[i] = np.nan
            try:
                y2_b[i] = C_ul[i,i-2,l]
            except IndexError:
                y2_b[i] = np.nan
            try:
                y2_c[i] = C_ul[i,i-3,l]
            except IndexError:
                y2_c[i] = np.nan
# --------------------------------
        ax2.plot(n_arr,y2_a,'-', label='$\Delta = 1$, C$_{n,n-1}$, '+lab)
        ax2.plot(n_arr,y2_b,'-.',label='$\Delta = 2$, C$_{n,n-2}$, '+lab)
        ax2.plot(n_arr,y2_c,'--',label='$\Delta = 3$, C$_{n,n-3}$, '+lab)
    
    ax2.set_title('Rate Coefficient of Collision de-excitation')
    ax2.set_xlabel('Primary level $n$')
    ax2.set_ylabel('Rate ($cm^{3} s^{-1}$)')
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.legend()
    
# =======================================================
    ax3 = plt.subplot(223)
# --------------------------------
    y3 = np.empty(n_arr.shape)
    
    for l, lab in enumerate(labels):
        for i,n in enumerate(n_arr):
            y3[i] = C_ni[i,l]
# --------------------------------
        ax3.plot(n_arr,y3,'-', label='C$_{n,i}$, '+lab)
    
    ax3.set_title('Rate Coefficient of Collision ionization')
    ax3.set_xlabel('Primary level $n$')
    ax3.set_ylabel('Rate ($cm^{3} s^{-1}$)')
    ax3.set_yscale('log')
    ax3.grid(True)
    ax3.legend()
   
# =======================================================

    ax4 = plt.subplot(224)
# --------------------------------
    y4 = np.empty(n_arr.shape)
    
    for l, lab in enumerate(labels):
        for i,n in enumerate(n_arr):
            y4[i] = C_in[i,l]
# --------------------------------
        ax4.plot(n_arr,y4,'-', label='C$_{i,n}$, '+lab)

    ax4.set_title('Rate Coefficient of Three-body recombination')
    ax4.set_xlabel('Primary level $n$')
    ax4.set_ylabel('Rate ($cm^{6} s^{-1}$)')
    ax4.set_yscale('log')
    ax4.grid(True)
    ax4.legend()

# =======================================================
    plt.tight_layout()
    if issave:
        plt.savefig('spont_emiss.png',format=format,dpi=300,bbox_inches='tight')
    if isshow:
        plt.show()
    plt.close(fig)
    return 0

#==============================================================================
def main(args):


    if args.Te is None:
        Te = [5000, 10000, 20000]
        labels = ['Te = 5e3 K', 'Te = 1e4 K', 'Te = 2e4 K']
    else:
        Te = [args.Te]
        labels = ['Te = {:e} K'.format(args.Te)]

    n_arr = np.arange(args.N_max) + 1 # ind: 0 -> N_max-1; n: 1 -> N_max
    C_ul = np.full([args.N_max,args.N_max,len(Te)],np.nan)
    C_lu = np.full([args.N_max,args.N_max,len(Te)],np.nan)
    C_ni = np.full([args.N_max,len(Te)],np.nan)
    C_in = np.full([args.N_max,len(Te)],np.nan)

    for i, n in enumerate(n_arr):
        for j, T in enumerate(Te):
            C_ni[i,j] = rrl.atom.collision_ionization_rate(n,T,method='V80')
            C_in[i,j] = rrl.atom.threebody_recombination_rate(n,T)

            n_u = n
            for k in range(i):
                n_l = n_arr[k]
                C_ul[i,k,j]= rrl.atom.collision_deexcitation_rate(n_u,n_l,T)
            n_l = n
            for k in range(i+1,args.N_max-1):
                n_u = n_arr[k]
                C_lu[i,k,j]= rrl.atom.collision_excitation_rate(n_l,n_u,T,method='Integ')

    plot(n_arr, C_lu, C_ul, C_ni, C_in, labels)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--N_max', type=int,   default=300, \
            help='upper limit Primary quantum level, default 300')
    parser.add_argument('--Te',   type=float,  \
            help='Electron Temperature in K, default 10000 K')
    args = parser.parse_args()
    start_time = time.time()
    main(args)
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
