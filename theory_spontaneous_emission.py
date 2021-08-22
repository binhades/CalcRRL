#!/usr/bin/python3
# Filename: bn.py
# Aim: to calculate the bn

import argparse, os, time
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl

def plot(n_arr,A,T,format='png',issave=False,isshow=True):
# =======================================================
    fig = plt.figure(figsize=(10,6))
# =======================================================

    y1_a = np.empty(n_arr.shape)
    y1_b = np.empty(n_arr.shape)
    y1_c = np.empty(n_arr.shape)
    
    for i,n in enumerate(n_arr):
        if i == 0:
            y1_a[i] = np.nan
            y1_b[i] = np.nan
            y1_c[i] = np.nan
        else:
            y1_a[i] = A[i,i-1]
            if i == 1:
                y1_b[i] = np.nan
                y1_c[i] = np.nan
            else:
                y1_b[i] = A[i,i-2]
                if i == 2:
                    y1_c[i] = np.nan
                else:
                    y1_c[i] = A[i,i-3]
        print('A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};'.format(i+1,i,y1_a[i],i+1,i-1,y1_b[i],i+1,i-2,y1_c[i]))
# --------------------------------
    ax1 = plt.subplot(221)
    ax1.plot(n_arr,y1_a,'-', label='$\Delta = 1$, A$_{n,n-1}$')
    ax1.plot(n_arr,y1_b,'-.',label='$\Delta = 2$, A$_{n,n-2}$')
    ax1.plot(n_arr,y1_c,'--',label='$\Delta = 3$, A$_{n,n-3}$')
    ax1.set_title('A$_{n,n-\Delta}$')
    ax1.set_xlabel('Primary level $n$')
    ax1.set_ylabel('Properbility ($s^{-1}$)')
    ax1.set_yscale('log')
    ax1.grid(True)
    ax1.legend()
    
# =======================================================
    y2_a = np.empty(n_arr.shape)
    y2_b = np.empty(n_arr.shape)
    y2_c = np.empty(n_arr.shape)
    for i,n in enumerate(n_arr):
        if i == 0:
            y2_a[i] = np.nan
            y2_b[i] = np.nan
            y2_c[i] = np.nan
        else:
            y2_a[i] = A[i,0]
            if i == 1:
                y2_b[i] = np.nan
                y2_c[i] = np.nan
            else:
                y2_b[i] = A[i,1]
                if i == 2:
                    y2_c[i] = np.nan
                else:
                    y2_c[i] = A[i,2]

        print('A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};'.format(i+1,1,y2_a[i],i+1,2,y2_b[i],i+1,3,y2_c[i]))
# --------------------------------
    ax2 = plt.subplot(222)
    ax2.plot(n_arr,y2_a,'-', label='n0 = 1, A$_{n,1}$')
    ax2.plot(n_arr,y2_b,'-.',label='n0 = 2, A$_{n,2}$')
    ax2.plot(n_arr,y2_c,'--',label='n0 = 3, A$_{n,3}$')
    ax2.set_title('A$_{n,n0}$')
    ax2.set_xlabel('Primary level $n$')
    ax2.set_ylabel('Properbility ($s^{-1}$)')
    ax2.set_yscale('log')
    ax2.grid(True)
    ax2.legend()
 # =======================================================
    y3_a = np.empty(n_arr.shape)
    y3_b = np.empty(n_arr.shape)
    y3_c = np.empty(n_arr.shape)
    for i,n in enumerate(n_arr):
        if i < 300-1:
            y3_a[i] = A[300-1,i]
            if i < 150 - 1:
                y3_b[i] = A[150-1,i]
                if i < 100 -1:
                    y3_c[i] = A[100-1,i]
                else:
                    y3_c[i] = np.nan
            else:
                y3_b[i] = np.nan
                y3_c[i] = np.nan
        else:
            y3_a[i] = np.nan
            y3_b[i] = np.nan
            y3_c[i] = np.nan

        print('A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};  A{:d}->{:d} = {:e};'.format(300,i+1,y3_a[i],150,i+1,y3_b[i],100,i+1,y3_c[i]))
# --------------------------------
    ax3 = plt.subplot(223)
    ax3.plot(n_arr,y3_a,'-', label='n0 = 300, A$_{300,n}$')
    ax3.plot(n_arr,y3_b,'-.',label='n0 = 150, A$_{150,n}$')
    ax3.plot(n_arr,y3_c,'--',label='n0 = 100, A$_{100,n}$')
    ax3.set_title('A$_{n0,n}$')
    ax3.set_xlabel('Primary level $n$')
    ax3.set_ylabel('Properbility ($s^{-1}$)')
    ax3.set_yscale('log')
    ax3.grid(True)
    ax3.legend()
 
# =======================================================
    ax4 = plt.subplot(224)
    ax4.plot(n_arr,T)
    ax4.set_title('Energy level life time')
    ax4.set_xlabel('Primary level $n$')
    ax4.set_ylabel('Lifetime ($s$)')
    ax4.set_yscale('log')
    ax4.grid(True)

# =======================================================
    plt.tight_layout()
    if issave:
        plt.savefig('spont_emiss.png',format=format,dpi=300,bbox_inches='tight')
    if isshow:
        plt.show()
    plt.close(fig)
    return 0

def main(args):

#==============================================================================
# Check functions from atom.py
# --------------------------------
#    n_u = args.N_u
#    n_l = args.N_l
#    f_ul = rrl.atom.oscillator_strength(n_u,n_l,method=args.method)
#    g_ul = rrl.atom.gaunt_factor(n_u,n_l,method=args.method)
#    A_ul0= rrl.atom.Einstein_coefficient_A(n_u,n_l,method=args.method)
#    A_ul1= rrl.atom.Einstein_coefficient_A(n_u,n_l)

# -----------------------------------------------------------------------------
# Print out results
# -----------------------------------------------------------------------------
#    print('N_u: {:d};\nN_l: {:d};\nMethod: {:d}'.format(n_u,n_l, args.method))
#    print('Oscillator Strength: f = {:f}'.format(f_ul))
#    print('Gaunt Factor:        g = {:f}'.format(g_ul))
#    print('Einstein A:          A = {:f}'.format(A_ul0))
#    print('Einstein A default:  A = {:f}'.format(A_ul1))
#    print('EnergyLeve Lifetime:     {:f}'.format(lt))
#==============================================================================
# TODO: n values and index confusing
    n_arr = np.arange(args.N_max) + 1 # ind: 0 -> N_max-1; n: 1 -> N_max

    A = np.zeros([args.N_max,args.N_max])
    T = np.zeros(args.N_max)
    for i, n_u in enumerate(n_arr):
        if i == 0: # n_u = 1, skip
            continue

        T[i] = rrl.atom.energy_level_lifetime(n_u)
        for j in range(i):
            n_l = n_arr[j]
            A[i,j]= rrl.atom.Einstein_coefficient_A(n_u,n_l)

    plot(n_arr,A,T)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--N_max', type=int,   default=300, \
            help='upper limit Primary quantum level, default 300')
    parser.add_argument('--method', type=int,   default=1,   help='method')

    args = parser.parse_args()

    start_time = time.time()
    main(args)
    print("--- {:.3f} seconds ---".format(time.time() - start_time))
