#!/usr/bin/python3

# given a source in a catalog and a set of datacubes of a serial of RRLs
# plot the line intensities v.s. RRLs (frequencies)

import argparse, os, glob, re
import numpy as np
import astropy 
import regions
import fastspec as fsp
import pyrrl as rrl
import astropy.units as u
from spectral_cube import SpectralCube

def plot_line_intensity(source, line_list, freq_list, flux_list, ferr_list):
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12,8))
    plt.errorbar(line_list,flux_list,yerr=ferr_list,marker='o')
    locs, labels = plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    ax1 = plt.gca()

    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(locs)
    ax2.set_xticklabels(freq_list, fontsize=14)
    ax2.tick_params(axis='x',direction='in',labelrotation=30)
    ax2.set_xlabel("Line Rest Frequency (MHz)", fontsize=14)

    ax1.set_xlabel("RRLs", fontsize=14)
    ax1.set_ylabel(r"T$_{MB}$ (K)", fontsize=14)
    ax1.tick_params(axis='x',direction='in',labelrotation=30)

#    plt.title(source,fontsize=14)
    plt.tight_layout()
    plt.savefig(source+'_RRL-flux'+'.png',format='png',dpi=300)
    plt.savefig(source+'_RRL-flux'+'.eps',format='eps',dpi=300)

    plt.show()

    return 0

def main(args):
    
    n_arr = np.arange(args.n_start,args.n_stop+1)
    N_arr = np.zeros(n_arr.shape,dtype=np.float64)

    N_1b = 1
    N_1s = rrl.calc.Saha_equation(1,Ne=args.Ne,Te=args.Te)
    print('n = {:3d} : Saha (LTE) Level population = {:e}, Boltzman ratio: {:e}'.format(1,N_1s,1))
    for i,n in enumerate(n_arr):
        r_N_u2l = rrl.calc.Boltzmann_equation(n,delta=n-1,T=args.Te)
        N_n = rrl.calc.Saha_equation(n,Ne=args.Ne,Te=args.Te)
        N_arr[i] = N_1b*r_N_u2l
        print('n = {:3d} : Saha (LTE) Level population = {:e}, Boltzman ratio: {:e}'.format(n,N_n, N_arr[i]))

    if args.plot:
        plot_line_intensity(source['gname'],line_list,freq_list,flux_list,ferr_list)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_start', type=int, default=40, help='the energy level to start')
    parser.add_argument('--n_stop', type=int, default=300, help='the energy level to stop')
    parser.add_argument('--Te', type=float, default=1e4, help='the Temperature')
    parser.add_argument('--Ne', type=float, default=1e3, help='the electron density')
    parser.add_argument('--plot', action='store_true', help='set to plot')

    args = parser.parse_args()
    main(args)

