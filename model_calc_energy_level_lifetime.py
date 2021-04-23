#!/home/bliu/anaconda3/bin/python

##!/usr/bin/python3
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
    #rate = rrl.calc.spontaneous_emission_rate_H(n_arr,delta=args.delta_n)
    A_ul,B_ul,B_lu = rrl.calc.Einstein_Coefficients_H(n_arr,delta=args.delta_n)
    life = 1./A_ul

    for i,n in enumerate(n_arr):
        print('n = {:3d} -> {:3d} : lifetime = {:e} s, A_ul = {:e}, B_ul = {:e}, B_lu = {:e}'.format(n,n-args.delta_n,life[i],A_ul[i],B_ul[i],B_lu[i]))

#    if args.plot:
#        plot_line_intensity(source['gname'],line_list,freq_list,flux_list,ferr_list)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_start', type=int, default=40, help='the energy level to start')
    parser.add_argument('--n_stop', type=int, default=300, help='the energy level to stop')
    parser.add_argument('--delta_n', type=int, default=1, help='the difference of energy levels to jump')
    parser.add_argument('--plot', action='store_true', help='set to plot')

    args = parser.parse_args()
    main(args)

