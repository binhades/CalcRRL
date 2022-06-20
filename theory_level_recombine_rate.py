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
    rate = rrl.calc.recombine_rate_H(n_arr,T=args.Tk)

    for i,n in enumerate(n_arr):
        print('n = {:d} -> An = {:e}'.format(n,rate[i]))

#    cube_list = sorted(glob.glob(args.file_prefix))
#    init = [0.1,55,15]
#    if len(cube_list) < 1:
#        print("cannot find file")
#        return 0
#
#    if args.catalog is None:
#        if args.l is None or args.b is None or args.size is None:
#            print('Nothing to do')
#            return 0
#        else:
#            glon = args.l
#            glat = args.b
#            size = args.size
#            gname='G{:07.3f}{:+.3f}'.format(glon,glat)
#        # source: {'gname:gname,'glon':glon,'glat':glat,'radius':radius}
#            source = {'gname':gname,'glon':glon,'glat':glat,'radius':size}
#    elif os.path.isfile(args.catalog):
#        name_list,source_list = load_catalog(args.catalog)
#        if args.source in name_list:
#            ind = name_list.index(args.source)
#            source = source_list[ind]
#        else:
#            print("Source not exists")
#            return 0
#    else:
#        print("Catalog not found:",args.catalog)
#        return 0
#
#    line_list=[]
#    freq_list=[]
#    flux_list=[]
#    ferr_list=[]
#    lines = fsp.rrl.get_lines(element='hydrogen',ltype='alpha')
#    for i, fcube in enumerate(cube_list):
#        line = parse_file(fcube)
#        freq = lines[line]
#        flux, ferr = rrl_from_cube(fcube, source,init=init)
#        print('{} ; Freq {:8.3f}'.format(line,freq))
#        print('TL: {:8.4f}, err: {:6.2f};'.format(flux,ferr))
#        line_list.append(line)
#        freq_list.append(freq)
#        flux_list.append(flux)
#        ferr_list.append(ferr)
#
#    if args.plot:
#        plot_line_intensity(source['gname'],line_list,freq_list,flux_list,ferr_list)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--Tk', type=float, default=10000, help='the kinematic temperature')
    parser.add_argument('--n_start', type=int, default=40, help='the energy level to start')
    parser.add_argument('--n_stop', type=int, default=300, help='the energy level to stop')
    parser.add_argument('--plot', action='store_true', help='set to plot')

    args = parser.parse_args()
    main(args)

