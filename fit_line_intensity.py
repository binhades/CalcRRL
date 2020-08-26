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

def parse_file(string):

    arr = re.search("_(H[0-9]{3}a)_",string)
    return arr.group(1)

def load_catalog(file_reg):
    reg_list = regions.read_ds9(file_reg)
    source_list = []
    name_list = []
    for reg in reg_list:
        # source: {'gname:gname,'glon':glon,'glat':glat,'radius':radius}
        source = fsp.db.region_to_source(reg)
        source_list.append(source)
        name_list.append(source['gname'])

    return name_list,source_list

def rrl_from_cube(file_cube,source,paras=None):

    cube = SpectralCube.read(file_cube)
    velo = cube.spectral_axis.value

    reg_str = "galactic; circle({},{},{})".format(source['glon'],source['glat'],1/60.)
    sub_cube = cube.subcube_from_ds9region(reg_str)
    spec = sub_cube.mean(axis=(1,2)).value
    print(type(spec))

    yfit,fpeak,vlsr,fwhm,e1,e2,e3 = rrl.spec.fit(velo,spec,paras=paras)

    err = np.sqrt(np.nanstd(spec-yfit)**2+e1**2)

    return fpeak,err
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

    cube_list = sorted(glob.glob(args.file_prefix))
    init = [0.1,55,15]
    if len(cube_list) < 1:
        print("cannot find file")
        return 0

    if args.catalog is None:
        if args.l is None or args.b is None or args.size is None:
            print('Nothing to do')
            return 0
        else:
            glon = args.l
            glat = args.b
            size = args.size
            gname='G{:07.3f}{:+.3f}'.format(glon,glat)
        # source: {'gname:gname,'glon':glon,'glat':glat,'radius':radius}
            source = {'gname':gname,'glon':glon,'glat':glat,'radius':size}
    elif os.path.isfile(args.catalog):
        name_list,source_list = load_catalog(args.catalog)
        if args.source in name_list:
            ind = name_list.index(args.source)
            source = source_list[ind]
        else:
            print("Source not exists")
            return 0
    else:
        print("Catalog not found:",args.catalog)
        return 0

    line_list=[]
    freq_list=[]
    flux_list=[]
    ferr_list=[]
    lines = fsp.rrl.get_lines(element='hydrogen',ltype='alpha')
    for i, fcube in enumerate(cube_list):
        line = parse_file(fcube)
        freq = lines[line]
        flux, ferr = rrl_from_cube(fcube, source,paras=init)
        print('{} ; Freq {:8.3f}'.format(line,freq))
        print('TL: {:8.4f}, err: {:6.2f};'.format(flux,ferr))
        line_list.append(line)
        freq_list.append(freq)
        flux_list.append(flux)
        ferr_list.append(ferr)

    if args.plot:
        plot_line_intensity(source['gname'],line_list,freq_list,flux_list,ferr_list)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_prefix', type=str, help='input RRL data cubes')
    parser.add_argument('--catalog', type=str, help='the source catalog')
    parser.add_argument('--source', type=str, help='the Source name to calculate')
    parser.add_argument('--l', type=float, help='the Source Glon')
    parser.add_argument('--b', type=float, help='the Source Glat')
    parser.add_argument('--size', type=float, help='the Source size')
    parser.add_argument('--plot', action='store_true', help='set to plot')

    args = parser.parse_args()
    main(args)

