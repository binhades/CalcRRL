#!/usr/bin/python3

# Filename : image_pixel_distri.py
# Aim: to convolve the image to the beam size we wanted.

# this method is adopted from the Spectra-Cube pacakage, which deals with 3D data.

import argparse, os
import numpy as np
import regions
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord


def plot_corr(x, y, label=None,ir=24,fsave=None):
    #fig=plt.figure(figsize=(8,6))
    fig, ax = plt.subplots(figsize=(8,6))
    fz=18 # fontsize
    IR = '{}$\mu$m'.format(ir)
    ylabel = 'IR '+IR+' intensity (mJy$\,$sr$^{-1}$)'
    if label is None:
        plt.plot(x,y,'o')
        print('All ')
    else:
        ind_o = np.where(label == 0)[0] # position is HII free
        ind_i = np.where(label == 1)[0] # position is inside HII
#        ind_r = np.where(label ==-1)[0]
        if len(ind_o) > 0:
            x_t = x[ind_o]
            y_t = y[ind_o]
            plt.scatter(x_t,y_t,marker='o',label='WIM',alpha=0.5)
            print('Out ')
        if len(ind_i) > 0:
            x_t = x[ind_i]
            y_t = y[ind_i]
            plt.scatter(x_t,y_t,marker='^',label='HII',alpha=0.5)
            print('In ')
    ax.legend(fontsize=fz,markerscale=3)
    plt.xlabel(r'RRL integral intensity (Jy$\,$beam$^{-1}$ km$\,$s$^{-1}$)',fontsize=fz)
    plt.ylabel(ylabel,fontsize=fz)
    plt.xticks(fontsize=fz)
    plt.yticks(fontsize=fz)
    plt.xlim([0,0.5])
    plt.ylim([0,200])
    plt.tight_layout()
    if fsave is not None:
        plt.savefig(fsave+'.png', format='png',dpi=300)
        #plt.savefig(fsave+'.eps', format='eps',dpi=300)
    plt.show()

def load_regions(file_reg):

    reg_list = regions.read_ds9(file_reg)

    return reg_list

def get_ir_flux(data, position,size,wcs):

    cut = Cutout2D(data,position,size,wcs=wcs)

    return np.nanmean(cut.data)

def get_label(position,reg_list):

    for reg in reg_list:
        diff = position.separation(reg.center)
        if reg.radius < 60*u.arcsec:
            radius = 60*u.arcsec
        else:
            radius = reg.radius
        if diff > radius:
            continue
        else:
            return 1

    return 0
def main(args):
    rrl_hdul = fits.open(args.file_rrl,'readonly')
    rrl_img = rrl_hdul[0].data
    rrl_wcs = wcs.WCS(rrl_hdul[0].header)

    ir_hdul = fits.open(args.file_ir,'readonly')
    ir_img = ir_hdul[0].data                # [Dec, RA]
    ir_wcs = wcs.WCS(ir_hdul[0].header)

    reg_list = load_regions(args.file_reg)

    x_arr = np.arange(rrl_hdul[0].header['NAXIS1'])
    y_arr = np.arange(rrl_hdul[0].header['NAXIS2'])

    coor = rrl_wcs.all_pix2world(x_arr,y_arr,0)
    label = np.empty(rrl_img.shape)
    ir_arr = np.empty(rrl_img.shape)
    size=(10*u.arcsec,10*u.arcsec)
    for i,y in enumerate(coor[1]): # Dec
        for j,x in enumerate(coor[0]): # RA
            position=SkyCoord(x,y,frame='fk5',unit='deg').transform_to('galactic')
            ir_arr[i,j] = get_ir_flux(ir_img,position,size,ir_wcs)
            label[i,j] = get_label(position, reg_list)
            if rrl_img[i,j] < args.noise:
                label[i,j] = -1
            if label[i,j] == 1:
                print('Ind:{:02d},{:02d} - Coor: {:6.3f},{:+6.3f} [{:+2.0f}] - rrl:{:7.4f} - ir:{:8.2f}'.format(i,j,position.l.value,position.b.value,label[i,j],rrl_img[i,j],ir_arr[i,j]))

    plot_corr(rrl_img.flatten(),ir_arr.flatten(),label=label.flatten(),fsave=args.fileout,ir=args.ir)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file_rrl', type=str, help='the input RRL file name')
    parser.add_argument('file_ir', type=str, help='the input IR file name')
    parser.add_argument('file_reg', type=str, help='the input region file name')
    parser.add_argument('--fileout', type=str, help='the output file name')
    parser.add_argument('--noise', type=float, default=0.02, help='the noise level of RRL map')
    parser.add_argument('--ir', type=float, default=24, help='the IR wavelength')
    args = parser.parse_args()
    main(args)
