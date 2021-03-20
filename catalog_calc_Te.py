#!/usr/bin/python3

# calc_te.py
# to Calculate electron temperature 
# 1. load source and line parameters from csv catalog
# 2. read continuum flux from image.

import argparse, os
import numpy as np
import astropy 
import regions
import fastspec as fsp
import pyrrl as rrl
import astropy.units as u

def load_catalog_csv(file_csv):
    import csv
    catalog = []
    with open(file_csv,mode='rt') as fcsv:
        reader = csv.DictReader(fcsv)
        for row in reader:
            catalog.append(row)
    return catalog 

def flux_from_image(file_img,l,b,r=10,method='fit'):

    hdul = astropy.io.fits.open(file_img,'readonly')
    wcs = astropy.wcs.WCS(hdul[0].header)
    img = np.squeeze(hdul[0].data)
    hdul.close()

    size = (r*u.arcmin,r*u.arcmin)
    pos = astropy.coordinates.SkyCoord(l,b,unit='deg',frame='galactic')

    img_cut = astropy.nddata.Cutout2D(img,position=pos,size=size,wcs=wcs)
    if method == 'fit':
        (Amp,xcenter,ycenter,FWHM_x,FWHM_y,offset0) = rrl.continuum.Gaussian2dFit(img_cut.data)
#    print("Peak Flux: {:8.4f}\nOffset: {:8.4f}\nX Size: {:6.2f} arcmin\nY Size: {:6.2f} arcmin".format(Amp,offset0,FWHM_x*np.abs(wcs.wcs.cdelt[1]*60), FWHM_y*np.abs(wcs.wcs.cdelt[1]*60)))
        return Amp
    elif method == 'peak':
        coor = np.round(img_cut.wcs.world_to_pixel(pos))
        return img_cut.data[round(coor[0]),round(coor[1])]
    elif method == 'mean':
        return np.mean(img_cut.data)


def main(args):

    if not os.path.isfile(args.file_csv):
        print("CSV file not found:",args.file_csv)
        return 0

    catalog = load_catalog_csv(args.file_csv)

    for item in catalog:
        #print(item)
        l = float(item['GLon'])
        b = float(item['GLat'])
        T_l  = float(item['Peak']  ) /1000. # mJy -> Jy
        d_V = float(item['FWHM']  )

        T_c = flux_from_image(args.cont_imag,l,b,r=10,method=args.method)
        T_e = rrl.calc.Te(T_c,T_l,d_V,freq=args.frequency,p_he=0.1)
        #print('Source:  {};\nPeak at: G{:07.3f}{:+07.3f};'.format(item['GName'],l,b))
        #print('TC: {:8.4f};\nTL: {:8.4f}, delta_V: {:6.2f} km/s;\nTe:{:6.0f} K'.format(T_c,T_l,d_V,T_e))
        print('{index} {gname} {tc:7.4f} {tl:7.4f} {te:5.0f}'.format(index=item['Index'],gname=item['GName'],tc=T_c,tl=T_l,te=T_e))

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=str, help='the source catalog')
    parser.add_argument('cont_imag', type=str, help='input continuum map')
    parser.add_argument('--method', type=str, help='method to get continuum: fit or peak')
    parser.add_argument('--frequency', type=float,default=1.4, help='the line rest frequency in GHz')
    parser.add_argument('--size', type=float,default=5, help='the size to average ')

    args = parser.parse_args()
    main(args)

