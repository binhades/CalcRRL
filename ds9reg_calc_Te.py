#!/usr/bin/python3

# calc_te.py
# to Calculate electron temperature 
# 1. load source from ds9-reg file
# 2. load source spectrum from data cube and then fit.
# 3. read continuum flux from image.

import argparse, os
import numpy as np
import astropy 
import regions
import fastspec as fsp
import pyrrl as rrl
import astropy.units as u
from spectral_cube import SpectralCube

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

def rrl_from_cube(file_cube,source):

    cube = SpectralCube.read(file_cube)
    velo = cube.spectral_axis.value

    reg_str = "galactic; circle({},{},{})".format(source['glon'],source['glat'],1/60.)
    sub_cube = cube.subcube_from_ds9region(reg_str)
    spec = sub_cube.mean(axis=(1,2))

    yfit,fpeak,vlsr,fwhm,e1,e2,e3 = rrl.spec.fit(velo,spec)

    return fpeak,fwhm

def flux_from_image(file_img,source):

    hdul = astropy.io.fits.open(file_img,'readonly')
    wcs = astropy.wcs.WCS(hdul[0].header)
    img = np.squeeze(hdul[0].data)
    hdul.close()

    l = source['glon']
    b = source['glat']
    r = source['radius']

    size = (25*u.arcmin,25*u.arcmin)
    pos = astropy.coordinates.SkyCoord(l,b,unit='deg',frame='galactic')
#    print(pos.transform_to('icrs'))

    img_cut = astropy.nddata.Cutout2D(img,position=pos,size=size,wcs=wcs)
    (Amp,xcenter,ycenter,FWHM_x,FWHM_y,offset0) = rrl.continuum.Gaussian2dFit(img_cut.data)
#    print("Peak Flux: {:8.4f}\nOffset: {:8.4f}\nX Size: {:6.2f} arcmin\nY Size: {:6.2f} arcmin".format(Amp,offset0,FWHM_x*np.abs(wcs.wcs.cdelt[1]*60), FWHM_y*np.abs(wcs.wcs.cdelt[1]*60)))

    coors = img_cut.wcs.wcs_pix2world([[xcenter,ycenter]],0)

    source['glon'] = coors[0][0]
    source['glat'] = coors[0][1]

    return Amp

def main(args):
    if not os.path.isfile(args.spec_cube):
        print("Cube file not found:",args.spec_cube)
        return 0
    if not os.path.isfile(args.cont_imag):
        print("Cube file not found:",args.cont_imag)
        return 0

    name_list,source_list = load_catalog(args.catalog)

    if args.source in name_list:
        ind = name_list.index(args.source)
        source = source_list[ind]
        tc = flux_from_image(args.cont_imag,source)
        tl, delta_v = rrl_from_cube(args.spec_cube, source)
        te = rrl.calc.Te(tc,tl,delta_v,freq=args.frequency,p_he=0.1)
        print('Source:  {};\nPeak at: G{:07.3f}{:+07.3f};'.format(args.source,source['glon'],source['glat']))
        print('TC: {:8.4f};\nTL: {:8.4f}, delta_V: {:6.2f} km/s;\nTe:{:6.0f} K'.format(tc,tl,delta_v,te))

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('spec_cube', type=str, help='input RRL data cube')
    parser.add_argument('cont_imag', type=str, help='input continuum map')
    parser.add_argument('catalog', type=str, help='the source catalog')
    parser.add_argument('--source', type=str, help='the Source name to calculate')
    parser.add_argument('--frequency', type=float,default=1.4, help='the line rest frequency in GHz')

    args = parser.parse_args()
    main(args)

