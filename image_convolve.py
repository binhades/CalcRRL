#!/home/bliu/anaconda3/bin/python

# Filename : tools_image_convolve.py
# Aim: to convolve the image to the beam size we wanted.

# this method is adopted from the Spectra-Cube pacakage, which deals with 3D data.

import argparse, os
from radio_beam import Beam
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astropy import convolution

def main(args):
    hdul = fits.open(args.file_in,'readonly')
    raw_img = hdul[0].data
    img_wcs = wcs.WCS(hdul[0].header)
    pixscale = wcs.utils.proj_plane_pixel_area(img_wcs.celestial)**0.5*u.deg

    raw_beam = Beam(args.beam*u.arcmin)
    new_beam = Beam(args.cbeam*u.arcmin)
    cov_beam = new_beam.deconvolve(raw_beam)
    cov_kernel = cov_beam.as_kernel(pixscale)

    new_img = convolution.convolve_fft(raw_img, cov_kernel, normalize_kernel=True)

    new_hdu = fits.PrimaryHDU(new_img)
    new_hdu.header = hdul[0].header.copy()

    new_hdu.writeto(args.fileout,overwrite=True)

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file_in', type=str, help='the input file name')
    parser.add_argument('beam', type=float, help='beam size of the raw image in arcmin')
    parser.add_argument('cbeam', type=float, help='beam size to convolve to in arcmin')
    parser.add_argument('--fileout', type=str, default='convolve.fits', help='the output file name')
    args = parser.parse_args()
    main(args)
