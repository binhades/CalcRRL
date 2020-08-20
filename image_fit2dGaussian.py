#!/home/bliu/anaconda3/bin/python

# Filename : db_update_continue.py
# Aim: Load fits image, print the flux of given ds9 regions.

import argparse
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import Angle, SkyCoord
import numpy as np
import scipy.optimize as opt
from astropy import units as u
from astropy.nddata import Cutout2D
import pyrrl as rrl
def Gaussian2d(xydata, xo, yo, sigma_x, sigma_y, amplitude, offset0, offset1, a, b):
#def Gaussian2d(xydata, xo, yo, sigma_x, sigma_y, amplitude, offset0):
    """Function to fit, returns 2D gaussian function as 1D array"""
    (x,y) = xydata
    xo = float(xo)
    yo = float(yo)    
    g = offset0 + offset1*(x+a)*(y+b) +  amplitude*np.exp( - (((x-xo)**2)/(2*sigma_x**2) + ((y-yo)**2)/(2*sigma_y**2)))
#    g = offset0 + amplitude*np.exp( - (((x-xo)**2)/(2*sigma_x**2) + ((y-yo)**2)/(2*sigma_y**2)))
    return g.ravel()

def Gaussian2dFit(img,radius = 10):
    """
    Parameter:
        img - image as numpy array
    """
    x = np.linspace(0, img.shape[1], img.shape[1])
    y = np.linspace(0, img.shape[0], img.shape[0])
    x, y = np.meshgrid(x, y)
    #Parameters: xpos, ypos, sigmaX, sigmaY, amp, baseline
    #initial_guess = (img.shape[1]/2,img.shape[0]/2,radius,radius,img.max(),0)
    initial_guess = (img.shape[1]/2,img.shape[0]/2,radius,radius,img.max(),0,0,0,0)
    # subtract background and rescale image into [0,1], with floor clipping
#    bg = np.percentile(img,5)
#    img_scaled = np.clip((img - bg) / (img.max() - bg),0,1)
    popt, pcov = opt.curve_fit(Gaussian2d, (x, y), 
                               img.ravel(), p0=initial_guess)
#                               bounds = ((img.shape[1]*0.4, img.shape[0]*0.4, 1, 1, 0.5, -0.1),
#                                     (img.shape[1]*0.6, img.shape[0]*0.6, img.shape[1]/2, img.shape[0]/2, 1.5, 0.5)))
    #xcenter, ycenter, sigmaX, sigmaY, amp, offset0, offset1, a, b = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]
    xcenter, ycenter, sigmaX, sigmaY, amp, offset0 = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]
    FWHM_x = np.abs(4*sigmaX*np.sqrt(-0.5*np.log(0.5)))
    FWHM_y = np.abs(4*sigmaY*np.sqrt(-0.5*np.log(0.5)))
#    print(initial_guess)
#    print(popt)
    return (amp,xcenter,ycenter,FWHM_x,FWHM_y,offset0)

def load_image(filein):

    hdul = fits.open(args.file_fits,'readonly')
    wcs = WCS(hdul[0].header)
    img = np.squeeze(hdul[0].data)
    hdul.close()
    return img, wcs

def load_source():
    pos = SkyCoord('18h53m20.71s +1d14m29.23s', frame='icrs') # 34.254 0.142
    size = (25*u.arcmin, 25*u.arcmin)
    return pos, size

def main(args):
    img, wcs = load_image(args.file_fits)

    pos, size = load_source()

    img_cut = Cutout2D(img,position=pos,size=size,wcs=wcs)

    (Amp,xcenter,ycenter,FWHM_x,FWHM_y,offset0) = Gaussian2dFit(img_cut.data)
    print(xcenter,ycenter)

    coors = img_cut.wcs.wcs_pix2world([[xcenter,ycenter]],0)
    print(coors)

    (Amp,xcenter,ycenter,FWHM_x,FWHM_y,offset0) = rrl.continuum.Gaussian2dFit(img_cut.data)

    print("Peak Flux: {:6.2f}\nOffset: {:6.2f}\nX Size: {:5.5f} arcmin\nY Size: {:5.5f} arcmin".format(Amp,offset0,FWHM_x*np.abs(wcs.wcs.cdelt[1]*60), FWHM_y*np.abs(wcs.wcs.cdelt[1]*60)))

    return 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('file_fits', type=str, help='the input file name')
    args = parser.parse_args()
    main(args)

