#! /usr/bin/python3

# Filename: db_fit_spec_mannual.py
# Aim: to load the data base and then to fit spectra.

import argparse, sqlite3
import numpy as np
from astropy.modeling import models 
from astropy.modeling.models import custom_model
from astropy.modeling.fitting import LevMarLSQFitter
import matplotlib.pyplot as plt
from fastspec import db


def gauss_fit_1p(velo,spec,toWrite=False,source='',paras=None):

    plt.gcf()
    if paras == None:
        plt.clf()
        paras = [0.1, 40., 10.]
        plt.plot(velo,spec)
        plt.xlim(-200,200)
        plt.title(source,size='x-large')
        plt.xlabel('V$_{LSR}$ (km/s)',size='x-large')
        plt.ylabel('Jy/beam',size='x-large')
    else:
        plt.gcf()

    r_s2f = 2.355 # sigma to fwhm: FWHM = 2.355*sigma
    gauss = models.Gaussian1D(amplitude=paras[0],mean=paras[1],stddev=paras[2]/r_s2f)
    gauss.amplitude.min = 0.
    gauss.mean.bounds=[0.,150.]
    gauss.stddev.bounds=[5./r_s2f,40./r_s2f]
    fit = LevMarLSQFitter()
    sp_fit = fit(gauss,velo,spec)
    fpeak = sp_fit.amplitude.value
    vlsr = sp_fit.mean.value
    fwhm = sp_fit.fwhm
    yfit = sp_fit(velo)
    try:
        if fit.fit_info['param_cov'] is not None:
            para_err = np.sqrt(np.diag(fit.fit_info['param_cov']))
        else:
            print('param_cov is None')
            para_err = np.zeros(len(paras))
    except:
        print('unknown fitting error')
        para_err = np.zeros(len(paras))

    print(source)
    print('Flux: {:6.4f}; Vlsr: {:5.1f}; FWHM: {:4.1f}'.format(fpeak, vlsr, fwhm))
    print(para_err)

    for old_line in plt.gca().lines + plt.gca().collections:
        old_line.remove()
    plt.plot(velo,spec,color='C0')
    plt.plot(velo,yfit,color='C3')

    if toWrite:
        plt.savefig(source+'.png',dpi=300,bbox_inches='tight')
    plt.show()

    return yfit, fpeak, vlsr, fwhm, para_err[0], para_err[1], para_err[2]

def gauss_fit_2p(velo,spec,toWrite=False,source='',paras=None):

    @custom_model
    def gaussian_2peak(x, amplitude1=1., mean1=-1., sigma1=1.,
            amplitude2=1., mean2=1.,  sigma2=1.):
        return (amplitude1 * np.exp(-0.5 * ((x - mean1) / sigma1)**2) + 
                amplitude2 * np.exp(-0.5 * ((x - mean2) / sigma2)**2))

    plt.gcf()
    if paras == None:
        plt.clf()
        paras = [0.1, 40., 10., 0.1, 100., 10.]
        plt.plot(velo,spec)
        plt.xlim(-200,200)
        plt.title(source,size='x-large')
        plt.xlabel('V$_{LSR}$ (km/s)',size='x-large')
        plt.ylabel('Jy/beam',size='x-large')
    else:
        plt.gcf()

    r_s2f = 2.35482 # sigma to fwhm: FWHM = 2.355*sigma
    gauss2 = gaussian_2peak(amplitude1=paras[0],mean1=paras[1],sigma1=paras[2]/r_s2f,amplitude2=paras[3],mean2=paras[4],sigma2=paras[5]/r_s2f)
    gauss2.amplitude1.min = 0.
    gauss2.amplitude2.min = 0.
    gauss2.mean1.bounds=[0.,150.]
    gauss2.mean2.bounds=[0.,150.]
    gauss2.sigma1.bounds=[5./r_s2f,40./r_s2f]
    gauss2.sigma2.bounds=[5./r_s2f,40./r_s2f]
    fit = LevMarLSQFitter()
    sp_fit = fit(gauss2,velo,spec,maxiter=50000)
    fpeak1 = sp_fit.amplitude1.value
    vlsr1 = sp_fit.mean1.value
    fwhm1 = sp_fit.sigma1.value*r_s2f
    fpeak2 = sp_fit.amplitude2.value
    vlsr2 = sp_fit.mean2.value
    fwhm2 = sp_fit.sigma2.value*r_s2f

    yfit = sp_fit(velo)
    try:
        if fit.fit_info['param_cov'] is not None:
            para_err = np.sqrt(np.diag(fit.fit_info['param_cov']))
        else:
            print('param_cov is None')
            print('fitting error: ',fit.fit_info['ierr'])
            print(fit.fit_info['message'])
            para_err = np.zeros(len(paras))
    except:
        print('unknown fitting error')
        para_err = np.zeros(len(paras))

    print(source)
    print('Flux1: {:6.4f}; Vlsr1: {:5.1f}; FWHM1: {:4.1f}'.format(fpeak1, vlsr1, fwhm1))
    print('Flux2: {:6.4f}; Vlsr2: {:5.1f}; FWHM2: {:4.1f}'.format(fpeak2, vlsr2, fwhm2))
    print(para_err)

    for old_line in plt.gca().lines + plt.gca().collections:
        old_line.remove()
    plt.plot(velo,spec,color='C0')
    plt.plot(velo,yfit,color='C3')
    if toWrite:
        plt.savefig(source+'.png',dpi=300,bbox_inches='tight')
    plt.show()

    return yfit,fpeak1,vlsr1,fwhm1,fpeak2,vlsr2,fwhm2,para_err[0], para_err[1], para_err[2]*r_s2f, para_err[3], para_err[4], para_err[5]*r_s2f

def main(args):

    db_conn = db.create_connection(args.file_db)

    with db_conn:

        plt.ion()
        plt.figure(figsize=(6,4))

        source_list = db.read_from_table(db_conn,'Source')
        spec_list = db.read_from_table(db_conn,'Spectrum')

        for source in spec_list:

            spec = np.frombuffer(source['Spec'],dtype=np.float32)
            velo = np.frombuffer(source['Velo'],dtype=np.float32)
            
            if spec.shape == velo.shape:
                paras = None
                while(1):
                    try:
                        if paras == None or len(paras) == 3:
                 
                            yfit,fpeak,vlsr,fwhm,fpeak_err,vlsr_err,fwhm_err = gauss_fit_1p(velo,spec,toWrite=args.savefig,source=source['GName'],paras=paras)
                            data1 = (fpeak,fpeak_err,vlsr,vlsr_err,fwhm,fwhm_err,source['GName'])
                            nline = 1
                        elif len(paras) == 6:
                            yfit,fpeak1,vlsr1,fwhm1,fpeak2,vlsr2,fwhm2,fpeak_err1,vlsr_err1,fwhm_err1,fpeak_err2,vlsr_err2,fwhm_err2 = gauss_fit_2p(velo,spec,toWrite=args.savefig,source=source['GName'],paras=paras)
                            data1 = (fpeak1,fpeak_err1,vlsr1,vlsr_err1,fwhm1,fwhm_err1,source['GName'])
                            data2 = (fpeak2,fpeak_err2,vlsr2,vlsr_err2,fwhm2,fwhm_err2,source['GName'])
                            nline = 2
                        paras = eval(input("Input initial paras: amp, vlsr, fwhm \n (1 to keep, 0 to quit)\n>: "))
                    except:
                        print("re-input:")
                        continue 
    
                    if paras == 0:
                        break
                    elif paras == 1:
                        db.update_specfit_yfit(db_conn,yfit,source['GName'])
                        db.update_specfit_para(db_conn,data1,table='SpFitHa',line_id=1)
                        if nline == 2:
                            db.update_specfit_para(db_conn,data2,table='SpFitHa',line_id=2)
                        break
        plt.close(fig='all')

    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_db',  type=str,help='The sqlite database file')
    parser.add_argument('--savefig',  action='store_true',help='save figure')
    args = parser.parse_args()
    main(args)
