#! /home/bliu/anaconda3/bin/python

# Filename: db_fit_spec.py
# Aim: to load the data base and then to fit spectra.

import argparse, sqlite3
import numpy as np
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import pyrrl as rrl

def spec_fit(velo,spec,toPlot=False,toWrite=False,source='',paras=None):

    if paras == None:
        paras = [0.1, 40., 10.]

    r_s2f = 2.355 # sigma to fwhm: FWHM = 2.355*sigma
    gauss = models.Gaussian1D(amplitude=paras[0],mean=paras[1],stddev=paras[2]/r_s2f)

    gauss.amplitude.min = 0.
    gauss.mean.bounds=[0.,150.]
    gauss.stddev.bounds=[5./r_s2f,40./r_s2f]

    fit = fitting.LevMarLSQFitter()

    sp_fit = fit(gauss,velo,spec)
    fpeak = sp_fit.amplitude.value
    vlsr = sp_fit.mean.value
    fwhm = sp_fit.fwhm
#    print(fit.fit_info['param_cov'])
#    print(fit.fit_info['cov_x'])

    para_err = np.sqrt(np.diag(fit.fit_info['param_cov']))
    print(para_err)
    yfit = sp_fit(velo)

    if toPlot:
        plt.figure(figsize=(6,4))
        plt.plot(velo,spec)
        plt.plot(velo,yfit)
        plt.xlim(-200,200)
        plt.title(source,size='x-large')
        plt.xlabel('V$_{LSR}$ (km/s)',size='x-large')
        plt.ylabel('Jy/beam',size='x-large')
        if toWrite:
            plt.savefig(source+'.png',dpi=300,bbox_inches='tight')
        else:
            plt.show()

        plt.close(fig='all')
    
    return yfit,fpeak, vlsr, fwhm, para_err[0], para_err[1], para_err[2]

def main(args):

    db_conn = rrl.db.create_connection(args.file_db)

    with db_conn:

        source_list = rrl.db.read_from_table(db_conn,'Source')
        spec_list = rrl.db.read_from_table(db_conn,'Spectrum')

        for source in spec_list:

            spec = np.frombuffer(source['Spec'],dtype=np.float32)
            velo = np.frombuffer(source['Velo'],dtype=np.float32)
            
            if spec.shape == velo.shape:
            
                yfit,fpeak,vlsr,fwhm,fpeak_err,vlsr_err,fwhm_err = spec_fit(velo,spec,toPlot=args.plot,toWrite=args.write,source=source['GName'])

                rrl.db.add_to_specfit(db_conn,(source['GName'],fpeak,fpeak_err,vlsr,vlsr_err,fwhm,fwhm_err),table='SpFitHa')
                rrl.db.update_specfit_yfit(db_conn,yfit,source['GName'])

    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_db',  type=str,help='The sqlite database file')
    parser.add_argument('--plot',  action='store_true',help='if set, plot the spectra with fitting')
    parser.add_argument('--write',  action='store_true',help='only works when "--plot" is set, write the plots into files')
    args = parser.parse_args()
    main(args)
