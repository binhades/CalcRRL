#! /home/bliu/anaconda3/bin/python

# Filename: db_fit_spec.py
# Aim: to load the data base and then to fit spectra.

import argparse
import numpy as np
import sqlite3
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt

def db_create_connection(db_file):

    conn = None
    try:
        conn = sqlite3.connect(db_file)
        conn.row_factory = sqlite3.Row
    except sqlite3.Error as err:
        print(err)

    return conn

def db_read_sources(conn):

    cur = conn.cursor()
    cur.execute('select * from catalog')

    return cur.fetchall()
   
def db_update_fitpara(conn,paras):

    sql = ''' UPDATE catalog
              SET Flux_peak = ? ,
                  Flux_peak_err = ? ,
                  Velo_lsr = ? ,
                  Velo_lsr_err = ? ,
                  FWHM = ? ,
                  FWHM_err = ?
              WHERE GName = ? '''

    cur = conn.cursor()
    cur.execute(sql,paras)
    conn.commit()

def spec_fit(velo,spec,toPlot=False,toWrite=False,source='',size=0,paras=None):

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

    if toPlot:
        plt.figure(figsize=(6,4))
        plt.plot(velo,spec)
        plt.plot(velo,sp_fit(velo))
        plt.xlim(-200,200)
        plt.title(source+' - Radius {:3.0f}"'.format(size),size='x-large')
        plt.xlabel('V$_{LSR}$ (km/s)',size='x-large')
        plt.ylabel('Jy/beam',size='x-large')
        if toWrite:
            plt.savefig(source+'.png',dpi=300,bbox_inches='tight')
        else:
            plt.show()

        plt.close(fig='all')
    
    return fpeak, vlsr, fwhm, para_err[0], para_err[1], para_err[2]

def main(args):

    db_conn = db_create_connection(args.file_db)

    with db_conn:

        source_list = db_read_sources(db_conn)

        for source in source_list:

            spec = np.frombuffer(source['Spectrum'],dtype=np.float32)
            velo = np.frombuffer(source['Velocity'],dtype=np.float32)
            
            if spec.shape == velo.shape:
            
                fpeak,vlsr,fwhm,fpeak_err,vlsr_err,fwhm_err = spec_fit(velo,spec,toPlot=args.plot,toWrite=args.write,source=source['GName'],size=source['Radius'])

                db_update_fitpara(db_conn,(fpeak,fpeak_err,vlsr,vlsr_err,fwhm,fwhm_err,source['GName']))

    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_db',  type=str,help='The sqlite database file')
    parser.add_argument('--plot',  action='store_true',help='if set, plot the spectra with fitting')
    parser.add_argument('--write',  action='store_true',help='only works when "--plot" is set, write the plots into files')
    args = parser.parse_args()
    main(args)
