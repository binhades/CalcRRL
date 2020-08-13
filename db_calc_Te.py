#! /home/bliu/anaconda3/bin/python

# Filename: db_fit_spec.py
# Aim: to load the data base and then to fit spectra.

import argparse, sqlite3
import numpy as np
from fastspec import db

def calc_te(T_c,T_l,d_v,freq=1.3,r_he2h=0.1):

    T_e = np.power(7103.3* np.power(freq,1.1)*(T_c/T_l)/d_v/(1+r_he2h),0.87)

    return T_e

def load_data(conn):
    try:
        cur = conn.cursor()
        cur.execute('''SELECT 
                            SpFitHa.GName,
                            Flux_peak1,
                            FWHM1,
                            Tb_vgps
                        FROM 
                            SpFitHa
                        LEFT JOIN PhyPara ON 
                            SpFitHa.GName = PhyPara.GName
                        ORDER BY
                            SpFitHa.GName ''')
        data_list = cur.fetchall()

    except sqlite3.Error as err:
        print(err)
        data_list = None

    return data_list

def main(args):

    db_conn = db.create_connection(args.file_db)

    with db_conn:

        data_list = load_data(db_conn)
        for source in data_list:

            t_c = source['Tb_vgps']
            t_l = source['Flux_peak1']*12/2.
            d_v = source['FWHM1'] 

            try:
                t_e = calc_te(t_c,t_l,d_v)
                if t_l > 0.5:
                    print(source['GName']+'  -  Te {:5.0f} K  -  Tc {:4.0f} K  -  Tl {:6.3f} K'.format(t_e,t_c,t_l))
            except:
                continue

    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_db',  type=str,help='The sqlite database file')
    args = parser.parse_args()
    main(args)
